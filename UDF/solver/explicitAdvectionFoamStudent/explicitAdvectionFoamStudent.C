/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    explicitAdvectionFoamStudent

Description
    Student implementation of the explicit finite-volume linear-advection
    solver, developed one numerical stage at a time.

    Stage 0 reads the mesh and fields.  Stage 1 computes the face volume flux.
    Stage 2 computes a CFL candidate time step.  Stage 3 computes the explicit
    convection residual.  Stage 4 performs a forward-Euler update.  Stage 5
    repeats the update in a complete time loop and writes scheduled results.
    Error analysis remains a later stage.
\*---------------------------------------------------------------------------*/

// argList.H 提供 OpenFOAM 程序的命令行解析能力。
// 所有 OpenFOAM solver 基本都会先经过 setRootCase.H，
// 由它读取 -case、-parallel、-help 等命令行参数。
#include "argList.H"

// dictionary.H 提供 dictionary 类及其查值接口，例如：
//
//     controlDict.lookupOrDefault<scalar>("maxCo", 0.2)
//
// 注意：maxCo 不是一个需要单独链接的 OpenFOAM 库；
// 它只是 system/controlDict 中的一个字典关键字。
#include "dictionary.H"

// dimensionedScalar.H 提供带物理量纲的标量类型。
// Stage 4 用它把普通数值 deltaT 标记成时间量，
// 从而让 deltaT*residual 的量纲检查能够通过。
#include "dimensionedScalar.H"

// volFields.H 提供体场类型：
//   volVectorField: 每个单元一个 vector，比如速度 U；
//   volScalarField: 每个单元一个 scalar，比如被输运标量 T。
#include "volFields.H"

// surfaceFields.H 提供定义在网格面上的场类型，
// 本阶段使用的 surfaceScalarField 就属于这一类。
#include "surfaceFields.H"

// fvcFlux.H 声明 fvc::flux(U)：
//
//     U (cell-centred velocity)
//         -> interpolate/evaluate U on faces
//         -> dot with face area vector Sf
//         -> face scalar flux
//
// 它对应有限体积公式：
//
//     F_f = U_f dot S_f
#include "fvcFlux.H"

// fvcSurfaceIntegrate.H 提供 fvc::surfaceSum(...)：
//
//     面场
//         -> 把每个 cell 周围的面值相加
//         -> 单元内部场
//
// Stage 2 用它计算：
//
//     sum_f(|F_cf|)
#include "fvcSurfaceIntegrate.H"

// fvcDiv.H 提供 fvc::div(...)：
//
//     面通量 + 单元场 + 离散格式名字
//         -> 体积归一化后的单元残差
//
// Stage 3 用它计算显式对流残差：
//
//     R_c = (1/V_c) sum_f F_cf T_f
#include "fvcDiv.H"

using namespace Foam;


int main(int argc, char *argv[])
{
    // setRootCase.H:
    //   1. 解析命令行参数；
    //   2. 定位当前 OpenFOAM case；
    //   3. 支持 `solver -case path/to/case` 这种调用方式。
    //
    // 注意：这不是普通函数调用，而是 OpenFOAM 的 include 式代码片段。
    // 它会在当前位置展开成若干行 C++ 代码。
    #include "setRootCase.H"

    // createTime.H:
    //   创建 Time runTime 对象。
    //
    // runTime 负责管理：
    //   - 当前时间目录，比如 0、0.05、1；
    //   - controlDict；
    //   - endTime、deltaT、writeInterval；
    //   - 什么时候写出结果。
    //
    // 在本机 OpenFOAM 14 中，当前时间目录名用 runTime.name() 取得；
    // Time::timeName(...) 主要是把 scalar 时间值格式化成目录名的静态函数。
    #include "createTime.H"

    // createMesh.H:
    //   从 constant/polyMesh 读取 fvMesh mesh。
    //
    // mesh 是有限体积法的几何和拓扑基础，里面有：
    //   - 单元数量 mesh.nCells()；
    //   - 单元体积 mesh.V()；
    //   - 面面积向量 mesh.Sf()；
    //   - owner/neighbour 拓扑；
    //   - 边界 patch 信息。
    #include "createMesh.H"

    // controlDict 是 system/controlDict 对应的字典对象。
    // 后续 maxCo、velocityField、advectedField 都从这里读。
    const dictionary& controlDict = runTime.controlDict();

    // velocityField 是一个可选关键字。
    // 如果 controlDict 里写了：
    //
    //     velocityField U;
    //
    // 就读取名为 U 的速度场；如果没写，就默认读 U。
    // 这样做的好处是：以后你想把速度场命名成 Umean、Uc 等，也不用改 C++。
    const word velocityName
    (
        controlDict.lookupOrDefault<word>("velocityField", "U")
    );

    // advectedField 同理，默认读取 T。
    // 在数学公式里它是 phi；在本工程里为了和 OpenFOAM 标量输运习惯接近，叫 T。
    const word advectedName
    (
        controlDict.lookupOrDefault<word>("advectedField", "T")
    );

    // Info 是 OpenFOAM 的标准输出流，类似 std::cout，
    // 但它会配合 OpenFOAM 的并行和日志机制。
    //
    // 这一段先确认程序真的进来了、读到了哪个 case、当前时间是多少、
    // 网格有多少个单元，以及准备读取哪些字段。
    Info<< "Student checkpoint: reading fields" << nl
        << "  case       = " << runTime.caseName() << nl
        << "  time       = " << runTime.name() << nl
        << "  cells      = " << mesh.nCells() << nl
        << "  velocity   = " << velocityName << nl
        << "  advected   = " << advectedName << nl
        << endl;

    // 读取速度场 U。
    //
    // volVectorField 表示“定义在单元上的向量场”。
    //
    // IOobject 这一段说明怎么从磁盘读写：
    //   velocityName        -> 字段文件名，默认是 U；
    //   runTime.name()      -> 当前时间目录，本阶段是 0；
    //   mesh                -> 这个字段属于哪个网格；
    //   MUST_READ           -> 文件必须存在，否则直接报错；
    //   NO_WRITE            -> 速度场在本题中给定，不由求解器写出。
    //
    // 这个对象对应数学公式里的速度：
    //
    //     boldsymbol{U}
    volVectorField U
    (
        IOobject
        (
            velocityName,
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // 读取被输运标量 T。
    //
    // volScalarField 表示“定义在单元上的标量场”。
    //
    // AUTO_WRITE 表示：当 runTime.write() 被调用时，T 会自动写到新的时间目录。
    // 后续真正实现时间推进后，输出目录中的 T 就是数值解。
    //
    // 这个对象对应题目数学公式中的标量未知量 phi。
    // 为了避免和 OpenFOAM 约定的面通量字段 phi 混淆，
    // 本工程把它命名为 T。也就是说：
    //
    //     mathematical phi  ->  code T
    //
    // 对应方程为：
    //
    //     partial_t T + div(U*T) = 0
    volScalarField T
    (
        IOobject
        (
            advectedName,
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // 打印量纲和范围。
    //
    // U 的量纲应该是 [0 1 -1 0 0 0 0]，即 m/s。
    // T 是无量纲标量，所以量纲是 [0 0 0 0 0 0 0]。
    //
    // min(T)、max(T) 用来检查初始场是否成功读入。
    // 对正弦波案例，范围应接近 [-1, 1]。
    Info<< "  U dimensions = " << U.dimensions() << nl
        << "  T dimensions = " << T.dimensions() << nl
        << "  T min       = " << min(T).value() << nl
        << "  T max       = " << max(T).value() << nl
        << endl;

    /*
     * Stage 1: face volume flux
     *
     * Finite-volume formula:
     *
     *     F_f = U_f dot S_f
     *
     * Here:
     *   - U is a volVectorField: one vector velocity per cell;
     *   - phi is a surfaceScalarField: one scalar flux per face;
     *   - fvc::flux(U) performs the face evaluation and dot product.
     *
     * The PDE unknown is called T in this project.  The name phi follows
     * the usual OpenFOAM convention for the face volume-flux field, so it
     * must not be confused with the mathematical scalar phi in the problem
     * statement.
     *
     * phi is calculated from U and mesh geometry at run time.  It is not an
     * initial field that must be read from 0/phi, hence NO_READ.  It is also
     * not written at this checkpoint because later stages have not started
     * a time loop yet, hence NO_WRITE.
     */
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux(U)
    );

    // min(phi) and max(phi) reduce all face values to two global scalars.
    //
    // For U = (1, 1, 0), the current mesh should produce both positive and
    // negative face fluxes because opposite face normals have opposite
    // orientations.  The magnitude is expected to be of order 1e-3 to 1e-2.
    Info<< "  phi min     = " << min(phi).value() << nl
        << "  phi max     = " << max(phi).value() << nl
        << endl;

    /*
     * Stage 2: CFL time step
     *
     * Mathematical definition:
     *
     *     Co_max
     *       = 0.5*deltaT*max_c
     *         (sum_f(|F_cf|)/V_c)
     *
     * We first calculate the cell-wise rate:
     *
     *     rate_c = sum_f(|F_cf|)/V_c
     *
     * and then choose deltaT so that the largest cell Courant number is
     * equal to the target maxCo:
     *
     *     deltaT = 2*maxCo/max(rate_c)
     *
     * This stage only calculates deltaT.  It does not update T and does not
     * advance runTime.
     */

    // Read the target Courant number from system/controlDict.
    //
    // If the case does not provide maxCo, use 0.2 as the default target.
    const scalar maxCo
    (
        controlDict.lookupOrDefault<scalar>("maxCo", 0.2)
    );

    if (maxCo <= 0)
    {
        FatalErrorInFunction
            << "maxCo must be positive, but received " << maxCo
            << exit(FatalError);
    }

    // `mag(phi)` is still a surfaceScalarField:
    //
    //     {F_f} -> {|F_f|}
    //
    // `fvc::surfaceSum(...)` changes the location of the result from faces
    // to cells by adding the values on the faces surrounding each cell:
    //
    //     sumPhi_c = sum_f(|F_cf|)
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    // Divide the face-flux sum by the cell volume:
    //
    //     rate_c = sum_f(|F_cf|)/V_c
    //
    // The resulting rate has units of 1/s.
    scalarField rate
    (
        sumPhi/mesh.V().primitiveField()
    );

    // Find the largest rate over the complete mesh.
    //
    // `gMax` means global maximum.  In a parallel calculation it also
    // performs the reduction across MPI processes.
    const scalar rateMax = gMax(rate);

    if (rateMax <= SMALL)
    {
        FatalErrorInFunction
            << "The maximum flux rate is zero or too small: "
            << rateMax << nl
            << "Cannot determine a positive CFL time step."
            << exit(FatalError);
    }

    // Rearrange Co_max = 0.5*deltaT*rateMax:
    //
    //     deltaT = 2*maxCo/rateMax
    const scalar deltaT = 2.0*maxCo/rateMax;

    Info<< "  rate max    = " << rateMax << nl
        << "  maxCo       = " << maxCo << nl
        << "  CFL deltaT  = " << deltaT << nl
        << endl;

    /*
     * Stage 3: explicit convection residual
     *
     * 这一阶段把旧时间层的 T 送入显式对流算子，得到体积归一化后的
     * 单元残差场。
     *
     * Formula:
     *
     *     R_c^n = (1/V_c) sum_f F_cf*T_f^n
     *
     * OpenFOAM interface:
     *
     *     fvc::div(phi, T, "div(phi,T)")
     *
     * The face interpolation is selected in fvSchemes:
     *
     *     div(phi,T) Gauss upwind;
     *
     * 中文提示：
     *   - `fvc` 表示显式 finite-volume calculus；
     *   - `fvm` 表示隐式 matrix 组装；
     *   - 本题要求显式，所以这里应使用 `fvc::div`；
     *   - `"div(phi,T)"` 这个字符串要和 system/fvSchemes 里的名字一致。
     *   - `fvc::div` 返回的已经是除过体积的残差。
     */
    tmp<volScalarField> tResidual
    (
        fvc::div(phi, T, "div(phi,T)")
    );

    const volScalarField& residual = tResidual();

    const scalar residualIntegral
    (
        gSum(residual.primitiveField()*mesh.V().primitiveField())
    );

    Info<< "  residual dimensions = " << residual.dimensions() << nl
        << "  residual min        = " << min(residual).value() << nl
        << "  residual max        = " << max(residual).value() << nl
        << "  residual integral   = " << residualIntegral << nl
        << endl;

    /*
     * Stage 5: complete explicit time loop
     *
     * The Stage 4 operation
     *
     *     T^(n+1) = T^n - deltaT*R^n
     *
     * is now repeated until endTime.
     *
     * The velocity field U is fixed in this teaching case, so phi is also
     * fixed and can be calculated once before the loop.  The residual must
     * still be recalculated every time step because T changes.
     *
     * We deliberately write the loop in an explicit educational form:
     *
     *     while (current time < end time)
     *         choose a safe step;
     *         set runTime deltaT;
     *         advance runTime to t^(n+1);
     *         calculate R from the old T;
     *         update T;
     *         correct boundary conditions;
     *         write according to controlDict.
     *
     * The last step uses min(deltaT, remainingTime), so the solver reaches
     * endTime without overshooting it.
     */
    const scalar endTime = runTime.endTime().value();
    const scalar timeTolerance
    (
        1.0e-12*max(1.0, mag(endTime))
    );
    const scalar startTime = runTime.value();
    label step = 0;

    Info<< "Starting Stage 5 time loop" << nl
        << "  start time   = " << runTime.value() << nl
        << "  end time     = " << endTime << nl
        << "  target deltaT = " << deltaT << nl
        << endl;

    while (endTime - runTime.value() > timeTolerance)
    {
        /*
         * Do not advance time by repeatedly adding deltaT to the previous
         * floating-point value.  That accumulation can eventually produce a
         * time name which OpenFOAM cannot distinguish from its rounded form.
         *
         * Instead, calculate the next target from the integer step number:
         *
         *     t_target = min(t_start + (step+1)*deltaT_CFL, endTime)
         *
         * The integer multiplication prevents long-run drift, and the final
         * min() makes the last step land exactly on endTime.
         */
        const scalar oldTime = runTime.value();
        const scalar targetTime
        (
            min
            (
                startTime + (step + 1)*deltaT,
                endTime
            )
        );

        // The actual step is the difference between the target and old time.
        const scalar stepDeltaT = targetTime - oldTime;

        if (stepDeltaT <= SMALL)
        {
            FatalErrorInFunction
                << "The remaining time is too small for another step: "
                << (endTime - oldTime)
                << exit(FatalError);
        }

        // Store the selected step in OpenFOAM's Time object.
        runTime.setDeltaT(stepDeltaT);

        // Advance to the precomputed target time.  T still contains T^n here.
        ++step;
        runTime.setTime(targetTime, step);

        Info<< "Time = " << runTime.name()
            << "  step = " << step
            << "  deltaT = " << stepDeltaT
            << "  maxCo = " << 0.5*stepDeltaT*rateMax
            << nl;

        /*
         * Stage 3 is repeated here because the advected field T has changed.
         *
         * The residual is still calculated from the old field value at the
         * beginning of this update:
         *
         *     R_c^n = (1/V_c) sum_f F_cf*T_f^n
         */
        tmp<volScalarField> tResidual
        (
            fvc::div(phi, T, "div(phi,T)")
        );

        const volScalarField& residual = tResidual();

        const scalar residualIntegral
        (
            gSum(residual.primitiveField()*mesh.V().primitiveField())
        );

        /*
         * Stage 4 is the forward-Euler update.
         *
         * deltaTDim supplies the time dimension, while residual has units
         * of 1/time.  Their product has the same dimensions as T.
         */
        const dimensionedScalar deltaTDim
        (
            "deltaT",
            dimTime,
            stepDeltaT
        );

        T = T - deltaTDim*residual;
        T.correctBoundaryConditions();

        Info<< "  residual min      = " << min(residual).value() << nl
            << "  residual max      = " << max(residual).value() << nl
            << "  residual integral = " << residualIntegral << nl
            << "  T min             = " << min(T).value() << nl
            << "  T max             = " << max(T).value() << nl
            << endl;

        // AUTO_WRITE on T and writeInterval in controlDict decide whether
        // this intermediate time level is written to disk.
        runTime.write();
    }

    /*
     * The regular runTime.write() call follows writeControl/writeInterval.
     * Therefore the last, shortened time step may reach endTime without
     * creating an output directory at that exact time.
     *
     * Force one final write after the loop:
     *
     *     t_final = endTime
     *
     * This is important for verification cases whose requested output time
     * is a mathematical value such as t = 2*pi.
     */
    const bool finalWriteOK = runTime.writeNow();

    if (!finalWriteOK)
    {
        FatalErrorInFunction
            << "Failed to write the final field at time "
            << runTime.name()
            << exit(FatalError);
    }

    Info<< "Stage 5 time loop completed." << nl
        << "  final time = " << runTime.value() << nl
        << "  time steps = " << step << nl
        << "  final field written at = " << runTime.name() << nl
        << "  T min      = " << min(T).value() << nl
        << "  T max      = " << max(T).value() << nl
        << "End" << nl
        << endl;

    return 0;
}

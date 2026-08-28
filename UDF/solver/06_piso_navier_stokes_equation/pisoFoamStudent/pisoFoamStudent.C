/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    pisoFoamStudent

Description
    第六个学生版求解器：不可压 Navier-Stokes 方程的 PISO 方法。

    本求解器用于第五题的第二种算法。它不是 projectionFoamStudent
    的改名版本，而是按 OpenFOAM 的 PISO 离散关系编写：

        A(U) U = H(U) - grad(p)

        rAU = 1/A(U)

        HbyA = rAU H(U)

        phiHbyA = flux(HbyA) + interpolate(rAU) ddtCorr(U, phi)

        div(phiHbyA) - laplacian(rAU, p) = 0

        phi = phiHbyA - flux(laplacian(rAU, p))

        U = HbyA - rAU grad(p)

    其中 p 是运动学压力 p/rho，nu 是运动黏度。空间离散格式不在
    本文件硬编码，而由 cases/<case>/system/fvSchemes 选择；线性求解器
    和 PISO 校正次数由 cases/<case>/system/fvSolution 选择。
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "dictionary.H"
#include "pisoControl.H"
#include "pressureReference.H"
#include "findRefCell.H"
#include "constrainPressure.H"
#include "constrainHbyA.H"
#include "adjustPhi.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"

#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcGrad.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

using namespace Foam;


int main(int argc, char *argv[])
{
    // 对应数学公式：
    //     无。这一步只确定当前 OpenFOAM case 的根路径。
    //
    // cases 入口：
    //     命令行：
    //         pisoFoamStudent -case <caseDir>
    //
    // OpenFOAM 接口：
    //     setRootCase.H
    //
    // 内部面/边界面：
    //     还没有读取网格，因此这里不区分内部面和边界面。
    //
    // 要改哪里：
    //     一般不用改。若换应用名，只改 Make/files 的 EXE 和
    //     cases/<case>/system/controlDict/application。
    #include "setRootCase.H"

    // 对应数学公式：
    //     t = t0, t1, ..., tn
    //
    //     Δt = t(n+1) - t(n)
    //
    //     计算从 startTime 推进到 endTime。
    //
    // cases 入口：
    //     system/controlDict:
    //         startTime
    //         endTime
    //         deltaT
    //         writeControl
    //         writeInterval
    //         maxCo
    //         steadyStateControl
    //         steadyVelocityTol
    //         steadyMassTol
    //         minimumSteadySteps
    //         requiredSteadySteps
    //
    // OpenFOAM 接口：
    //     createTime.H -> Time runTime
    //
    // 内部面/边界面：
    //     时间对象只控制步进，不区分面类型。
    //
    // 要改哪里：
    //     改时间步、终止时间、写出频率和稳态判据时，改 system/controlDict。
    #include "createTime.H"

    // 对应数学公式：
    //     对不可压流动，OpenFOAM 用面体积通量 Φ_f = U_f · S_f
    //     计算最大 Courant 数：
    //
    //         Co_max = max_c [ Δt / (2 V_c) Σ_f |Φ_f| ]
    //
    //     若开启自适应时间步，则下一步时间步取：
    //
    //         Δt_new = min(1.2 Δt_old, maxCo / Co_max Δt_old, maxDeltaT)
    //
    // cases 入口：
    //     system/controlDict:
    //         adjustTimeStep
    //         maxCo
    //         maxDeltaT
    //         deltaT
    //
    // OpenFOAM 接口：
    //     createTimeControls.H 读取 adjustTimeStep、maxCo、maxDeltaT。
    //     CourantNo.H 根据 phi、mesh.V() 和当前 Δt 计算 Co。
    //     本求解器用同样的限制公式调用 runTime.setDeltaT(...) 更新 Δt。
    //
    // 内部面/边界面：
    //     Σ_f |Φ_f| 同时累加单元周围的内部面和边界面通量；
    //     对封闭腔流，壁面法向边界通量理论上为 0，
    //     但内部面通量仍决定速度场中的 CFL 限制。
    //
    // 要改哪里：
    //     想固定时间步：system/controlDict/adjustTimeStep false。
    //     想放宽或收紧稳定限制：system/controlDict/maxCo。
    //     想限制最大增长幅度：system/controlDict/maxDeltaT。
    #include "createTimeControls.H"

    // 对应数学公式：
    //     计算区域：
    //         Ω = union of all cells Ω_c
    //
    //     单元体积：
    //         V_c = volume(Ω_c)
    //
    //     面面积向量：
    //         S_f = |S_f| n_f
    //
    //     内部面方向：
    //         f connects owner O_f and neighbour N_f
    //
    //     对某个单元 c 的外法向符号：
    //         sigma_cf = +1 if c = O_f
    //         sigma_cf = -1 if c = N_f
    //
    // cases 入口：
    //     constant/polyMesh:
    //         points
    //         faces
    //         owner
    //         neighbour
    //         boundary
    //
    //     规则网格入口：
    //         system/blockMeshDict
    //
    //     Gmsh/混合非结构网格入口：
    //         mesh/mesh.msh
    //         system/createPatchDict
    //
    // OpenFOAM 接口：
    //     createMesh.H -> fvMesh mesh
    //
    // 内部面/边界面：
    //     内部面由 owner/neighbour 同时存在来定义；
    //     边界面由 constant/polyMesh/boundary 中的 patch 定义。
    //     本求解器后续调用 fvm/fvc 算子时，OpenFOAM 会自动遍历这两类面。
    //
    // 要改哪里：
    //     改区域、网格类型、边界 patch 名称时，改网格生成脚本或
    //     system/blockMeshDict、system/createPatchDict。
    #include "createMesh.H"

    // 对应数学公式：
    //     PISO 在同一个时间步内做多次压力--速度校正：
    //
    //         U* -> p(1) -> U(1) -> p(2) -> U(2) -> ...
    //
    //     非正交循环处理压力 Laplacian 中的非正交修正：
    //
    //         grad(p)_f · S_f
    //       = orthogonal part + non-orthogonal correction
    //
    // cases 入口：
    //     system/fvSolution/PISO:
    //         nCorrectors
    //         nNonOrthogonalCorrectors
    //         momentumPredictor
    //         pRefCell
    //         pRefValue
    //
    // OpenFOAM 接口：
    //     pisoControl piso(mesh)
    //     piso.correct()
    //     piso.correctNonOrthogonal()
    //     piso.finalNonOrthogonalIter()
    //
    // 内部面/边界面：
    //     PISO 控制器不直接遍历面；它控制后面的压力方程何时装配、
    //     何时用最后一次非正交校正更新守恒通量。
    //
    // 要改哪里：
    //     改 PISO 校正次数时，改 system/fvSolution/PISO/nCorrectors。
    //     改非正交修正次数时，改
    //     system/fvSolution/PISO/nNonOrthogonalCorrectors。
    pisoControl piso(mesh);

    // 对应数学公式：
    //     ν = μ/rho
    //
    //     不可压常黏度 Navier-Stokes 方程：
    //
    //         div(U) = 0
    //
    //         dU/dt + div(U U) = -grad(p) + div(nu grad(U))
    //
    //     p 是运动学压力，即真实压力除以密度：
    //
    //         p = p_physical / rho
    //
    // cases 入口：
    //     constant/physicalProperties:
    //         nu <value>;
    //
    //     0/U:
    //         internalField
    //         boundaryField
    //
    //     0/p:
    //         internalField
    //         boundaryField
    //
    // OpenFOAM 接口：
    //     IOdictionary
    //     dimensionedScalar
    //     volVectorField U
    //     volScalarField p
    //     createPhi.H -> surfaceScalarField phi
    //
    // 内部面/边界面：
    //     U 和 p 是单元中心体场，同时带有 boundaryField。
    //     phi 是面通量场，内部面保存 owner 到 neighbour 的通量，
    //     边界面保存 patch 外法向通量。
    //
    // 要改哪里：
    //     改 Reynolds 数时改 constant/physicalProperties/nu。
    //     改初始速度时改 0/U/internalField。
    //     改顶盖速度和无滑移壁面时改 0/U/boundaryField。
    //     改压力边界条件时改 0/p/boundaryField。
    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        "nu",
        dimKinematicViscosity,
        physicalProperties.lookup("nu")
    );

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"
    mesh.schemes().setFluxRequired(p.name());

    // 对应数学公式：
    //     压力只通过 grad(p) 进入速度方程；
    //     如果所有压力边界都是零法向梯度，则压力可以整体加常数。
    //
    //     因此需要固定：
    //         p(c_ref) = p_ref
    //
    // cases 入口：
    //     system/fvSolution/PISO:
    //         pRefCell
    //         pRefValue
    //
    // OpenFOAM 接口：
    //     setRefCell(p, piso.dict(), pRefCell, pRefValue)
    //     fvScalarMatrix::setReference(pRefCell, pRefValue)
    //
    // 内部面/边界面：
    //     压力参考值作用在一个单元中心，不是边界条件。
    //
    // 要改哪里：
    //     若压力有 fixedValue 边界，一般不需要参考单元；
    //     若都是 zeroGradient，则保留 pRefCell/pRefValue。
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, piso.dict(), pRefCell, pRefValue);

    // 对应数学公式：
    //     稳态监测量：
    //
    //         max_c |U_c(n+1) - U_c(n)| <= eps_U
    //
    //         max_c |div(phi)_c| <= eps_mass
    //
    // cases 入口：
    //     system/controlDict:
    //         steadyStateControl
    //         steadyVelocityTol
    //         steadyMassTol
    //         minimumSteadySteps
    //         requiredSteadySteps
    //
    // OpenFOAM 接口：
    //     runTime.controlDict().lookupOrDefault
    //
    // 内部面/边界面：
    //     稳态速度变化量来自单元中心；
    //     质量残差来自所有面通量对每个单元的散度。
    //
    // 要改哪里：
    //     改是否提前停止、稳态阈值和连续满足步数时，改 system/controlDict。
    volVectorField UPrevious
    (
        IOobject
        (
            "UPrevious",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );

    const dictionary& controlDict = runTime.controlDict();
    const Switch steadyStateControl
    (
        controlDict.lookupOrDefault<Switch>("steadyStateControl", true)
    );
    const scalar steadyVelocityTol
    (
        controlDict.lookupOrDefault<scalar>("steadyVelocityTol", 1e-6)
    );
    const scalar steadyMassTol
    (
        controlDict.lookupOrDefault<scalar>("steadyMassTol", 1e-8)
    );
    const label minimumSteadySteps
    (
        controlDict.lookupOrDefault<label>("minimumSteadySteps", 1000)
    );
    const label requiredSteadySteps
    (
        controlDict.lookupOrDefault<label>("requiredSteadySteps", 20)
    );
    label steadyStepCount = 0;

    Info<< "\nStarting pisoFoamStudent\n"
        << "  case = " << runTime.caseName() << nl
        << "  cells = " << mesh.nCells() << nl
        << "  nu = " << nu.value() << nl
        << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        UPrevious = U;

        // 对应数学公式：
        //     本时间步开始时，用上一时间步守恒面通量 Φ_f 估计：
        //
        //         Co_max = max_c [ Δt / (2 V_c) Σ_f |Φ_f| ]
        //
        //     若 adjustTimeStep = true，则按 maxCo 自动更新 Δt。
        //
        // cases 入口：
        //     system/controlDict/adjustTimeStep
        //     system/controlDict/maxCo
        //     system/controlDict/maxDeltaT
        //
        // OpenFOAM 接口：
        //     readTimeControls.H 支持 runTimeModifiable 时重读控制参数；
        //     CourantNo.H 输出 meanCo 和 maxCo；
        //     runTime.setDeltaT(...) 更新本时间步使用的 Δt。
        //
        // 内部面/边界面：
        //     CourantNo.H 调用 fvc::surfaceSum(mag(phi))，
        //     OpenFOAM 会按 mesh 的 owner/neighbour/boundary 自动汇总。
        //
        // 要改哪里：
        //     速度过快导致发散时降低 maxCo；
        //     运行太慢且 Co 很小时提高 maxDeltaT 或初始 deltaT。
        #include "readTimeControls.H"
        #include "CourantNo.H"

        if (adjustTimeStep)
        {
            scalar newDeltaT = 1.2*runTime.deltaTValue();

            if (CoNum > small)
            {
                newDeltaT = min
                (
                    newDeltaT,
                    maxCo/CoNum*runTime.deltaTValue()
                );
            }

            newDeltaT = min(newDeltaT, maxDeltaT);
            runTime.setDeltaT(newDeltaT);

            Info<< "deltaT = " << runTime.deltaTValue() << endl;
        }

        // 对应数学公式：
        //     动量预测方程：
        //
        //         dU/dt + div(phi, U) - div(nu grad(U)) = -grad(p)
        //
        //     有限体积离散后写成矩阵：
        //
        //         A(U) U = H(U) - grad(p)
        //
        //     对控制体 Ω_c：
        //
        //         V_c (U_c(n+1) - U_c(n))/Δt
        //       + sum_f phi_f U_f
        //       - sum_f (nu grad(U))_f · S_cf
        //       = - V_c grad(p)_c
        //
        // cases 入口：
        //     system/fvSchemes:
        //         ddtSchemes/default
        //         divSchemes/div(phi,U)
        //         laplacianSchemes/laplacian(nu,U)
        //         gradSchemes/grad(p)
        //         interpolationSchemes/default
        //         snGradSchemes/default
        //
        //     system/fvSolution/solvers/U:
        //         solver
        //         smoother
        //         tolerance
        //         relTol
        //
        //     system/fvSolution/PISO:
        //         momentumPredictor
        //
        // OpenFOAM 接口：
        //     fvm::ddt(U)
        //     fvm::div(phi, U)
        //     fvm::laplacian(nu, U)
        //     fvc::grad(p)
        //     fvVectorMatrix UEqn
        //
        // 内部面/边界面：
        //     内部面对流和扩散由 owner/neighbour 两侧单元装配；
        //     边界面由 0/U/boundaryField 和 0/p/boundaryField 提供面值或梯度。
        //
        // 要改哪里：
        //     改对流格式时，改 system/fvSchemes/divSchemes/div(phi,U)。
        //     改扩散非正交处理时，改 system/fvSchemes/laplacianSchemes
        //     和 snGradSchemes。
        //     改动量方程线性求解器时，改 system/fvSolution/solvers/U。
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // 对应数学公式：
        //     PISO corrector 循环：
        //
        //         for k = 1, 2, ..., nCorrectors
        //
        //     每次校正从同一个离散动量矩阵关系出发：
        //
        //         A(U) U = H(U) - grad(p)
        //
        //     令：
        //
        //         rAU = 1/A(U)
        //
        //         HbyA = rAU H(U)
        //
        //     得到速度的压力响应形式：
        //
        //         U = HbyA - rAU grad(p)
        //
        // cases 入口：
        //     system/fvSolution/PISO:
        //         nCorrectors
        //
        // OpenFOAM 接口：
        //     piso.correct()
        //
        // 内部面/边界面：
        //     corrector 本身控制循环次数；
        //     真正区分内部面/边界面发生在后面的 flux、laplacian 和边界约束。
        //
        // 要改哪里：
        //     若要至少两次 PISO 校正，设置 nCorrectors 2 或更大。
        while (piso.correct())
        {
            // 对应数学公式：
            //     动量矩阵对角响应系数：
            //
            //         rAU_c = 1/A_c
            //
            //     显式动量贡献除以对角系数：
            //
            //         HbyA_c = rAU_c H_c
            //
            //     此时速度可写为：
            //
            //         U_c = HbyA_c - rAU_c grad(p)_c
            //
            // cases 入口：
            //     system/fvSchemes 中的动量方程离散格式会改变 UEqn.A()
            //     和 UEqn.H()。
            //
            // OpenFOAM 接口：
            //     UEqn.A()
            //     UEqn.H()
            //     constrainHbyA(...)
            //
            // 内部面/边界面：
            //     UEqn.A() 是单元中心对角系数；
            //     UEqn.H() 已经包含内部面 neighbour 贡献和边界贡献；
            //     constrainHbyA 会根据 0/U 和 0/p 的边界条件修正边界一致性。
            //
            // 要改哪里：
            //     改 rAU/HbyA 的数学含义，本质上要改 UEqn 的装配方式；
            //     常规算例只改 fvSchemes 和 fvSolution。
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

            // 对应数学公式：
            //     动量一致的预测面通量：
            //
            //         phiHbyA_f = flux(HbyA)_f
            //                   + interpolate(rAU)_f ddtCorr(U, phi)_f
            //
            //     第一项来自 HbyA 的面插值；
            //     第二项是瞬态 Rhie-Chow 型时间修正，用来保持单元速度、
            //     面通量和离散时间项的一致性。
            //
            // cases 入口：
            //     system/fvSchemes:
            //         interpolationSchemes/default
            //         ddtSchemes/default
            //
            //     0/U/boundaryField:
            //         movingTop
            //         no-slip walls
            //         empty front/back
            //
            // OpenFOAM 接口：
            //     fvc::flux(HbyA)
            //     fvc::interpolate(rAU)
            //     fvc::ddtCorr(U, phi)
            //
            // 内部面/边界面：
            //     内部面 phiHbyA_f 位于 owner/neighbour 之间；
            //     边界面 phiHbyA_f 位于 patch 上，壁面不可穿透条件通过
            //     0/U 的 fixedValue 和 constrainPressure/adjustPhi 保持一致。
            //
            // 要改哪里：
            //     一般不直接改这行；若中心线速度出现压力速度解耦，
            //     优先检查 ddtSchemes、interpolationSchemes 和边界条件。
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            // 对应数学公式：
            //     对封闭不可压区域，边界净通量必须满足：
            //
            //         sum_boundary phi_b = 0
            //
            //     顶盖驱动方腔中固壁不可穿透：
            //
            //         U · n = 0 on all walls
            //
            // cases 入口：
            //     0/U/boundaryField:
            //         movingTop fixedValue (1 0 0)
            //         leftWall/rightWall/bottomWall fixedValue (0 0 0)
            //
            //     0/p/boundaryField:
            //         walls zeroGradient
            //
            // OpenFOAM 接口：
            //     if (p.needReference())
            //         adjustPhi(phiHbyA, U, p)
            //
            // 内部面/边界面：
            //     adjustPhi 主要修正边界 patch 通量，使封闭域通量平衡；
            //     内部面通量不作为边界净通量入口。
            //
            // 要改哪里：
            //     若从封闭腔流改为入口出口流，要重新检查 U/p 边界条件，
            //     并确认是否仍需要 adjustPhi。
            if (p.needReference())
            {
                adjustPhi(phiHbyA, U, p);
            }

            // 对应数学公式：
            //     压力边界要与预测通量和速度修正相容：
            //
            //         U = HbyA - rAU grad(p)
            //
            //     在壁面不可穿透边界上：
            //
            //         U · n = 0
            //
            //     因此压力法向梯度需要使修正速度不穿透壁面。
            //
            // cases 入口：
            //     0/p/boundaryField
            //     0/U/boundaryField
            //     constant/polyMesh/boundary
            //
            // OpenFOAM 接口：
            //     constrainPressure(p, U, phiHbyA, rAU)
            //
            // 内部面/边界面：
            //     这一步主要作用在边界面；
            //     内部面的压力跳跃由后面的 laplacian(rAU,p) 矩阵处理。
            //
            // 要改哪里：
            //     若壁面、入口、出口或周期边界类型改变，先改 0/U 和 0/p；
            //     通常不改求解器这一行。
            constrainPressure(p, U, phiHbyA, rAU);

            // 对应数学公式：
            //     PISO 压力方程来自连续性：
            //
            //         div(U) = 0
            //
            //     代入：
            //
            //         U = HbyA - rAU grad(p)
            //
            //     得：
            //
            //         div(HbyA) - div(rAU grad(p)) = 0
            //
            //     即：
            //
            //         laplacian(rAU, p) = div(phiHbyA)
            //
            //     对控制体 Ω_c：
            //
            //         sum_f (rAU grad(p))_f · S_cf
            //       = sum_f phiHbyA_f
            //
            // cases 入口：
            //     system/fvSchemes:
            //         laplacianSchemes/laplacian(rAU,p)
            //         snGradSchemes/default
            //
            //     system/fvSolution/solvers/p:
            //         solver
            //         tolerance
            //         relTol
            //
            //     system/fvSolution/PISO:
            //         nNonOrthogonalCorrectors
            //         pRefCell
            //         pRefValue
            //
            // OpenFOAM 接口：
            //     fvm::laplacian(rAU, p)
            //     fvc::div(phiHbyA)
            //     fvScalarMatrix pEqn
            //     pEqn.flux()
            //
            // 内部面/边界面：
            //     内部面：laplacian 根据 owner/neighbour 组装压力系数；
            //     边界面：laplacian 根据 0/p/boundaryField 贡献边界通量；
            //     非正交修正由 correctNonOrthogonal 循环和 fvSchemes 共同决定。
            //
            // 要改哪里：
            //     改压力求解器，改 system/fvSolution/solvers/p。
            //     改非正交修正格式，改 system/fvSchemes/laplacianSchemes
            //     和 snGradSchemes。
            while (piso.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                // 对应数学公式：
                //     用压力方程的离散面通量修正守恒体积通量：
                //
                //         phi_f = phiHbyA_f - flux(laplacian(rAU,p))_f
                //
                //     这里的 flux(laplacian(rAU,p))_f 就是 pEqn.flux()。
                //
                // cases 入口：
                //     system/fvSchemes/laplacianSchemes/laplacian(rAU,p)
                //     system/fvSchemes/snGradSchemes/default
                //     constant/polyMesh/boundary
                //
                // OpenFOAM 接口：
                //     piso.finalNonOrthogonalIter()
                //     fvScalarMatrix::flux()
                //
                // 内部面/边界面：
                //     pEqn.flux() 同时包含内部面的压力通量和边界面的压力通量；
                //     只在最后一次非正交循环更新 phi，避免中间修正覆盖最终守恒通量。
                //
                // 要改哪里：
                //     不要把这里改成 phi = fvc::flux(U)；
                //     若要研究守恒通量误差，应围绕 pEqn.flux() 和 fvc::div(phi) 检查。
                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            // 对应数学公式：
            //     压力校正后的单元中心速度：
            //
            //         U_c = HbyA_c - rAU_c grad(p)_c
            //
            //     校正目标：
            //
            //         div(phi) = 0
            //
            // cases 入口：
            //     system/fvSchemes/gradSchemes/grad(p)
            //     0/U/boundaryField
            //     0/p/boundaryField
            //
            // OpenFOAM 接口：
            //     fvc::grad(p)
            //     volVectorField::correctBoundaryConditions()
            //
            // 内部面/边界面：
            //     U 是单元中心速度；
            //     边界速度在 correctBoundaryConditions 中按 patch 类型重新约束；
            //     真正用于连续性的仍是上一段得到的守恒面通量 phi。
            //
            // 要改哪里：
            //     改压力梯度格式，改 system/fvSchemes/gradSchemes。
            //     改壁面速度、顶盖速度或 empty 面，改 0/U/boundaryField。
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        // 对应数学公式：
        //     质量残差：
        //
        //         R_mass,c = div(phi)_c
        //
        //     离散守恒要求：
        //
        //         R_mass,c = 0
        //
        //     稳态速度变化：
        //
        //         ΔU_max = max_c |U_c(n+1) - U_c(n)|
        //
        // cases 入口：
        //     system/controlDict 的稳态阈值；
        //     constant/polyMesh 的内部面和边界面；
        //     0/U、0/p 的边界条件。
        //
        // OpenFOAM 接口：
        //     fvc::div(phi)
        //     mag(...)
        //     gMax(...)
        //
        // 内部面/边界面：
        //     fvc::div(phi) 对每个单元累加内部面通量和边界面通量；
        //     对封闭方腔，边界法向通量应接近零。
        //
        // 要改哪里：
        //     如果质量残差偏大，优先检查 pressure corrector、
        //     pEqn.flux() 更新、压力边界条件和 mesh patch 类型。
        const volScalarField continuity(fvc::div(phi));
        tmp<volScalarField> tContinuityMag(mag(continuity));
        const volScalarField& continuityMag = tContinuityMag();

        tmp<volScalarField> tVelocityChange(mag(U - UPrevious));
        const volScalarField& velocityChange = tVelocityChange();

        const scalar maxMassResidual =
            gMax(continuityMag.primitiveField());
        const scalar maxVelocityChange =
            gMax(velocityChange.primitiveField());

        Info<< "  max |div(phi)| = "
            << maxMassResidual
            << nl
            << "  max |U - Uprevious| = "
            << maxVelocityChange
            << nl << endl;

        runTime.write();

        // 对应数学公式：
        //     若：
        //
        //         max_c |U_c(n+1) - U_c(n)| <= eps_U
        //
        //     且：
        //
        //         max_c |div(phi)_c| <= eps_mass
        //
        //     并且已经至少推进 minimumSteadySteps 步，
        //     连续满足 requiredSteadySteps 次，则认为稳态。
        //
        // cases 入口：
        //     system/controlDict:
        //         steadyStateControl
        //         steadyVelocityTol
        //         steadyMassTol
        //         minimumSteadySteps
        //         requiredSteadySteps
        //
        // OpenFOAM 接口：
        //     Time::timeIndex()
        //     Time::writeAndEnd()
        //
        // 内部面/边界面：
        //     稳态判断使用单元速度和单元散度；
        //     散度来自内部面与边界面通量共同贡献。
        //
        // 要改哪里：
        //     报告稳态太早或太晚时，先改 controlDict 中的稳态阈值。
        if
        (
            steadyStateControl
         && runTime.timeIndex() >= minimumSteadySteps
         && maxVelocityChange <= steadyVelocityTol
         && maxMassResidual <= steadyMassTol
        )
        {
            ++steadyStepCount;
        }
        else
        {
            steadyStepCount = 0;
        }

        Info<< "  steady consecutive steps = "
            << steadyStepCount << " / " << requiredSteadySteps
            << nl << endl;

        if
        (
            steadyStateControl
         && steadyStepCount >= requiredSteadySteps
        )
        {
            Info<< "Steady state reached at Time = "
                << runTime.userTimeName() << nl << endl;
            runTime.writeAndEnd();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //

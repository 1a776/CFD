/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    explicitAdvectionDiffusionFoamStudent

Description
    第三题对流扩散方程全显式有限体积求解器。

    本文件按第二题求解器的注释标准书写。每个关键代码片段前都说明：

        对应数学公式
        OpenFOAM 接口
        cases 中相关文件和字段
        字段含义或可选写法

    主方程：

        ∂φ/∂t + ∇·(U φ) - ∇·(μ∇φ) = 0

    显式右端残差：

        R(c,n) = [-∇·(U φ(n)) + ∇·(μ∇φ(n))] in cell c

    显式更新：

        φ(c,n+1) = φ(c,n) + Δt R(c,n)

    OpenFOAM 核心接口：

        faceFlux = fvc::flux(U)
        Rphi = -fvc::div(faceFlux, phi) + fvc::laplacian(mu, phi)
\*---------------------------------------------------------------------------*/

// 对应数学公式：
//     无。这里是 OpenFOAM application 入口。
//
// OpenFOAM 接口：
//     argList
//
// cases 中相关文件和字段：
//     explicitAdvectionDiffusionFoamStudent -case path/to/case
#include "argList.H"

// 对应数学公式：
//     无。这里用于读取 OpenFOAM 字典。
//
// OpenFOAM 接口：
//     dictionary
//     IOdictionary
//
// cases 中相关文件和字段：
//     system/controlDict
//     constant/transportProperties
#include "dictionary.H"
#include "IOdictionary.H"

// 对应数学公式：
//     [μ] = L²/T
//     [Δt] = T
//
// OpenFOAM 接口：
//     dimensionedScalar
//
// cases 中相关文件和字段：
//     constant/transportProperties
//         mu [0 2 -1 0 0 0 0] value;
#include "dimensionedScalar.H"

// 对应数学公式：
//     U(c)
//     φ(c)
//     R(c)
//
// OpenFOAM 接口：
//     volVectorField
//     volScalarField
//
// cases 中相关文件和字段：
//     0/U
//     0/phi
//     constant/polyMesh
#include "volFields.H"

// 对应数学公式：
//     F(f) = U(f) · S(f)
//
// OpenFOAM 接口：
//     surfaceScalarField
//
// cases 中相关文件和字段：
//     constant/polyMesh/faces
//     constant/polyMesh/owner
//     constant/polyMesh/neighbour
//     constant/polyMesh/boundary
#include "surfaceFields.H"

// 对应数学公式：
//     F(f) = U(f) · S(f)
//
// OpenFOAM 接口：
//     fvc::flux(U)
//
// cases 中相关文件和字段：
//     0/U
//     system/fvSchemes/interpolationSchemes
#include "fvcFlux.H"

// 对应数学公式：
//     [∇·(Uφ)](c)
//
// OpenFOAM 接口：
//     fvc::div(faceFlux, phi, "div(faceFlux,phi)")
//
// cases 中相关文件和字段：
//     system/fvSchemes
//         divSchemes
//         {
//             div(faceFlux,phi) Gauss upwind;
//         }
#include "fvcDiv.H"

// 对应数学公式：
//     [∇·(μ∇φ)](c)
//
// OpenFOAM 接口：
//     fvc::laplacian(mu, phi)
//
// cases 中相关文件和字段：
//     system/fvSchemes
//         laplacianSchemes
//         {
//             laplacian(mu,phi) Gauss linear corrected;
//         }
#include "fvcLaplacian.H"

// 对应数学公式：
//     Σ over faces f of cell c |F(f)|
//     Σ over faces f of cell c D(f)
//
// OpenFOAM 接口：
//     fvc::surfaceSum(surfaceField)
//
// cases 中相关文件和字段：
//     constant/polyMesh/owner
//     constant/polyMesh/neighbour
//     constant/polyMesh/boundary
#include "fvcSurfaceIntegrate.H"

using namespace Foam;


int main(int argc, char *argv[])
{
    // 对应数学公式：
    //     无。这里确定当前运行的 case 路径。
    //
    // OpenFOAM 接口：
    //     setRootCase.H
    //
    // cases 中相关文件和字段：
    //     由命令行 -case path/to/case 指定
    #include "setRootCase.H"

    // 对应数学公式：
    //     t
    //     Δt
    //     t(end)
    //
    // OpenFOAM 接口：
    //     createTime.H
    //
    // cases 中相关文件和字段：
    //     system/controlDict
    //         startTime
    //         endTime
    //         writeControl
    //         writeInterval
    #include "createTime.H"

    // 对应数学公式：
    //     Ω = union of all cells
    //     V(c) = volume of cell c
    //     S(f) = face area vector
    //
    // OpenFOAM 接口：
    //     createMesh.H
    //
    // cases 中相关文件和字段：
    //     求解器实际读取：
    //
    //     constant/polyMesh/points
    //         -> 网格点坐标，决定计算区域 Ω
    //
    //     constant/polyMesh/faces
    //         -> 面拓扑，决定每个面由哪些点组成
    //
    //     constant/polyMesh/owner
    //         -> 每个面的 owner 单元
    //
    //     constant/polyMesh/neighbour
    //         -> 内部面的 neighbour 单元；有 neighbour 的面是内部面
    //
    //     constant/polyMesh/boundary
    //         -> 边界 patch 名称、patch 类型、边界面范围
    //
    //     网格生成入口：
    //
    //     system/blockMeshDict
    //         vertices -> 计算区域顶点
    //         blocks   -> 四边形/六面体结构网格和 N 的划分
    //         boundary -> patch 名称和 patch 类型
    //
    //     或者 Gmsh 网格生成脚本和 .geo/.msh 文件：
    //         physical groups -> OpenFOAM 边界 patch
    //         recombine false -> 三角形网格
    //         recombine true  -> 四边形网格
    //
    // 字段含义或可选写法：
    //     求解器不直接读取 blockMeshDict 或 .geo。
    //     blockMesh/gmshToFoam 先生成 constant/polyMesh，求解器再读取 polyMesh。
    #include "createMesh.H"

    // 对应数学公式：
    //     Δt 相关控制参数
    //
    // OpenFOAM 接口：
    //     runTime.controlDict()
    //
    // cases 中相关文件和字段：
    //     system/controlDict
    //         scalarField
    //         velocityField
    //         advectionDiffusionCo
    //         maxDeltaT
    const dictionary& controlDict = runTime.controlDict();

    /*
     * 第一部分：显式离散总框架
     *
     * 对应数学公式：
     *
     *     ∂φ/∂t + ∇·(Uφ) - ∇·(μ∇φ) = 0
     *
     *     dφ(c)/dt = R(c)
     *
     *     R(c) = [-∇·(Uφ) + ∇·(μ∇φ)] in cell c
     *
     *     φ(c,n+1) = φ(c,n) + Δt R(c,n)
     *
     * OpenFOAM 接口：
     *
     *     volVectorField U
     *     volScalarField phi
     *     surfaceScalarField faceFlux
     *     fvc::div(faceFlux, phi, "div(faceFlux,phi)")
     *     fvc::laplacian(mu, phi)
     *
     * cases 中相关文件和字段：
     *
     *     0/U
     *         internalField  -> U(c)
     *         boundaryField  -> 速度边界条件
     *
     *     0/phi
     *         internalField  -> φ(c,0)
     *         boundaryField  -> φ 的边界条件
     *
     *     constant/transportProperties
     *         mu             -> μ
     *
     *     system/fvSchemes
     *         divSchemes       -> 对流格式
     *         laplacianSchemes -> 扩散格式
     *         snGradSchemes    -> 非正交修正
     *
     *     system/controlDict
     *         advectionDiffusionCo -> 显式对流扩散稳定系数
     *         maxDeltaT            -> 最大时间步
     */

    // 对应数学公式：
    //     U(c,n)
    //
    // OpenFOAM 接口：
    //     controlDict.lookupOrDefault<word>("velocityField", "U")
    //
    // cases 中相关文件和字段：
    //     system/controlDict
    //         velocityField U;
    //
    // 字段含义或可选写法：
    //     velocityField U;     -> 读取 0/U
    //     velocityField Urot;  -> 读取 0/Urot
    const word velocityName
    (
        controlDict.lookupOrDefault<word>("velocityField", "U")
    );

    // 对应数学公式：
    //     φ(c,n)
    //
    // OpenFOAM 接口：
    //     controlDict.lookupOrDefault<word>("scalarField", "phi")
    //
    // cases 中相关文件和字段：
    //     system/controlDict
    //         scalarField phi;
    //
    // 字段含义或可选写法：
    //     scalarField phi; -> 读取 0/phi
    //     scalarField T;   -> 读取 0/T
    const word scalarName
    (
        controlDict.lookupOrDefault<word>("scalarField", "phi")
    );

    // 对应数学公式：
    //     U(c,n)
    //
    // OpenFOAM 接口：
    //     volVectorField
    //     IOobject::MUST_READ
    //     IOobject::NO_WRITE
    //
    // cases 中相关文件和字段：
    //     0/U
    //         internalField -> 单元中心速度
    //         boundaryField -> 速度边界条件
    //
    // 字段含义或可选写法：
    //     正弦波算例：
    //         U = (a,a,0)
    //
    //     旋转尖峰算例：
    //         U = (-y,x,0)
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

    // 对应数学公式：
    //     φ(c,n)
    //
    // OpenFOAM 接口：
    //     volScalarField
    //     IOobject::MUST_READ
    //     IOobject::AUTO_WRITE
    //
    // cases 中相关文件和字段：
    //     0/phi
    //         internalField -> 初始条件 φ(c,0) 或 φ(c,t0)
    //         boundaryField -> 边界条件
    //
    // 字段含义或可选写法：
    //     internalField nonuniform List<scalar> (...) -> 解析函数写到每个单元
    //
    //     type cyclic;        -> 周期边界，正弦波算例使用
    //     type fixedValue;    -> φ = φ_exact(x_b,t)，旋转尖峰严格边界使用
    //     type zeroGradient;  -> ∂φ/∂n = 0，若用于旋转尖峰需说明是近似
    //     type fixedGradient; -> ∂φ/∂n = given gradient
    volScalarField phi
    (
        IOobject
        (
            scalarName,
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // 对应数学公式：
    //     μ
    //
    // OpenFOAM 接口：
    //     IOdictionary
    //
    // cases 中相关文件和字段：
    //     constant/transportProperties
    //         mu [0 2 -1 0 0 0 0] value;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // 对应数学公式：
    //     [μ] = L²/T
    //
    // OpenFOAM 接口：
    //     dimensionedScalar("mu", dimensionSet(0,2,-1,0,0), dictionary)
    //
    // cases 中相关文件和字段：
    //     constant/transportProperties
    //         mu [0 2 -1 0 0 0 0] 1;
    //
    // 字段含义或可选写法：
    //     正弦波对流扩散：mu = 1
    //     旋转尖峰：mu = epsilon = 1e-3
    dimensionedScalar mu
    (
        "mu",
        dimensionSet(0, 2, -1, 0, 0),
        transportProperties
    );

    // 对应数学公式：
    //     Δt = α min_c V(c) / [0.5 Σ_f |F_cf| + Σ_f D_cf]
    //
    // OpenFOAM 接口：
    //     controlDict.lookupOrDefault<scalar>("advectionDiffusionCo", 0.5)
    //
    // cases 中相关文件和字段：
    //     system/controlDict
    //         advectionDiffusionCo 0.5;
    //
    // 字段含义或可选写法：
    //     α 越小越稳定，计算越慢。
    const scalar advectionDiffusionCo
    (
        controlDict.lookupOrDefault<scalar>("advectionDiffusionCo", 0.5)
    );

    if (advectionDiffusionCo <= 0)
    {
        FatalErrorInFunction
            << "advectionDiffusionCo must be positive, but received "
            << advectionDiffusionCo
            << exit(FatalError);
    }

    // 对应数学公式：
    //     Δt = min(Δt, maxDeltaT)
    //
    // OpenFOAM 接口：
    //     controlDict.lookupOrDefault<scalar>("maxDeltaT", great)
    //
    // cases 中相关文件和字段：
    //     system/controlDict
    //         maxDeltaT 1e-4;
    //
    // 字段含义或可选写法：
    //     maxDeltaT -> 用户给定最大时间步；不写则不额外限制
    const scalar maxDeltaT
    (
        controlDict.lookupOrDefault<scalar>("maxDeltaT", great)
    );

    Info<< "Student explicit advection-diffusion solver" << nl
        << "  case                 = " << runTime.caseName() << nl
        << "  time                 = " << runTime.name() << nl
        << "  cells                = " << mesh.nCells() << nl
        << "  velocityField        = " << velocityName << nl
        << "  scalarField          = " << scalarName << nl
        << "  mu                   = " << mu.value() << nl
        << "  advectionDiffusionCo = " << advectionDiffusionCo << nl
        << "  phi min              = " << min(phi).value() << nl
        << "  phi max              = " << max(phi).value() << nl
        << endl;

    /*
     * 第二部分：空间离散
     *
     * 对应数学公式：
     *
     *     R(c,n) = [-∇·(Uφ(n)) + ∇·(μ∇φ(n))] in cell c
     *
     *     对流内部面：
     *         F(f) = U(f) · S(f)
     *         A(f) = F(f) φ_up(f)
     *
     *     扩散内部面：
     *         H(f) = J_D(f) · S(f)
     *         J_D = -μ∇φ
     *
     * OpenFOAM 接口：
     *
     *     fvc::flux(U)
     *     fvc::div(faceFlux, phi, "div(faceFlux,phi)")
     *     fvc::laplacian(mu, phi)
     *
     * cases 中相关文件和字段：
     *
     *     system/fvSchemes
     *         divSchemes
     *             div(faceFlux,phi) Gauss upwind;
     *
     *         laplacianSchemes
     *             laplacian(mu,phi) Gauss linear corrected;
     *
     *         snGradSchemes
     *             default corrected;
     *
     *     0/phi
     *         boundaryField
     */

    // 对应数学公式：
    //     F(f) = U(f) · S(f)
    //
    // OpenFOAM 接口：
    //     fvc::flux(U)
    //
    // 内部面/边界面在哪里分辨：
    //     fvc::flux(U) 生成所有面的 surfaceScalarField。
    //     内部面由 constant/polyMesh/owner 和 neighbour 描述。
    //     边界面由 constant/polyMesh/boundary 和 0/U/boundaryField 描述。
    //
    // cases 中相关文件和字段：
    //     0/U
    //         internalField
    //         boundaryField
    //
    //     system/fvSchemes
    //         interpolationSchemes
    //
    // 字段含义或可选写法：
    //     faceFlux 只是运行时由 U 和 mesh 计算出的面通量，不需要在 0/ 里手写。
    surfaceScalarField faceFlux
    (
        IOobject
        (
            "faceFlux",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux(U)
    );

    // 对应数学公式：
    //     D(f) = μ |S(f)| / |d(f)|
    //     A_diff(c) = Σ over faces f of cell c D(f)
    //
    // OpenFOAM 接口：
    //     mesh.magSf()
    //     mesh.deltaCoeffs()
    //     fvc::surfaceSum(faceDiffCoeff)
    //
    // cases 中相关文件和字段：
    //     constant/polyMesh/faces
    //     constant/polyMesh/owner
    //     constant/polyMesh/neighbour
    //
    // 字段含义或可选写法：
    //     这里仅用于估计显式扩散时间步。
    //     真正的扩散残差仍由 fvc::laplacian(mu, phi) 按 fvSchemes 计算。
    surfaceScalarField faceDiffCoeff
    (
        IOobject
        (
            "faceDiffCoeff",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mu*mesh.magSf()*mesh.deltaCoeffs()
    );

    // 对应数学公式：
    //     advectiveDen(c) = 0.5 Σ over faces f of cell c |F_cf|
    //
    // OpenFOAM 接口：
    //     fvc::surfaceSum(mag(faceFlux))
    //
    // cases 中相关文件和字段：
    //     constant/polyMesh/owner
    //     constant/polyMesh/neighbour
    //
    // 字段含义或可选写法：
    //     0.5 系数对应闭合控制体中流出通量约为绝对通量和的一半，
    //     与第一题 CFL 估计 Co = 0.5 Δt Σ|F|/V 保持一致。
    scalarField advectiveDen
    (
        0.5*fvc::surfaceSum(mag(faceFlux))().primitiveField()
    );

    // 对应数学公式：
    //     diffusiveDen(c) = Σ over faces f of cell c D(f)
    //
    // OpenFOAM 接口：
    //     fvc::surfaceSum(faceDiffCoeff)
    //
    // cases 中相关文件和字段：
    //     constant/polyMesh
    //
    // 字段含义或可选写法：
    //     这一项给出显式扩散的 O(h²) 时间步限制。
    scalarField diffusiveDen
    (
        fvc::surfaceSum(faceDiffCoeff)().primitiveField()
    );

    // 对应数学公式：
    //     denom(c) = advectiveDen(c) + diffusiveDen(c)
    //
    // OpenFOAM 接口：
    //     scalarField
    //
    // cases 中相关文件和字段：
    //     无直接字段；由 U、mu 和 mesh 共同决定
    scalarField combinedDen
    (
        advectiveDen + diffusiveDen
    );

    const scalar combinedDenMax
    (
        gMax(combinedDen/(mesh.V().primitiveField() + VSMALL))
    );

    if (combinedDenMax <= SMALL)
    {
        FatalErrorInFunction
            << "The combined advection-diffusion rate is zero or too small: "
            << combinedDenMax << nl
            << "Cannot determine a positive explicit time step."
            << exit(FatalError);
    }

    /*
     * 第三部分：时间离散
     *
     * 对应数学公式：
     *
     *     Δt = α min_c V(c) / [0.5 Σ_f |F_cf| + Σ_f D_cf]
     *
     *     φ(c,n+1) = φ(c,n) + Δt R(c,n)
     *
     * OpenFOAM 接口：
     *
     *     runTime.setDeltaT(Δt)
     *     runTime.setTime(t, step)
     *     runTime.write()
     *
     * cases 中相关文件和字段：
     *
     *     system/controlDict
     *         advectionDiffusionCo
     *         maxDeltaT
     *         endTime
     *         writeControl
     *         writeInterval
     */

    const scalar endTime = runTime.endTime().value();
    const scalar startTime = runTime.value();
    const scalar timeTolerance = 1e-12*max(1.0, mag(endTime));
    const scalar stableDeltaT
    (
        min(advectionDiffusionCo/combinedDenMax, maxDeltaT)
    );

    if (stableDeltaT <= SMALL)
    {
        FatalErrorInFunction
            << "Calculated advection-diffusion deltaT is too small: "
            << stableDeltaT
            << exit(FatalError);
    }

    Info<< "  combined rate max   = " << combinedDenMax << nl
        << "  stable deltaT       = " << stableDeltaT << nl
        << "  faceFlux min        = " << min(faceFlux).value() << nl
        << "  faceFlux max        = " << max(faceFlux).value() << nl
        << endl;

    label step = 0;

    while (endTime - runTime.value() > timeTolerance)
    {
        // 对应数学公式：
        //     R(c,n) = [-∇·(Uφ(n)) + ∇·(μ∇φ(n))] in cell c
        //
        // OpenFOAM 接口：
        //     -fvc::div(faceFlux, phi, "div(faceFlux,phi)")
        //     +fvc::laplacian(mu, phi)
        //
        // 内部面/边界面在哪里分辨：
        //     本求解器这里不手写 forAll internal faces / boundary faces。
        //
        //     fvc::div 和 fvc::laplacian 会根据 mesh 和 boundaryField 自动处理：
        //
        //         内部面：
        //             constant/polyMesh/owner
        //             constant/polyMesh/neighbour
        //
        //         边界面：
        //             constant/polyMesh/boundary
        //             0/phi/boundaryField
        //             0/U/boundaryField
        //
        // cases 中相关文件和字段：
        //     system/fvSchemes
        //         div(faceFlux,phi) Gauss upwind;
        //         laplacian(mu,phi) Gauss linear corrected;
        //
        //     0/phi
        //         boundaryField
        //
        // 字段含义或可选写法：
        //     对流格式可改为：
        //         div(faceFlux,phi) Gauss linearUpwind grad(phi);
        //
        //     扩散格式可改为：
        //         laplacian(mu,phi) Gauss linear uncorrected;
        //         laplacian(mu,phi) Gauss linear limited 0.5;
        tmp<volScalarField> tRphi
        (
            -fvc::div(faceFlux, phi, "div(faceFlux,phi)")
          +  fvc::laplacian(mu, phi)
        );

        const volScalarField& Rphi = tRphi();

        // 对应数学公式：
        //     t(target) = min(t(start) + (step + 1) Δt, t(end))
        //     Δt(step) = t(target) - t(old)
        //
        // OpenFOAM 接口：
        //     runTime.value()
        //
        // cases 中相关文件和字段：
        //     system/controlDict
        //         endTime
        //
        // 字段含义或可选写法：
        //     最后一步自动缩短 Δt，使最终输出严格落在 endTime
        const scalar oldTime = runTime.value();
        const scalar targetTime
        (
            min(startTime + (step + 1)*stableDeltaT, endTime)
        );
        const scalar stepDeltaT = targetTime - oldTime;

        if (stepDeltaT <= SMALL)
        {
            FatalErrorInFunction
                << "The remaining time is too small for another step: "
                << (endTime - oldTime)
                << exit(FatalError);
        }

        // 对应数学公式：
        //     Δt = Δt(step)
        //
        // OpenFOAM 接口：
        //     runTime.setDeltaT(Δt)
        //
        // cases 中相关文件和字段：
        //     system/controlDict
        //         writeControl
        //         writeInterval
        runTime.setDeltaT(stepDeltaT);

        // 对应数学公式：
        //     t(n+1) = t(n) + Δt
        //
        // OpenFOAM 接口：
        //     runTime.setTime(t, step)
        //
        // cases 中相关文件和字段：
        //     system/controlDict
        //         startTime
        //         endTime
        ++step;
        runTime.setTime(targetTime, step);

        Info<< "Time = " << runTime.name()
            << "  step = " << step
            << "  deltaT = " << stepDeltaT
            << "  advectionDiffusionCo = "
            << stepDeltaT*combinedDenMax
            << nl;

        // 对应数学公式：
        //     [Δt R] = [φ]
        //
        // OpenFOAM 接口：
        //     dimensionedScalar("deltaT", dimTime, stepDeltaT)
        //
        // cases 中相关文件和字段：
        //     无直接字段；stepDeltaT 来自上面的稳定时间步计算
        const dimensionedScalar deltaTDim
        (
            "deltaT",
            dimTime,
            stepDeltaT
        );

        // 对应数学公式：
        //     φ(c,n+1) = φ(c,n) + Δt R(c,n)
        //
        // OpenFOAM 接口：
        //     volScalarField = volScalarField + dimensionedScalar*volScalarField
        //
        // cases 中相关文件和字段：
        //     0/phi
        //         internalField
        //
        // 字段含义或可选写法：
        //     这一行就是显式 Euler。
        //     system/fvSchemes/ddtSchemes 不控制这一行源码更新。
        phi = phi + deltaTDim*Rphi;

        // 对应数学公式：
        //     边界面项 A(cb) + H(cb)
        //
        //     周期边界：
        //         φ left = φ right
        //
        //     Dirichlet 边界：
        //         φ = φ_b
        //
        //     Neumann 边界：
        //         ∇φ · n = g_b
        //
        // OpenFOAM 接口：
        //     phi.correctBoundaryConditions()
        //
        // cases 中相关文件和字段：
        //     0/phi
        //         boundaryField
        //
        // 字段含义或可选写法：
        //     正弦波算例：
        //         cyclic patches
        //
        //     旋转尖峰严格复现：
        //         time-dependent fixedValue from exact solution
        //
        //     若旋转尖峰临时使用 zeroGradient：
        //         这是工程近似，需要在报告中说明
        phi.correctBoundaryConditions();

        // 对应数学公式：
        //     M = Σ over cells c φ(c) V(c)
        //
        // OpenFOAM 接口：
        //     gSum(phi.primitiveField()*mesh.V().primitiveField())
        //
        // cases 中相关文件和字段：
        //     constant/polyMesh
        //     0/phi
        //
        // 字段含义或可选写法：
        //     M 用于监测总量变化；是否守恒取决于边界条件、速度场和扩散通量。
        const scalar mass
        (
            gSum(phi.primitiveField()*mesh.V().primitiveField())
        );

        Info<< "  Rphi min = " << min(Rphi).value() << nl
            << "  Rphi max = " << max(Rphi).value() << nl
            << "  phi min  = " << min(phi).value() << nl
            << "  phi max  = " << max(phi).value() << nl
            << "  mass     = " << mass << nl
            << endl;

        // 对应数学公式：
        //     φ(c,n+1)
        //
        // OpenFOAM 接口：
        //     runTime.write()
        //
        // cases 中相关文件和字段：
        //     system/controlDict
        //         writeControl
        //         writeInterval
        runTime.write();
    }

    // 最终输出：保证 endTime 有结果文件。
    //
    // 第三题验收时刻：
    //     正弦波算例：t = 1
    //     旋转尖峰算例：一圈旋转后的目标时刻
    //
    // 浮点数累加可能让循环在 |t(end)-t| 小于 timeTolerance 时结束，
    // 例如 t=0.9999999999999848。该差值不再做一次数值更新，但这里把
    // 时间标签精确设为 endTime，保证最终目录和后处理严格对应题目时刻。
    if (mag(endTime - runTime.value()) <= timeTolerance)
    {
        runTime.setTime(endTime, step);
    }

    const bool finalWriteOK = runTime.writeNow();

    if (!finalWriteOK)
    {
        FatalErrorInFunction
            << "Failed to write the final field at time "
            << runTime.name()
            << exit(FatalError);
    }

    Info<< "End" << nl << endl;

    return 0;
}

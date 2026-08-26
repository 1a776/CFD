/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    explicitDiffusionFoamStudent

Description
    第二题扩散方程全显式有限体积求解器框架。

    注释按固定格式写：

        对应数学公式
        OpenFOAM 接口
        cases 中相关文件和字段
        字段含义或可选写法

    源码按数值方法顺序组织：

        第一部分：显式离散总框架
        第二部分：空间离散
        第三部分：时间离散

    主方程：

        ∂φ/∂t - ∇·(μ∇φ) = 0

    显式残差：

        R(c,n) = [∇·(μ∇φ(n))] in cell c

    显式更新：

        φ(c,n+1) = φ(c,n) + Δt R(c,n)
\*---------------------------------------------------------------------------*/

// 对应数学公式：
//     无。这里是 OpenFOAM application 入口。
//
// OpenFOAM 接口：
//     argList
//
// cases 中相关文件和字段：
//     explicitDiffusionFoamStudent -case path/to/case
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
//     φ(c)
//     R(c)
//     D(f)
//
// OpenFOAM 接口：
//     volScalarField
//     surfaceScalarField
//
// cases 中相关文件和字段：
//     0/phi
//     constant/polyMesh
#include "volFields.H"
#include "surfaceFields.H"

// 对应数学公式：
//     R(c) = [∇·(μ∇φ)] in cell c
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
//     A(c) = Σ over faces f of cell c D(f)
//
// OpenFOAM 接口：
//     fvc::surfaceSum(faceDiffCoeff)
//
// cases 中相关文件和字段：
//     constant/polyMesh
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
    //         -> 网格点坐标，决定计算区域 Ω 的几何位置
    //
    //     constant/polyMesh/faces
    //         -> 面由哪些点组成，决定四边形面/三角形面等拓扑
    //
    //     constant/polyMesh/owner
    //         -> 每个面的 owner 单元
    //
    //     constant/polyMesh/neighbour
    //         -> 内部面的 neighbour 单元；有 neighbour 的面是内部面
    //
    //     constant/polyMesh/boundary
    //         -> 边界 patch 名称、类型、边界面范围
    //
    //     网格生成入口：
    //
    //     system/blockMeshDict
    //         vertices -> 计算区域顶点，例如单位正方形 [0,1]×[0,1]
    //         blocks   -> 四边形/六面体结构网格和 N 的划分
    //         boundary -> patch 名称和 patch 类型
    //
    //     或者 Gmsh 网格生成脚本和 .geo/.msh 文件：
    //         几何点、线、面 -> 计算区域 Ω
    //         recombine false -> 三角形网格
    //         recombine true  -> 四边形网格
    //         physical groups -> OpenFOAM 边界 patch
    //
    // 字段含义或可选写法：
    //     求解器不直接读取 blockMeshDict 或 .geo。
    //     blockMesh/gmshToFoam 先生成 constant/polyMesh，求解器再读取 polyMesh。
    #include "createMesh.H"

    // 对应数学公式：
    //     Δt = σ min over cells c [ V(c) / A(c) ]
    //
    // OpenFOAM 接口：
    //     runTime.controlDict()
    //
    // cases 中相关文件和字段：
    //     system/controlDict
    //         scalarField  -> 选择未知量文件名，例如 phi 或 T
    //         diffusionCo  -> 显式扩散稳定系数 σ
    //         maxDeltaT    -> 最大时间步限制
    const dictionary& controlDict = runTime.controlDict();


    /*
     * 第一部分：显式离散总框架
     *
     * 对应数学公式：
     *
     *     ∂φ/∂t - ∇·(μ∇φ) = 0
     *     dφ(c)/dt = R(c)
     *     R(c) = [∇·(μ∇φ)] in cell c
     *     φ(c,n+1) = φ(c,n) + Δt R(c,n)
     *
     * OpenFOAM 接口：
     *
     *     volScalarField phi
     *     dimensionedScalar mu
     *     fvc::laplacian(mu, phi)
     *     runTime.setDeltaT(Δt)
     *
     * cases 中相关文件和字段：
     *
     *     0/phi
     *         internalField  -> φ(c,0)
     *         boundaryField  -> 边界面条件
     *
     *     constant/polyMesh
     *         points         -> 计算区域 Ω 的几何坐标
     *         faces          -> 网格面
     *         owner          -> 面 owner 单元
     *         neighbour      -> 内部面 neighbour 单元
     *         boundary       -> 边界 patch
     *
     *     system/blockMeshDict 或 Gmsh 网格文件
     *         -> 生成 constant/polyMesh 的前处理入口
     *         -> 用来选择计算区域、四边形网格、三角形网格和分辨率 N
     *
     *     constant/transportProperties
     *         mu             -> μ
     *
     *     system/controlDict
     *         scalarField    -> 选择 φ 字段文件名
     *         diffusionCo    -> σ
     *         maxDeltaT      -> 最大 Δt
     *         endTime        -> 终止时间
     *
     *     system/fvSchemes
     *         laplacianSchemes -> 空间扩散离散格式
     */

    // 对应数学公式：
    //     φ
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
    //     φ(c,n)
    //
    // OpenFOAM 接口：
    //     volScalarField
    //     IOobject::MUST_READ  -> 必须从当前时间目录读取字段
    //     IOobject::AUTO_WRITE -> 写时间步时自动输出字段
    //
    // cases 中相关文件和字段：
    //     0/phi
    //         internalField -> 单元中心初值 φ(c,0)
    //         boundaryField -> 边界面条件
    //
    // 字段含义或可选写法：
    //     internalField uniform value;
    //         -> 所有单元使用同一个初值
    //
    //     internalField nonuniform List<scalar> (...);
    //         -> 每个单元中心给一个初值，适合解析函数离散到网格中心
    //
    //     type fixedValue;    -> φ = given value
    //     type zeroGradient;  -> ∂φ/∂n = 0
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
    //     ∂φ/∂t = ∇·(μ∇φ)
    //
    // OpenFOAM 接口：
    //     IOdictionary
    //
    // cases 中相关文件和字段：
    //     constant/transportProperties
    //         mu [0 2 -1 0 0 0 0] value;
    //
    // 字段含义或可选写法：
    //     value -> 扩散系数 μ
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
    //     ∂φ/∂t = ∇·(μ∇φ)
    //
    // OpenFOAM 接口：
    //     dimensionedScalar("mu", dimViscosity, transportProperties)
    //
    // cases 中相关文件和字段：
    //     constant/transportProperties
    //         mu [0 2 -1 0 0 0 0] value;
    //
    // 字段含义或可选写法：
    //     [0 2 -1 0 0 0 0] -> μ 的量纲 L²/T
    //     value             -> μ 的数值
    dimensionedScalar mu
    (
        "mu",
        dimensionSet(0, 2, -1, 0, 0),
        transportProperties
    );

    // 对应数学公式：
    //     Δt = σ min over cells c [ V(c) / A(c) ]
    //
    // OpenFOAM 接口：
    //     controlDict.lookupOrDefault<scalar>("diffusionCo", 0.5)
    //
    // cases 中相关文件和字段：
    //     system/controlDict
    //         diffusionCo 0.5;
    //
    // 字段含义或可选写法：
    //     diffusionCo -> 显式扩散稳定系数 σ
    const scalar diffusionCo
    (
        controlDict.lookupOrDefault<scalar>("diffusionCo", 0.5)
    );

    if (diffusionCo <= 0)
    {
        FatalErrorInFunction
            << "diffusionCo must be positive, but received "
            << diffusionCo
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
    //     maxDeltaT -> 用户给定的最大时间步；不写则不额外限制
    const scalar maxDeltaT
    (
        controlDict.lookupOrDefault<scalar>("maxDeltaT", great)
    );

    Info<< "Student explicit diffusion solver" << nl
        << "  case        = " << runTime.caseName() << nl
        << "  time        = " << runTime.name() << nl
        << "  cells       = " << mesh.nCells() << nl
        << "  scalarField = " << scalarName << nl
        << "  mu          = " << mu.value() << nl
        << "  diffusionCo = " << diffusionCo << nl
        << "  phi min     = " << min(phi).value() << nl
        << "  phi max     = " << max(phi).value() << nl
        << endl;


    /*
     * 第二部分：空间离散
     *
     * 对应数学公式：
     *     ∂φ/∂t = ∇·(μ∇φ)
     *     R(c,n) = 1/V(c) Σ over faces f of cell c [ μ(f) ∇φ(f,n) · S(f) ]
     *
     * OpenFOAM 接口：
     *     fvc::laplacian(mu, phi)
     *
     * cases 中相关文件和字段：
     *     system/fvSchemes
     *         laplacianSchemes
     *             laplacian(mu,phi) Gauss linear corrected;
     *
     *     system/fvSchemes
     *         snGradSchemes
     *             default corrected;
     *
     *     0/phi
     *         boundaryField
     *
     * 字段含义或可选写法：
     *     Gauss       -> 有限体积面通量求和
     *     linear      -> 面值线性插值
     *     corrected   -> 非正交修正
     *     uncorrected -> 不做非正交修正
     *     limited ψ   -> 限制非正交修正强度
     */

    // 对应数学公式：
    //     R(c,n) = 1/V(c) Σ over faces f of cell c [ μ(f) ∇φ(f,n) · S(f) ]
    //
    //     内部面 f：
    //         f connects owner cell O(f) and neighbour cell N(f)
    //         one face flux enters O(f) and N(f) with opposite signs
    //
    //     边界面 f：
    //         f lies on one boundary patch
    //         φ(f) or ∂φ/∂n on f is provided by 0/phi boundaryField
    //
    // OpenFOAM 接口：
    //     fvc::laplacian(mu, phi)
    //
    // 内部面/边界面在哪里分辨：
    //     本求解器这里不写 if internalFace / if boundaryFace。
    //
    //     fvc::laplacian(mu, phi) 会根据 mesh 拓扑自动分辨：
    //
    //         内部面：
    //             constant/polyMesh/owner
    //             constant/polyMesh/neighbour
    //
    //         边界面：
    //             constant/polyMesh/boundary
    //             0/phi/boundaryField
    //
    //     有 neighbour cell 的面是内部面。
    //     没有 neighbour cell、属于某个 patch 的面是边界面。
    //
    // cases 中相关文件和字段：
    //     system/fvSchemes
    //         laplacianSchemes
    //             laplacian(mu,phi) Gauss linear corrected;
    //
    //     system/fvSchemes
    //         snGradSchemes
    //             default corrected;
    //
    //     0/phi
    //         boundaryField
    //
    // 字段含义或可选写法：
    //     laplacian(mu,phi) Gauss linear corrected;
    //     laplacian(mu,phi) Gauss linear uncorrected;
    //     laplacian(mu,phi) Gauss linear limited 0.5;
    //
    //     zeroGradient   in 0/phi -> ∂φ/∂n = 0
    //     fixedValue     in 0/phi -> φ = given value
    //     fixedGradient  in 0/phi -> ∂φ/∂n = given gradient
    auto computeDiffusionResidual = [&]() -> tmp<volScalarField>
    {
        return fvc::laplacian(mu, phi);
    };

    // 对应数学公式：
    //     D(f) = μ |S(f)| / |d(f)|
    //     A(c) = Σ over faces f of cell c D(f)
    //
    // OpenFOAM 接口：
    //     mesh.magSf()       -> |S(f)|
    //     mesh.deltaCoeffs() -> 1 / |d(f)|
    //
    // cases 中相关文件和字段：
    //     constant/polyMesh/faces
    //     constant/polyMesh/points
    //     constant/polyMesh/owner
    //     constant/polyMesh/neighbour
    //
    // 字段含义或可选写法：
    //     这些文件由 blockMesh、gmshToFoam 等网格工具生成
    //     quad/tri 网格的差别最终体现在这些 polyMesh 文件里
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
    //     A(c) = Σ over faces f of cell c D(f)
    //
    // OpenFOAM 接口：
    //     fvc::surfaceSum(faceDiffCoeff)
    //
    // cases 中相关文件和字段：
    //     constant/polyMesh/owner
    //     constant/polyMesh/neighbour
    //
    // 字段含义或可选写法：
    //     owner/neighbour 决定一个面属于哪些单元，进而决定 D(f) 加到哪个 A(c)
    scalarField aDiff
    (
        fvc::surfaceSum(faceDiffCoeff)().primitiveField()
    );

    /*
     * 第三部分：时间离散
     *
     * 对应数学公式：
     *     Δt = σ min over cells c [ V(c) / A(c) ]
     *     φ(c,n+1) = φ(c,n) + Δt R(c,n)
     *
     * OpenFOAM 接口：
     *     runTime.setDeltaT(Δt)
     *     runTime.setTime(t, step)
     *     runTime.write()
     *
     * cases 中相关文件和字段：
     *     system/controlDict
     *         diffusionCo
     *         maxDeltaT
     *         endTime
     *         writeControl
     *         writeInterval
     *
     * 字段含义或可选写法：
     *     diffusionCo  -> σ
     *     maxDeltaT    -> 最大 Δt
     *     endTime      -> 最终计算时间
     *
     *     本求解器的显式 Euler 时间推进写在源码中。
     *     system/fvSchemes/ddtSchemes 不控制这一行显式更新。
     */

    // 对应数学公式：
    //     t(start)
    //     t(end)
    //
    // OpenFOAM 接口：
    //     runTime.endTime()
    //     runTime.value()
    //
    // cases 中相关文件和字段：
    //     system/controlDict
    //         startTime
    //         endTime
    //
    // 字段含义或可选写法：
    //     startTime -> 初始计算时间
    //     endTime   -> 终止计算时间
    const scalar endTime = runTime.endTime().value();
    const scalar startTime = runTime.value();
    const scalar timeTolerance = 1e-12*max(1.0, mag(endTime));
    label step = 0;

    while (endTime - runTime.value() > timeTolerance)
    {
        // 对应数学公式：
        //     R(c,n) = [∇·(μ∇φ(n))] in cell c
        //
        // OpenFOAM 接口：
        //     computeDiffusionResidual()
        //     fvc::laplacian(mu, phi)
        //
        // cases 中相关文件和字段：
        //     system/fvSchemes/laplacianSchemes
        //     system/fvSchemes/snGradSchemes
        //     0/phi/boundaryField
        //
        // 字段含义或可选写法：
        //     这里使用第二部分定义的空间离散入口，读取当前 φ(c,n)
        tmp<volScalarField> tRphi
        (
            computeDiffusionResidual()
        );

        const volScalarField& Rphi = tRphi();

        // 对应数学公式：
        //     Δt(c) = V(c) / A(c)
        //     Δt = σ min over cells c Δt(c)
        //
        // OpenFOAM 接口：
        //     mesh.V() -> V(c)
        //     gMin(candidateDeltaT)
        //
        // cases 中相关文件和字段：
        //     system/controlDict
        //         diffusionCo
        //         maxDeltaT
        //
        // 字段含义或可选写法：
        //     diffusionCo -> σ，越小越稳定，时间步越小
        //     maxDeltaT   -> 用户额外限制的最大时间步
        scalarField candidateDeltaT
        (
            mesh.V().primitiveField()/(aDiff + VSMALL)
        );

        const scalar stableDeltaT
        (
            min(diffusionCo*gMin(candidateDeltaT), maxDeltaT)
        );

        if (stableDeltaT <= SMALL)
        {
            FatalErrorInFunction
                << "Calculated diffusion deltaT is too small: "
                << stableDeltaT
                << exit(FatalError);
        }

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
        //
        // 字段含义或可选写法：
        //     runTime.write() 后续会按 writeControl/writeInterval 判断是否写出
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
        //
        // 字段含义或可选写法：
        //     step 是当前显式推进步数
        ++step;
        runTime.setTime(targetTime, step);

        Info<< "Time = " << runTime.name()
            << "  step = " << step
            << "  deltaT = " << stepDeltaT
            << nl;

        // 对应数学公式：
        //     [Δt R] = [φ]
        //
        // OpenFOAM 接口：
        //     dimensionedScalar("deltaT", dimTime, stepDeltaT)
        //
        // cases 中相关文件和字段：
        //     无直接字段；stepDeltaT 来自上面的稳定时间步计算
        //
        // 字段含义或可选写法：
        //     dimTime 保证 OpenFOAM 场运算量纲正确
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
        //     system/fvSchemes/ddtSchemes 不控制这一行。
        phi = phi + deltaTDim*Rphi;

        // 对应数学公式：
        //     fixedValue    -> φ = φ(b)
        //     zeroGradient  -> ∂φ/∂n = 0
        //     fixedGradient -> ∂φ/∂n = g(b)
        //
        // OpenFOAM 接口：
        //     phi.correctBoundaryConditions()
        //
        // cases 中相关文件和字段：
        //     0/phi
        //         boundaryField
        //
        // 字段含义或可选写法：
        //     每次更新内部场后，重新根据 boundaryField 更新边界场
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
        //     M 用于检查总量变化；是否守恒取决于边界条件和方程设置
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
        //
        // 字段含义或可选写法：
        //     writeControl timeStep;    -> 按步数输出
        //     writeControl runTime;     -> 按物理时间输出
        //     writeInterval value;      -> 输出间隔
        runTime.write();
    }

    // 最终输出：保证 endTime 有结果文件。
    //
    // 第二题验收时刻：
    //
    //     t = 0.2
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

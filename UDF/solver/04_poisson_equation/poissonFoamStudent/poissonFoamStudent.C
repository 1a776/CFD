/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    poissonFoamStudent

Description
    第四题 Poisson 方程有限体积求解器。

    本文件按前三题求解器的注释标准书写。每个关键代码片段前都说明：

        对应数学公式
        OpenFOAM 接口
        cases 中相关文件和字段
        字段含义或可选写法

    主方程：

        ∇²φ = ω

    控制体有限体积离散骨架：

        Σ over faces f of cell c [ ∇φ(f) · S(c,f) ] = V(c) ω(c)

    OpenFOAM 核心接口：

        fvm::laplacian(phi) == omega

    与前三题的主要区别：

        前三题是瞬态显式求解器：

            φ(c,n+1) = φ(c,n) + Δt R(c,n)

        第四题是稳态 Poisson 边值问题：

            A φ = b

        所以这里不计算 CFL，不计算显式时间步，也不手动更新 φ。
        fvm::laplacian(phi) 会根据网格、边界条件和 fvSchemes 自动组装
        fvScalarMatrix；phiEqn.solve() 再根据 fvSolution 求解线性方程组。
\*---------------------------------------------------------------------------*/

// 对应数学公式：
//     无。这里是 OpenFOAM application 入口。
//
// OpenFOAM 接口：
//     argList
//
// cases 中相关文件和字段：
//     poissonFoamStudent -case path/to/case
#include "argList.H"

// 对应数学公式：
//     无。这里用于读取 OpenFOAM 字典。
//
// OpenFOAM 接口：
//     dictionary
//
// cases 中相关文件和字段：
//     system/controlDict
#include "dictionary.H"

// 对应数学公式：
//     φ(c)
//     ω(c)
//
// OpenFOAM 接口：
//     volScalarField
//
// cases 中相关文件和字段：
//     0/phi
//     0/omega
//     constant/polyMesh
#include "volFields.H"

// 对应数学公式：
//     A φ = b
//
// OpenFOAM 接口：
//     fvScalarMatrix
//
// cases 中相关文件和字段：
//     system/fvSolution
//         solvers
//         {
//             phi
//             {
//                 solver          PCG;
//                 preconditioner  DIC;
//                 tolerance       1e-12;
//                 relTol          0;
//             }
//         }
#include "fvMatrices.H"

// 对应数学公式：
//     ∇²φ = ω
//
// OpenFOAM 接口：
//     fvm::laplacian(phi)
//
// cases 中相关文件和字段：
//     system/fvSchemes
//         laplacianSchemes
//         {
//             laplacian(phi) Gauss linear corrected;
//         }
#include "fvmLaplacian.H"

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
    //     稳态问题没有真实时间推进。
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
    //
    // 字段含义或可选写法：
    //     本求解器只求一个稳态 Poisson 方程。
    //     runTime 仍然需要存在，因为 OpenFOAM 用它管理 0/、constant/、
    //     system/、输出时间目录和字段写出。
    #include "createTime.H"

    // 对应数学公式：
    //     Ω = union of all cells
    //     V(c) = volume of cell c
    //     S(c,f) = outward face area vector of cell c
    //
    // OpenFOAM 接口：
    //     createMesh.H
    //
    // cases 中相关文件和字段：
    //     求解器实际读取：
    //
    //     constant/polyMesh/points
    //         -> 网格点坐标，决定计算区域 Ω，例如 [0,1]×[0,1]
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
    //     无。这里读取本求解器的字段名和非正交修正次数。
    //
    // OpenFOAM 接口：
    //     runTime.controlDict()
    //
    // cases 中相关文件和字段：
    //     system/controlDict
    //         solutionField phi;
    //         sourceField   omega;
    //         nNonOrthogonalCorrectors 0;
    //
    // 字段含义或可选写法：
    //     solutionField -> 未知量字段名，默认 phi
    //     sourceField   -> 源项字段名，默认 omega
    //     nNonOrthogonalCorrectors -> 非正交修正次数，默认 0
    const dictionary& controlDict = runTime.controlDict();

    const word solutionName
    (
        controlDict.lookupOrDefault<word>("solutionField", "phi")
    );

    const word sourceName
    (
        controlDict.lookupOrDefault<word>("sourceField", "omega")
    );

    const label nNonOrthogonalCorrectors
    (
        controlDict.lookupOrDefault<label>("nNonOrthogonalCorrectors", 0)
    );

    if (nNonOrthogonalCorrectors < 0)
    {
        FatalErrorInFunction
            << "nNonOrthogonalCorrectors must be non-negative, but received "
            << nNonOrthogonalCorrectors
            << exit(FatalError);
    }

    /*
     * 第一部分：读入未知量、源项和边界条件
     *
     * 对应数学公式：
     *
     *     ∇²φ = ω
     *
     *     φ is unknown
     *     ω is known source
     *
     *     第四题解析验证案例：
     *
     *         φ_exact(x,y) = cos(pi x) cos(pi y)
     *
     *         ω(x,y) = -2 pi^2 cos(pi x) cos(pi y)
     *
     *         φ on boundary = φ_exact on boundary
     *
     * OpenFOAM 接口：
     *
     *     volScalarField phi
     *     volScalarField omega
     *
     * cases 中相关文件和字段：
     *
     *     0/phi
     *         internalField  -> 初始猜测 φ(c)，例如 uniform 0
     *         boundaryField  -> Dirichlet 边界 φ = φ_exact
     *
     *     0/omega
     *         internalField  -> 源项 ω(c)
     *         boundaryField  -> 通常和内部场一致，或按 patch 给出
     *
     * 字段含义或可选写法：
     *
     *     对 Poisson 方程，0/phi 的 internalField 只是线性求解器初始猜测，
     *     不像前三题那样代表显式时间推进的初始条件。
     *
     *     0/phi/boundaryField 中：
     *
     *         type fixedValue
     *
     *     对应数学上的 Dirichlet 条件：
     *
     *         φ|boundary = φ_exact|boundary
     *
     *     如果以后换成 Neumann 边界，则可能使用：
     *
     *         type fixedGradient
     *         type zeroGradient
     *
     *     但第四题原题给的是 Dirichlet 边界，因此本验证案例应使用 fixedValue。
     */

    // 对应数学公式：
    //     φ(c)
    //
    // OpenFOAM 接口：
    //     volScalarField
    //     IOobject::MUST_READ  -> 必须从当前时间目录读取字段
    //     IOobject::AUTO_WRITE -> 求解后自动写出字段
    //
    // cases 中相关文件和字段：
    //     0/phi
    //         dimensions    -> 通常为 [0 0 0 0 0 0 0]
    //         internalField -> φ 的初始猜测
    //         boundaryField -> φ 的边界条件
    volScalarField phi
    (
        IOobject
        (
            solutionName,
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // 对应数学公式：
    //     ω(c)
    //
    // OpenFOAM 接口：
    //     volScalarField
    //     IOobject::MUST_READ -> 必须从当前时间目录读取源项
    //
    // cases 中相关文件和字段：
    //     0/omega
    //         dimensions    -> 如果 φ 无量纲，则 ω 的量纲是 [0 -2 0 0 0 0 0]
    //         internalField -> ω(c) = -2 pi^2 cos(pi x_c) cos(pi y_c)
    //
    // 字段含义或可选写法：
    //     为了让求解器通用，这里不把第四题解析源项硬编码在 C++ 中。
    //     解析函数由 prepare 脚本或手写 0/omega 离散到单元中心。
    volScalarField omega
    (
        IOobject
        (
            sourceName,
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Student Poisson solver" << nl
        << "  case                     = " << runTime.caseName() << nl
        << "  time                     = " << runTime.name() << nl
        << "  cells                    = " << mesh.nCells() << nl
        << "  solutionField            = " << solutionName << nl
        << "  sourceField              = " << sourceName << nl
        << "  nNonOrthogonalCorrectors = "
        << nNonOrthogonalCorrectors << nl
        << "  phi dimensions           = " << phi.dimensions() << nl
        << "  omega dimensions         = " << omega.dimensions() << nl
        << "  phi min                  = " << min(phi).value() << nl
        << "  phi max                  = " << max(phi).value() << nl
        << "  omega min                = " << min(omega).value() << nl
        << "  omega max                = " << max(omega).value() << nl
        << endl;

    /*
     * 第二部分：空间离散和方程组装
     *
     * 对应数学公式：
     *
     *     ∇²φ = ω
     *
     * 控制体 c 上的有限体积形式：
     *
     *     Σ over faces f of cell c [ ∇φ(f) · S(c,f) ] = V(c) ω(c)
     *
     * 内部面 f：
     *
     *     f connects owner cell O(f) and neighbour cell N(f)
     *
     *     正交主项大致对应：
     *
     *         ∇φ(f) · S(f) ≈ D(f) [ φ(N(f)) - φ(O(f)) ]
     *
     *     非正交修正对应：
     *
     *         ∇φ(f) · S(f)
     *             ≈ D(f) [ φ(N(f)) - φ(O(f)) ] + C(f)
     *
     * 边界面 f：
     *
     *     f belongs to one boundary patch
     *
     *     Dirichlet 边界给出：
     *
     *         φ(f) = φ_boundary(f)
     *
     * OpenFOAM 接口：
     *
     *     fvm::laplacian(phi)
     *
     *     这个接口返回 fvScalarMatrix，不是返回一个 volScalarField。
     *     它会把上面的面通量离散关系组装为：
     *
     *         A φ = b
     *
     * cases 中相关文件和字段：
     *
     *     system/fvSchemes
     *         gradSchemes
     *         {
     *             default Gauss linear;
     *         }
     *
     *         laplacianSchemes
     *         {
     *             laplacian(phi) Gauss linear corrected;
     *         }
     *
     *         snGradSchemes
     *         {
     *             default corrected;
     *         }
     *
     * 字段含义或可选写法：
     *
     *     Gauss
     *         -> 有限体积面通量求和
     *
     *     linear
     *         -> 面值或系数线性插值
     *
     *     corrected
     *         -> 根据网格非正交性加入修正项 C(f)
     *
     *     uncorrected
     *         -> 不做非正交修正
     *
     * 内部面/边界面在哪里分辨：
     *
     *     本求解器不手写 if internalFace / if boundaryFace。
     *
     *     fvm::laplacian(phi) 会根据 mesh 拓扑自动分辨：
     *
     *         内部面：
     *             constant/polyMesh/owner
     *             constant/polyMesh/neighbour
     *
     *         边界面：
     *             constant/polyMesh/boundary
     *             0/phi/boundaryField
     *
     *     有 neighbour cell 的面是内部面。
     *     没有 neighbour cell、属于某个 patch 的面是边界面。
     */

    /*
     * 第三部分：线性方程求解
     *
     * 对应数学公式：
     *
     *     A φ = b
     *
     *     residual = b - A φ
     *
     * OpenFOAM 接口：
     *
     *     fvScalarMatrix phiEqn
     *     phiEqn.solve()
     *     phi.correctBoundaryConditions()
     *
     * cases 中相关文件和字段：
     *
     *     system/fvSolution
     *         solvers
     *         {
     *             phi
     *             {
     *                 solver          PCG;
     *                 preconditioner  DIC;
     *                 tolerance       1e-12;
     *                 relTol          0;
     *             }
     *         }
     *
     * 字段含义或可选写法：
     *
     *     solver
     *         -> 线性方程组求解器，例如 PCG、PBiCGStab、smoothSolver
     *
     *     preconditioner
     *         -> 预条件器，例如 DIC、DILU
     *
     *     tolerance / relTol
     *         -> 线性求解器停止准则
     *
     *     这些参数不写在 C++ 里，而由 cases 中的 fvSolution 决定。
     */

    for
    (
        label nonOrth = 0;
        nonOrth <= nNonOrthogonalCorrectors;
        ++nonOrth
    )
    {
        Info<< "Solving Poisson equation, non-orthogonal correction "
            << nonOrth << " of " << nNonOrthogonalCorrectors << endl;

        // 对应数学公式：
        //     ∇²φ = ω
        //
        //     有限体积离散：
        //
        //         Σ over faces f of cell c [ ∇φ(f) · S(c,f) ] = V(c) ω(c)
        //
        //     线性代数形式：
        //
        //         A φ = b
        //
        // OpenFOAM 接口：
        //     fvm::laplacian(phi) == omega
        //
        // cases 中相关文件和字段：
        //     system/fvSchemes/laplacianSchemes/laplacian(phi)
        //     system/fvSchemes/snGradSchemes/default
        //     0/phi/boundaryField
        //     0/omega/internalField
        //
        // 字段含义或可选写法：
        //     这行是第四题求解器真正的核心。
        //     它不是显式残差公式，也不是手动更新 φ。
        //     OpenFOAM 会根据 cases 中的 mesh、scheme、boundary condition
        //     自动组装 fvScalarMatrix。
        fvScalarMatrix phiEqn
        (
            fvm::laplacian(phi) == omega
        );

        // 对应数学公式：
        //     求解 A φ = b
        //
        // OpenFOAM 接口：
        //     phiEqn.solve()
        //
        // cases 中相关文件和字段：
        //     system/fvSolution/solvers/phi
        //
        // 字段含义或可选写法：
        //     线性求解器输出的 residual 是代数残差，
        //     和前三题中手动定义的显式空间残差不是同一个用法。
        phiEqn.solve();

        // 对应数学公式：
        //     φ on boundary = φ_boundary
        //
        // OpenFOAM 接口：
        //     phi.correctBoundaryConditions()
        //
        // cases 中相关文件和字段：
        //     0/phi/boundaryField
        //
        // 字段含义或可选写法：
        //     求解后刷新 patch 上的 fixedValue、zeroGradient、
        //     fixedGradient 等边界字段。
        phi.correctBoundaryConditions();
    }

    Info<< "Solved phi min = " << min(phi).value() << nl
        << "Solved phi max = " << max(phi).value() << nl
        << endl;

    // 对应数学公式：
    //     无。这里写出数值解 φ(c)。
    //
    // OpenFOAM 接口：
    //     phi.write()
    //
    // cases 中相关文件和字段：
    //     system/controlDict
    //         writeFormat
    //         writePrecision
    //
    // 字段含义或可选写法：
    //     第四题是稳态问题，求解完成后直接写当前时间目录下的 phi。
    phi.write();

    Info<< "End" << nl << endl;

    return 0;
}

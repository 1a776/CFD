/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    projectionFoamStudent

Description
    第五题：不可压 Navier-Stokes 方程的传统压力投影法。

    本文件中的每一个主要代码块都按以下顺序说明：

        对应数学公式
        OpenFOAM 接口
        cases 中的入口文件和字段
        代码实际作用

    本实现是独立的 projection application，不是 PISO 的别名。
    时间步内的核心流程是：

        预测速度：

            (U* − Uⁿ) / Δt
          + ∇·(Uⁿ U*)
          − ∇·(ν∇U*) = 0

        预测通量：
            Φ* = U* · S

        压力投影：
            ∇·(Δt ∇p) = ∇·Φ*

        通量修正：
            Φⁿ⁺¹ = Φ* − flux(Δt∇p)

        速度修正：
            Uⁿ⁺¹ = U* − Δt∇p

    动量方程和压力方程由 fvm 接口组装为线性方程组；
    离散格式由 system/fvSchemes 选择；
    线性求解器由 system/fvSolution 选择。
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "dictionary.H"
#include "pisoControl.H"
#include "pressureReference.H"
#include "findRefCell.H"
#include "constrainPressure.H"
#include "adjustPhi.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"

#include "fvcFlux.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvcDdt.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

using namespace Foam;


int main(int argc, char *argv[])
{
    // 对应数学公式：
    //     无；解析命令行和 -case 路径。
    //
    // OpenFOAM 接口：
    //     setRootCase.H
    //
    // cases 入口：
    //     命令行参数 `projectionFoamStudent -case <caseDir>`。
    #include "setRootCase.H"

    // 对应数学公式：
    //     tⁿ、Δt、t_end。
    //
    // OpenFOAM 接口：
    //     createTime.H
    //
    // cases 入口：
    //     system/controlDict 的 startTime、deltaT、endTime、
    //     writeControl 和 writeInterval。
    #include "createTime.H"

    // 对应数学公式：
    //     Ω、V₍c₎、S₍cf₎、owner、neighbour。
    //
    // OpenFOAM 接口：
    //     createMesh.H -> fvMesh mesh
    //
    // cases 入口：
    //     constant/polyMesh/points、faces、owner、neighbour、boundary。
    //     blockMeshDict 或 Gmsh 只负责生成这些最终网格文件。
    #include "createMesh.H"

    // 对应数学公式：
    //     投影法可以设置多个压力投影校正：
    //         pᵏ⁺¹ → Uᵏ⁺¹ → pᵏ⁺²
    //
    // OpenFOAM 接口：
    //     pisoControl
    //
    // cases 入口：
    //     system/fvSolution 的 PISO 子字典，至少可设置：
    //         nCorrectors
    //         nNonOrthogonalCorrectors
    //
    // 代码作用：
    //     这里只借用 OpenFOAM 的循环控制接口；
    //     压力系数仍明确使用 Δt，而不是 PISO 动量矩阵的 rAU。
    pisoControl projection(mesh);

    // 对应数学公式：
    //     ν：运动黏度；
    //     p：压力；
    //     U：速度；
    //     Φ：面体积通量。
    //
    // OpenFOAM 接口：
    //     physicalProperties、volScalarField、volVectorField、createPhi.H
    //
    // cases 入口：
    //     constant/physicalProperties:
    //         nu  <value>;
    //     0/U:
    //         internalField、boundaryField；
    //     0/p:
    //         internalField、boundaryField；
    //     0/phi:
    //         可由 createPhi.H 创建，通常不要求手写。
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
    //     ΔUⁿ = ||Uⁿ - Uⁿ⁻¹||∞
    //
    // OpenFOAM 接口：
    //     volVectorField
    //
    // cases 入口：
    //     system/controlDict 中的稳态判据参数控制阈值；
    //     U 的当前值来自上一个时间步。
    //
    // 代码作用：
    //     保存当前时间步开始前的速度场。
    //     这与 max |div(U)| 不同：前者判断是否稳态，后者判断是否无散。
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

    // 对应数学公式：
    //     压力只由梯度决定，因此纯 Neumann 压力存在常数零空间。
    //
    // OpenFOAM 接口：
    //     setRefCell、pisoControl::dict、fvScalarMatrix::setReference
    //
    // cases 入口：
    //     system/fvSolution/PISO：
    //         pRefCell  0;
    //         pRefValue 0;
    //
    // 代码作用：
    //     给压力矩阵固定一个参考值，避免线性系统奇异。
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, projection.dict(), pRefCell, pRefValue);

    Info<< "\nStarting projectionFoamStudent\n"
        << "  case = " << runTime.caseName() << nl
        << "  cells = " << mesh.nCells() << nl
        << "  nu = " << nu.value() << nl
        << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        UPrevious = U;

        // 对应数学公式：
        //     (U* − Uⁿ) / Δt
        //   + ∇·(Uⁿ U*)
        //   − ∇·(ν∇U*) = 0
        //
        // 对一个控制体 Ω_c 的积分形式：
        //
        //     V₍c₎ (U*₍c₎ − Uⁿ₍c₎) / Δt
        //   + Σ_f Φⁿ₍f₎ U*₍f₎
        //   − Σ_f (ν∇U*)₍f₎ · S₍cf₎ = 0
        //
        // OpenFOAM 接口：
        //     fvm::ddt、fvm::div、fvm::laplacian、fvVectorMatrix
        //
        // cases 入口：
        //     system/fvSchemes：
        //         ddtSchemes
        //         divSchemes/div(phi,U)
        //         laplacianSchemes/laplacian(nu,U)
        //         gradSchemes、snGradSchemes
        //     system/fvSolution/solvers/U：
        //         预测动量方程的线性求解器和容差。
        //
        // 代码作用：
        //     先复制 U^n 得到预测场 U*，再由 fvm 接口直接组装
        //     A_U U* = b_U
        //
        //     当前压力不放入预测方程，
        //     压力由下一步投影恢复。
        volVectorField UStar
        (
            IOobject
            (
                "UStar",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U
        );

        fvVectorMatrix UStarEqn
        (
            fvm::ddt(UStar)
          + fvm::div(phi, UStar)
          - fvm::laplacian(nu, UStar)
        );

        // UStarEqn 的源项没有额外显式右端，因此直接求解：
        //
        //     A_U U* = b_U
        //
        // `solve(fvMatrix)` 是 OpenFOAM 对矩阵方程的标准接口。
        solve(UStarEqn);

        // 对应数学公式：
        //     Φ*₍f₎ = U*₍f₎ · S₍f₎。
        //
        //     内部面：
        //         f ∈ internalFaces，面两侧由 owner/neighbour 连接；
        //
        //     边界面：
        //         f ∈ boundary patch，面值由 boundaryField 和 patch 类型给出。
        //
        // OpenFOAM 接口：
        //     fvc::flux(UStar)
        //
        // cases 入口：
        //     system/fvSchemes/interpolationSchemes；
        //     0/U 的 boundaryField；
        //     constant/polyMesh/boundary。
        //
        // 代码作用：
        //     把预测单元速度转换为预测面通量；
        //     内部面和边界面由 fvMesh 的 owner/neighbour/boundary 自动区分。
        surfaceScalarField phiStar
        (
            IOobject
            (
                "phiStar",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::flux(UStar)
        );

        // 对应数学公式：
        //     α = Δt。
        //
        // OpenFOAM 接口：
        //     volScalarField
        //
        // cases 入口：
        //     system/controlDict/deltaT；
        //     runTime.deltaT() 是当前时间步实际使用的 Δt。
        //
        // 代码作用：
        //     用带量纲的体场作为压力方程的扩散系数，
        //     使 fvm::laplacian(alpha,p) 生成 Δt∇²p。
        volScalarField dtCoeff
        (
            IOobject
            (
                "dtCoeff",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "dtCoeff",
                dimTime,
                runTime.deltaTValue()
            )
        );

        // 对应数学公式：
        //     ∇·(Δt∇p) = ∇·Φ*。
        //
        // 对一个控制体 Ω_c 的离散形式：
        //
        //     Σ_f (Δt∇p)_f · S_cf
        //       = Σ_f Φ*_f
        //
        // OpenFOAM 接口：
        //     fvm::laplacian、fvc::div、fvScalarMatrix
        //
        // cases 入口：
        //     system/fvSchemes：
        //         laplacianSchemes/laplacian(dtCoeff,p)
        //         snGradSchemes
        //     system/fvSolution/solvers/p：
        //         压力投影方程的线性求解器、tolerance、relTol。
        //     0/p/boundaryField：
        //         压力边界条件。
        //
        // 代码作用：
        //     压力方程的右端是预测通量的离散散度；
        //     非正交修正由 OpenFOAM 的 fvSchemes 和 projection 控制器处理。
        while (projection.correct())
        {
            while (projection.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(dtCoeff, p)
                 ==
                    fvc::div(phiStar)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                // 对应数学公式：
                //     Φⁿ⁺¹₍f₎ = Φ*₍f₎ − (Δt∇p)₍f₎ · S₍f₎。
                //
                // OpenFOAM 接口：
                //     fvScalarMatrix::flux()
                //
                // cases 入口：
                //     constant/polyMesh/boundary；
                //     system/fvSchemes 中的 laplacian(dtCoeff,p)；
                //     0/p 的边界条件。
                //
                // 代码作用：
                //     pEqn.flux() 是本次压力方程对应的面通量；
                //     只在最后一次非正交校正后更新 phi，避免中间迭代覆盖。
                if (projection.finalNonOrthogonalIter())
                {
                    phi = phiStar - pEqn.flux();
                }
            }

            // 对应数学公式：
            //     Uⁿ⁺¹ = U* − Δt∇p。
            //
            // OpenFOAM 接口：
            //     fvc::grad、volVectorField::correctBoundaryConditions
            //
            // cases 入口：
            //     system/fvSchemes/gradSchemes；
            //     0/p/boundaryField；
            //     0/U/boundaryField。
            //
            // 代码作用：
            //     用压力梯度去除预测速度的离散散度，
            //     并重新应用速度边界条件。
            U = UStar - dtCoeff*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        // 对应数学公式：
        //     ∇·Uⁿ⁺¹ ≈ 0。
        //
        // OpenFOAM 接口：
        //     fvc::div(phi)
        //
        // cases 入口：
        //     constant/polyMesh/boundary 和 0/U、0/p 的边界条件。
        //
        // 代码作用：
        //     输出投影后的连续性检查量，便于判断压力校正是否生效。
        const volScalarField continuity
        (
            fvc::div(phi)
        );

        tmp<volScalarField> tContinuityMag(mag(continuity));
        const volScalarField& continuityMag = tContinuityMag();

        tmp<volScalarField> tVelocityChange(mag(U - UPrevious));
        const volScalarField& velocityChange = tVelocityChange();
        const scalar maxMassResidual =
            gMax(continuityMag.primitiveField());
        const scalar maxVelocityChange =
            gMax(velocityChange.primitiveField());

        Info<< "  max |div(U)| = "
            << maxMassResidual
            << nl
            << "  max |U - Uprevious| = "
            << maxVelocityChange
            << nl << endl;

        runTime.write();

        // 对应数学公式：
        //     max |Uⁿ - Uⁿ⁻¹| <= ε_U
        //     max |div(Uⁿ)|   <= ε_mass
        //
        //     两个条件连续满足 requiredSteadySteps 次，
        //     且当前步数不小于 minimumSteadySteps，才判定为稳态。
        //
        // OpenFOAM 接口：
        //     Time::timeIndex()
        //     Time::writeAndEnd()
        //
        // cases 入口：
        //     system/controlDict：
        //         steadyStateControl
        //         steadyVelocityTol
        //         steadyMassTol
        //         minimumSteadySteps
        //         requiredSteadySteps
        //
        // 代码作用：
        //     将“满足不可压缩”与“达到稳态”分开判断。
        //     writeAndEnd() 会保留当前稳态场后结束时间循环。
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
    }

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //

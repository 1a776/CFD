# UDF 求解器说明

本目录存放本项目自行开发的 OpenFOAM 求解器，不修改 `/opt/openfoam14/src` 中的 OpenFOAM 源码。

## 求解器族

```text
UDF/
└── solver/
    ├── 01_advection_equation/
    │   └── explicitAdvectionFoamStudent/
    │       ├── explicitAdvectionFoamStudent.C
    │       └── Make/
    │           ├── files
    │           └── options
    ├── 02_diffusion_equation/
    │   └── explicitDiffusionFoamStudent/
    │       ├── explicitDiffusionFoamStudent.C
    │       └── Make/
    │           ├── files
    │           └── options
    └── 03_advection_diffusion_equation/
        └── explicitAdvectionDiffusionFoamStudent/
            ├── explicitAdvectionDiffusionFoamStudent.C
            └── Make/
                ├── files
                └── options
    └── 04_poisson_equation/
        └── poissonFoamStudent/
            ├── poissonFoamStudent.C
            └── Make/
                ├── files
                └── options
```

## 当前求解器

### `01_advection_equation`

二维线性对流方程显式有限体积求解器：

```text
∂T/∂t + ∇·(U T) = 0
```

核心 OpenFOAM 接口：

```text
fvc::flux(U)
fvc::div(phi, T)
```

编译产物：

```text
build/01_advection_equation/bin/explicitAdvectionFoamStudent
```

### `02_diffusion_equation`

二维扩散方程显式有限体积求解器：

```text
∂φ/∂t - ∇·(μ∇φ) = 0
```

核心 OpenFOAM 接口：

```text
fvc::laplacian(mu, phi)
```

编译产物：

```text
build/02_diffusion_equation/bin/explicitDiffusionFoamStudent
```

### `03_advection_diffusion_equation`

第三题对流扩散方程显式有限体积求解器：

```text
∂φ/∂t + ∇·(Uφ) - ∇·(μ∇φ) = 0
```

显式右端残差：

```text
R = -∇·(Uφ) + ∇·(μ∇φ)
```

核心 OpenFOAM 接口：

```text
fvc::flux(U)
-fvc::div(faceFlux, phi)
fvc::laplacian(mu, phi)
```

源码注释按第二题求解器标准书写：每个关键代码块都说明对应数学公式、OpenFOAM 接口、
case 文件入口和字段含义。边界面不在源码里手工遍历，而由 `fvc::div`、
`fvc::laplacian` 根据 `constant/polyMesh/boundary`、`0/U` 和 `0/phi` 自动处理。

编译产物：

```text
build/03_advection_diffusion_equation/bin/explicitAdvectionDiffusionFoamStudent
```

### `04_poisson_equation`

第四题 Poisson 方程有限体积求解器：

```text
∇²φ = ω
```

核心 OpenFOAM 接口：

```text
fvm::laplacian(phi) == omega
phiEqn.solve()
```

注释风格：

- 每个关键代码块前都说明对应数学公式；
- 每个关键代码块前都说明 OpenFOAM 接口；
- 每个关键代码块前都说明 cases 中相关文件和字段；
- 每个关键代码块前都说明字段含义或可选写法。

编译产物：

```text
build/04_poisson_equation/bin/poissonFoamStudent
```

## 编译

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
source /opt/openfoam14/etc/bashrc
sh scripts/build_student_solver.sh
```

编译脚本会依次编译：

```text
build/01_advection_equation/bin/explicitAdvectionFoamStudent
build/02_diffusion_equation/bin/explicitDiffusionFoamStudent
build/03_advection_diffusion_equation/bin/explicitAdvectionDiffusionFoamStudent
build/04_poisson_equation/bin/poissonFoamStudent
```

## 与案例的关系

求解器源码、案例、配置、数据和图片按求解器族命名空间组织。例如：

```text
cases/01_advection_equation/
scripts/configs/01_advection_equation/
data/01_advection_equation/
figures/01_advection_equation/
```

后续第三题案例建议放在：

```text
cases/03_advection_diffusion_equation/
scripts/configs/03_advection_diffusion_equation/
data/03_advection_diffusion_equation/
figures/03_advection_diffusion_equation/
```

第三题第一个已集成的验证案例为周期正弦波：

```text
scripts/configs/03_advection_diffusion_equation/01_sine_wave_quad_upwind.json
scripts/configs/03_advection_diffusion_equation/02_sine_wave_tri_upwind.json
```

第三题第二个已集成的验证案例为旋转尖峰：

```text
scripts/configs/03_advection_diffusion_equation/03_rotating_peak_quad_upwind.json
scripts/configs/03_advection_diffusion_equation/04_rotating_peak_tri_upwind.json
scripts/configs/03_advection_diffusion_equation/05_rotating_peak_quad_analyticDirichlet_upwind.json
scripts/configs/03_advection_diffusion_equation/06_rotating_peak_tri_analyticDirichlet_upwind.json
```

两组配置分别使用结构化四边形网格和 Gmsh 生成的三角形棱柱网格；均使用
$U=(1,1,0)$、$\mu=1$、x/y 周期边界、一阶迎风对流格式和
`Gauss linear corrected` 扩散格式，在 $t=1$ 与解析解

$$\phi(x,y,t)=\exp(-8\pi^2\mu t)\sin(2\pi[x+y-(u+v)t])$$

比较归一化 $L_1$、$L_2$、$L_\infty$ 误差。生成后的 `0.orig/phi`、`0.orig/U`、
`constant/transportProperties`、`system/fvSchemes`、`system/controlDict` 和
`system/fvSolution` 都在文件内写明了公式、求解器接口和对应字段。

### 已完成验证

已在 2026-08-26 对该案例的四边形 $N=10,20,40,80$ 全部完成；并在
2026-08-27 对三角形 $N=10,20,40,80$ 全部完成。三角形版本采用同一配置接口，
其中 `Allrun` 在 Gmsh 网格转换、`checkMesh` 和单元中心写出完成后，再按真实三角形
单元中心写入 `0.orig/phi`。

- 编译：`wmake UDF/solver/03_advection_diffusion_equation/explicitAdvectionDiffusionFoamStudent`；
- 可执行文件检查：`explicitAdvectionDiffusionFoamStudent -help`；
- 网格检查：四边形和三角形的四个分辨率均报告 `Mesh OK.`；
- 时间推进：八个案例均严格写出 `t=1` 的最终场，且求解日志以 `End` 结束；
- 三角形单元数：$N=10,20,40,80$ 分别为 $200,800,3200,12800$；
- 后处理：两类网格均已生成逐网格场对比、对角线剖面、振幅历史、稳定性历史，以及跨网格收敛表和图；
- 三角形绝对 $L_1$ 残余：$7.54\times10^{-14}$、$5.40\times10^{-14}$、
  $4.68\times10^{-14}$、$5.49\times10^{-14}$，与最终解析场已经低于双精度有效
  判别尺度的预期一致。

该原题规定 $\mu=1$、$t=1$，所以精确正弦波振幅为
$\exp(-8\pi^2)\approx5.12\times10^{-35}$，远低于双精度舍入量级。求解稳定、
绝对场值接近零，但按题目指定的“以精确解范数归一化”的误差会被极小分母放大；
因此这组归一化误差及其观察收敛阶只应作为原题定义的复现记录，不能据此判断格式的
实际收敛性。单案例 `summary.json` 同时保存了绝对误差，供后续报告说明。

第三题第二个旋转尖峰案例已在 2026-08-27 完成四边形和三角形 $N=20,40,80$
验证。该案例使用区域 $[-1,1]^2$、速度 $U=(-y,x,0)$、扩散系数
$\epsilon=10^{-3}$，OpenFOAM 计算时间 $\tau$ 从 0 推进到 $2\pi$；解析解中的扩散
年龄为 $t_0+\tau$，其中 $t_0=\pi/2$，初始峰中心为 $(0,0.5)$。原题说明初始条件
和边界条件均由上述解析尖峰函数给出，因此严格复现应使用随时间变化的解析
Dirichlet 边界；当前案例先采用 `fixedValue 0` 近似，因为尖峰始终位于区域内部，
边界解析值为指数小量。

为保留原近似案例并给出更严格的原题版本，2026-08-27 新增两组解析 Dirichlet
边界案例：

- `05_rotating_peak_quad_analyticDirichlet_upwind`：四边形网格，侧边界
  `xMin/xMax/yMin/yMax` 使用 `codedFixedValue` 按解析尖峰函数实时计算；
- `06_rotating_peak_tri_analyticDirichlet_upwind`：三角形棱柱网格，侧边界同样使用
  `codedFixedValue`，三角形预处理拆分出 `Allrun.pre`，可单独完成 Gmsh、
  `gmshToFoam`、`createPatch`、`checkMesh`、单元中心/体积写出和初始场刷新。

解析边界对应公式为

$$\phi_b(x,y,\tau)=\frac{1}{4\pi\epsilon(t_0+\tau)}\exp\left(-\frac{(x-\hat{x})^2+(y-\hat{y})^2}{4\epsilon(t_0+\tau)}\right)$$

其中 $\hat{x}=x_0\cos\tau-y_0\sin\tau$，$\hat{y}=x_0\sin\tau+y_0\cos\tau$。
OpenFOAM 实现入口在生成的 `0.orig/phi` 的 `boundaryField`，关键字段为
`type codedFixedValue` 和 `name rotatingPeak_<patch>`。

- 四边形 $N=20,40,80$ 的归一化 $L_1$：`1.406084`、`1.184313`、`0.925823`；
- 三角形 $N=20,40,80$ 的归一化 $L_1$：`1.195397`、`0.964530`、`0.698101`；
- 两类网格的 `checkMesh` 均通过，求解日志均以 `End` 结束；
- 最终时间均严格为 `6.283185307179586`；
- 后处理已生成 `field_comparison.png`、`peak_profile.png`、
  `advection_diffusion_stability_history.png` 和原题要求的最终等值线图
  `contour_final.png`；
- 结果已写入 `data/03_advection_diffusion_equation/analysis/03_rotating_peak_quad_upwind/`
  和 `data/03_advection_diffusion_equation/analysis/04_rotating_peak_tri_upwind/`。

解析 Dirichlet 新案例已完成如下轻量验证：

- 四边形 `N20` 完整运行到 `t=2*pi`，四个 `codedFixedValue` 动态边界库均成功编译，
  求解日志以 `End` 结束，`contour_final.png` 已生成；
- 三角形 `N20` 完整运行到 `t=2*pi`，Gmsh 网格、`createPatch`、`checkMesh`、
  `codedFixedValue` 动态边界库和后处理均通过；
- 三角形 `N40/N80` 已执行 `Allrun.pre`，完成网格生成、`checkMesh` 和解析边界
  `0.orig/phi` 刷新，尚未完整求解。

2026-08-27 又对第三题全部脚本和产物做了一次一致性复查：

- `scripts/configs/03_advection_diffusion_equation/` 下 4 个 JSON 均能被解析，且均指向
  `solverFamily=03_advection_diffusion_equation` 和
  `solver=explicitAdvectionDiffusionFoamStudent`；
- `scripts/common/advection_diffusion_tools.py`、`case_config.py`、`foam_case.py`、
  `postprocess_case.py`、`study_analysis.py` 以及顶层运行脚本均通过 `py_compile`；
- 14 个已生成 case 的 `metadata.json` 与 `system/controlDict`、`system/fvSchemes`、
  `constant/transportProperties` 中的 `application`、`scalarField`、`endTime`、
  `div(faceFlux,phi)`、`laplacian(mu,phi)` 和 `mu` 一致；
- 临时生成的第三题 OpenFOAM 字典能够被 `foamDictionary` 读取到
  `divSchemes`、`laplacianSchemes`、`snGradSchemes`、`application` 和 `mu`；
- `sh scripts/build_student_solver.sh` 可重新编译第三题求解器；
- `git diff --check` 未发现空白或补丁格式问题。

典型第三题 case 中需要关注：

```text
0/phi
    internalField  -> 初始条件
    boundaryField  -> 标量边界条件

0/U
    internalField  -> 速度场
    boundaryField  -> 速度边界条件

constant/transportProperties
    mu             -> 扩散系数

system/fvSchemes
    div(faceFlux,phi)       -> 对流格式
    laplacian(mu,phi)       -> 扩散格式
    snGradSchemes           -> 非正交修正

system/controlDict
    application             -> explicitAdvectionDiffusionFoamStudent
    scalarField             -> phi
    velocityField           -> U
    advectionDiffusionCo    -> 显式稳定系数
    maxDeltaT               -> 最大时间步
```

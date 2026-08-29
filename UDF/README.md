# UDF 求解器说明

本目录存放本项目自行开发的 OpenFOAM 求解器。每个求解器只负责对应控制方程的离散和时间推进；具体实验案例、网格、边界条件、数值格式、分辨率、后处理数据和图片归档，均在下表中直接链接到 [`cases/README.md`](../cases/README.md) 的对应小节。


## 求解器与 cases/README 具体位置

| 对应题目 | 求解器可执行程序 | 对应案例的位置 |
|---|---|---|
| [1. 对流方程](../cases/README.md#advection-equation) | [`explicitAdvectionFoamStudent`](solver/01_advection_equation/explicitAdvectionFoamStudent/explicitAdvectionFoamStudent.C) | [正弦波平移](../cases/README.md#advection-sine-wave-translation)：[实验配置](../cases/README.md#advection-sine-wave-translation-config)、[数据归档](../cases/README.md#advection-sine-wave-translation-archive)<br>[复杂轮廓固体旋转](../cases/README.md#advection-solid-rotation)：[实验配置](../cases/README.md#advection-solid-rotation-config)、[数据归档](../cases/README.md#advection-solid-rotation-archive) |
| [2. 扩散方程](../cases/README.md#diffusion-equation) | [`explicitDiffusionFoamStudent`](solver/02_diffusion_equation/explicitDiffusionFoamStudent/explicitDiffusionFoamStudent.C) | [间断初值扩散](../cases/README.md#diffusion-discontinuous-initial)：[实验配置](../cases/README.md#diffusion-discontinuous-initial-config)、[数据归档](../cases/README.md#diffusion-discontinuous-initial-archive)<br>[Gaussian 扩散](../cases/README.md#diffusion-gaussian)：[实验配置](../cases/README.md#diffusion-gaussian-config)、[数据归档](../cases/README.md#diffusion-gaussian-archive) |
| [3. 对流-扩散方程](../cases/README.md#advection-diffusion-equation) | [`explicitAdvectionDiffusionFoamStudent`](solver/03_advection_diffusion_equation/explicitAdvectionDiffusionFoamStudent/explicitAdvectionDiffusionFoamStudent.C) | [正弦波平移](../cases/README.md#advection-diffusion-sine-wave-translation)：[实验配置](../cases/README.md#advection-diffusion-sine-wave-translation-config)、[数据归档](../cases/README.md#advection-diffusion-sine-wave-translation-archive)<br>[旋转尖峰平流扩散](../cases/README.md#advection-diffusion-rotating-peak)：[实验配置](../cases/README.md#advection-diffusion-rotating-peak-config)、[数据归档](../cases/README.md#advection-diffusion-rotating-peak-archive) |
| [4. Poisson 方程](../cases/README.md#poisson-equation) | [`poissonFoamStudent`](solver/04_poisson_equation/poissonFoamStudent/poissonFoamStudent.C) | [制造解 Poisson 方程](../cases/README.md#poisson-manufactured-solution)：[实验配置](../cases/README.md#poisson-manufactured-solution-config)、[数据归档](../cases/README.md#poisson-manufactured-solution-archive) |
| [5. Navier-Stokes 方程](../cases/README.md#navier-stokes-equation) | [`projectionFoamStudent`](solver/05_navier_stokes_equation/projectionFoamStudent/projectionFoamStudent.C) | 压力投影法对应 [方腔顶盖驱动流](../cases/README.md#navier-stokes-lid-driven-cavity)：[实验配置](../cases/README.md#navier-stokes-lid-driven-cavity-config)、[数据归档](../cases/README.md#navier-stokes-lid-driven-cavity-archive)<br>压力投影法对应 [等边三角腔顶盖驱动流](../cases/README.md#navier-stokes-triangular-cavity)：[实验配置](../cases/README.md#navier-stokes-triangular-cavity-config)、[数据归档](../cases/README.md#navier-stokes-triangular-cavity-archive) |
| [5. Navier-Stokes 方程](../cases/README.md#navier-stokes-equation) | [`pisoFoamStudent`](solver/06_piso_navier_stokes_equation/pisoFoamStudent/pisoFoamStudent.C) | PISO 法对应 [方腔顶盖驱动流](../cases/README.md#navier-stokes-lid-driven-cavity)：[实验配置](../cases/README.md#navier-stokes-lid-driven-cavity-config)、[数据归档](../cases/README.md#navier-stokes-lid-driven-cavity-archive)<br>PISO 法对应 [等边三角腔顶盖驱动流](../cases/README.md#navier-stokes-triangular-cavity)：[实验配置](../cases/README.md#navier-stokes-triangular-cavity-config)、[数据归档](../cases/README.md#navier-stokes-triangular-cavity-archive) |

## 题目说明

<a id="udf-advection-equation"></a>
### 1. 对流方程

`explicitAdvectionFoamStudent` 用于二维线性对流方程：

```text
∂T/∂t + ∇·(U T) = 0
```

对应题目包括 [正弦波平移](../cases/README.md#advection-sine-wave-translation) 和 [复杂轮廓固体旋转](../cases/README.md#advection-solid-rotation)。具体设置分别见正弦波平移的[实验配置](../cases/README.md#advection-sine-wave-translation-config)和[数据归档](../cases/README.md#advection-sine-wave-translation-archive)，以及复杂轮廓固体旋转的[实验配置](../cases/README.md#advection-solid-rotation-config)和[数据归档](../cases/README.md#advection-solid-rotation-archive)。

<a id="udf-diffusion-equation"></a>
### 2. 扩散方程

`explicitDiffusionFoamStudent` 用于二维扩散方程：

```text
∂φ/∂t - ∇·(μ∇φ) = 0
```

对应题目包括 [间断初值扩散](../cases/README.md#diffusion-discontinuous-initial) 和 [Gaussian 扩散](../cases/README.md#diffusion-gaussian)。具体设置分别见间断初值扩散的[实验配置](../cases/README.md#diffusion-discontinuous-initial-config)和[数据归档](../cases/README.md#diffusion-discontinuous-initial-archive)，以及 Gaussian 扩散的[实验配置](../cases/README.md#diffusion-gaussian-config)和[数据归档](../cases/README.md#diffusion-gaussian-archive)。

<a id="udf-advection-diffusion-equation"></a>
### 3. 对流-扩散方程

`explicitAdvectionDiffusionFoamStudent` 用于二维对流-扩散方程：

```text
∂φ/∂t + ∇·(Uφ) - ∇·(μ∇φ) = 0
```

对应题目包括 [正弦波平移](../cases/README.md#advection-diffusion-sine-wave-translation) 和 [旋转尖峰平流扩散](../cases/README.md#advection-diffusion-rotating-peak)。具体设置分别见正弦波平移的[实验配置](../cases/README.md#advection-diffusion-sine-wave-translation-config)和[数据归档](../cases/README.md#advection-diffusion-sine-wave-translation-archive)，以及旋转尖峰平流扩散的[实验配置](../cases/README.md#advection-diffusion-rotating-peak-config)和[数据归档](../cases/README.md#advection-diffusion-rotating-peak-archive)。

<a id="udf-poisson-equation"></a>
### 4. Poisson 方程

`poissonFoamStudent` 用于二维稳态 Poisson 方程：

```text
∇²φ = ω
```

对应题目为 [制造解 Poisson 方程](../cases/README.md#poisson-manufactured-solution)。具体解析解、源项、Dirichlet 边界、网格加密设置、误差范数和数据归档见该案例的[实验配置](../cases/README.md#poisson-manufactured-solution-config)和[数据归档](../cases/README.md#poisson-manufactured-solution-archive)。

<a id="udf-navier-stokes-equation"></a>
### 5. Navier-Stokes 方程

`projectionFoamStudent` 和 `pisoFoamStudent` 都用于二维不可压 Navier-Stokes 方程案例，区别在于压力-速度耦合算法：

| 求解器 | 算法 | 对应案例 |
|---|---|---|
| `projectionFoamStudent` | 压力投影法 | [方腔顶盖驱动流](../cases/README.md#navier-stokes-lid-driven-cavity)：[实验配置](../cases/README.md#navier-stokes-lid-driven-cavity-config)、[数据归档](../cases/README.md#navier-stokes-lid-driven-cavity-archive)<br>[等边三角腔顶盖驱动流](../cases/README.md#navier-stokes-triangular-cavity)：[实验配置](../cases/README.md#navier-stokes-triangular-cavity-config)、[数据归档](../cases/README.md#navier-stokes-triangular-cavity-archive) |
| `pisoFoamStudent` | PISO 法 | [方腔顶盖驱动流](../cases/README.md#navier-stokes-lid-driven-cavity)：[实验配置](../cases/README.md#navier-stokes-lid-driven-cavity-config)、[数据归档](../cases/README.md#navier-stokes-lid-driven-cavity-archive)<br>[等边三角腔顶盖驱动流](../cases/README.md#navier-stokes-triangular-cavity)：[实验配置](../cases/README.md#navier-stokes-triangular-cavity-config)、[数据归档](../cases/README.md#navier-stokes-triangular-cavity-archive) |

具体 Reynolds 数、混合网格设置、时间步长、稳态判据、中心线数据、Ghia 参考对比、主涡数据和图片归档，分别见方腔顶盖驱动流的[实验配置](../cases/README.md#navier-stokes-lid-driven-cavity-config)和[数据归档](../cases/README.md#navier-stokes-lid-driven-cavity-archive)，以及等边三角腔顶盖驱动流的[实验配置](../cases/README.md#navier-stokes-triangular-cavity-config)和[数据归档](../cases/README.md#navier-stokes-triangular-cavity-archive)。

## 编译

在项目根目录执行：

```bash
export OPENFOAM_BASHRC="${OPENFOAM_BASHRC:-/opt/openfoam14/etc/bashrc}"
. "$OPENFOAM_BASHRC"
sh scripts/build_student_solver.sh
```

编译脚本会依次生成以下可执行文件：

```text
build/01_advection_equation/bin/explicitAdvectionFoamStudent
build/02_diffusion_equation/bin/explicitDiffusionFoamStudent
build/03_advection_diffusion_equation/bin/explicitAdvectionDiffusionFoamStudent
build/04_poisson_equation/bin/poissonFoamStudent
build/05_navier_stokes_equation/bin/projectionFoamStudent
build/06_piso_navier_stokes_equation/bin/pisoFoamStudent
```

## 配置与结果索引

本 README 只说明求解器和题目的对应关系。各题目的完整实验配置、OpenFOAM 案例目录、分析数据、图片和报告链接以上表列出的 `cases/README.md` 具体小节为准。

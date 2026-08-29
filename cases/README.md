## `cases/` 案例目录说明
- `cases/` 是 OpenFOAM 原始案例目录；真正用于报告分析的数据和图片放在 `data/`、`figures/` 和 `report/` 中。

## 具体题目的相关设计的实验案例以及实验配置

### 1. 对流方程

### 案例一：正弦波平移

本案例对应二维对流方程的正弦波平移算例。所有实验统一放在同一物理域、同一速度场和同一时间控制下，只比较网格类型、空间格式和分辨率的变化。所有算例都采用单位正方形周期区域，速度场取 `U=(1,1,0)`，标量场为 `T`，终止时间为 `1.0`，最大 CFL 为 `0.2`。实验代号按罗马数字 I-IV 编排，其中 I、II 对应四边形网格，III、IV 对应三角形网格；I、III 是一阶迎风，II、IV 是线性迎风。

该组生成的图片主要包括数值场与精确解的对比图、对角线剖面对比图、振幅随时间变化图和 CFL 历史图，重点观察波形平移后是否保持幅值、是否出现数值耗散以及是否有非物理过冲。对应的分析数据主要包括 `L1`、`L2`、`Linf` 误差、`L1` 观察收敛阶、最终振幅、最大 CFL 和质量误差，其中 `L1` 与收敛阶最能直接反映网格加密后的精度变化，`L2` 和 `Linf` 则分别补充平均误差和峰值偏差。

#### 相关实验配置

| 实验代号 | 案例 | 网格类型 | 空间格式 | 分辨率 N | 固定设置 |
|---|---|---|---|---|---|
| I | `01_sine_wave_quad_upwind` | 四边形 | `Gauss upwind` | `10, 20, 40, 80` | 单位正方形周期区域、周期边界、`U=(1,1,0)`、`T`、`endTime=1.0`、`maxCo=0.2` |
| II | `02_sine_wave_quad_linearUpwind` | 四边形 | `Gauss linearUpwind grad(T)` | `10, 20, 40, 80` | 单位正方形周期区域、周期边界、`U=(1,1,0)`、`T`、`endTime=1.0`、`maxCo=0.2` |
| III | `03_sine_wave_tri_upwind` | 三角形 | `Gauss upwind` | `10, 20, 40, 80` | 单位正方形周期区域、周期边界、`U=(1,1,0)`、`T`、`endTime=1.0`、`maxCo=0.2` |
| IV | `03_sine_wave_tri_linearUpwind` | 三角形 | `Gauss linearUpwind grad(T)` | `10, 20, 40, 80` | 单位正方形周期区域、周期边界、`U=(1,1,0)`、`T`、`endTime=1.0`、`maxCo=0.2` |

#### 相关实验的数据归档

- I `01_sine_wave_quad_upwind`
  - OpenFOAM目录: [`cases/01_advection_equation/01_sine_wave_quad/`](./01_advection_equation/01_sine_wave_quad/)
  - 数据: [`data/01_advection_equation/cases/01_sine_wave_quad/`](../data/01_advection_equation/cases/01_sine_wave_quad/)
  - 图片: [`figures/01_advection_equation/cases/01_sine_wave_quad/`](../figures/01_advection_equation/cases/01_sine_wave_quad/)
  - 数据汇总: [`data/01_advection_equation/analysis/01_sine_wave_quad/analysis.md`](../data/01_advection_equation/analysis/01_sine_wave_quad/analysis.md)
  - 报告: [`report/01_advection_equation/report.md`](../report/01_advection_equation/report.md)
  - 证据索引: [`report/01_advection_equation/evidence_index.md`](../report/01_advection_equation/evidence_index.md)
- II `02_sine_wave_quad_linearUpwind`
  - OpenFOAM目录: [`cases/01_advection_equation/02_sine_wave_quad_linearUpwind/`](./01_advection_equation/02_sine_wave_quad_linearUpwind/)
  - 数据: [`data/01_advection_equation/cases/02_sine_wave_quad_linearUpwind/`](../data/01_advection_equation/cases/02_sine_wave_quad_linearUpwind/)
  - 图片: [`figures/01_advection_equation/cases/02_sine_wave_quad_linearUpwind/`](../figures/01_advection_equation/cases/02_sine_wave_quad_linearUpwind/)
  - 数据汇总: [`data/01_advection_equation/analysis/02_sine_wave_quad_linearUpwind/analysis.md`](../data/01_advection_equation/analysis/02_sine_wave_quad_linearUpwind/analysis.md)
  - 报告: [`report/01_advection_equation/report.md`](../report/01_advection_equation/report.md)
  - 证据索引: [`report/01_advection_equation/evidence_index.md`](../report/01_advection_equation/evidence_index.md)
- III `03_sine_wave_tri_upwind`
  - OpenFOAM目录: [`cases/01_advection_equation/03_sine_wave_tri_upwind/`](./01_advection_equation/03_sine_wave_tri_upwind/)
  - 数据: [`data/01_advection_equation/cases/03_sine_wave_tri_upwind/`](../data/01_advection_equation/cases/03_sine_wave_tri_upwind/)
  - 图片: [`figures/01_advection_equation/cases/03_sine_wave_tri_upwind/`](../figures/01_advection_equation/cases/03_sine_wave_tri_upwind/)
  - 数据汇总: [`data/01_advection_equation/analysis/03_sine_wave_tri_upwind/analysis.md`](../data/01_advection_equation/analysis/03_sine_wave_tri_upwind/analysis.md)
  - 报告: [`report/01_advection_equation/report.md`](../report/01_advection_equation/report.md)
  - 证据索引: [`report/01_advection_equation/evidence_index.md`](../report/01_advection_equation/evidence_index.md)
- IV `03_sine_wave_tri_linearUpwind`
  - OpenFOAM目录: [`cases/01_advection_equation/03_sine_wave_tri_linearUpwind/`](./01_advection_equation/03_sine_wave_tri_linearUpwind/)
  - 数据: [`data/01_advection_equation/cases/03_sine_wave_tri_linearUpwind/`](../data/01_advection_equation/cases/03_sine_wave_tri_linearUpwind/)
  - 图片: [`figures/01_advection_equation/cases/03_sine_wave_tri_linearUpwind/`](../figures/01_advection_equation/cases/03_sine_wave_tri_linearUpwind/)
  - 数据汇总: [`data/01_advection_equation/analysis/03_sine_wave_tri_linearUpwind/analysis.md`](../data/01_advection_equation/analysis/03_sine_wave_tri_linearUpwind/analysis.md)
  - 报告: [`report/01_advection_equation/report.md`](../report/01_advection_equation/report.md)
  - 证据索引: [`report/01_advection_equation/evidence_index.md`](../report/01_advection_equation/evidence_index.md)

### 案例二：复杂轮廓固体旋转

本案例对应二维对流方程的复杂轮廓固体旋转算例。所有实验同样放在单位正方形周期域内，速度场取绕 `(0.5,0.5)` 的刚体旋转，角速度为 `1`，终止时间为 `2π`，最大 CFL 为 `0.2`。初始场由切口圆盘、圆锥和光滑峰值三部分构成。实验代号按罗马数字 V-VI 编排，其中 V 对应四边形网格，VI 对应三角形网格；两组都采用一阶迎风格式，通过 `N=50,100,200` 的结果比较轮廓保持、数值耗散和质量守恒。

该组生成的图片主要包括初始场与最终场的对比图、最终等值线图和 CFL 历史图，重点观察切口是否闭合、圆盘边界是否变钝、峰值是否被抹平，以及一圈旋转后整体轮廓是否回到原位。对应的分析数据主要包括一圈旋转后的 `cycleL1AgainstInitial`、最大 CFL、归一化质量误差、最终最大值和单元数；其中 `cycleL1AgainstInitial` 反映整体形状偏差，质量误差反映守恒性，最终最大值则反映数值耗散对峰值的削弱程度。

#### 相关实验配置

| 实验代号 | 案例 | 网格类型 | 空间格式 | 分辨率 N | 固定设置 |
|---|---|---|---|---|---|
| V | `04_solid_rotation_quad_upwind` | 四边形 | `Gauss upwind` | `50, 100, 200` | 单位正方形区域、刚体旋转速度场 `u=(0.5-y, x-0.5, 0)`、旋转中心 `(0.5,0.5)`、角速度 `1`、`t=2π`、`maxCo=0.2`、初始场为切口圆盘/圆锥/光滑峰值、外边界标量为零 |
| VI | `04_solid_rotation_tri_upwind` | 三角形 | `Gauss upwind` | `50, 100, 200` | 单位正方形区域、刚体旋转速度场 `u=(0.5-y, x-0.5, 0)`、旋转中心 `(0.5,0.5)`、角速度 `1`、`t=2π`、`maxCo=0.2`、初始场为切口圆盘/圆锥/光滑峰值、外边界标量为零 |

#### 相关实验的数据归档

- V `04_solid_rotation_quad_upwind`
  - OpenFOAM目录: [`cases/01_advection_equation/04_solid_rotation_quad_upwind/`](./01_advection_equation/04_solid_rotation_quad_upwind/)
  - 数据: [`data/01_advection_equation/cases/04_solid_rotation_quad_upwind/`](../data/01_advection_equation/cases/04_solid_rotation_quad_upwind/)
  - 图片: [`figures/01_advection_equation/cases/04_solid_rotation_quad_upwind/`](../figures/01_advection_equation/cases/04_solid_rotation_quad_upwind/)
  - 报告: [`report/01_advection_equation/report.md`](../report/01_advection_equation/report.md)
  - 证据索引: [`report/01_advection_equation/evidence_index.md`](../report/01_advection_equation/evidence_index.md)
- VI `04_solid_rotation_tri_upwind`
  - OpenFOAM目录: [`cases/01_advection_equation/04_solid_rotation_tri_upwind/`](./01_advection_equation/04_solid_rotation_tri_upwind/)
  - 数据: [`data/01_advection_equation/cases/04_solid_rotation_tri_upwind/`](../data/01_advection_equation/cases/04_solid_rotation_tri_upwind/)
  - 图片: [`figures/01_advection_equation/cases/04_solid_rotation_tri_upwind/`](../figures/01_advection_equation/cases/04_solid_rotation_tri_upwind/)
  - 报告: [`report/01_advection_equation/report.md`](../report/01_advection_equation/report.md)
  - 证据索引: [`report/01_advection_equation/evidence_index.md`](../report/01_advection_equation/evidence_index.md)

### 2. 扩散方程

### 案例一：间断初值扩散

本案例对应二维扩散方程的间断初值验证。计算域为 `[-5,5] \times [-5,5]`，初始场为中心方块指示函数，外边界采用齐次 Neumann 条件，扩散系数取 `\mu=1`，终止时间为 `0.2`。实验代号按罗马数字 I-II 编排，其中 I 对应四边形网格，II 对应三角形网格；两组都在 `N=10,20,40,80` 下比较扩散平滑过程、误差收敛和守恒性。

该组生成的图片主要包括数值场与解析解的对比图、典型剖面图、时间演化图，以及网格加密下的误差收敛曲线和观察收敛阶图。前两类图主要看间断是否被正确平滑、峰值是否衰减、剖面是否与解析解一致；后两类图主要看误差是否随网格加密而下降，以及收敛阶是否稳定。对应的分析数据主要包括各分辨率的 `L1`、`L2` 误差和由相邻网格计算得到的观察收敛阶；其中 `L1` 更直接反映整体偏差，`L2` 更能体现平均意义下的误差水平，观察收敛阶则用来判断网格加密后的精度提升趋势。

#### 相关实验配置

| 实验代号 | 案例 | 网格类型 | 空间格式 | 分辨率 N | 固定设置 |
|---|---|---|---|---|---|
| I | `01_discontinuous_quad` | 四边形 | `Gauss linear corrected` | `10, 20, 40, 80` | `solver=explicitDiffusionFoamStudent`，`phi`，`mu=1`，`endTime=0.2`，`diffusionCo=0.45`，`maxDeltaT=0.001`，中心方块初值，齐次 Neumann 边界 |
| II | `02_discontinuous_tri` | 三角形 | `Gauss linear corrected` | `10, 20, 40, 80` | `solver=explicitDiffusionFoamStudent`，`phi`，`mu=1`，`endTime=0.2`，`diffusionCo=0.45`，`maxDeltaT=0.001`，中心方块初值，齐次 Neumann 边界 |

#### 相关实验的数据归档

- I `01_discontinuous_quad`
  - OpenFOAM目录: [`cases/02_diffusion_equation/01_discontinuous_quad/`](./02_diffusion_equation/01_discontinuous_quad/)
  - 数据: [`data/02_diffusion_equation/cases/01_discontinuous_quad/`](../data/02_diffusion_equation/cases/01_discontinuous_quad/)
  - 图片: [`figures/02_diffusion_equation/cases/01_discontinuous_quad/`](../figures/02_diffusion_equation/cases/01_discontinuous_quad/)
  - 数据汇总: [`data/02_diffusion_equation/analysis/01_discontinuous_quad/analysis.md`](../data/02_diffusion_equation/analysis/01_discontinuous_quad/analysis.md)
  - 报告: [`report/02_diffusion_equation/report.md`](../report/02_diffusion_equation/report.md)
  - 证据索引: [`report/02_diffusion_equation/evidence_index.md`](../report/02_diffusion_equation/evidence_index.md)
- II `02_discontinuous_tri`
  - OpenFOAM目录: [`cases/02_diffusion_equation/02_discontinuous_tri/`](./02_diffusion_equation/02_discontinuous_tri/)
  - 数据: [`data/02_diffusion_equation/cases/02_discontinuous_tri/`](../data/02_diffusion_equation/cases/02_discontinuous_tri/)
  - 图片: [`figures/02_diffusion_equation/cases/02_discontinuous_tri/`](../figures/02_diffusion_equation/cases/02_discontinuous_tri/)
  - 数据汇总: [`data/02_diffusion_equation/analysis/02_discontinuous_tri/analysis.md`](../data/02_diffusion_equation/analysis/02_discontinuous_tri/analysis.md)
  - 报告: [`report/02_diffusion_equation/report.md`](../report/02_diffusion_equation/report.md)
  - 证据索引: [`report/02_diffusion_equation/evidence_index.md`](../report/02_diffusion_equation/evidence_index.md)

### 案例二：Gaussian 扩散

本案例对应二维扩散方程的 Gaussian 初值验证。计算域同样为 `[-5,5] \times [-5,5]`，初始场为光滑 Gaussian 分布，外边界采用与解析解一致的时间变化 Dirichlet 条件，扩散系数仍取 `\mu=1`，终止时间为 `0.2`。实验代号按罗马数字 III-IV 编排，其中 III 对应四边形网格，IV 对应三角形网格；两组都在 `N=10,20,40,80` 下比较光滑场的扩散形态、边界一致性和误差收敛。

该组生成的图片主要包括数值场与解析解的对比图、典型剖面图，以及时间演化图，重点观察峰值衰减、分布扩散是否与解析解一致、边界是否与理论值贴合，以及剖面是否随网格加密逐步逼近精确解。对应的分析数据主要包括各分辨率 `summary.json` 中的 `L1`、`L2` 误差和由相邻网格计算得到的观察收敛阶；其中 `L1` 主要看整体误差，`L2` 用来补充平均误差水平，观察收敛阶则用来判断网格加密后的收敛趋势是否稳定。

#### 相关实验配置

| 实验代号 | 案例 | 网格类型 | 空间格式 | 分辨率 N | 固定设置 |
|---|---|---|---|---|---|
| III | `03_gaussian_quad` | 四边形 | `Gauss linear corrected` | `10, 20, 40, 80` | `solver=explicitDiffusionFoamStudent`，`phi`，`mu=1`，`endTime=0.2`，`diffusionCo=0.45`，`maxDeltaT=0.001`，Gaussian 初值，解析解 Dirichlet 边界 |
| IV | `04_gaussian_tri` | 三角形 | `Gauss linear corrected` | `10, 20, 40, 80` | `solver=explicitDiffusionFoamStudent`，`phi`，`mu=1`，`endTime=0.2`，`diffusionCo=0.45`，`maxDeltaT=0.001`，Gaussian 初值，解析解 Dirichlet 边界 |

#### 相关实验的数据归档

- III `03_gaussian_quad`
  - OpenFOAM目录: [`cases/02_diffusion_equation/03_gaussian_quad/`](./02_diffusion_equation/03_gaussian_quad/)
  - 数据: [`data/02_diffusion_equation/cases/03_gaussian_quad/`](../data/02_diffusion_equation/cases/03_gaussian_quad/)
  - 图片: [`figures/02_diffusion_equation/cases/03_gaussian_quad/`](../figures/02_diffusion_equation/cases/03_gaussian_quad/)
  - 数据汇总: 无独立 `analysis.md`，汇总见 [`report/02_diffusion_equation/report.md`](../report/02_diffusion_equation/report.md) 第 7 节
  - 报告: [`report/02_diffusion_equation/report.md`](../report/02_diffusion_equation/report.md)
  - 证据索引: [`report/02_diffusion_equation/evidence_index.md`](../report/02_diffusion_equation/evidence_index.md)
- IV `04_gaussian_tri`
  - OpenFOAM目录: [`cases/02_diffusion_equation/04_gaussian_tri/`](./02_diffusion_equation/04_gaussian_tri/)
  - 数据: [`data/02_diffusion_equation/cases/04_gaussian_tri/`](../data/02_diffusion_equation/cases/04_gaussian_tri/)
  - 图片: [`figures/02_diffusion_equation/cases/04_gaussian_tri/`](../figures/02_diffusion_equation/cases/04_gaussian_tri/)
  - 数据汇总: 无独立 `analysis.md`，汇总见 [`report/02_diffusion_equation/report.md`](../report/02_diffusion_equation/report.md) 第 7 节
  - 报告: [`report/02_diffusion_equation/report.md`](../report/02_diffusion_equation/report.md)
  - 证据索引: [`report/02_diffusion_equation/evidence_index.md`](../report/02_diffusion_equation/evidence_index.md)

### 3. 对流-扩散方程

### 案例一：正弦波平移

本案例对应二维对流-扩散方程的正弦波平移验证。计算域为单位正方形周期区域，初始场取正弦波，速度场为常量平移，扩散系数取题面给定值，终止时间为 `1.0`。实验代号按罗马数字 I-II 编排，其中 I 对应四边形网格，II 对应三角形网格；两组都在 `N=10,20,40,80` 下比较平移后的波形保持、数值耗散和误差收敛。

该组生成的图片主要包括数值场与精确解的对比图、典型剖面图、振幅变化图，以及误差收敛曲线和观察收敛阶图。前两类图主要看波形平移后是否还能保持峰谷位置和幅值，后两类图主要看误差是否随网格加密而下降，以及收敛阶是否稳定。对应的分析数据主要包括各分辨率的 `L1`、`L2` 误差和由相邻网格计算得到的观察收敛阶；其中 `L1` 更直接反映整体偏差，`L2` 用来补充平均误差水平，观察收敛阶则用来判断网格加密后的精度提升趋势。

#### 相关实验配置

| 实验代号 | 案例 | 网格类型 | 空间格式 | 分辨率 N | 固定设置 |
|---|---|---|---|---|---|
| I | `01_sine_wave_quad_upwind` | 四边形 | `Gauss upwind` | `10, 20, 40, 80` | `solver=explicitAdvectionDiffusionFoamStudent`，`phi`，`U=(1,1,0)`，`mu`，`endTime=1.0`，`advectionDiffusionCo=0.45`，`maxDeltaT=0.001`，周期边界 |
| II | `02_sine_wave_tri_upwind` | 三角形 | `Gauss upwind` | `10, 20, 40, 80` | `solver=explicitAdvectionDiffusionFoamStudent`，`phi`，`U=(1,1,0)`，`mu`，`endTime=1.0`，`advectionDiffusionCo=0.45`，`maxDeltaT=0.001`，周期边界 |

#### 相关实验的数据归档

- I `01_sine_wave_quad_upwind`
  - OpenFOAM目录: [`cases/03_advection_diffusion_equation/01_sine_wave_quad_upwind/`](./03_advection_diffusion_equation/01_sine_wave_quad_upwind/)
  - 数据: [`data/03_advection_diffusion_equation/cases/01_sine_wave_quad_upwind/`](../data/03_advection_diffusion_equation/cases/01_sine_wave_quad_upwind/)
  - 图片: [`figures/03_advection_diffusion_equation/cases/01_sine_wave_quad_upwind/`](../figures/03_advection_diffusion_equation/cases/01_sine_wave_quad_upwind/)
  - 数据汇总: [`data/03_advection_diffusion_equation/analysis/01_sine_wave_quad_upwind/analysis.md`](../data/03_advection_diffusion_equation/analysis/01_sine_wave_quad_upwind/analysis.md)
  - 报告: [`report/03_advection_diffusion_equation/report.md`](../report/03_advection_diffusion_equation/report.md)
  - 证据索引: [`report/03_advection_diffusion_equation/evidence_index.md`](../report/03_advection_diffusion_equation/evidence_index.md)
- II `02_sine_wave_tri_upwind`
  - OpenFOAM目录: [`cases/03_advection_diffusion_equation/02_sine_wave_tri_upwind/`](./03_advection_diffusion_equation/02_sine_wave_tri_upwind/)
  - 数据: [`data/03_advection_diffusion_equation/cases/02_sine_wave_tri_upwind/`](../data/03_advection_diffusion_equation/cases/02_sine_wave_tri_upwind/)
  - 图片: [`figures/03_advection_diffusion_equation/cases/02_sine_wave_tri_upwind/`](../figures/03_advection_diffusion_equation/cases/02_sine_wave_tri_upwind/)
  - 数据汇总: [`data/03_advection_diffusion_equation/analysis/02_sine_wave_tri_upwind/analysis.md`](../data/03_advection_diffusion_equation/analysis/02_sine_wave_tri_upwind/analysis.md)
  - 报告: [`report/03_advection_diffusion_equation/report.md`](../report/03_advection_diffusion_equation/report.md)
  - 证据索引: [`report/03_advection_diffusion_equation/evidence_index.md`](../report/03_advection_diffusion_equation/evidence_index.md)

### 案例二：旋转尖峰平流扩散

本案例对应二维对流-扩散方程的旋转尖峰验证。计算域为 `[-1,1]\times[-1,1]`，初始场为旋转尖峰，速度场为刚体旋转，扩散系数取较小常数，终止时间为 `2\pi`。实验代号按罗马数字 III-VI 编排，其中 III、V 对应四边形网格，IV、VI 对应三角形网格；两类边界设置分别比较零 Dirichlet 近似与解析 Dirichlet，在 `N=20,40,80` 下观察轮廓保持、边界影响和误差收敛。

该组生成的图片主要包括初始场与最终场对比图、最终等值线图、典型剖面图，以及误差收敛曲线和观察收敛阶图。前两类图主要看尖峰轮廓是否被抹平、峰值是否衰减、切口区域是否保持合理，后两类图主要看误差是否随网格加密而下降，以及不同边界处理下的差异是否明显。对应的分析数据主要包括各分辨率的 `L1`、`L2` 误差和由相邻网格计算得到的观察收敛阶；其中 `L1` 更直接反映整体形状偏差，`L2` 用来补充平均误差水平，观察收敛阶则用来判断网格加密后的精度提升趋势。

#### 相关实验配置

| 实验代号 | 案例 | 网格类型 | 空间格式 | 分辨率 N | 固定设置 |
|---|---|---|---|---|---|
| III | `03_rotating_peak_quad_upwind` | 四边形 | `Gauss upwind` | `20, 40, 80` | `solver=explicitAdvectionDiffusionFoamStudent`，`phi`，`U=(-y,x,0)`，`mu=1e-3`，`endTime=2π`，`advectionDiffusionCo=0.45`，`maxDeltaT=0.001`，零 Dirichlet 近似 |
| IV | `04_rotating_peak_tri_upwind` | 三角形 | `Gauss upwind` | `20, 40, 80` | `solver=explicitAdvectionDiffusionFoamStudent`，`phi`，`U=(-y,x,0)`，`mu=1e-3`，`endTime=2π`，`advectionDiffusionCo=0.45`，`maxDeltaT=0.001`，零 Dirichlet 近似 |
| V | `05_rotating_peak_quad_analyticDirichlet_upwind` | 四边形 | `Gauss upwind` | `20, 40, 80` | `solver=explicitAdvectionDiffusionFoamStudent`，`phi`，`U=(-y,x,0)`，`mu=1e-3`，`endTime=2π`，`advectionDiffusionCo=0.45`，`maxDeltaT=0.001`，解析 Dirichlet |
| VI | `06_rotating_peak_tri_analyticDirichlet_upwind` | 三角形 | `Gauss upwind` | `20, 40, 80` | `solver=explicitAdvectionDiffusionFoamStudent`，`phi`，`U=(-y,x,0)`，`mu=1e-3`，`endTime=2π`，`advectionDiffusionCo=0.45`，`maxDeltaT=0.001`，解析 Dirichlet |

#### 相关实验的数据归档

- III `03_rotating_peak_quad_upwind`
  - OpenFOAM目录: [`cases/03_advection_diffusion_equation/03_rotating_peak_quad_upwind/`](./03_advection_diffusion_equation/03_rotating_peak_quad_upwind/)
  - 数据: [`data/03_advection_diffusion_equation/cases/03_rotating_peak_quad_upwind/`](../data/03_advection_diffusion_equation/cases/03_rotating_peak_quad_upwind/)
  - 图片: [`figures/03_advection_diffusion_equation/cases/03_rotating_peak_quad_upwind/`](../figures/03_advection_diffusion_equation/cases/03_rotating_peak_quad_upwind/)
  - 数据汇总: [`data/03_advection_diffusion_equation/analysis/03_rotating_peak_quad_upwind/analysis.md`](../data/03_advection_diffusion_equation/analysis/03_rotating_peak_quad_upwind/analysis.md)
  - 报告: [`report/03_advection_diffusion_equation/report.md`](../report/03_advection_diffusion_equation/report.md)
  - 证据索引: [`report/03_advection_diffusion_equation/evidence_index.md`](../report/03_advection_diffusion_equation/evidence_index.md)
- IV `04_rotating_peak_tri_upwind`
  - OpenFOAM目录: [`cases/03_advection_diffusion_equation/04_rotating_peak_tri_upwind/`](./03_advection_diffusion_equation/04_rotating_peak_tri_upwind/)
  - 数据: [`data/03_advection_diffusion_equation/cases/04_rotating_peak_tri_upwind/`](../data/03_advection_diffusion_equation/cases/04_rotating_peak_tri_upwind/)
  - 图片: [`figures/03_advection_diffusion_equation/cases/04_rotating_peak_tri_upwind/`](../figures/03_advection_diffusion_equation/cases/04_rotating_peak_tri_upwind/)
  - 数据汇总: [`data/03_advection_diffusion_equation/analysis/04_rotating_peak_tri_upwind/analysis.md`](../data/03_advection_diffusion_equation/analysis/04_rotating_peak_tri_upwind/analysis.md)
  - 报告: [`report/03_advection_diffusion_equation/report.md`](../report/03_advection_diffusion_equation/report.md)
  - 证据索引: [`report/03_advection_diffusion_equation/evidence_index.md`](../report/03_advection_diffusion_equation/evidence_index.md)
- V `05_rotating_peak_quad_analyticDirichlet_upwind`
  - OpenFOAM目录: [`cases/03_advection_diffusion_equation/05_rotating_peak_quad_analyticDirichlet_upwind/`](./03_advection_diffusion_equation/05_rotating_peak_quad_analyticDirichlet_upwind/)
  - 数据: [`data/03_advection_diffusion_equation/cases/05_rotating_peak_quad_analyticDirichlet_upwind/`](../data/03_advection_diffusion_equation/cases/05_rotating_peak_quad_analyticDirichlet_upwind/)
  - 图片: [`figures/03_advection_diffusion_equation/cases/05_rotating_peak_quad_analyticDirichlet_upwind/`](../figures/03_advection_diffusion_equation/cases/05_rotating_peak_quad_analyticDirichlet_upwind/)
  - 数据汇总: [`data/03_advection_diffusion_equation/analysis/05_rotating_peak_quad_analyticDirichlet_upwind/analysis.md`](../data/03_advection_diffusion_equation/analysis/05_rotating_peak_quad_analyticDirichlet_upwind/analysis.md)
  - 报告: [`report/03_advection_diffusion_equation/report.md`](../report/03_advection_diffusion_equation/report.md)
  - 证据索引: [`report/03_advection_diffusion_equation/evidence_index.md`](../report/03_advection_diffusion_equation/evidence_index.md)
- VI `06_rotating_peak_tri_analyticDirichlet_upwind`
  - OpenFOAM目录: [`cases/03_advection_diffusion_equation/06_rotating_peak_tri_analyticDirichlet_upwind/`](./03_advection_diffusion_equation/06_rotating_peak_tri_analyticDirichlet_upwind/)
  - 数据: [`data/03_advection_diffusion_equation/cases/06_rotating_peak_tri_analyticDirichlet_upwind/`](../data/03_advection_diffusion_equation/cases/06_rotating_peak_tri_analyticDirichlet_upwind/)
  - 图片: [`figures/03_advection_diffusion_equation/cases/06_rotating_peak_tri_analyticDirichlet_upwind/`](../figures/03_advection_diffusion_equation/cases/06_rotating_peak_tri_analyticDirichlet_upwind/)
  - 数据汇总: [`data/03_advection_diffusion_equation/analysis/06_rotating_peak_tri_analyticDirichlet_upwind/analysis.md`](../data/03_advection_diffusion_equation/analysis/06_rotating_peak_tri_analyticDirichlet_upwind/analysis.md)
  - 报告: [`report/03_advection_diffusion_equation/report.md`](../report/03_advection_diffusion_equation/report.md)
  - 证据索引: [`report/03_advection_diffusion_equation/evidence_index.md`](../report/03_advection_diffusion_equation/evidence_index.md)

### 4. Poisson 方程

### 案例一：制造解 Poisson 方程

本案例对应二维 Poisson 方程的制造解验证。计算域为单位正方形，解析解与源项按题面给定，四条边均采用制造解 Dirichlet 边界，目标是检查离散一致性和网格收敛。实验代号按罗马数字 I-II 编排，其中 I 对应四边形网格，II 对应三角形网格；两组都在 `N=10,20,40,80` 下比较误差下降和收敛阶。

该组生成的图片主要包括数值解与解析解对比图、误差分布图、全分辨率对比图，以及观察收敛阶图。前两类图主要看数值解是否贴合解析解、误差是否均匀分布，后两类图主要看误差是否随网格加密而下降，以及收敛阶是否接近理论值。对应的分析数据主要包括各分辨率的 `normalizedL1`、`normalizedL2` 误差和由相邻网格计算得到的观察收敛阶；其中 `normalizedL1` 更直接反映整体偏差，`normalizedL2` 用来补充平均误差水平，观察收敛阶则用来判断网格加密后的精度提升趋势。

#### 相关实验配置

| 实验代号 | 案例 | 网格类型 | 空间格式 | 分辨率 N | 固定设置 |
|---|---|---|---|---|---|
| I | `01_poisson_manufactured_quad` | 四边形 | `Gauss linear corrected` | `10, 20, 40, 80` | `solver=poissonFoamStudent`，`phi`，`omega`，制造解 Dirichlet 边界，稳态求解 |
| II | `02_poisson_manufactured_tri` | 三角形 | `Gauss linear corrected` | `10, 20, 40, 80` | `solver=poissonFoamStudent`，`phi`，`omega`，制造解 Dirichlet 边界，稳态求解 |

#### 相关实验的数据归档

- I `01_poisson_manufactured_quad`
  - OpenFOAM目录: [`cases/04_poisson_equation/01_poisson_manufactured_quad/`](./04_poisson_equation/01_poisson_manufactured_quad/)
  - 数据: [`data/04_poisson_equation/cases/01_poisson_manufactured_quad/`](../data/04_poisson_equation/cases/01_poisson_manufactured_quad/)
  - 图片: [`figures/04_poisson_equation/cases/01_poisson_manufactured_quad/`](../figures/04_poisson_equation/cases/01_poisson_manufactured_quad/)
  - 数据汇总: [`data/04_poisson_equation/analysis/01_poisson_manufactured_quad/analysis.md`](../data/04_poisson_equation/analysis/01_poisson_manufactured_quad/analysis.md)
  - 报告: [`report/04_poisson_equation/report.md`](../report/04_poisson_equation/report.md)
  - 证据索引: [`report/04_poisson_equation/evidence_index.md`](../report/04_poisson_equation/evidence_index.md)
- II `02_poisson_manufactured_tri`
  - OpenFOAM目录: [`cases/04_poisson_equation/02_poisson_manufactured_tri/`](./04_poisson_equation/02_poisson_manufactured_tri/)
  - 数据: [`data/04_poisson_equation/cases/02_poisson_manufactured_tri/`](../data/04_poisson_equation/cases/02_poisson_manufactured_tri/)
  - 图片: [`figures/04_poisson_equation/cases/02_poisson_manufactured_tri/`](../figures/04_poisson_equation/cases/02_poisson_manufactured_tri/)
  - 数据汇总: [`data/04_poisson_equation/analysis/02_poisson_manufactured_tri/analysis.md`](../data/04_poisson_equation/analysis/02_poisson_manufactured_tri/analysis.md)
  - 报告: [`report/04_poisson_equation/report.md`](../report/04_poisson_equation/report.md)
  - 证据索引: [`report/04_poisson_equation/evidence_index.md`](../report/04_poisson_equation/evidence_index.md)

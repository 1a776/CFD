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

本案例对应二维稳态 Poisson 方程的制造解验证。计算域为单位正方形
`[0,1] × [0,1]`，解析解为
$\phi_{\mathrm{exact}}=\cos(\pi x)\cos(\pi y)$，源项为
$\omega=-2\pi^2\cos(\pi x)\cos(\pi y)$，四条边均采用与解析解一致的
Dirichlet 边界条件。两组实验使用相同的控制方程、源项、边界条件、
拉普拉斯离散格式和线性求解设置，仅改变网格拓扑和网格生成后端：
I 采用结构化四边形网格，II 采用三角形棱柱网格。每组均在
`N=10,20,40,80` 下进行网格加密，用于分析离散误差和网格收敛行为。

该案例生成的图片主要包括各分辨率下数值解、解析解与误差场的对比图，
不同分辨率的场分布汇总图，以及 `normalizedL1`、`normalizedL2` 和
`normalizedLinf` 随网格加密变化的误差曲线和观察收敛阶图。对应的分析数据
主要包括每个分辨率的单元数、三种误差范数、相邻网格之间计算得到的
`L1`、`L2` 和 `Linf` 观察收敛阶、最终场幅值、归一化质量误差、网格检查状态
和求解结束状态。其中，`normalizedL1` 用于衡量整体场误差，`normalizedL2`
用于补充均方意义下的误差水平，`normalizedLinf` 用于检查局部最大偏差；
观察收敛阶用于判断误差是否随网格加密按预期下降，并比较四边形和三角形
网格在相同名义分辨率下的收敛表现。

#### 相关实验配置

| 实验代号 | 案例 | 网格类型 | 网格后端 | 分辨率 N | 固定设置 |
|---|---|---|---|---|---|
| I | `01_poisson_manufactured_quad` | 四边形 | `blockMesh` | `10, 20, 40, 80` | 单位正方形区域、$\phi_{\mathrm{exact}}=\cos(\pi x)\cos(\pi y)$、$\omega=-2\pi^2\cos(\pi x)\cos(\pi y)$、四边 Dirichlet、`Gauss linear corrected`、`GAMG`、线性求解容差 `1e-12`、非正交修正 `2` 次、稳态计算 |
| II | `02_poisson_manufactured_tri` | 三角形棱柱 | `gmsh` | `10, 20, 40, 80` | 单位正方形区域、$\phi_{\mathrm{exact}}=\cos(\pi x)\cos(\pi y)$、$\omega=-2\pi^2\cos(\pi x)\cos(\pi y)$、四边 Dirichlet、`Gauss linear corrected`、`GAMG`、线性求解容差 `1e-12`、非正交修正 `2` 次、稳态计算 |

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

### 5. Navier-Stokes 方程

### 案例一：方腔顶盖驱动流

本案例对应二维不可压 Navier-Stokes 方程的方腔顶盖驱动流。计算域为
单位正方形 `[0,1] × [0,1]`，左、右、下壁面采用无滑移条件，上壁面以
`U=(1,0,0)` 匀速运动，运动黏度按 $\nu=1/Re$ 设置。实验比较压力投影法
和 PISO 法在不同 Reynolds 数及混合网格分辨率下的计算结果，当前归档结果
覆盖 `Re=1000`、`Re=3200`，以及 `hybrid40` 和 `hybrid80` 两个网格等级。
混合网格在壁面附近布置结构化网格层，在内部保留非结构化区域，用于同时
分辨壁面速度梯度和腔体内部主循环。

该案例生成的图片主要包括速度场与流线图，以及方腔两条中心线
`u(0.5,y)`、`v(x,0.5)` 与 Ghia 参考数据的对比图。速度场和流线图用于
观察顶盖剪切层、主涡和整体循环结构；中心线图用于检查不同 Reynolds 数、
网格等级和压力-速度耦合算法对速度分布的影响。对应的分析数据主要包括
`u_centerline.csv`、`v_centerline.csv`、单元数、实际结束时间、速度分量的
最大值和最小值、相对于 Ghia 数据的 RMSE 和最大绝对误差，以及网格检查和
稳态求解日志。该案例没有解析解，因此不使用 Poisson 案例中的 `L1`、`L2`
收敛阶，而是以中心线 RMSE、最大绝对误差、主循环形态和稳态判据作为主要
比较依据。

#### 相关实验配置

| 实验代号 | 案例 | 算法 | 网格 | Reynolds 数 | 分辨率 / 单元数 | 固定设置 |
|---|---|---|---|---:|---|---|
| I | `07_lid_driven_cavity_projection_Re1000_hybrid40` | 压力投影法 | 混合网格 `hybrid40` | `1000` | `N=40` / `2646` | 单位正方形、顶盖 `U=(1,0,0)`、其余壁面无滑移、$\nu=0.001$、`Δt=0.001`、`endTime=100`、`maxCo=0.2`、对流项 `Gauss linearUpwind grad(U)`、扩散项 `Gauss linear corrected`、`GAMG`、线性容差 `1e-8`、压力校正 `2` 次、非正交修正 `1` 次、稳态速度容差 `1e-6`、质量容差 `1e-8` |
| II | `08_lid_driven_cavity_projection_Re1000_hybrid80` | 压力投影法 | 混合网格 `hybrid80` | `1000` | `N=80` / `9282` | 单位正方形、顶盖 `U=(1,0,0)`、其余壁面无滑移、$\nu=0.001$、`Δt=0.0005`、`endTime=100`、`maxCo=0.2`、对流项 `Gauss linearUpwind grad(U)`、扩散项 `Gauss linear corrected`、`GAMG`、线性容差 `1e-8`、压力校正 `2` 次、非正交修正 `1` 次、稳态速度容差 `1e-6`、质量容差 `1e-8` |
| III | `11_lid_driven_cavity_projection_Re3200_hybrid80` | 压力投影法 | 混合网格 `hybrid80` | `3200` | `N=80` / `9282` | 单位正方形、顶盖 `U=(1,0,0)`、其余壁面无滑移、$\nu=0.0003125$、`Δt=0.0005`、`endTime=150`、`maxCo=0.2`、对流项 `Gauss linearUpwind grad(U)`、扩散项 `Gauss linear corrected`、`GAMG`、线性容差 `1e-8`、压力校正 `2` 次、非正交修正 `1` 次、稳态速度容差 `1e-6`、质量容差 `1e-8` |
| IV | `16_lid_driven_cavity_piso_Re1000_hybrid40` | PISO 法 | 混合网格 `hybrid40` | `1000` | `N=40` / `2646` | 单位正方形、顶盖 `U=(1,0,0)`、其余壁面无滑移、$\nu=0.001$、`Δt=0.001`、`endTime=100`、`maxCo=0.2`、`maxDeltaT=0.01`、对流项 `Gauss linearUpwind grad(U)`、扩散项 `Gauss linear corrected`、`GAMG`、线性容差 `1e-8`、压力校正 `2` 次、非正交修正 `1` 次、稳态速度容差 `1e-6`、质量容差 `1e-8` |
| V | `17_lid_driven_cavity_piso_Re1000_hybrid80` | PISO 法 | 混合网格 `hybrid80` | `1000` | `N=80` / `9282` | 单位正方形、顶盖 `U=(1,0,0)`、其余壁面无滑移、$\nu=0.001$、`Δt=0.0005`、`endTime=100`、`maxCo=0.2`、`maxDeltaT=0.01`、对流项 `Gauss linearUpwind grad(U)`、扩散项 `Gauss linear corrected`、`GAMG`、线性容差 `1e-8`、压力校正 `2` 次、非正交修正 `1` 次、稳态速度容差 `1e-6`、质量容差 `1e-8` |
| VI | `18_lid_driven_cavity_piso_Re3200_hybrid40` | PISO 法 | 混合网格 `hybrid40` | `3200` | `N=40` / `2646` | 单位正方形、顶盖 `U=(1,0,0)`、其余壁面无滑移、$\nu=0.0003125$、`Δt=0.001`、`endTime=150`、`maxCo=0.2`、`maxDeltaT=0.01`、对流项 `Gauss linearUpwind grad(U)`、扩散项 `Gauss linear corrected`、`GAMG`、线性容差 `1e-8`、压力校正 `2` 次、非正交修正 `1` 次、稳态速度容差 `1e-6`、质量容差 `1e-8` |
| VII | `19_lid_driven_cavity_piso_Re3200_hybrid80` | PISO 法 | 混合网格 `hybrid80` | `3200` | `N=80` / `9282` | 单位正方形、顶盖 `U=(1,0,0)`、其余壁面无滑移、$\nu=0.0003125$、`Δt=0.0005`、`endTime=150`、`maxCo=0.2`、`maxDeltaT=0.01`、对流项 `Gauss linearUpwind grad(U)`、扩散项 `Gauss linear corrected`、`GAMG`、线性容差 `1e-8`、压力校正 `2` 次、非正交修正 `1` 次、稳态速度容差 `1e-6`、质量容差 `1e-8` |

#### 相关实验的数据归档

- I `07_lid_driven_cavity_projection_Re1000_hybrid40`
  - OpenFOAM目录: [`cases/05_navier_stokes_equation/07_lid_driven_cavity_projection_Re1000_hybrid40/`](./05_navier_stokes_equation/07_lid_driven_cavity_projection_Re1000_hybrid40/)
  - 数据: [`data/05_navier_stokes_equation/cases/07_lid_driven_cavity_projection_Re1000_hybrid40/`](../data/05_navier_stokes_equation/cases/07_lid_driven_cavity_projection_Re1000_hybrid40/)
  - 图片: [`figures/05_navier_stokes_equation/cases/07_lid_driven_cavity_projection_Re1000_hybrid40/`](../figures/05_navier_stokes_equation/cases/07_lid_driven_cavity_projection_Re1000_hybrid40/)
  - 数据汇总: [`data/05_navier_stokes_equation/cases/07_lid_driven_cavity_projection_Re1000_hybrid40/summary.json`](../data/05_navier_stokes_equation/cases/07_lid_driven_cavity_projection_Re1000_hybrid40/summary.json)
  - 报告: [`report/05_navier_stokes_equation/projection/report.md`](../report/05_navier_stokes_equation/projection/report.md)
  - 证据索引: [`report/05_navier_stokes_equation/projection/evidence_index.md`](../report/05_navier_stokes_equation/projection/evidence_index.md)
- II `08_lid_driven_cavity_projection_Re1000_hybrid80`
  - OpenFOAM目录: [`cases/05_navier_stokes_equation/08_lid_driven_cavity_projection_Re1000_hybrid80/`](./05_navier_stokes_equation/08_lid_driven_cavity_projection_Re1000_hybrid80/)
  - 数据: [`data/05_navier_stokes_equation/cases/08_lid_driven_cavity_projection_Re1000_hybrid80/`](../data/05_navier_stokes_equation/cases/08_lid_driven_cavity_projection_Re1000_hybrid80/)
  - 图片: [`figures/05_navier_stokes_equation/cases/08_lid_driven_cavity_projection_Re1000_hybrid80/`](../figures/05_navier_stokes_equation/cases/08_lid_driven_cavity_projection_Re1000_hybrid80/)
  - 数据汇总: [`data/05_navier_stokes_equation/cases/08_lid_driven_cavity_projection_Re1000_hybrid80/summary.json`](../data/05_navier_stokes_equation/cases/08_lid_driven_cavity_projection_Re1000_hybrid80/summary.json)
  - 报告: [`report/05_navier_stokes_equation/projection/report.md`](../report/05_navier_stokes_equation/projection/report.md)
  - 证据索引: [`report/05_navier_stokes_equation/projection/evidence_index.md`](../report/05_navier_stokes_equation/projection/evidence_index.md)
- III `11_lid_driven_cavity_projection_Re3200_hybrid80`
  - OpenFOAM目录: [`cases/05_navier_stokes_equation/11_lid_driven_cavity_projection_Re3200_hybrid80/`](./05_navier_stokes_equation/11_lid_driven_cavity_projection_Re3200_hybrid80/)
  - 数据: [`data/05_navier_stokes_equation/cases/11_lid_driven_cavity_projection_Re3200_hybrid80/`](../data/05_navier_stokes_equation/cases/11_lid_driven_cavity_projection_Re3200_hybrid80/)
  - 图片: [`figures/05_navier_stokes_equation/cases/11_lid_driven_cavity_projection_Re3200_hybrid80/`](../figures/05_navier_stokes_equation/cases/11_lid_driven_cavity_projection_Re3200_hybrid80/)
  - 数据汇总: [`data/05_navier_stokes_equation/cases/11_lid_driven_cavity_projection_Re3200_hybrid80/summary.json`](../data/05_navier_stokes_equation/cases/11_lid_driven_cavity_projection_Re3200_hybrid80/summary.json)
  - 报告: [`report/05_navier_stokes_equation/projection/report.md`](../report/05_navier_stokes_equation/projection/report.md)
  - 证据索引: [`report/05_navier_stokes_equation/projection/evidence_index.md`](../report/05_navier_stokes_equation/projection/evidence_index.md)
- IV `16_lid_driven_cavity_piso_Re1000_hybrid40`
  - OpenFOAM目录: [`cases/05_navier_stokes_equation/16_lid_driven_cavity_piso_Re1000_hybrid40/`](./05_navier_stokes_equation/16_lid_driven_cavity_piso_Re1000_hybrid40/)
  - 数据: [`data/05_navier_stokes_equation/cases/16_lid_driven_cavity_piso_Re1000_hybrid40/`](../data/05_navier_stokes_equation/cases/16_lid_driven_cavity_piso_Re1000_hybrid40/)
  - 图片: [`figures/05_navier_stokes_equation/cases/16_lid_driven_cavity_piso_Re1000_hybrid40/`](../figures/05_navier_stokes_equation/cases/16_lid_driven_cavity_piso_Re1000_hybrid40/)
  - 数据汇总: [`data/05_navier_stokes_equation/cases/16_lid_driven_cavity_piso_Re1000_hybrid40/summary.json`](../data/05_navier_stokes_equation/cases/16_lid_driven_cavity_piso_Re1000_hybrid40/summary.json)
  - 报告: [`report/05_navier_stokes_equation/piso/report.md`](../report/05_navier_stokes_equation/piso/report.md)
  - 证据索引: [`report/05_navier_stokes_equation/piso/evidence_index.md`](../report/05_navier_stokes_equation/piso/evidence_index.md)
- V `17_lid_driven_cavity_piso_Re1000_hybrid80`
  - OpenFOAM目录: [`cases/05_navier_stokes_equation/17_lid_driven_cavity_piso_Re1000_hybrid80/`](./05_navier_stokes_equation/17_lid_driven_cavity_piso_Re1000_hybrid80/)
  - 数据: [`data/05_navier_stokes_equation/cases/17_lid_driven_cavity_piso_Re1000_hybrid80/`](../data/05_navier_stokes_equation/cases/17_lid_driven_cavity_piso_Re1000_hybrid80/)
  - 图片: [`figures/05_navier_stokes_equation/cases/17_lid_driven_cavity_piso_Re1000_hybrid80/`](../figures/05_navier_stokes_equation/cases/17_lid_driven_cavity_piso_Re1000_hybrid80/)
  - 数据汇总: [`data/05_navier_stokes_equation/cases/17_lid_driven_cavity_piso_Re1000_hybrid80/summary.json`](../data/05_navier_stokes_equation/cases/17_lid_driven_cavity_piso_Re1000_hybrid80/summary.json)
  - 报告: [`report/05_navier_stokes_equation/piso/report.md`](../report/05_navier_stokes_equation/piso/report.md)
  - 证据索引: [`report/05_navier_stokes_equation/piso/evidence_index.md`](../report/05_navier_stokes_equation/piso/evidence_index.md)
- VI `18_lid_driven_cavity_piso_Re3200_hybrid40`
  - OpenFOAM目录: [`cases/05_navier_stokes_equation/18_lid_driven_cavity_piso_Re3200_hybrid40/N40/`](./05_navier_stokes_equation/18_lid_driven_cavity_piso_Re3200_hybrid40/N40/)
  - 数据: [`data/05_navier_stokes_equation/cases/18_lid_driven_cavity_piso_Re3200_hybrid40/`](../data/05_navier_stokes_equation/cases/18_lid_driven_cavity_piso_Re3200_hybrid40/)
  - 图片: [`figures/05_navier_stokes_equation/cases/18_lid_driven_cavity_piso_Re3200_hybrid40/`](../figures/05_navier_stokes_equation/cases/18_lid_driven_cavity_piso_Re3200_hybrid40/)
  - 数据汇总: [`data/05_navier_stokes_equation/cases/18_lid_driven_cavity_piso_Re3200_hybrid40/summary.json`](../data/05_navier_stokes_equation/cases/18_lid_driven_cavity_piso_Re3200_hybrid40/summary.json)
  - 报告: [`report/05_navier_stokes_equation/piso/report.md`](../report/05_navier_stokes_equation/piso/report.md)
  - 证据索引: [`report/05_navier_stokes_equation/piso/evidence_index.md`](../report/05_navier_stokes_equation/piso/evidence_index.md)
- VII `19_lid_driven_cavity_piso_Re3200_hybrid80`
  - OpenFOAM目录: [`cases/05_navier_stokes_equation/19_lid_driven_cavity_piso_Re3200_hybrid80/N80/`](./05_navier_stokes_equation/19_lid_driven_cavity_piso_Re3200_hybrid80/N80/)
  - 数据: [`data/05_navier_stokes_equation/cases/19_lid_driven_cavity_piso_Re3200_hybrid80/`](../data/05_navier_stokes_equation/cases/19_lid_driven_cavity_piso_Re3200_hybrid80/)
  - 图片: [`figures/05_navier_stokes_equation/cases/19_lid_driven_cavity_piso_Re3200_hybrid80/`](../figures/05_navier_stokes_equation/cases/19_lid_driven_cavity_piso_Re3200_hybrid80/)
  - 数据汇总: [`data/05_navier_stokes_equation/cases/19_lid_driven_cavity_piso_Re3200_hybrid80/summary.json`](../data/05_navier_stokes_equation/cases/19_lid_driven_cavity_piso_Re3200_hybrid80/summary.json)
  - 报告: [`report/05_navier_stokes_equation/piso/report.md`](../report/05_navier_stokes_equation/piso/report.md)
  - 证据索引: [`report/05_navier_stokes_equation/piso/evidence_index.md`](../report/05_navier_stokes_equation/piso/evidence_index.md)

### 案例二：等边三角腔顶盖驱动流

本案例对应二维等边三角腔顶盖驱动流。三角腔顶点取
$A=(-\sqrt{3},0)$、$B=(\sqrt{3},0)$、$C=(0,-3)$，上边界以
`U=(1,0,0)` 匀速运动，左右边界采用无滑移条件，特征长度取 `a=1`，
Reynolds 数设置为 `100`、`200` 和 `500`。所有当前归档实验采用
`hybrid80` 混合网格，壁面附近设置结构化层，内部为非结构化三角形区域，
再沿厚度方向挤出为三维棱柱网格。压力投影法和 PISO 法分别覆盖三个
Reynolds 数，用于比较压力-速度耦合算法对三角腔主循环、角区流动和涡结构
的影响。

该案例生成的图片主要包括速度场与流线图、中心线和水平剖面与参考数据的
对比图，以及流函数、涡量和主涡结构图。速度场和流线图用于观察顶盖剪切层
及腔体主循环，中心线和水平剖面用于比较速度分布，流函数和涡量图用于定位
主涡及可能的局部涡结构。对应的分析数据主要包括 `u_centerline.csv`、
`v_horizontal.csv`、主涡坐标 `(x,y)`、主涡流函数值 `psi`、主涡涡量
`omega`、与文献参考主涡的 `|dx|`、`|dy|` 和 `|dpsi|`、速度分量范围、
实际结束时间和网格单元数。该案例不存在统一的解析解，因此主要采用主涡
位置、流函数、中心线剖面和流线形态进行验证，不使用 `L1`、`L2` 观察收敛阶。
其中，压力投影法的 Re=100 结果只到约 `t=0.6`，不能作为稳态结论；压力
投影法 Re=200、500 的当前归档也缺少完整求解日志和最终正时间场，只能作为
摘要级结果使用，PISO 三个 Reynolds 数的稳态日志和结果较完整。

#### 相关实验配置

| 实验代号 | 案例 | 算法 | 网格 | Reynolds 数 | 分辨率 / 单元数 | 固定设置 |
|---|---|---|---|---:|---|---|
| VIII | `26_triangular_cavity_projection_Re100_hybrid80` | 压力投影法 | 混合网格 `hybrid80` | `100` | `N=80` / `7308` | 等边三角腔、顶盖 `U=(1,0,0)`、左右壁面无滑移、$\nu=0.01$、`Δt=0.001`、`endTime=80`、`maxCo=0.2`、对流项 `Gauss linearUpwind grad(U)`、扩散项 `Gauss linear corrected`、`GAMG`、线性容差 `1e-8`、压力校正 `2` 次、非正交修正 `1` 次、稳态速度容差 `1e-6`、质量容差 `1e-8` |
| IX | `27_triangular_cavity_projection_Re200_hybrid80` | 压力投影法 | 混合网格 `hybrid80` | `200` | `N=80` / `7308` | 等边三角腔、顶盖 `U=(1,0,0)`、左右壁面无滑移、$\nu=0.005$、`Δt=0.0005`、`endTime=100`、`maxCo=0.2`、对流项 `Gauss linearUpwind grad(U)`、扩散项 `Gauss linear corrected`、`GAMG`、线性容差 `1e-8`、压力校正 `2` 次、非正交修正 `1` 次、稳态速度容差 `1e-6`、质量容差 `1e-8` |
| X | `28_triangular_cavity_projection_Re500_hybrid80` | 压力投影法 | 混合网格 `hybrid80` | `500` | `N=80` / `7308` | 等边三角腔、顶盖 `U=(1,0,0)`、左右壁面无滑移、$\nu=0.002$、`Δt=0.00025`、`endTime=140`、`maxCo=0.2`、对流项 `Gauss linearUpwind grad(U)`、扩散项 `Gauss linear corrected`、`GAMG`、线性容差 `1e-8`、压力校正 `2` 次、非正交修正 `1` 次、稳态速度容差 `1e-6`、质量容差 `1e-8` |
| XI | `29_triangular_cavity_piso_Re100_hybrid80` | PISO 法 | 混合网格 `hybrid80` | `100` | `N=80` / `7308` | 等边三角腔、顶盖 `U=(1,0,0)`、左右壁面无滑移、$\nu=0.01$、`Δt=0.001`、`endTime=80`、`maxCo=0.2`、`maxDeltaT=0.01`、对流项 `Gauss linearUpwind grad(U)`、扩散项 `Gauss linear corrected`、`GAMG`、线性容差 `1e-8`、压力校正 `2` 次、非正交修正 `1` 次、稳态速度容差 `1e-6`、质量容差 `1e-8` |
| XII | `30_triangular_cavity_piso_Re200_hybrid80` | PISO 法 | 混合网格 `hybrid80` | `200` | `N=80` / `7308` | 等边三角腔、顶盖 `U=(1,0,0)`、左右壁面无滑移、$\nu=0.005$、`Δt=0.0005`、`endTime=100`、`maxCo=0.2`、`maxDeltaT=0.01`、对流项 `Gauss linearUpwind grad(U)`、扩散项 `Gauss linear corrected`、`GAMG`、线性容差 `1e-8`、压力校正 `2` 次、非正交修正 `1` 次、稳态速度容差 `1e-6`、质量容差 `1e-8` |
| XIII | `31_triangular_cavity_piso_Re500_hybrid80` | PISO 法 | 混合网格 `hybrid80` | `500` | `N=80` / `7308` | 等边三角腔、顶盖 `U=(1,0,0)`、左右壁面无滑移、$\nu=0.002$、`Δt=0.00025`、`endTime=140`、`maxCo=0.2`、`maxDeltaT=0.01`、对流项 `Gauss linearUpwind grad(U)`、扩散项 `Gauss linear corrected`、`GAMG`、线性容差 `1e-8`、压力校正 `2` 次、非正交修正 `1` 次、稳态速度容差 `1e-6`、质量容差 `1e-8` |

#### 相关实验的数据归档

- VIII `26_triangular_cavity_projection_Re100_hybrid80`
  - OpenFOAM目录: 当前未保留该案例的 OpenFOAM 目录，仅保留结果摘要和图片
  - 数据: [`data/05_navier_stokes_equation/cases/26_triangular_cavity_projection_Re100_hybrid80/`](../data/05_navier_stokes_equation/cases/26_triangular_cavity_projection_Re100_hybrid80/)
  - 图片: [`figures/05_navier_stokes_equation/cases/26_triangular_cavity_projection_Re100_hybrid80/`](../figures/05_navier_stokes_equation/cases/26_triangular_cavity_projection_Re100_hybrid80/)
  - 数据汇总: [`data/05_navier_stokes_equation/cases/26_triangular_cavity_projection_Re100_hybrid80/summary.json`](../data/05_navier_stokes_equation/cases/26_triangular_cavity_projection_Re100_hybrid80/summary.json)
  - 报告: [`report/05_navier_stokes_equation/projection/report.md`](../report/05_navier_stokes_equation/projection/report.md)
  - 证据索引: [`report/05_navier_stokes_equation/projection/evidence_index.md`](../report/05_navier_stokes_equation/projection/evidence_index.md)
- IX `27_triangular_cavity_projection_Re200_hybrid80`
  - OpenFOAM目录: [`cases/05_navier_stokes_equation/27_triangular_cavity_projection_Re200_hybrid80/N80/`](./05_navier_stokes_equation/27_triangular_cavity_projection_Re200_hybrid80/N80/)
  - 数据: [`data/05_navier_stokes_equation/cases/27_triangular_cavity_projection_Re200_hybrid80/`](../data/05_navier_stokes_equation/cases/27_triangular_cavity_projection_Re200_hybrid80/)
  - 图片: [`figures/05_navier_stokes_equation/cases/27_triangular_cavity_projection_Re200_hybrid80/`](../figures/05_navier_stokes_equation/cases/27_triangular_cavity_projection_Re200_hybrid80/)
  - 数据汇总: [`data/05_navier_stokes_equation/cases/27_triangular_cavity_projection_Re200_hybrid80/summary.json`](../data/05_navier_stokes_equation/cases/27_triangular_cavity_projection_Re200_hybrid80/summary.json)
  - 报告: [`report/05_navier_stokes_equation/projection/report.md`](../report/05_navier_stokes_equation/projection/report.md)
  - 证据索引: [`report/05_navier_stokes_equation/projection/evidence_index.md`](../report/05_navier_stokes_equation/projection/evidence_index.md)
- X `28_triangular_cavity_projection_Re500_hybrid80`
  - OpenFOAM目录: [`cases/05_navier_stokes_equation/28_triangular_cavity_projection_Re500_hybrid80/N80/`](./05_navier_stokes_equation/28_triangular_cavity_projection_Re500_hybrid80/N80/)
  - 数据: [`data/05_navier_stokes_equation/cases/28_triangular_cavity_projection_Re500_hybrid80/`](../data/05_navier_stokes_equation/cases/28_triangular_cavity_projection_Re500_hybrid80/)
  - 图片: [`figures/05_navier_stokes_equation/cases/28_triangular_cavity_projection_Re500_hybrid80/`](../figures/05_navier_stokes_equation/cases/28_triangular_cavity_projection_Re500_hybrid80/)
  - 数据汇总: [`data/05_navier_stokes_equation/cases/28_triangular_cavity_projection_Re500_hybrid80/summary.json`](../data/05_navier_stokes_equation/cases/28_triangular_cavity_projection_Re500_hybrid80/summary.json)
  - 报告: [`report/05_navier_stokes_equation/projection/report.md`](../report/05_navier_stokes_equation/projection/report.md)
  - 证据索引: [`report/05_navier_stokes_equation/projection/evidence_index.md`](../report/05_navier_stokes_equation/projection/evidence_index.md)
- XI `29_triangular_cavity_piso_Re100_hybrid80`
  - OpenFOAM目录: [`cases/05_navier_stokes_equation/29_triangular_cavity_piso_Re100_hybrid80/N80/`](./05_navier_stokes_equation/29_triangular_cavity_piso_Re100_hybrid80/N80/)
  - 数据: [`data/05_navier_stokes_equation/cases/29_triangular_cavity_piso_Re100_hybrid80/`](../data/05_navier_stokes_equation/cases/29_triangular_cavity_piso_Re100_hybrid80/)
  - 图片: [`figures/05_navier_stokes_equation/cases/29_triangular_cavity_piso_Re100_hybrid80/`](../figures/05_navier_stokes_equation/cases/29_triangular_cavity_piso_Re100_hybrid80/)
  - 数据汇总: [`data/05_navier_stokes_equation/cases/29_triangular_cavity_piso_Re100_hybrid80/summary.json`](../data/05_navier_stokes_equation/cases/29_triangular_cavity_piso_Re100_hybrid80/summary.json)
  - 报告: [`report/05_navier_stokes_equation/piso/report.md`](../report/05_navier_stokes_equation/piso/report.md)
  - 证据索引: [`report/05_navier_stokes_equation/piso/evidence_index.md`](../report/05_navier_stokes_equation/piso/evidence_index.md)
- XII `30_triangular_cavity_piso_Re200_hybrid80`
  - OpenFOAM目录: [`cases/05_navier_stokes_equation/30_triangular_cavity_piso_Re200_hybrid80/N80/`](./05_navier_stokes_equation/30_triangular_cavity_piso_Re200_hybrid80/N80/)
  - 数据: [`data/05_navier_stokes_equation/cases/30_triangular_cavity_piso_Re200_hybrid80/`](../data/05_navier_stokes_equation/cases/30_triangular_cavity_piso_Re200_hybrid80/)
  - 图片: [`figures/05_navier_stokes_equation/cases/30_triangular_cavity_piso_Re200_hybrid80/`](../figures/05_navier_stokes_equation/cases/30_triangular_cavity_piso_Re200_hybrid80/)
  - 数据汇总: [`data/05_navier_stokes_equation/cases/30_triangular_cavity_piso_Re200_hybrid80/summary.json`](../data/05_navier_stokes_equation/cases/30_triangular_cavity_piso_Re200_hybrid80/summary.json)
  - 报告: [`report/05_navier_stokes_equation/piso/report.md`](../report/05_navier_stokes_equation/piso/report.md)
  - 证据索引: [`report/05_navier_stokes_equation/piso/evidence_index.md`](../report/05_navier_stokes_equation/piso/evidence_index.md)
- XIII `31_triangular_cavity_piso_Re500_hybrid80`
  - OpenFOAM目录: [`cases/05_navier_stokes_equation/31_triangular_cavity_piso_Re500_hybrid80/N80/`](./05_navier_stokes_equation/31_triangular_cavity_piso_Re500_hybrid80/N80/)
  - 数据: [`data/05_navier_stokes_equation/cases/31_triangular_cavity_piso_Re500_hybrid80/`](../data/05_navier_stokes_equation/cases/31_triangular_cavity_piso_Re500_hybrid80/)
  - 图片: [`figures/05_navier_stokes_equation/cases/31_triangular_cavity_piso_Re500_hybrid80/`](../figures/05_navier_stokes_equation/cases/31_triangular_cavity_piso_Re500_hybrid80/)
  - 数据汇总: [`data/05_navier_stokes_equation/cases/31_triangular_cavity_piso_Re500_hybrid80/summary.json`](../data/05_navier_stokes_equation/cases/31_triangular_cavity_piso_Re500_hybrid80/summary.json)
  - 报告: [`report/05_navier_stokes_equation/piso/report.md`](../report/05_navier_stokes_equation/piso/report.md)
  - 证据索引: [`report/05_navier_stokes_equation/piso/evidence_index.md`](../report/05_navier_stokes_equation/piso/evidence_index.md)

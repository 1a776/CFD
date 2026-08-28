## `cases/` 案例目录说明
- `cases/` 是 OpenFOAM 原始案例目录；真正用于报告分析的数据和图片放在`data/`、`figures/` 和 `report/` 中。

## 具体题目的相关设计的的实验案例以及实验配置

### 1.对流方程
  ### 案例一：正弦波平移
  ####  相关的实验设计及其配置
  - 在三角形网格和四边形网格上分别采用一阶迎风和线性插值配置N=10,20,40,80,Conmax=0.2,并分析L1 L2 误差及收敛阶，图像包括...
  - 相关实验配置具体说明
    （用表格说明）
  - 相关实验的数据
    用表格说明 OpenFOAM的cases目录，后处理数据，表格，report报告，证据索引




  
  - 第一题总报告和证据索引
  - 第一题报告系统分析二维线性对流方程中正弦波平移和复杂轮廓固体旋转两个算例在不同网格与格式下的数值表现。
  - 汇总报告：[`report/01_advection_equation/report.md`](../report/01_advection_equation/report.md)，分析各案例的题目定义、网格设置、格式选择、运行结果和主要结论。
  - 证据索引：[`report/01_advection_equation/evidence_index.md`](../report/01_advection_equation/evidence_index.md)，把报告结论对应到 JSON 配置、OpenFOAM 运行目录、日志文件、数据表和图片。
  - 覆盖的案例目录：
  - `01_advection_equation/01_sine_wave_quad/`
  - `01_advection_equation/02_sine_wave_quad_linearUpwind/`
  - `01_advection_equation/03_sine_wave_tri_upwind/`
  - `01_advection_equation/03_sine_wave_tri_linearUpwind/`
  - `01_advection_equation/04_solid_rotation_quad_upwind/`
  - `01_advection_equation/04_solid_rotation_tri_upwind/`

- `01_sine_wave_quad/`
  - 题目/算例：第一题，正弦波平移。
  - 网格类型：结构化四边形网格。
  - 空间插值格式：`Gauss upwind`，一阶迎风。
  - `N`：`10, 20, 40, 80`。
  - 速度或速度模型：常速度 `U=(1,1,0)`。
  - 终止时间：`1.0`。
  - 主要检查内容：`L1` 误差、网格收敛阶、数值耗散。
  - JSON 配置：[`01_sine_wave_quad_upwind.json`](../scripts/configs/01_advection_equation/01_sine_wave_quad_upwind.json)。
  - 结果数据：[`data/01_advection_equation/cases/01_sine_wave_quad/`](../data/01_advection_equation/cases/01_sine_wave_quad/)，[`data/01_advection_equation/analysis/01_sine_wave_quad/`](../data/01_advection_equation/analysis/01_sine_wave_quad/)。
  - 结果图片：[`figures/01_advection_equation/cases/01_sine_wave_quad/`](../figures/01_advection_equation/cases/01_sine_wave_quad/)，[`figures/01_advection_equation/analysis/01_sine_wave_quad/`](../figures/01_advection_equation/analysis/01_sine_wave_quad/)。
  - 分析说明：[`analysis.md`](../data/01_advection_equation/analysis/01_sine_wave_quad/analysis.md) 主要记录正弦波一阶迎风四边形网格的 `L1/L2/Linf` 误差、质量守恒、CFL 历史、振幅衰减和观察收敛阶。

- `02_sine_wave_quad_linearUpwind/`
  - 题目/算例：第一题，正弦波平移。
  - 网格类型：结构化四边形网格。
  - 空间插值格式：`Gauss linearUpwind grad(T)`，线性迎风。
  - `N`：`10, 20, 40, 80`。
  - 速度或速度模型：常速度 `U=(1,1,0)`。
  - 终止时间：`1.0`。
  - 主要检查内容：与一阶迎风比较误差、耗散和过冲。
  - JSON 配置：[`02_sine_wave_quad_linearUpwind.json`](../scripts/configs/01_advection_equation/02_sine_wave_quad_linearUpwind.json)。
  - 结果数据：[`data/01_advection_equation/cases/02_sine_wave_quad_linearUpwind/`](../data/01_advection_equation/cases/02_sine_wave_quad_linearUpwind/)，[`data/01_advection_equation/analysis/02_sine_wave_quad_linearUpwind/`](../data/01_advection_equation/analysis/02_sine_wave_quad_linearUpwind/)。
  - 结果图片：[`figures/01_advection_equation/cases/02_sine_wave_quad_linearUpwind/`](../figures/01_advection_equation/cases/02_sine_wave_quad_linearUpwind/)，[`figures/01_advection_equation/analysis/02_sine_wave_quad_linearUpwind/`](../figures/01_advection_equation/analysis/02_sine_wave_quad_linearUpwind/)。
  - 分析说明：[`analysis.md`](../data/01_advection_equation/analysis/02_sine_wave_quad_linearUpwind/analysis.md) 主要记录正弦波线性迎风四边形网格的误差、收敛阶、耗散变化和与一阶迎风相比的精度表现。

- `03_sine_wave_tri_upwind/`
  - 题目/算例：第一题，正弦波平移。
  - 网格类型：Gmsh 三角形棱柱网格。
  - 空间插值格式：`Gauss upwind`，一阶迎风。
  - `N`：`10, 20, 40, 80`。
  - 速度或速度模型：常速度 `U=(1,1,0)`。
  - 终止时间：`1.0`。
  - 主要检查内容：三角形网格上的误差和收敛阶。
  - JSON 配置：[`03_sine_wave_tri_upwind.json`](../scripts/configs/01_advection_equation/03_sine_wave_tri_upwind.json)。
  - 结果数据：[`data/01_advection_equation/cases/03_sine_wave_tri_upwind/`](../data/01_advection_equation/cases/03_sine_wave_tri_upwind/)，[`data/01_advection_equation/analysis/03_sine_wave_tri_upwind/`](../data/01_advection_equation/analysis/03_sine_wave_tri_upwind/)。
  - 结果图片：[`figures/01_advection_equation/cases/03_sine_wave_tri_upwind/`](../figures/01_advection_equation/cases/03_sine_wave_tri_upwind/)，[`figures/01_advection_equation/analysis/03_sine_wave_tri_upwind/`](../figures/01_advection_equation/analysis/03_sine_wave_tri_upwind/)。
  - 分析说明：[`analysis.md`](../data/01_advection_equation/analysis/03_sine_wave_tri_upwind/analysis.md) 主要记录正弦波一阶迎风三角形棱柱网格的误差、收敛阶、真实 cell centre/cell volume 积分和网格类型影响。

- `03_sine_wave_tri_linearUpwind/`
  - 题目/算例：第一题，正弦波平移。
  - 网格类型：Gmsh 三角形棱柱网格。
  - 空间插值格式：`Gauss linearUpwind grad(T)`。
  - `N`：`10, 20, 40, 80`。
  - 速度或速度模型：常速度 `U=(1,1,0)`。
  - 终止时间：`1.0`。
  - 主要检查内容：三角形网格上线性迎风的精度和稳定性。
  - JSON 配置：[`03_sine_wave_tri_linearUpwind.json`](../scripts/configs/01_advection_equation/03_sine_wave_tri_linearUpwind.json)。
  - 结果数据：[`data/01_advection_equation/cases/03_sine_wave_tri_linearUpwind/`](../data/01_advection_equation/cases/03_sine_wave_tri_linearUpwind/)，[`data/01_advection_equation/analysis/03_sine_wave_tri_linearUpwind/`](../data/01_advection_equation/analysis/03_sine_wave_tri_linearUpwind/)。
  - 结果图片：[`figures/01_advection_equation/cases/03_sine_wave_tri_linearUpwind/`](../figures/01_advection_equation/cases/03_sine_wave_tri_linearUpwind/)，[`figures/01_advection_equation/analysis/03_sine_wave_tri_linearUpwind/`](../figures/01_advection_equation/analysis/03_sine_wave_tri_linearUpwind/)。
  - 分析说明：[`analysis.md`](../data/01_advection_equation/analysis/03_sine_wave_tri_linearUpwind/analysis.md) 主要记录正弦波线性迎风三角形棱柱网格的误差、收敛阶、真实 cell centre/cell volume 积分和稳定性表现。

- `04_solid_rotation_quad_upwind/`
  - 题目/算例：第一题，复杂轮廓固体旋转。
  - 网格类型：结构化四边形网格。
  - 空间插值格式：`Gauss upwind`，一阶迎风。
  - `N`：`50, 100, 200`。
  - 速度或速度模型：绕 `(0.5,0.5)` 旋转，角速度 `1`。
  - 终止时间：`2π`。
  - 主要检查内容：最终轮廓、振幅衰减、数值耗散。
  - JSON 配置：[`04_solid_rotation_quad_upwind.json`](../scripts/configs/01_advection_equation/04_solid_rotation_quad_upwind.json)。
  - 结果数据：[`data/01_advection_equation/cases/04_solid_rotation_quad_upwind/`](../data/01_advection_equation/cases/04_solid_rotation_quad_upwind/)。
  - 结果图片：[`figures/01_advection_equation/cases/04_solid_rotation_quad_upwind/`](../figures/01_advection_equation/cases/04_solid_rotation_quad_upwind/)。

- `04_solid_rotation_tri_upwind/`
  - 题目/算例：第一题，复杂轮廓固体旋转。
  - 网格类型：Gmsh 三角形棱柱网格。
  - 空间插值格式：`Gauss upwind`，一阶迎风。
  - `N`：`50, 100, 200`。
  - 速度或速度模型：绕 `(0.5,0.5)` 旋转，角速度 `1`。
  - 终止时间：`2π`。
  - 主要检查内容：三角形和四边形网格的旋转结果比较。
  - JSON 配置：[`04_solid_rotation_tri_upwind.json`](../scripts/configs/01_advection_equation/04_solid_rotation_tri_upwind.json)。
  - 结果数据：[`data/01_advection_equation/cases/04_solid_rotation_tri_upwind/`](../data/01_advection_equation/cases/04_solid_rotation_tri_upwind/)。
  - 结果图片：[`figures/01_advection_equation/cases/04_solid_rotation_tri_upwind/`](../figures/01_advection_equation/cases/04_solid_rotation_tri_upwind/)。

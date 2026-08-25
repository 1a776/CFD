# 画图、误差和收敛性分析命令

本项目不再使用 `studies/`。所有案例脚本放在 `scripts/<caseName>/`，每个网格的 OpenFOAM 数据放在 `cases/<caseName>/Nxx/`。

## 1. 一键运行一个案例的全部 N

以一阶迎风案例为例：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
source /opt/openfoam14/etc/bashrc
sh scripts/build_student_solver.sh
python3 scripts/01_sine_wave_quad/run_all_N.py --overwrite
```

默认运行：

```text
N10, N20, N40, N80
```

二阶迎风案例：

```bash
python3 scripts/02_sine_wave_quad_linearUpwind/run_all_N.py --overwrite
```

脚本对每一个 N 依次执行：

```text
生成/修改 blockMeshDict
    -> 生成 0.orig/T
    -> blockMesh
    -> checkMesh
    -> explicitAdvectionFoamStudent
    -> 解析日志和最终场
    -> 生成单案例数据和图片
```

## 2. 单个 N 的后处理

例如后处理 `N20`：

```bash
python3 scripts/01_sine_wave_quad/plot_results.py \
    --case-dir cases/01_advection_equation/01_sine_wave_quad/N20
```

快速打印误差：

```bash
python3 scripts/01_sine_wave_quad/measure_error.py \
    --case-dir cases/01_advection_equation/01_sine_wave_quad/N20
```

单个 N 的输出：

```text
data/01_advection_equation/cases/01_sine_wave_quad/N20/
├── summary.json
├── time_history.csv
├── field_data.csv
└── error_field.csv

figures/01_advection_equation/cases/01_sine_wave_quad/N20/
├── field_comparison.png
├── diagonal_profile.png
├── amplitude_history.png
└── cfl_history.png
```

## 3. 所有 N 的收敛性分析

完整运行脚本会自动完成收敛性分析。如果每个 N 已经有 `summary.json`，可以只重新分析：

```bash
python3 scripts/01_sine_wave_quad/collect_results.py
python3 scripts/01_sine_wave_quad/analyze_results.py
```

也可以单独重画收敛图：

```bash
python3 scripts/01_sine_wave_quad/plot_convergence.py
```

收敛分析数据：

```text
data/01_advection_equation/analysis/01_sine_wave_quad/
├── raw_results.csv
├── convergence_summary.csv
├── run_manifest.json
└── analysis.md
```

收敛分析图片：

```text
figures/01_advection_equation/analysis/01_sine_wave_quad/
├── convergence_errors.png
├── convergence_order.png
└── all_N_comparison.png
```

## 4. 每张单案例图的含义

### `field_comparison.png`

四个面板分别是初始场、数值最终场、解析最终场和数值误差场：

```text
error = numerical - exact
```

它主要检查波形位置、周期边界传播、数值耗散和局部误差分布。

### `diagonal_profile.png`

提取近似 $x=y$ 的单元中心值，把数值解与解析解画在同一条线上，用于观察波形相位和振幅的差异。

### `amplitude_history.png`

纵轴为：

$$A(t)=\frac{T_{max}(t)-T_{min}(t)}{2}$$

虚线 $A=1$ 是解析正弦波振幅。曲线下降表示数值耗散；曲线超过 1 表示出现过冲。

### `cfl_history.png`

显示每一步实际最大 CFL，并与目标 `CFL=0.2` 对比，用于检查显式时间推进是否按目标 CFL 工作。

## 5. 收敛性图的含义

### `convergence_errors.png`

横轴是每条边的网格数 $N$，纵轴是 $t=1$ 时的归一化误差。若网格加密后误差下降，说明解在趋近解析解。

图中还会绘制一条一阶斜率参考线，用于比较当前格式的误差下降速度。

### `convergence_order.png`

相邻网格的观察收敛阶使用：

$$p=\frac{\log(E_N/E_{2N})}{\log(2)}$$

其中 $E_N$ 是分辨率为 $N$ 的误差。对于一阶迎风格式，目标是观察到接近 1 的收敛阶。

### `all_N_comparison.png`

左图比较所有 N 的 L1 误差，右图比较最终波幅，帮助同时观察精度和数值耗散。

## 6. 数据字段

`summary.json` 和 `raw_results.csv` 中的主要字段：

| 字段 | 含义 |
|---|---|
| `resolution` | 每条边的单元数 N |
| `nCells` | 总单元数 |
| `normalizedL1` | 归一化 L1 误差 |
| `normalizedL2` | 归一化 L2 误差 |
| `normalizedLinf` | 归一化 L∞ 误差 |
| `maxCo` | 实际最大 CFL |
| `finalAmplitude` | 最终数值波幅 |
| `normalizedMassError` | 归一化质量误差 |
| `meshOK` | `checkMesh` 是否通过 |
| `solverEnded` | 求解器是否正常结束 |
| `solverFatal` | 是否出现 fatal error |

## 7. 检查是否产生伪第 201 步

当前 `N20` 的目标时间步数为 200，因为：

```text
endTime = 1.0
deltaT  = 0.005
1.0 / 0.005 = 200
```

查看时间历史行数：

```bash
wc -l data/01_advection_equation/cases/01_sine_wave_quad/N20/time_history.csv
```

结果应为：

```text
201
```

其中 1 行是 CSV 表头，剩下 200 行是实际时间步，不代表存在第 201 个时间步。

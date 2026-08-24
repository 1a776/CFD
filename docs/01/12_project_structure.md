# 工程目录结构和结果归档规则

本项目采用“案例数据、案例脚本、解析数据、图片”分离的结构。目标是：每一个 `N` 都能独立运行，每次运行后都能找到对应的原始场、日志、数据和图片。

## 1. 总体结构

```text
student_project/
├── README.md
├── UDF/
├── build/bin/explicitAdvectionFoamStudent
├── cases/
│   ├── 01_sine_wave_quad/
│   │   ├── N10/
│   │   ├── N20/
│   │   ├── N40/
│   │   └── N80/
│   ├── 02_sine_wave_quad_linearUpwind/
│   │   ├── N10/
│   │   ├── N20/
│   │   ├── N40/
│   │   └── N80/
│   └── 后续案例/
├── scripts/
│   ├── build_student_solver.sh
│   ├── 01_sine_wave_quad/
│   ├── 02_sine_wave_quad_linearUpwind/
│   └── common/
├── data/
│   ├── cases/
│   └── analysis/
├── figures/
│   ├── cases/
│   └── analysis/
└── docs/
    ├── 01/
    └── compare/
```

项目中不再使用 `studies/`，也不建立 `archive/`。旧的外部收敛数据不作为当前分析输入。

## 2. `cases/` 放什么

`cases/<caseName>/Nxx/` 是完整的 OpenFOAM 算例。以 `N20` 为例：

```text
cases/01_sine_wave_quad/N20/
├── Allrun
├── Allclean
├── case.foam
├── 0.orig/
│   ├── U
│   └── T
├── 0/
├── constant/
│   └── polyMesh/
├── system/
│   ├── blockMeshDict
│   ├── controlDict
│   ├── fvSchemes
│   └── fvSolution
├── 0.05/
├── ...
├── 1/
├── log.blockMesh
├── log.checkMesh
├── log.explicitAdvectionFoamStudent
└── run.batch.log
```

其中：

| 内容 | 含义 |
|---|---|
| `0.orig/` | 初始场源文件，长期保留 |
| `0/` | 运行时由 `0.orig/` 恢复的目录 |
| `system/` | 当前 N 的 OpenFOAM 字典 |
| `constant/polyMesh/` | `blockMesh` 生成的运行网格 |
| 时间目录 | 求解器写出的场文件 |
| `log.*` | 各 OpenFOAM 程序的日志 |
| `Allrun` | 该 N 的一键运行脚本 |
| `Allclean` | 清理该 N 的可再生成产物 |

Python 脚本不放在 `cases/` 中。

## 3. `scripts/` 放什么

脚本按案例名称集中保存：

```text
scripts/
├── build_student_solver.sh
├── 01_sine_wave_quad/
│   ├── create_initial_fields.py
│   ├── measure_error.py
│   ├── plot_results.py
│   ├── run_all_N.py
│   ├── collect_results.py
│   ├── plot_convergence.py
│   ├── analyze_results.py
│   └── README.md
├── 02_sine_wave_quad_linearUpwind/
│   ├── create_initial_fields.py
│   ├── measure_error.py
│   ├── plot_results.py
│   ├── run_all_N.py
│   ├── collect_results.py
│   ├── analyze_results.py
│   ├── plot_convergence.py
│   └── README.md
└── common/
    ├── advection_tools.py
    ├── postprocess_case.py
    ├── run_suite.py
    └── study_analysis.py
```

脚本职责如下：

| 脚本 | 作用 |
|---|---|
| `create_initial_fields.py` | 根据指定 N 生成 `sin(2*pi*(x+y))` |
| `measure_error.py` | 快速打印一个 N 的最终误差 |
| `plot_results.py` | 生成一个 N 的数据和四张图 |
| `run_all_N.py` | 生成并运行 N10、N20、N40、N80 |
| `collect_results.py` | 收集每个 N 的 `summary.json` |
| `analyze_results.py` | 计算观察收敛阶并写分析报告 |
| `plot_convergence.py` | 重新绘制所有 N 的分析图 |
| `common/` | 两个案例共享的实现，不作为案例入口使用 |

## 4. `data/` 放什么

单个 N 的数据按案例和分辨率归档：

```text
data/cases/01_sine_wave_quad/N20/
├── summary.json
├── time_history.csv
├── field_data.csv
└── error_field.csv
```

含义是：

- `summary.json`：误差、质量守恒、CFL、振幅、网格和求解器状态；
- `time_history.csv`：每个时间步的时间、步长、CFL、残差、极值和振幅；
- `field_data.csv`：每个单元的坐标、初值、数值解、精确解和误差；
- `error_field.csv`：只保留最终误差相关列，方便后续读取。

所有 N 的研究数据放在：

```text
data/analysis/01_sine_wave_quad/
├── raw_results.csv
├── convergence_summary.csv
├── run_manifest.json
└── analysis.md
```

## 5. `figures/` 放什么

单个 N 的图片：

```text
figures/cases/01_sine_wave_quad/N20/
├── field_comparison.png
├── diagonal_profile.png
├── amplitude_history.png
└── cfl_history.png
```

全部 N 的分析图片：

```text
figures/analysis/01_sine_wave_quad/
├── convergence_errors.png
├── convergence_order.png
└── all_N_comparison.png
```

## 6. 运行命令

编译：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
sh scripts/build_student_solver.sh
```

运行 01 案例的全部 N：

```bash
source /opt/openfoam14/etc/bashrc
python3 scripts/01_sine_wave_quad/run_all_N.py --overwrite
```

运行 02 案例的全部 N：

```bash
source /opt/openfoam14/etc/bashrc
python3 scripts/02_sine_wave_quad_linearUpwind/run_all_N.py --overwrite
```

单独运行一个 N：

```bash
sh cases/01_sine_wave_quad/N20/Allrun
python3 scripts/01_sine_wave_quad/plot_results.py \
    --case-dir cases/01_sine_wave_quad/N20
```

## 7. 清理规则

可以重新生成的内容由各 N 目录下的 `Allclean` 清理：

```text
cases/<caseName>/Nxx/0/
cases/<caseName>/Nxx/<timeDirectory>/
cases/<caseName>/Nxx/constant/polyMesh/
cases/<caseName>/Nxx/log.*
cases/<caseName>/Nxx/postProcessing/
```

`Allclean` 不删除 `0.orig/`、`system/`、`UDF/` 和项目级 Python 脚本。

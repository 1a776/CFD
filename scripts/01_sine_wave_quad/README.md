# 01_sine_wave_quad 脚本

本目录管理 `01_sine_wave_quad` 案例的单网格后处理和多网格分析。案例数据本身只放在：

```text
cases/01_sine_wave_quad/N10/
cases/01_sine_wave_quad/N20/
cases/01_sine_wave_quad/N40/
cases/01_sine_wave_quad/N80/
```

## 一键运行全部 N

先编译求解器并加载 OpenFOAM 环境：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
sh scripts/build_student_solver.sh
source /opt/openfoam14/etc/bashrc
```

运行四个网格：

```bash
python3 scripts/01_sine_wave_quad/run_all_N.py --overwrite
```

脚本会依次完成：

1. 以 `N20` 的字典为模板生成 `N10/N40/N80`；
2. 修改 `blockMeshDict` 的网格数；
3. 按单元中心生成 `0.orig/T`；
4. 运行 `blockMesh`、`checkMesh` 和学生求解器；
5. 对每个 N 写出完整数据和四张单案例图；
6. 计算 L1、L2、Linf 误差和观察收敛阶；
7. 写出收敛图和 `analysis.md`。

## 单个 N

例如只处理已有的 `N20`：

```bash
python3 scripts/01_sine_wave_quad/create_initial_fields.py \
    --case-dir cases/01_sine_wave_quad/N20

python3 scripts/01_sine_wave_quad/measure_error.py \
    --case-dir cases/01_sine_wave_quad/N20

python3 scripts/01_sine_wave_quad/plot_results.py \
    --case-dir cases/01_sine_wave_quad/N20
```

## 输出位置

单个 N 的数据：

```text
data/cases/01_sine_wave_quad/N20/
├── summary.json
├── time_history.csv
├── field_data.csv
└── error_field.csv
```

单个 N 的图：

```text
figures/cases/01_sine_wave_quad/N20/
├── field_comparison.png
├── diagonal_profile.png
├── amplitude_history.png
└── cfl_history.png
```

所有 N 的汇总数据：

```text
data/analysis/01_sine_wave_quad/
├── raw_results.csv
├── convergence_summary.csv
├── run_manifest.json
└── analysis.md
```

所有 N 的对比图：

```text
figures/analysis/01_sine_wave_quad/
├── convergence_errors.png
├── convergence_order.png
└── all_N_comparison.png
```

## 只重画分析结果

已有每个 N 的 `summary.json` 后，可以不重新运行求解器：

```bash
python3 scripts/01_sine_wave_quad/collect_results.py
python3 scripts/01_sine_wave_quad/analyze_results.py
```

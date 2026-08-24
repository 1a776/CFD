# 02_sine_wave_quad_linearUpwind 脚本

本目录管理二阶迎风 `linearUpwind` 案例。它与一阶迎风案例使用相同的网格、初始函数、速度、周期边界和 CFL 设置，区别在于 `system/fvSchemes` 中的 `div(phi,T)`。

运行全部网格：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
sh scripts/build_student_solver.sh
source /opt/openfoam14/etc/bashrc
python3 scripts/02_sine_wave_quad_linearUpwind/run_all_N.py --overwrite
```

默认网格为：

```text
N10, N20, N40, N80
```

单个网格的案例、数据和图片分别位于：

```text
cases/02_sine_wave_quad_linearUpwind/N20/
data/cases/02_sine_wave_quad_linearUpwind/N20/
figures/cases/02_sine_wave_quad_linearUpwind/N20/
```

批量汇总结果位于：

```text
data/analysis/02_sine_wave_quad_linearUpwind/
figures/analysis/02_sine_wave_quad_linearUpwind/
```

单个 N 的后处理命令：

```bash
python3 scripts/02_sine_wave_quad_linearUpwind/plot_results.py \
    --case-dir cases/02_sine_wave_quad_linearUpwind/N20
```

已有全部 N 的单案例结果后，可以重新收集数据、计算观察收敛阶并重画总体图：

```bash
python3 scripts/02_sine_wave_quad_linearUpwind/collect_results.py
python3 scripts/02_sine_wave_quad_linearUpwind/analyze_results.py
python3 scripts/02_sine_wave_quad_linearUpwind/plot_convergence.py
```

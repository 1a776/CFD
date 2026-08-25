# 分阶段检查点

## Checkpoint 0：工程能编译

命令：

```bash
source /opt/openfoam14/etc/bashrc
cd /home/a776/workdocuments/上交船舶/slover/student_project
sh scripts/build_student_solver.sh
build/01_advection_equation/bin/explicitAdvectionFoamStudent -help
```

预期：

- `wmake` 返回成功；
- 帮助信息能正常显示；
- 没有修改 `/opt/openfoam14/src`。

## Checkpoint 1：案例能生成网格

命令：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project/cases/01_advection_equation/01_sine_wave_quad
source /opt/openfoam14/etc/bashrc
sh N20/Allclean
python3 ../../scripts/01_sine_wave_quad/create_initial_fields.py \
  --case-dir N20
cp -R N20/0.orig N20/0
cd N20
blockMesh > log.blockMesh 2>&1
checkMesh > log.checkMesh 2>&1
```

检查：

- `constant/polyMesh` 存在；
- `0/U` 和 `0/T` 存在；
- `checkMesh` 没有致命错误；
- 网格是 `20 x 20 x 1`。

## Checkpoint 2：求解器能读入字段

命令：

```bash
sh Allrun
```

最初 Stage 0 骨架阶段的预期日志包括：

```text
Student checkpoint: reading fields
cells      = 400
velocity   = U
advected   = T
No numerical update has been implemented yet.
```

这段日志只用于回看最初的 Stage 0。当前版本已经继续完成 Stage 1--Stage 5，
不应再把“没有时间推进”当作当前版本的验收标准。

历史上在 Stage 1 和 Stage 2 刚写入时，新增的日志还会包括：

```text
phi min     = ...
phi max     = ...
rate max    = ...
maxCo       = ...
CFL deltaT  = ...
```

在 Stage 2 单独教学阶段，`deltaT` 只是根据 CFL 公式计算出的候选时间步。
当前完整求解器已经把它接入 Stage 5 时间循环。Stage 2 历史日志曾写成：

```text
Stage 2 CFL estimate completed.
No time update has been implemented yet.
```

## Checkpoint 3：你完成 Stage 1

要求：

- 增加 `fvcFlux.H`；
- 创建 `surfaceScalarField phi`；
- 打印 `min(phi)` 和 `max(phi)`；
- 不加入时间循环。

验证问题：

1. 为什么 `phi` 是面场？
2. 为什么 `fvc::flux(U)` 的结果不是 `volScalarField`？
3. `phi` 的量纲是什么？

## Checkpoint 4：Stage 2 CFL 时间步

Stage 2 源码已经写入，下一次编译运行时检查：

- 用 `mag(phi)` 计算面通量绝对值；
- 用 `fvc::surfaceSum` 得到单元面通量总和；
- 除以 `mesh.V()`；
- 用 `gMax` 取最大值；
- 根据 `maxCo` 计算 `deltaT`。

当前代码中这段计算的逐行解释见：

```text
docs/05_stage2_cfl.md
```

最近一次实际验证结果：

```text
rate max    = 80
maxCo       = 0.2
CFL deltaT  = 0.005
```

同时确认：

- `checkMesh` 历史结果为 `Mesh OK`；
- `0/T` 与 `0.orig/T` 一致；
- 没有产生新的时间目录；
- `Stage 2 CFL estimate completed.` 正常输出。

这说明 Stage 2 的 CFL 计算正确执行，但时间推进尚未开始。

特别要区分：

```text
fvc::surfaceSum(mag(phi))
    -> CFL 需要的面通量绝对值总和

fvc::div(phi,T,"div(phi,T)")
    -> Stage 3 才使用的有向对流散度
```

当前源码还增加了：

- `maxCo` 非正时终止并报错；
- `rateMax` 过小或为零时终止并报错。

验证问题：

1. 为什么 CFL 里要有 `1/2`？
2. 为什么分母是单元体积？
3. 最后一小步为什么要取 `min(deltaT, remaining)`？

## Checkpoint 5：Stage 3 显式对流残差

Stage 3 的教学任务是加入：

```cpp
#include "fvcDiv.H"

tmp<volScalarField> tResidual
(
    fvc::div(phi, T, "div(phi,T)")
);
```

它对应：

$$
R_c^n=\frac{1}{V_c}\sum_fF_{cf}T_f^n.
$$

如果手动写 owner/neighbour 面循环，未除体积的面通量和只作为中间量：

$$
\mathcal Q_c^n=\sum_fF_{cf}T_f^n.
$$

完成面通量装配后，$R_c^n=\mathcal Q_c^n/V_c$。
因此 `fvc::div` 返回的就是已经除以体积的最终残差 $R_c$。

详细教学说明见：

```text
docs/06_stage3_convection_residual.md
```

验证命令：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project/cases/01_advection_equation/01_sine_wave_quad
grep -E "residual dimensions|residual min|residual max|residual integral|Stage 3|No time update" \
    log.explicitAdvectionFoamStudent
```

你亲自写完后，要求：

- 残差量纲是 $1/\mathrm{s}$；
- 残差最小值和最大值是有限数；
- 周期案例中 `residual integral` 接近零；
- `T` 没有被修改；
- 没有生成新的时间目录。

## Checkpoint 6：Stage 4 前向 Euler 更新

Stage 4 的前向 Euler 更新已经接入求解器，并在 Stage 5 的每一个时间步中执行。
如果要单独学习它，重点对照：

要求：

1. 用 `fvc::div(phi,T,"div(phi,T)")` 计算残差；
2. 用前向 Euler 更新 `T`；
3. 调用 `T.correctBoundaryConditions()`；
4. 输出更新后的最小值和最大值；
5. 每个时间步的 `T min` 和 `T max` 都是有限数；
6. 日志中没有 `nan` 或 `inf`；
7. 输出时间目录中的 `T` 随时间变化。

如果一步之后场值变成 `nan` 或极大值，先停下来检查量纲、时间步和符号。

## Checkpoint 7：Stage 5 完整时间循环

Stage 5 已经写入当前学生版求解器。验收命令：

```bash
source /opt/openfoam14/etc/bashrc
cd /home/a776/workdocuments/上交船舶/slover/student_project
sh scripts/build_student_solver.sh
cd cases/01_advection_equation/01_sine_wave_quad/N20
sh Allrun
```

要求：

- 从 `runTime.value()` 循环到 `endTime`；
- 每一步设置 `deltaT`；
- 每一步推进 `runTime`；
- 每一步计算残差并更新；
- 按 `writeControl` 输出；
- 最终时间达到 `1`；
- 不因浮点舍入多出一个伪时间步。

## Checkpoint 8：Stage 6 结果后处理

Stage 6 的代码和文档已经加入。单案例图像验收命令：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
python3 scripts/01_sine_wave_quad/plot_results.py
```

要求：

- 生成四联场图、对角线剖面图、振幅历史图和 CFL 历史图；
- 生成 `data/01_advection_equation/cases/01_sine_wave_quad/N20/summary.json`；
- 生成 `data/01_advection_equation/cases/01_sine_wave_quad/N20/time_history.csv`；
- 检查 `normalizedL1`、守恒误差和最大 CFL。

多网格收敛图验收命令：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
python3 scripts/01_sine_wave_quad/collect_results.py
python3 scripts/01_sine_wave_quad/plot_convergence.py
```

详细解释见：

```text
docs/01_advection_equation/09_stage6_visualization_and_convergence.md
docs/compare/命令.md
```

第一版不要急着追求误差很小。先确认误差随网格加密下降，并且 CFL 没有超标。

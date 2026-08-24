# 学习路线

本工程采用“一个公式、一个代码片段、一次验证”的节奏。

## 总路线

| 阶段 | 数学目标 | 代码目标 | 验证方式 |
|---|---|---|---|
| Stage 0 | 明确网格和场 | 读 `mesh`、`U`、`T` | 编译并运行 |
| Stage 1 | $F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_f$ | 创建面通量 `phi` | 打印通量范围 |
| Stage 2 | CFL 时间步 | 计算 `rate` 和 `deltaT` | 打印 `Co` |
| Stage 3 | $R_c^n=\sum_f F_{cf}T_f^n/V_c$ | 调用显式散度 | 打印残差范围 |
| Stage 4 | $T^{n+1}=T^n-\Delta tR^n$ | 更新内部场 | 检查单步量纲和范围 |
| Stage 5 | $t^n\to t^{n+1}$ | 加入时间循环和写出 | 推进到 `t=1` |
| Stage 6 | 离散守恒和误差 | 图像、质量积分、L1 误差 | 多网格收敛 |

Stage 0 的详细逐行讲解见：

```text
docs/02_stage0_code_walkthrough.md
```

Stage 1 面通量的详细讲解见：

```text
docs/04_stage1_face_flux.md
```

Stage 2 CFL 时间步的详细讲解见：

```text
docs/05_stage2_cfl.md
```

Stage 3 显式对流残差的详细讲解见：

```text
docs/06_stage3_convection_residual.md
```

## 每个阶段的工作习惯

先在纸上写出：

1. 当前要实现的公式；
2. 公式中每个量的物理意义；
3. 这个量是 cell、face 还是 global scalar；
4. OpenFOAM 对应的数据类型；
5. 更新发生在旧值、临时值还是新值上。

然后只改一小段代码。改完以后立即：

```bash
sh scripts/build_student_solver.sh
```

编译通过后再运行：

```bash
sh cases/01_sine_wave_quad/N20/Allrun
```

## 当前状态

Stage 1 到 Stage 5 的源码链条已经写入，Stage 6 的后处理脚本也已经加入。
当前完整链条是：

$$
\mathrm{Co}_{\max}
=\frac{\Delta t}{2}
\max_c\left(
\frac{\sum_f|F_{cf}|}{V_c}
\right).
$$

Stage 2 的源码已经在 `explicitAdvectionFoamStudent.C` 中完成：

1. 读取 `controlDict` 中的 `maxCo`；
2. 用 `fvc::surfaceSum(mag(phi))` 得到每个 cell 的通量绝对值总和；
3. 除以 `mesh.V()` 得到每个 cell 的 `rate`；
4. 用 `gMax(rate)` 得到全局最大值；
5. 计算并输出 `deltaT`。

Stage 3 和 Stage 4 现在已经接在一起：

```cpp
fvc::div(phi, T, "div(phi,T)")
```

计算：

$$
R_c^n=\frac{1}{V_c}\sum_fF_{cf}T_f^n.
$$

这里的 $R_c^n$ 从定义上已经除以控制体体积，正好是 OpenFOAM
`fvc::div(phi,T,"div(phi,T)")` 返回的体积归一化离散散度。

随后 Stage 4 用前向 Euler 更新：

$$
T_c^{n+1}=T_c^n-\Delta tR_c^n.
$$

Stage 5 在循环中重复计算残差、更新 `T`、推进 `runTime` 并按
`controlDict` 写出时间层。

Stage 6 使用：

```text
plot_results.py
    -> 场图、剖面图、振幅图、CFL 图

plot_convergence.py
    -> 多网格误差图和观察收敛阶
```

详细说明见：

```text
docs/09_stage6_visualization_and_convergence.md
```

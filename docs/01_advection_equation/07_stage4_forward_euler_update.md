# Stage 4：前向 Euler 单步更新

本文档讲解如何把 Stage 3 得到的体积归一化残差 `residual` 接进显式时间推进。

本阶段只做一个更新步骤，不进入完整时间循环，也不写出新的时间目录。

## 1. 本阶段要实现的公式

Stage 3 已经计算出：

$$R_c^n=\frac{1}{V_c}\sum_{f\in\partial\Omega_c}F_{cf}^nT_f^n.$$

守恒型对流方程半离散以后是：

$$\frac{\mathrm dT_c}{\mathrm dt}=-R_c.$$

前向 Euler 把时间导数近似为：

$$\left.\frac{\mathrm dT_c}{\mathrm dt}\right|_{t^n}
\approx\frac{T_c^{n+1}-T_c^n}{\Delta t}.$$

代入半离散方程，得到：

$$\frac{T_c^{n+1}-T_c^n}{\Delta t}=-R_c^n.$$

解出新时间层：

$$T_c^{n+1}=T_c^n-\Delta tR_c^n.$$

注意：这里没有再除以 $V_c$，因为 `residual` 已经是：

$$R_c^n=\frac{1}{V_c}\sum_fF_{cf}^nT_f^n.$$

如果你使用的是未除体积的 `fluxSum`，即 $\mathcal Q_c^n$，更新式才是：

$$T_c^{n+1}=T_c^n-\frac{\Delta t}{V_c}\mathcal Q_c^n.$$

## 2. 本阶段写入的核心代码

本阶段新增的头文件是：

```cpp
#include "dimensionedScalar.H"
```

它来自 OpenFOAM 的核心类型系统，用于声明 `dimensionedScalar`。
它不是一个需要你额外安装或单独链接的库；`Make/options` 仍然使用已有的
`finiteVolume` 和 `meshTools` 依赖。

源码中新增：

```cpp
const dimensionedScalar deltaTDim("deltaT", dimTime, deltaT);

T = T - deltaTDim*residual;
T.correctBoundaryConditions();

Info<< "  updated T min = " << min(T).value() << nl
    << "  updated T max = " << max(T).value() << nl
    << endl;
```

对应关系是：

| 数学对象 | 代码对象 |
|---|---|
| $\Delta t$ | `deltaTDim` |
| $R^n$ | `residual` |
| $T^n$ | 更新前的 `T` |
| $T^{n+1}$ | 更新后的 `T` |

## 3. 为什么需要 `dimensionedScalar`

Stage 2 得到的 `deltaT` 是普通 `scalar`：

```cpp
const scalar deltaT = 2.0*maxCo/rateMax;
```

它只是一个数。可是 OpenFOAM 的场运算会检查物理量纲。`residual` 的量纲是
时间倒数：

$$[R]=T_{\mathrm{time}}^{-1}.$$

为了让 OpenFOAM 知道 $\Delta t$ 是时间量纲，要写成：

```cpp
const dimensionedScalar deltaTDim("deltaT", dimTime, deltaT);
```

这样：

$$[\Delta tR]=T_{\mathrm{time}}\cdot T_{\mathrm{time}}^{-1}=1.$$

`deltaTDim*residual` 与 `T` 同为无量纲标量场，可以相减。

## 4. 为什么直接写 `T = T - deltaTDim*residual`

这一行就是前向 Euler：

```cpp
T = T - deltaTDim*residual;
```

数学上对应：

$$T_c^{n+1}=T_c^n-\Delta tR_c^n.$$

右端的 `residual` 已经由 Stage 3 使用旧的 `T` 算出，所以这仍然是显式更新。
本阶段不要改成：

```cpp
fvm::ddt(T) + fvm::div(phi, T)
```

那会进入隐式矩阵路线，不是本题正在训练的显式求解器。

## 5. 为什么要调用 `T.correctBoundaryConditions()`

更新内部场以后，边界上的 `T` 也要同步。代码是：

```cpp
T.correctBoundaryConditions();
```

对于当前周期正弦波案例，`cyclic` 边界需要把配对边界上的值重新对应起来；
对于 `empty` 边界，也需要保持二维算例的边界语义正确。

简单说，内部 cell 更新以后，要让边界 patch 跟着恢复到一致状态。

## 6. 本阶段仍然没有做什么

Stage 4 只做一次更新。它仍然没有：

- `runTime++`；
- `runTime.write()`；
- `while` 时间循环；
- 推进到 `t = 1`；
- 计算最终误差。

所以运行后你应该看到 `updated T min/max`，但不应该看到新的时间目录。

## 7. 运行后的验收标准

由你执行以下命令：

```bash
source /opt/openfoam14/etc/bashrc
cd /home/a776/workdocuments/上交船舶/slover/student_project
sh scripts/build_student_solver.sh
build/01_advection_equation/bin/explicitAdvectionFoamStudent -help
cd cases/01_advection_equation/01_sine_wave_quad/N20
sh Allrun
```

然后检查日志：

```bash
grep -E "cells|phi min|phi max|rate max|maxCo|CFL deltaT|residual dimensions|residual min|residual max|residual integral|updated T min|updated T max|Stage 4|No time loop" \
    log.explicitAdvectionFoamStudent
```

### 7.1 编译和启动验收

1. `sh scripts/build_student_solver.sh` 成功返回；
2. `build/01_advection_equation/bin/explicitAdvectionFoamStudent -help` 能输出帮助信息；
3. `Allrun` 成功完成 `blockMesh`、`checkMesh` 和学生版求解器；
4. 日志中应看到 `cells = 400` 或等价的 `cells      = 400`；
5. `checkMesh` 日志不出现致命错误，并包含 `Mesh OK`。

### 7.2 Stage 2--Stage 3 数值验收

1. 对当前 `20 x 20` 网格，预期：

```text
rate max    = 80
maxCo       = 0.2
CFL deltaT  = 0.005
```

2. `residual dimensions` 应为时间倒数量纲，通常类似：

```text
[0 0 -1 0 0 0 0]
```

3. `residual min`、`residual max` 必须是有限数；
4. 周期正弦波案例中，`residual integral` 应接近 `0`。

### 7.3 Stage 4 单步更新验收

1. `updated T min` 和 `updated T max` 是有限数；
2. 它们应保持在初始场范围 `[-1,1]` 附近，不应出现明显爆炸；
3. 日志中出现 `Stage 4 forward Euler update completed.`；
4. 日志中出现 `No time loop has been implemented yet.`；
5. 当前源码没有 `runTime++`、`runTime.write()` 和 `while` 时间循环；
6. 磁盘上不应生成 `0.005`、`0.01`、`1` 等新的时间目录，只应保留
   `0` 和 `0.orig` 这两个初始目录；
7. 因为当前没有调用 `runTime.write()`，可以用下面的命令确认磁盘上的
   `0/T` 没有被单步内存更新覆盖：

```bash
cmp 0/T 0.orig/T
```

如果这些都成立，说明 Stage 4 的单步显式更新已经接通。

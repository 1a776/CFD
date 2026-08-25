# Stage 5：完整时间循环与结果输出

本文档把 Stage 4 的单步前向 Euler 更新扩展成完整的时间推进。
本阶段开始以后，求解器会从初始时间推进到 `controlDict` 中的 `endTime`，
并按照 `writeControl` 和 `writeInterval` 输出结果。

本阶段仍然不做误差收敛率分析。误差分析属于 Stage 6。

## 1. Stage 5 要解决什么问题

Stage 4 只有一次更新：

$$T_c^{n+1}=T_c^n-\Delta tR_c^n.$$

它只能验证一件事：一个时间步的公式、量纲和边界更新是否接通。

真正的瞬态计算需要不断重复：

$$
T^0\longrightarrow T^1\longrightarrow T^2
\longrightarrow\cdots\longrightarrow T^N.
$$

最终要求是：

$$
t^N=t_{\mathrm{end}}.
$$

当前算例中：

```foam
startTime   0;
endTime     1;
maxCo       0.2;
```

所以求解器应从 $t=0$ 推进到 $t=1$。

## 2. 一个时间步的数学顺序

设当前时间为 $t^n$，当前标量场为 $T^n$。

### 2.1 计算当前剩余时间

$$
t_{\mathrm{remain}}=t_{\mathrm{end}}-t^n.
$$

代码：

```cpp
const scalar remainingTime = endTime - runTime.value();
```

`runTime.value()` 是当前物理时间 $t^n$，`endTime` 是终止时间。

### 2.2 选择本步时间步

Stage 2 已经得到 CFL 候选时间步：

$$
\Delta t_{\mathrm{CFL}}
=
\frac{2\,\mathrm{Co}_{\mathrm{target}}}
{\max_c\left(\sum_f|F_{cf}|/V_c\right)}.
$$

正常情况下：

$$
\Delta t^n=\Delta t_{\mathrm{CFL}}.
$$

最后一步不能超过终止时间，因此取：

$$
\Delta t^n
=
\min\left(\Delta t_{\mathrm{CFL}},t_{\mathrm{end}}-t^n\right).
$$

代码：

```cpp
const scalar stepDeltaT = min(deltaT, remainingTime);
```

这里的 `deltaT` 是 Stage 2 得到的 CFL 候选值，`stepDeltaT` 是本次循环真正采用的时间步。

### 2.3 把时间步交给 OpenFOAM

```cpp
runTime.setDeltaT(stepDeltaT);
```

这一步只修改 `Time` 对象保存的时间步长，还没有让物理时间前进。

### 2.4 推进时间标签

```cpp
runTime++;
```

数学上：

$$
t^n\longrightarrow t^{n+1}=t^n+\Delta t^n.
$$

注意：此时 `T` 仍然保存更新前的 $T^n$。这正好符合 OpenFOAM 瞬态求解器
常用的组织方式：先进入新时间层，再用旧场计算方程右端项并更新场。

### 2.5 重新计算空间残差

```cpp
tmp<volScalarField> tResidual
(
    fvc::div(phi, T, "div(phi,T)")
);
```

它对应：

$$
R_c^n=\frac{1}{V_c}\sum_fF_{cf}^nT_f^n.
$$

为什么必须每一步重算？因为 $T$ 在上一步更新后已经变化：

$$
T^n\neq T^{n-1}.
$$

虽然本算例速度 $\boldsymbol U$ 固定，所以面通量 `phi` 可以复用；
但对流标量 `T` 变了，残差一定要重新计算。

### 2.6 执行前向 Euler 更新

```cpp
const dimensionedScalar deltaTDim("deltaT", dimTime, stepDeltaT);

T = T - deltaTDim*residual;
T.correctBoundaryConditions();
```

数学上：

$$
T_c^{n+1}=T_c^n-\Delta t^nR_c^n.
$$

`residual` 已经除以控制体体积，所以这里不能再写一遍 `/mesh.V()`。

### 2.7 输出当前时间层

```cpp
runTime.write();
```

这并不意味着每一步一定生成目录。是否真正写出由
`system/controlDict` 控制：

```foam
writeControl    timeStep;
writeInterval   10;
```

含义是每 10 个时间步写一次。当前 `deltaT=0.005` 时，理论上会输出：

```text
0.05, 0.10, 0.15, ..., 1.00
```

`T` 的 `IOobject` 使用了 `AUTO_WRITE`，所以 `runTime.write()` 会把它写出。
速度场 `U` 是 `NO_WRITE`，不会被重复输出。

### 2.8 为什么循环结束后还要强制写一次

`runTime.write()` 服从 `writeControl` 和 `writeInterval`。因此最后一个时间步虽然
已经把求解器推进到

$$
t_{\mathrm{end}},
$$

但如果这个时间步不是第 `writeInterval` 个步，最后场可能不会自动生成对应的时间目录。
这会导致后处理读取到上一个写出时刻，而不是题目要求的最终时刻。

因此时间循环结束后，学生版求解器额外调用：

```cpp
const bool finalWriteOK = runTime.writeNow();
```

它强制把当前时间层的 `T` 写出。对于刚体旋转案例，当前时间层就是：

$$
t_{\mathrm{end}}=2\pi.
$$

同时，案例生成脚本使用高精度把 JSON 中的 `endTime` 写入 `controlDict`，并将
`timePrecision` 提高到 17，避免把
$2\pi=6.283185307179586\ldots$ 截断为 `6.28319` 或在时间目录名中丢失有效数字。

## 3. 为什么 `phi` 在循环外，而 `residual` 在循环内

当前代码的结构是：

```text
U
  -> fvc::flux(U)
  -> phi                  # 循环外计算一次

T^n
  -> fvc::div(phi,T)
  -> residual
  -> T^(n+1)               # 每一步重复
```

原因是本算例假设：

$$
\boldsymbol U(\boldsymbol x,t)=\boldsymbol U(\boldsymbol x).
$$

也就是速度场不随时间变化。因此：

$$
\phi^n=\phi.
$$

如果以后速度场也需要求解或随时间变化，就不能简单复用循环外的 `phi`，
而要在每个时间步重新执行：

```cpp
phi = fvc::flux(U);
```

这属于后续扩展，不是当前固定速度线性对流案例的要求。

## 4. `runTime`、`T` 和时间目录的关系

### 4.1 时间对象

```cpp
Time runTime;
```

由 `createTime.H` 创建。它管理：

- 当前时间数值；
- 当前时间目录名；
- `endTime`；
- `deltaT`；
- 输出控制；
- 时间推进索引。

### 4.2 内存中的更新

```cpp
T = T - deltaTDim*residual;
```

只改变内存里的 `T`。

### 4.3 磁盘上的输出

```cpp
runTime.write();
```

才会根据输出策略把 `T` 写入当前时间目录。

所以这三件事要区分：

| 操作 | 数学/工程含义 |
|---|---|
| `runTime.setDeltaT(...)` | 选择本步 $\Delta t$ |
| `runTime++` | 进入 $t^{n+1}$ |
| `runTime.write()` | 按规则写出当前字段 |

## 5. 当前代码中每个 Stage 5 位置

在 `explicitAdvectionFoamStudent.C` 中：

| 代码 | 对应内容 |
|---|---|
| `runTime.endTime().value()` | 读取 $t_{\mathrm{end}}$ |
| `while (runTime.value() < endTime - SMALL)` | 时间循环条件 |
| `remainingTime` | 计算剩余时间 |
| `stepDeltaT` | 处理最后一步的实际时间步 |
| `runTime.setDeltaT(stepDeltaT)` | 设置 OpenFOAM 时间步 |
| `runTime++` | 时间从 $t^n$ 进入 $t^{n+1}$ |
| `fvc::div(phi,T,...)` | 计算当前残差 |
| `T = T - deltaTDim*residual` | 前向 Euler 更新 |
| `T.correctBoundaryConditions()` | 修正边界条件 |
| `runTime.write()` | 输出当前时间层 |

完整链条是：

```text
current time
    -> remaining time
    -> stepDeltaT
    -> setDeltaT
    -> runTime++
    -> fvc::div(phi,T)
    -> T = T - deltaT*residual
    -> correctBoundaryConditions
    -> runTime.write
    -> next time step
```

## 6. 本阶段的日志

每个时间步会输出：

```text
Time = ...
step = ...
deltaT = ...
maxCo = ...
residual min = ...
residual max = ...
residual integral = ...
T min = ...
T max = ...
```

最后会输出：

```text
Stage 5 time loop completed.
final time = 1
time steps = ...
```

其中 `maxCo` 的计算是：

$$
\mathrm{Co}_{\max}^{n}
=
\frac{\Delta t^n}{2}
\max_c\left(\frac{\sum_f|F_{cf}|}{V_c}\right).
$$

正常时间步应接近目标值 $0.2$；最后一步如果被缩短，实际值可以小于 $0.2$。

## 7. 运行与验收

本轮不由助手运行。请你执行：

```bash
source /opt/openfoam14/etc/bashrc
cd /home/a776/workdocuments/上交船舶/slover/student_project
sh scripts/build_student_solver.sh
build/bin/explicitAdvectionFoamStudent -help
cd cases/01_sine_wave_quad/N20
sh Allrun
```

检查求解器日志：

```bash
grep -E "Starting Stage 5|Time =|step =|deltaT =|maxCo =|residual integral|T min|T max|Stage 5 time loop completed|final time|time steps" \
    log.explicitAdvectionFoamStudent
```

检查网格：

```bash
grep -E "Mesh OK|cells:" log.checkMesh
```

检查输出时间目录：

```bash
find . -maxdepth 1 -type d -printf '%f\n' | sort -V
```

### 7.1 编译验收

- `sh scripts/build_student_solver.sh` 成功；
- `explicitAdvectionFoamStudent -help` 能正常输出；
- `/opt/openfoam14/src` 没有被修改。

### 7.2 时间推进验收

- 日志包含 `Starting Stage 5 time loop`；
- 日志包含多个不同的 `Time =`；
- 日志包含 `Stage 5 time loop completed.`；
- `final time` 应为 `1` 或数值误差范围内接近 `1`；
- `time steps` 对当前网格和 `maxCo=0.2` 应接近 `200`；
- 正常步的 `maxCo` 应接近 `0.2`，最后一步可以略小；
- 日志中不应再出现 `No time loop has been implemented yet.`。

### 7.3 输出验收

- `0.05`、`0.10`、...、`1` 等写出目录应按 `writeInterval 10` 出现；
- `0/T` 保留初始场；
- 输出时间目录中的 `T` 文件存在；
- `T min`、`T max` 始终是有限数；
- 不应出现 `nan`、`inf` 或明显爆炸。

### 7.4 守恒与边界验收

- 周期案例中每一步 `residual integral` 应接近 `0`；
- `xMin/xMax` 和 `yMin/yMax` 仍是周期 patch；
- `zMin/zMax` 仍是 `empty`；
- 每一步都执行 `T.correctBoundaryConditions()`；
- 场值总体应保持在合理范围内。

## 8. Stage 5 还没有完成什么

本阶段已经实现完整时间循环，但仍然没有：

- 解析精确解并计算最终 L1 误差；
- 自动运行多个网格分辨率；
- 输出收敛率；
- 手写三角形网格或非结构网格；
- 处理随时间变化的速度场；
- 并行运行验证。

这些属于 Stage 6 或后续扩展。

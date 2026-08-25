# Stage 2：CFL 稳定时间步

本文档讲解如何根据 Stage 1 得到的面通量 `phi` 计算显式方法允许的时间步长。

本阶段只计算和检查时间步，不更新 `T`，也不进入时间循环。

## 0. Stage 2 使用的头文件和源码

学生版求解器新增：

```cpp
#include "fvcSurfaceIntegrate.H"
```

这个头文件的职责是声明：

```cpp
fvc::surfaceSum(...)
```

OpenFOAM 14 本地相关文件是：

```text
/opt/openfoam14/src/finiteVolume/finiteVolume/fvc/fvcSurfaceIntegrate.H
/opt/openfoam14/src/finiteVolume/finiteVolume/fvc/fvcSurfaceIntegrate.C
```

其中 `.H` 主要提供函数声明，`.C` 包含模板实现。由于当前 `wmake` 使用
`-DNoRepository`，头文件末尾会包含对应的模板实现：

```cpp
#ifdef NoRepository
    #include "fvcSurfaceIntegrate.C"
#endif
```

所以你的求解器可以直接调用：

```cpp
fvc::surfaceSum(mag(phi))
```

而不需要自己重新写 owner/neighbour 遍历。

## 1. 为什么不能随便选择 `deltaT`

显式 Euler 更新的形式是：

$$T_c^{n+1}=T_c^n-\Delta tR_c^n.$$

信息通过速度从一个控制体传播到相邻控制体。如果 $\Delta t$ 太大，
一次更新可能跨过很多个网格单元，数值解就可能出现：

- 振荡；
- 负值或超出初始范围；
- 误差快速放大；
- 数值发散。

因此显式方法必须限制时间步长。CFL 数就是“一个时间步内，流动跨过多少网格尺度”的无量纲衡量。

本工程使用的定义是：

$$\mathrm{Co}_c
=\frac{\Delta t}{2V_c}
\sum_{f\in\partial\Omega_c}|F_{cf}|.$$

全网格的最大 Courant 数是：

$$\mathrm{Co}_{\max}
=\max_c\left[
\frac{\Delta t}{2V_c}
\sum_f|F_{cf}|
\right].$$

代码中的 `maxCo` 是目标上限：

```foam
maxCo 0.2;
```

因此需要选择 $\Delta t$，使：

$$\mathrm{Co}_{\max}\leq\mathrm{maxCo}.$$

## 2. 为什么使用通量绝对值

对于一个封闭控制体，带符号的面通量可能互相抵消：

$$\sum_fF_{cf}=0.$$

例如均匀速度通过一个矩形控制体时，左侧和右侧通量符号相反，
上侧和下侧通量也可能符号相反。

但是稳定性关心的不是“净流出多少”，而是“总共有多少流动穿过面”。
因此需要：

$$\sum_f|F_{cf}|.$$

`abs` 的作用是避免正负通量抵消后错误地得到一个很小的稳定性指标。

## 3. `mag(phi)`：对所有面通量取绝对值

Stage 1 中：

```cpp
surfaceScalarField phi;
```

表示每个面一个有符号通量：

$$\{F_f\}_{f=1}^{N_f}.$$

调用：

```cpp
mag(phi)
```

后得到：

$$\{|F_f|\}_{f=1}^{N_f}.$$

它仍然是面场，类型仍然是面标量场，只是每个值都变成了非负数。

## 4. `fvc::surfaceSum(mag(phi))`

相关接口声明在：

```text
/opt/openfoam14/src/finiteVolume/finiteVolume/fvc/fvcSurfaceIntegrate.H
```

接口大致是：

```cpp
tmp<VolInternalField<Type>> surfaceSum
(
    const SurfaceField<Type>& ssf
);
```

输入是面场，输出是单元内部场：

```text
surface field
    -> 对每个 cell 周围的面求和
        -> cell field
```

所以：

```cpp
fvc::surfaceSum(mag(phi))
```

对应数学量：

$$
\left\{
\sum_{f\in\partial\Omega_c}|F_{cf}|
\right\}_{c=1}^{N_c}.
$$

它的结果是每个 cell 一个数，而不是每个 face 一个数。

OpenFOAM 源码中，内部面处理大致是：

```cpp
vf[owner[facei]] += ssf[facei];
vf[neighbour[facei]] += ssf[facei];
```

意思是：一个内部面属于两个相邻单元，因此该面的绝对通量要同时加到
owner 单元和 neighbour 单元的总和中。

边界面则通过边界 patch 的 `faceCells()` 加回所属单元。

### 4.1 为什么这里可以使用 `surfaceSum`

这里的输入是：

$$
|F_f|\geq0.
$$

因此一个内部面上的同一个绝对通量可以同时计入相邻的两个控制体：

$$
\sum_{f\in\partial\Omega_c}|F_{cf}|.
$$

这和计算有方向的净通量不同。`surfaceSum` 在这里的任务只是累加“通量大小”，
服务于 CFL 稳定性估计；它不是 Stage 3 中带符号的对流散度。

需要区分：

| 操作 | 数学意义 | 代码 |
|---|---|---|
| 面通量绝对值求和 | $\sum_f|F_{cf}|$ | `fvc::surfaceSum(mag(phi))` |
| 有方向的对流散度 | $\frac{1}{V_c}\sum_fF_{cf}T_f$ | `fvc::div(phi,T,"div(phi,T)")` |

前者用于 CFL，后者用于 PDE 更新。两者都涉及面，但数学任务不同。

## 5. 为什么要使用 `.primitiveField()`

`fvc::surfaceSum(mag(phi))` 返回的是一个 OpenFOAM 临时场对象：

```cpp
tmp<VolInternalField<scalar>>
```

为了把它保存成普通的单元标量数组，可以写成：

```cpp
scalarField sumPhi
(
    fvc::surfaceSum(mag(phi))().primitiveField()
);
```

这里分成三步理解：

| 表达式 | 作用 |
|---|---|
| `fvc::surfaceSum(mag(phi))` | 返回临时单元场 |
| `()` | 访问临时对象中的实际场 |
| `.primitiveField()` | 取得单元内部的普通数值数组 |

结果是：

```text
sumPhi[celli]
```

表示第 `celli` 个控制体周围的：

$$\sum_f|F_{cf}|.$$

### 5.1 `tmp<T>` 是什么

OpenFOAM 的很多显式算子返回：

```cpp
tmp<T>
```

它可以理解为：

```text
一个由 OpenFOAM 管理生命周期的临时 T 对象
```

本例中：

```cpp
fvc::surfaceSum(mag(phi))
```

的返回类型近似为：

```cpp
tmp<VolInternalField<scalar>>
```

它不是最终的 `scalarField`，而是暂时包着一个“每个 cell 一个 scalar”的场。
使用 `tmp` 的好处是：复杂表达式中产生的中间场可以自动管理，减少不必要的复制和手动内存释放。

### 5.2 后面的 `()` 不是函数参数

表达式：

```cpp
fvc::surfaceSum(mag(phi))()
```

看起来像连续调用两个函数，但含义不同：

```text
fvc::surfaceSum(mag(phi))
    -> 调用 surfaceSum 函数，返回 tmp<...>

()
    -> 调用 tmp 对象的 operator()()
    -> 取出它内部包着的实际 VolInternalField
```

OpenFOAM 14 本地 `tmp` 实现中的核心逻辑是：

```cpp
template<class T>
inline const T& Foam::tmp<T>::operator()() const
{
    return *ptr_;
}
```

所以这一步的结果是对临时内部场的常量引用：

```cpp
const VolInternalField<scalar>&
```

它没有重新计算 `surfaceSum`，也没有把面场变成普通数组；它只是打开 `tmp` 包装。

### 5.3 `.primitiveField()` 是什么

OpenFOAM 的几何场不仅有内部 cell 值，还可能有边界 patch 值。对于一个体场，可以概念性地写成：

```text
GeometricField
├── primitiveField()   -> cell 内部值
└── boundaryField()    -> 边界 patch 值
```

`primitiveField()` 的本地接口返回内部场的常量引用：

```cpp
const typename Internal::FieldType&
primitiveField() const;
```

因此：

```cpp
fvc::surfaceSum(mag(phi))().primitiveField()
```

最终得到的是：

```cpp
const scalarField&
```

它的元素可以写成：

```cpp
sumPhi[celli]
```

对应：

$$
\left(\sum_f|F_{cf}|\right)_{celli}.
$$

这里不是把边界贡献丢掉了。边界面的贡献在 `surfaceSum` 计算阶段已经加到了相应的 owner cell；
`.primitiveField()` 只是取出已经算好的 cell 结果。

### 5.4 为什么还要构造一个 `scalarField`

源码写成：

```cpp
scalarField sumPhi
(
    fvc::surfaceSum(mag(phi))().primitiveField()
);
```

左边的 `scalarField sumPhi` 是一个独立的普通内部数组，保存：

$$
\{\,\sum_f|F_{cf}|\,\}_{c=1}^{N_c}.
$$

这样后面可以直接进行逐元素除法：

```cpp
sumPhi/mesh.V().primitiveField()
```

数学上就是对每个 cell 做：

$$
\frac{\sum_f|F_{cf}|}{V_c}.
$$

## 6. `mesh.V()`：每个控制体的体积

有限体积公式中需要除以控制体体积：

$$
r_c=\frac{1}{V_c}\sum_f|F_{cf}|.
$$

OpenFOAM 中：

```cpp
mesh.V()
```

返回每个 cell 的体积场。

本案例中：

$$V_c=0.05\times0.05\times0.1=0.00025.$$

因此可以计算：

```cpp
scalarField rate
(
    sumPhi/mesh.V().primitiveField()
);
```

数学上对应：

$$
r_c=
\frac{\sum_f|F_{cf}|}{V_c}.
$$

由于通量的量纲是 $L^3/T$，体积的量纲是 $L^3$，所以：

$$[r_c]=\frac{L^3/T}{L^3}=\frac{1}{T}.$$

## 7. `gMax`：找全局最大值

我们需要的是：

$$r_{\max}=\max_c r_c.$$

OpenFOAM 中可以使用：

```cpp
scalar rateMax = gMax(rate);
```

`gMax` 的含义是 global maximum：

- 串行运行时：求当前进程中所有 cell 的最大值；
- 并行运行时：在所有进程之间进行全局归约，得到整个网格的最大值。

这就是为什么不应只用某个局部 cell 的值决定时间步。

## 8. 根据目标 `maxCo` 反推时间步

根据：

$$
\mathrm{Co}_{\max}
=\frac{\Delta t}{2}r_{\max},
$$

得到：

$$
\Delta t
=\frac{2\,\mathrm{Co}_{\mathrm{target}}}{r_{\max}}.
$$

在代码中，`maxCo` 可以从 `controlDict` 读取：

```cpp
const scalar maxCo
(
    controlDict.lookupOrDefault<scalar>("maxCo", 0.2)
);
```

随后根据 `rateMax` 计算：

```cpp
const scalar deltaT = 2.0*maxCo/rateMax;
```

本阶段还可以输出：

```cpp
Info<< "  rate max    = " << rateMax << nl
    << "  maxCo       = " << maxCo << nl
    << "  deltaT      = " << deltaT << nl
    << endl;
```

## 9. `maxCo`、`rate`、`rateMax` 和 `deltaT` 的类型

### 9.1 `maxCo` 为什么是 `scalar`

`maxCo` 是整个计算案例共享的一个目标数，不是每个 cell 一个值，也不是每个 face 一个值：

$$\mathrm{maxCo}=0.2.$$

所以使用：

```cpp
const scalar maxCo
(
    controlDict.lookupOrDefault<scalar>("maxCo", 0.2)
);
```

逐项解释：

| 代码 | 含义 |
|---|---|
| `const` | 读取后不在本阶段修改 |
| `scalar` | 一个普通浮点数 |
| `controlDict` | `system/controlDict` 对应的字典 |
| `lookupOrDefault<scalar>` | 读取 scalar；缺失时使用默认值 |
| `"maxCo"` | 字典关键字 |
| `0.2` | 默认目标值 |

当前案例的字典中有：

```foam
maxCo 0.2;
```

因此运行时会使用案例参数，而不是默认值。

### 9.2 `rate` 为什么是 `scalarField`

每个 cell 都有一个：

$$
r_c=\frac{\sum_f|F_{cf}|}{V_c}.
$$

所以 `rate` 不是一个普通 `scalar`，而是：

```cpp
scalarField rate;
```

其大小应等于 cell 数：

```text
rate.size() = mesh.nCells()
```

### 9.3 `rateMax` 为什么又变回 `scalar`

`gMax(rate)` 把所有 cell 的 rate 归约成一个数：

$$
r_{\max}=\max_c r_c.
$$

所以：

```cpp
const scalar rateMax = gMax(rate);
```

返回的是全局标量。

### 9.4 `deltaT` 为什么是 `scalar`

当前实现中：

```cpp
const scalar deltaT = 2.0*maxCo/rateMax;
```

它表示根据 CFL 公式计算出的全局时间步长：

$$[\Delta t]=T.$$

后续 Stage 5 会把它交给：

```cpp
runTime.setDeltaT(deltaT);
```

当前 Stage 2 还没有调用这句，因此现在只是计算时间步，不推进时间。

## 10. 两个安全检查为什么存在

代码中有：

```cpp
if (maxCo <= 0)
{
    FatalErrorInFunction
        << "maxCo must be positive, but received " << maxCo
        << exit(FatalError);
}
```

目标 Courant 数必须是正数。若为零或负数：

- `deltaT` 会变成零或负数；
- 时间推进失去意义；
- 用户输入错误不应被静默接受。

`FatalErrorInFunction` 会报告当前函数和错误信息，
`exit(FatalError)` 会终止程序，而不是继续用错误参数计算。

第二个检查是：

```cpp
if (rateMax <= SMALL)
{
    FatalErrorInFunction
        << "The maximum flux rate is zero or too small: "
        << rateMax << nl
        << "Cannot determine a positive CFL time step."
        << exit(FatalError);
}
```

因为公式：

$$
\Delta t=\frac{2\,\mathrm{maxCo}}{r_{\max}}
$$

要求 $r_{\max}>0$。如果速度场为零，所有面通量也可能为零，此时：

$$r_{\max}=0.$$

直接相除会产生除零问题。当前代码选择明确报错，提醒后续需要决定“静止场的时间步策略”，
而不是让错误的无穷大时间步悄悄进入时间循环。

`SMALL` 是 OpenFOAM 提供的小正数阈值，用于避免把接近零的浮点数当成安全分母。

## 11. 当前 `.C` 代码的逐行数学映射

当前源码中的核心段落是：

```cpp
const scalar maxCo
(
    controlDict.lookupOrDefault<scalar>("maxCo", 0.2)
);

scalarField sumPhi
(
    fvc::surfaceSum(mag(phi))().primitiveField()
);

scalarField rate
(
    sumPhi/mesh.V().primitiveField()
);

const scalar rateMax = gMax(rate);

const scalar deltaT = 2.0*maxCo/rateMax;
```

逐行映射为：

| C++ 代码 | 数学对象 | 所在位置 |
|---|---|---|
| `maxCo` | $\mathrm{Co}_{target}$ | 全局 scalar |
| `mag(phi)` | $\{|F_f|\}$ | face |
| `surfaceSum(...)` | $\{\sum_f|F_{cf}|\}$ | cell |
| `mesh.V()` | $\{V_c\}$ | cell |
| `sumPhi/mesh.V()` | $\{r_c\}$ | cell |
| `gMax(rate)` | $r_{\max}$ | 全局 scalar |
| `2.0*maxCo/rateMax` | $\Delta t$ | 全局 scalar |

整个 Stage 2 的数据位置变化是：

```text
surfaceScalarField phi
    -> surfaceScalarField mag(phi)
        -> tmp<VolInternalField<scalar>>
            -> const scalarField& primitiveField()
                -> scalarField sumPhi
                    -> scalarField rate
                        -> scalar rateMax
                            -> scalar deltaT
```

这条链展示了有限体积代码中“位置”的变化：

```text
face 数据
    -> cell 数据
        -> 全局数据
```

## 12. 当前案例的手算结果

当前速度为：

$$\boldsymbol U=(1,1,0).$$

对于一个 cell：

- 两个 x 方向面各有 $|F_f|=0.005$；
- 两个 y 方向面各有 $|F_f|=0.005$；
- 两个 z 方向面通量为 $0$。

因此：

$$
\sum_f|F_{cf}|
=0.005+0.005+0.005+0.005
=0.02.
$$

又因为：

$$V_c=0.00025,$$

所以：

$$
r_c=\frac{0.02}{0.00025}=80\ \mathrm{s}^{-1}.
$$

当前所有单元相同，因此：

$$r_{\max}=80\ \mathrm{s}^{-1}.$$

当：

$$\mathrm{maxCo}=0.2,$$

允许的时间步为：

$$
\Delta t
=\frac{2\times0.2}{80}
=0.005\ \mathrm{s}.
$$

这说明 `controlDict` 中原来的：

```foam
deltaT 0.01;
```

对于目标 `maxCo 0.2` 来说偏大。Stage 2 需要根据 CFL 计算得到 `0.005`，
后续时间推进时使用这个值。

## 13. 本阶段的代码边界

Stage 2 只完成：

```text
phi
    -> mag(phi)
        -> surfaceSum
            -> 除以 mesh.V()
                -> gMax
                    -> deltaT
```

本阶段暂时不做：

- 更新 `T`；
- 调用 `fvc::div`；
- `runTime++`；
- 时间循环；
- 写出新的时间目录。

## 14. Stage 2 当前状态

Stage 2 代码已经写入：

```text
UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/explicitAdvectionFoamStudent.C
```

本轮只完成了源码编辑和静态对应关系检查，尚未编译和运行。下一次运行时，
建议检查结果是：

```text
rate max  approximately 80
maxCo     = 0.2
deltaT    approximately 0.005
```

Stage 2 仍然没有更新 `T`。如果这些结果确认正确，下一阶段才进入：

```cpp
fvc::div(phi, T, "div(phi,T)")
```

## 15. 按代码执行顺序重新读一遍

前面的章节分别解释了公式和 API。下面把当前 `.C` 中的 Stage 2
按“程序真正执行的顺序”串起来。建议你打开：

```text
UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/explicitAdvectionFoamStudent.C
```

对照下面的代码阅读。

### 15.1 第一步：从案例参数中读取目标 CFL 数

当前代码是：

```cpp
const scalar maxCo
(
    controlDict.lookupOrDefault<scalar>("maxCo", 0.2)
);
```

这不是在计算 CFL，而是在读取“希望 CFL 不超过多少”的用户参数。

这一段直接依赖的头文件是：

```cpp
#include "dictionary.H"
```

这里要区分“参数名”和“代码接口”：

| 对象 | 是什么 |
|---|---|
| `maxCo` | `system/controlDict` 中的关键字 |
| `controlDict` | 一个 OpenFOAM 字典对象 |
| `lookupOrDefault` | `dictionary` 类提供的查值成员函数 |
| `<scalar>` | 告诉模板按浮点标量读取 |

因此不存在一个单独的：

```cpp
#include "maxCo.H"
```

`maxCo` 不是库类，也不是函数名，而是用户输入字典里的一个字符串键。
真正的代码依赖是 `dictionary.H` 中声明的：

```cpp
template<class T>
T lookupOrDefault(const word&, const T&) const;
```

调用：

```cpp
controlDict.lookupOrDefault<scalar>("maxCo", 0.2)
```

可以理解为：

1. 在 `controlDict` 中查找名为 `"maxCo"` 的条目；
2. 按 `scalar` 类型解析它；
3. 如果找不到，就返回默认值 `0.2`。

之前没有显式写出 `dictionary.H` 时，`dictionary` 可能通过 `Time.H` 等头文件
被间接包含，因此代码有机会编译通过；但这依赖头文件的传递关系。现在显式
写出 `dictionary.H`，代码的直接依赖更清楚，也更适合教学和后续维护。

对应的文件是：

```text
cases/01_advection_equation/01_sine_wave_quad/N20/system/controlDict
```

其中写有：

```foam
maxCo 0.2;
```

可以把这行 C++ 展开成四个问题：

| 代码片段 | 要回答的问题 | 本例答案 |
|---|---|---|
| `controlDict` | 从哪里取参数？ | `system/controlDict` |
| `"maxCo"` | 取哪个名字？ | 目标 Courant 数 |
| `<scalar>` | 按什么类型解析？ | 一个浮点数 |
| `0.2` | 如果没写怎么办？ | 默认使用 `0.2` |

`const` 表示本阶段读出来以后不再改这个目标值。它不是说整个程序永远不能
有别的 CFL 数，而是说这个 C++ 变量在当前作用域内只读。

数学上：

$$
\mathrm{Co}_{\mathrm{target}}=0.2.
$$

这里的 `maxCo` 是一个全局数，所以是 `scalar`，不是 `scalarField`。
它没有“第几个 cell 的 maxCo”，整个案例只设一个目标上限。

### 15.2 第二步：把每个面通量变成非负数

前一阶段已经得到：

```cpp
surfaceScalarField phi;
```

其每个面上的值是带方向的体积通量：

$$
F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_f.
$$

因此：

```cpp
mag(phi)
```

对应：

$$
\lvert F_f\rvert.
$$

这一行仍然是“面数据”。`mag` 只改变数值，不改变数据位置：

```text
surfaceScalarField phi
    -> surfaceScalarField mag(phi)
```

为什么要去掉符号？因为一个 cell 的左、右两个面可能有相反符号。若直接求
带符号和，均匀流场可能得到接近零的净通量；CFL 需要的是所有穿过面的“流动规模”，
所以必须求：

$$
\sum_{f\in\partial\Omega_c}\lvert F_{cf}\rvert.
$$

### 15.3 第三步：从面位置变成 cell 位置

当前代码的完整表达式是：

```cpp
fvc::surfaceSum(mag(phi))
```

它的输入和输出不是同一个位置：

```text
输入：每个 face 一个 |F_f|
输出：每个 cell 一个 sum_f |F_cf|
```

OpenFOAM 14 的声明是：

```cpp
template<class Type>
tmp<VolInternalField<Type>> surfaceSum
(
    const SurfaceField<Type>&
);
```

当 `Type = scalar` 时，可以把它读成：

```cpp
tmp<VolInternalField<scalar>>
```

也就是“一个暂时拥有 cell 内部标量结果的对象”。

它做的不是散度，而是纯加法。对内部面，源码逻辑可以概念性地写成：

```cpp
sum[owner[facei]]    += absFlux[facei];
sum[neighbour[facei]] += absFlux[facei];
```

这里的 `owner` 和 `neighbour` 是网格拓扑中的两个相邻 cell 编号。由于输入已经
是非负的 `absFlux`，同一内部面要计入它两侧 cell 的总和。

边界面则通过该 patch 的 `faceCells()` 加到它所属的 cell。

注意：`surfaceSum` 的结果还没有除以体积，因此它代表：

$$
\sum_f\lvert F_{cf}\rvert,
$$

而不是：

$$
\frac{1}{V_c}\sum_f\lvert F_{cf}\rvert.
$$

### 15.4 第四步：理解后面的 `()`

源码写成：

```cpp
fvc::surfaceSum(mag(phi))()
```

这里确实出现了两个括号，但它们是两种不同的 C++ 语法：

```text
fvc::surfaceSum(mag(phi))
    -> 调用 surfaceSum 函数
    -> 得到 tmp<VolInternalField<scalar>>

()
    -> 调用这个 tmp 对象的 operator()()
    -> 取得它包裹的实际场
```

OpenFOAM 14 的 `tmp` 实现中，常量版本的核心含义是：

```cpp
const T& tmp<T>::operator()() const
{
    return *ptr_;
}
```

因此第二个 `()`：

- 不会再次计算 `surfaceSum`；
- 不会访问某个 cell；
- 不表示“调用一个带参数的数学函数”；
- 只是把临时包装打开，得到其中的场对象。

把它类比成一个盒子：

```text
tmp<VolInternalField<scalar>>
    -> 盒子
operator()()
    -> 打开盒子
VolInternalField<scalar>
    -> 盒子里的实际 cell 数据
```

### 15.5 第五步：理解 `.primitiveField()`

继续看：

```cpp
fvc::surfaceSum(mag(phi))().primitiveField()
```

`VolInternalField<scalar>` 是 OpenFOAM 对“体场内部部分”的抽象。调用
`.primitiveField()` 后，取出它底层的内部数值数组：

```cpp
const scalarField&
```

它的第 `celli` 个元素满足：

$$
\texttt{sumPhi[celli]}
=
\sum_{f\in\partial\Omega_{\texttt{celli}}}\lvert F_{cf}\rvert.
$$

这里有一个容易误解的地方：`.primitiveField()` 不是重新遍历网格，也不是删掉
边界贡献。边界贡献在 `surfaceSum` 内部已经加到相应 cell；此处只是把最终的
内部 cell 数值取出来。

为什么这里只取内部数组？因为 CFL 的公式需要的是每个控制体的：

$$
V_c,\qquad \sum_f\lvert F_{cf}\rvert.
$$

它不需要一个独立的边界 patch 场。边界面已经属于某个控制体的边界，已经参与
了前面的求和。

### 15.6 第六步：构造 `sumPhi`

源码是：

```cpp
scalarField sumPhi
(
    fvc::surfaceSum(mag(phi))().primitiveField()
);
```

这一步把前面得到的内部数组复制或构造成一个独立的 `scalarField` 变量。

`scalarField` 可以暂时理解为：

```text
长度等于 cell 数量的一维浮点数组
```

在当前网格上：

```text
sumPhi.size() = mesh.nCells() = 400
```

所以：

```cpp
sumPhi[0], sumPhi[1], ..., sumPhi[399]
```

分别对应 400 个控制体的面通量绝对值总和。

这一步的数据位置链为：

```text
face values
    -> surfaceSum
    -> cell values
    -> scalarField sumPhi
```

### 15.7 第七步：除以每个 cell 的体积

源码是：

```cpp
scalarField rate
(
    sumPhi/mesh.V().primitiveField()
);
```

数学上对每个 cell 同时做：

$$
r_c
=
\frac{\sum_f\lvert F_{cf}\rvert}{V_c}.
$$

这里不是把 `sumPhi` 除以一个统一的总面积，而是逐 cell 除以该 cell 自己的
体积。因此这是一个逐元素运算：

```text
rate[0] = sumPhi[0] / V[0]
rate[1] = sumPhi[1] / V[1]
...
```

`mesh.V()` 返回网格中每个 cell 的体积。当前网格是：

$$
\Delta x=\Delta y=0.05,\qquad \Delta z=0.1,
$$

因此：

$$
V_c=0.05\times0.05\times0.1=0.00025.
$$

在 OpenFOAM 14 的 `fvMesh.H` 中，`V()` 的返回类型是：

```cpp
const DimensionedField<scalar, fvMesh>& V() const;
```

这句话可以拆成：

| 部分 | 含义 |
|---|---|
| `DimensionedField` | 带有量纲信息的场 |
| `scalar` | 每个 cell 的值是一个实数 |
| `fvMesh` | 这个场和有限体积网格的 cell 关联 |
| `const ...&` | 只读引用，不复制整个体积场 |

因此 `mesh.V()` 比普通 `scalarField` 更丰富：它不仅有一组体积数值，还知道
这些数值的量纲是体积。当前代码继续调用：

```cpp
mesh.V().primitiveField()
```

是把这个带场包装的对象取成底层 `scalarField`，用于和同样已经取出底层数值的
`sumPhi` 做逐元素除法。也就是说，当前 Stage 2 的计算链是：

```text
带量纲的体积场 mesh.V()
    -> primitiveField()
    -> 纯数值数组 V_c
```

这一步适合教学上观察数组和公式的对应关系，但也意味着这段除法是在底层数值
数组上进行的，量纲检查不再由这两个 `scalarField` 自动表达。后续如果把时间推进
写成 OpenFOAM 的完整场表达式，可以再使用带量纲的 `dimensionedScalar`，让
`deltaT * R` 的量纲关系显式保留下来。

面通量的量纲是：

$$
[F_f]
=[\boldsymbol U_f][\boldsymbol S_f]
=\frac{L}{T}L^2
=\frac{L^3}{T}.
$$

于是：

$$
[r_c]
=\frac{L^3/T}{L^3}
=T^{-1}.
$$

所以 `rate` 不是速度，也不是通量；它是一个“单位时间内相对于控制体体积的
通量规模”，量纲为 `1/s`。

### 15.8 第八步：取所有 cell 中最大的 rate

源码是：

```cpp
const scalar rateMax = gMax(rate);
```

前面 `rate` 有 400 个数，现在 `gMax` 把它们压缩成一个数：

$$
r_{\max}=\max_{c=1,\ldots,N_c}r_c.
$$

为什么必须取最大值？因为显式稳定性必须对最不利的 cell 负责。只使用平均
值或某个 cell 的值，不能保证整个网格满足 CFL 约束。

`gMax` 中的 `g` 是 global 的意思。串行时它是整个当前网格的最大值；并行时还
需要在进程之间做归约，最终得到全局最大值。

类型因此发生变化：

```text
scalarField rate
    -> gMax(rate)
    -> scalar rateMax
```

### 15.9 第九步：由目标 CFL 反解时间步

源码是：

```cpp
const scalar deltaT = 2.0*maxCo/rateMax;
```

它直接对应代数变形：

$$
\mathrm{Co}_{\max}
=\frac{\Delta t}{2}r_{\max}
$$

因此：

$$
\Delta t
=\frac{2\,\mathrm{Co}_{\mathrm{target}}}{r_{\max}}.
$$

当前算例中：

$$
\Delta t
=\frac{2\times0.2}{80}
=0.005\ \mathrm{s}.
$$

到这里，程序只知道“下一步最多建议用多长时间”。它还没有把这个数字交给
`runTime`，所以：

- 当前物理时间仍然是 `0`；
- `T` 的内部值仍然是初始正弦波；
- 没有产生 `0.005` 时间目录；
- 没有执行一次 Euler 更新。

## 16. 为什么现在不调用 `runTime.setDeltaT`

你可能会问：既然已经算出了 `deltaT`，为什么不立刻写：

```cpp
runTime.setDeltaT(deltaT);
```

原因是当前教学阶段刻意把三个概念分开：

| 概念 | 当前 Stage 2 是否完成 |
|---|---|
| 计算一个稳定候选时间步 | 是 |
| 把它设置为时间管理器的下一步 | 否 |
| 用它推进 PDE | 否 |

如果现在调用 `setDeltaT`，它只会改变 `runTime` 保存的时间步参数；只要不调用
`runTime++`，物理时间仍不前进。但从教学上看，当前先验证 CFL 公式更清楚，
避免把“算出时间步”和“完成时间推进”混成一件事。

进入 Stage 5 后，才会形成：

```cpp
runTime.setDeltaT(deltaT);
runTime++;
```

并且在每一个时间步中执行：

```text
旧 T
    -> fvc::div(phi,T)
    -> 残差 R^n
    -> T = T - deltaT*R^n
    -> 新 T
```

## 17. 两个错误检查：为什么要主动终止

### 17.1 `maxCo <= 0`

代码：

```cpp
if (maxCo <= 0)
{
    FatalErrorInFunction
        << "maxCo must be positive, but received " << maxCo
        << exit(FatalError);
}
```

如果目标 CFL 是零或负数，得到的 `deltaT` 将不是一个有效的正时间步。
`FatalErrorInFunction` 会把错误定位到当前函数，`exit(FatalError)` 终止程序。

这体现了数值程序的基本原则：

```text
用户参数不满足数学前提
    -> 立即报告
    -> 不继续生成看似正常的结果
```

### 17.2 `rateMax <= SMALL`

代码：

```cpp
if (rateMax <= SMALL)
{
    FatalErrorInFunction
        << "The maximum flux rate is zero or too small: "
        << rateMax << nl
        << "Cannot determine a positive CFL time step."
        << exit(FatalError);
}
```

公式的分母要求：

$$
r_{\max}>0.
$$

如果速度场为零，可能有：

$$
\phi=0,\qquad r_{\max}=0.
$$

此时直接计算 `2.0*maxCo/rateMax` 会除零。`SMALL` 是 OpenFOAM 提供的小正阈值，
用来把“数值上已经接近零”的分母也视为不安全。

当前代码选择报错，而不是擅自把静止场的时间步设成无穷大。以后如果要支持
静止场，需要单独设计策略，例如使用用户给定的 `deltaT`，不能在这里默默掩盖。

## 18. 本阶段的总数据流

把所有对象放在一张图中：

```text
system/controlDict
    |
    | lookupOrDefault<scalar>("maxCo", 0.2)
    v
scalar maxCo
    |
    | target Co
    |
volVectorField U + fvMesh mesh
    |
    | fvc::flux(U)
    v
surfaceScalarField phi
    |
    | mag(phi)
    v
surfaceScalarField |phi|
    |
    | fvc::surfaceSum(...)
    v
tmp<VolInternalField<scalar>>
    |
    | () + primitiveField()
    v
scalarField sumPhi
    |
    | / mesh.V().primitiveField()
    v
scalarField rate
    |
    | gMax(rate)
    v
scalar rateMax
    |
    | 2.0*maxCo/rateMax
    v
scalar deltaT
```

这条链可以用一句话概括：

```text
速度和几何 -> 面通量 -> 每个 cell 的通量规模 -> 最坏 cell -> 全局时间步
```

它还没有进入：

```text
面通量 + T -> 对流残差 -> Euler 更新 -> 新时间层
```

## 19. 读懂 Stage 2 后的自测题

在编译运行之前，你应该能口头回答：

1. `maxCo` 为什么是一个 `scalar`，而 `rate` 是一个 `scalarField`？
2. `mag(phi)` 为什么仍然在 face 上？
3. `surfaceSum` 为什么把结果放到 cell 上？
4. `()` 是重新计算，还是打开 `tmp`？
5. `.primitiveField()` 为什么不等于“删除边界贡献”？
6. 为什么 `rate` 要除以 `mesh.V()`？
7. 为什么使用 `gMax(rate)` 而不是 `average(rate)`？
8. 为什么当前已经得到 `deltaT`，却还没有 `T^{n+1}`？

如果这八个问题都能用“位置、类型、公式、程序阶段”解释清楚，Stage 2
就已经真正理解了，而不是只记住一行 API。

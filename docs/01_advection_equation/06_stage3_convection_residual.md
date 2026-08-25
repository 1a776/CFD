# Stage 3：显式对流残差

本文档讲解如何从 Stage 1 的面通量 `phi` 和 Stage 0 的单元标量场 `T`
计算守恒型线性对流方程的体积归一化对流残差。这个残差就是
OpenFOAM 的 `fvc::div(phi,T,...)` 返回的离散散度场。

本阶段只计算残差，不更新 `T`，不推进 `runTime`。

## 1. 本阶段要实现的公式

第 1 题的守恒型线性对流方程是：

$$\frac{\partial T}{\partial t}+\nabla\cdot(\boldsymbol{U}T)=0.$$

对控制体 $\Omega_c$ 积分，并使用高斯散度定理：

$$V_c\frac{\mathrm dT_c}{\mathrm dt}+\sum_{f\in\partial\Omega_c}F_{cf}T_f=0.$$

其中：

$$F_{cf}=\boldsymbol{U}_f\cdot\boldsymbol{S}_{cf}.$$

本文统一采用更贴近 OpenFOAM 的残差定义：残差已经除以控制体体积，
也就是单元上的体积归一化离散散度：

$$R_c^n=\frac{1}{V_c}\sum_{f\in\partial\Omega_c}F_{cf}^nT_f^n.$$

如果手动写 owner/neighbour 面循环，可以先把未除体积的面通量和记为中间量：

$$\mathcal Q_c^n=\sum_{f\in\partial\Omega_c}F_{cf}^nT_f^n.$$

然后再定义：

$$R_c^n=\frac{\mathcal Q_c^n}{V_c}.$$

注意：$\mathcal Q_c^n$ 只是装配缓存，不叫最终残差。最终残差一律指
$R_c^n$。

因此 PDE 半离散形式是：

$$\frac{\mathrm dT_c}{\mathrm dt}=-R_c.$$

注意符号：Stage 3 只计算空间项 $R_c^n$，此时还没有执行时间更新：

$$T_c^{n+1}=T_c^n-\Delta tR_c^n.$$

减号属于 Stage 4 的前向 Euler 更新，不属于 Stage 3。

## 2. 本阶段新增的头文件

源码新增：

```cpp
#include "fvcDiv.H"
```

它声明 OpenFOAM 的显式散度接口：

```cpp
fvc::div(...)
```

OpenFOAM 14 本地相关文件是：

```text
/opt/openfoam14/src/finiteVolume/finiteVolume/fvc/fvcDiv.H
/opt/openfoam14/src/finiteVolume/finiteVolume/fvc/fvcDiv.C
```

这条调用链里，`fvcDiv.C` 先根据 `fvSchemes` 选离散格式，再把工作交给具体的
convection scheme。对于本题的 `Gauss upwind`，会继续走到：

```text
/opt/openfoam14/src/finiteVolume/finiteVolume/convectionSchemes/gaussConvectionScheme/gaussConvectionScheme.C
/opt/openfoam14/src/finiteVolume/finiteVolume/fvc/fvcSurfaceIntegrate.C
```

其中 `gaussConvectionScheme.C` 负责把面通量和面值组装成单元结果，真正的除体积
动作则在 `fvcSurfaceIntegrate.C` 末尾完成：

```cpp
ivf /= mesh.Vsc();
```

`fvcDiv.H` 中和本阶段直接相关的接口形式是：

```cpp
template<class Type>
tmp<VolField<Type>> div
(
    const surfaceScalarField&,
    const VolField<Type>&,
    const word& name
);
```

当 `Type = scalar` 时，它可以理解为：

```cpp
surfaceScalarField + volScalarField + word
    -> tmp<volScalarField>
```

也就是：

```text
面通量 phi + 单元标量 T + 离散格式名字
    -> 单元标量残差 R
```

## 3. 本阶段写入的核心代码

源码中新增：

```cpp
tmp<volScalarField> tResidual
(
    fvc::div(phi, T, "div(phi,T)")
);

const volScalarField& residual = tResidual();

const scalar residualIntegral
(
    gSum(residual.primitiveField()*mesh.V().primitiveField())
);

Info<< "  residual dimensions = " << residual.dimensions() << nl
    << "  residual min        = " << min(residual).value() << nl
    << "  residual max        = " << max(residual).value() << nl
    << "  residual integral   = " << residualIntegral << nl
    << endl;
```

这段代码没有修改 `T`。它只是把旧时间层的 `T` 放进对流算子，计算一个临时的
体积归一化残差场 $R$。

### 3.1 如果你要写未除体积的原始版本

如果你不是直接调用 `fvc::div`，而是想亲手做面循环，那么代码更像这样：

```cpp
scalarField fluxSum(mesh.nCells(), 0.0);

// internal faces
forAll(owner, facei)
{
    const label own = owner[facei];
    const label nei = neighbour[facei];
    const scalar F = phi[facei];

    const scalar Tf = (F >= 0 ? T[own] : T[nei]);
    const scalar H  = F * Tf;

    fluxSum[own] += H;
    fluxSum[nei] -= H;
}
```

这时 `fluxSum` 对应的是未除体积的原始面通量和：

$$
\mathcal Q_c^n=\sum_{f\in\partial\Omega_c}F_{cf}^nT_f^n.
$$

如果要得到最终残差，仍然要再除一次体积：

$$
R_c^n=\frac{\mathcal Q_c^n}{V_c}.
$$

所以两种写法的区别是：

| 写法 | 结果 | 是否已除体积 |
|---|---|---|
| `fvc::div(phi, T, "div(phi,T)")` | `residual` | 是 |
| 手写 `fluxSum` | `\mathcal Q_c` | 否 |

如果你想用 `fluxSum` 直接更新，那么要写成：

$$
T_c^{n+1}=T_c^n-\frac{\Delta t}{V_c}\mathcal Q_c^n.
$$

如果你先做归一化，再更新，就回到当前主文档的写法：

$$
T_c^{n+1}=T_c^n-\Delta tR_c^n.
$$

## 4. `fvc::div(phi,T,"div(phi,T)")` 逐项解释

### 4.1 `fvc`

`fvc` 是 finite volume calculus 的缩写。它表示显式有限体积计算：

```text
输入已有字段
    -> 直接计算出一个新字段
```

本阶段需要旧时间层空间项，所以用 `fvc`。如果调用 `fvc::div`，
得到的是已经除以体积的 $R^n$。

不要写成：

```cpp
fvm::div(phi, T)
```

`fvm` 会组装矩阵，属于隐式离散路线，不符合本题“显式方法”的训练目标。

### 4.2 第一个参数 `phi`

`phi` 的类型是：

```cpp
surfaceScalarField
```

数学上是每个面上的体积通量：

$$F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_f.$$

它告诉 OpenFOAM 两件事：

1. 每个面有多少流量穿过；
2. 流动方向从面通量的符号判断。

如果 $F_f>0$，通量沿 owner 到 neighbour 的方向；如果 $F_f<0$，方向相反。
这正是迎风格式判断上游单元所需要的信息。

### 4.3 第二个参数 `T`

`T` 的类型是：

```cpp
volScalarField
```

数学上是单元中心的标量：

$$\{T_c^n\}_{c=1}^{N_c}.$$

`fvc::div` 需要根据 `fvSchemes` 中的插值格式把 cell 中心的 `T_c`
变成面上的 $T_f$。

本算例中选择：

```foam
div(phi,T) Gauss upwind;
```

所以面值是上风值：

$$
T_f=
\begin{cases}
T_{\mathrm{owner}},&F_f\geq0,\\
T_{\mathrm{neighbour}},&F_f<0.
\end{cases}
$$

### 4.4 第三个参数 `"div(phi,T)"`

这个字符串不是随便写的注释，它是 OpenFOAM 查找离散格式的名字。

源码中写：

```cpp
fvc::div(phi, T, "div(phi,T)")
```

字典中必须有匹配项：

```foam
divSchemes
{
    default         none;
    div(phi,T)      Gauss upwind;
}
```

对应关系是：

```text
C++ 字符串 "div(phi,T)"
    -> fvSchemes 中的 div(phi,T)
        -> Gauss upwind
```

如果 C++ 写成 `"div(phi, T)"`，多了一个空格，就可能找不到字典项。
由于 `default none;`，找不到时 OpenFOAM 会报错，而不是悄悄使用默认格式。

这对教学很好：拼写错会暴露出来。

## 5. `Gauss upwind` 到底做了什么

`Gauss` 指的是用高斯散度定理把体积分散度变成面积分通量：

$$\int_{\Omega_c}\nabla\cdot(\boldsymbol{U}T)\,\mathrm d\Omega
=
\int_{\partial\Omega_c}\boldsymbol{U}T\cdot\boldsymbol{n}\,\mathrm dS.$$

离散后就是：

$$\sum_fF_{cf}T_f.$$

`upwind` 指的是面值 $T_f$ 的选择方式。由于 `T` 原本在 cell 中心，
而通量需要面上的值，所以必须进行面插值。

一阶迎风格式的思想是：

```text
流从哪里来，就使用哪里的 T
```

这比中心插值更耗散，但更稳健，是初学显式对流方程时最合适的第一版。

## 6. 返回类型为什么是 `tmp<volScalarField>`

调用：

```cpp
fvc::div(phi, T, "div(phi,T)")
```

返回：

```cpp
tmp<volScalarField>
```

分开看：

| 部分 | 含义 |
|---|---|
| `tmp<...>` | OpenFOAM 管理的临时对象 |
| `vol` | 结果定义在 cell 上 |
| `Scalar` | 每个 cell 一个标量 |
| `Field` | 一整个场 |

为什么结果回到 cell 上？因为散度是控制体上的体积平均量：

$$R_c=\frac{1}{V_c}\sum_fF_{cf}T_f.$$

每个 cell 得到一个 $R_c$，所以它是 `volScalarField`。

## 7. 为什么写 `tResidual` 和 `residual`

代码先保存临时对象：

```cpp
tmp<volScalarField> tResidual
(
    fvc::div(phi, T, "div(phi,T)")
);
```

这里的 `t` 可以读作 temporary。它表示：

```text
这个场是中间结果，用完可以释放
```

随后：

```cpp
const volScalarField& residual = tResidual();
```

这里的 `()` 与 Stage 2 中 `tmp` 的解释相同：它打开 `tmp`，取得内部实际场。

为什么不直接写成一行？因为教学阶段需要多次访问这个空间散度场：

- 打印量纲；
- 打印最小值；
- 打印最大值；
- 计算体积分。

保存成 `residual` 引用后，后面的代码更容易读，并且它和公式中的
$R_c$ 一一对应。

## 8. 残差的量纲

`T` 是无量纲：

$$[T]=1.$$

`phi` 是体积通量：

$$[\phi]=\frac{L^3}{T_{\mathrm{time}}}.$$

散度计算包含除以 cell 体积：

$$[R]=\frac{L^3/T_{\mathrm{time}}}{L^3}=\frac{1}{T_{\mathrm{time}}}.$$

而未归一化的面通量和 $\mathcal Q_c$ 的量纲是：

$$[\mathcal Q]=\frac{L^3}{T_{\mathrm{time}}}.$$

因此日志中的：

```text
residual dimensions
```

应该显示时间倒数的量纲。由于 OpenFOAM 的时间维指数为 `-1`，
你应预期它类似：

```text
[0 0 -1 0 0 0 0]
```

具体输出以 OpenFOAM 14 的格式为准。

## 9. 为什么计算 `residualIntegral`

源码中：

```cpp
const scalar residualIntegral
(
    gSum(residual.primitiveField()*mesh.V().primitiveField())
);
```

数学上对应体积归一化残差的体积分：

$$\sum_cR_cV_c=\sum_c\mathcal Q_c.$$

把 $\mathcal Q_c=\sum_fF_{cf}T_f$ 代入：

$$\sum_cR_cV_c=\sum_c\mathcal Q_c=\sum_c\sum_fF_{cf}T_f.$$

在周期边界、守恒通量和封闭拓扑下，内部面贡献应该成对抵消，总体积分应接近零。

因此它是一个很好的 Stage 3 检查量：

```text
residual integral approximately 0
```

这里的 `gSum` 和 Stage 2 的 `gMax` 类似，都是全局归约：

- 串行时，对当前进程所有 cell 求和；
- 并行时，还会跨进程求和。

## 10. Stage 2 和 Stage 3 的关键区别

这两个阶段都使用了面通量 `phi`，但目的完全不同。

| 阶段 | 代码 | 数学意义 | 是否带方向 | 用途 |
|---|---|---|---|---|
| Stage 2 | `fvc::surfaceSum(mag(phi))` | $\sum_f|F_{cf}|$ | 不带方向 | CFL 稳定性 |
| Stage 3 中间装配 | 手动累加面通量 | $\mathcal Q_c=\sum_fF_{cf}T_f$ | 带方向 | 面通量守恒装配 |
| Stage 3 残差 | `fvc::div(phi,T,"div(phi,T)")` | $R_c=\mathcal Q_c/V_c$ | 带方向 | PDE 空间散度 |

Stage 2 用绝对值，因为它关心总穿越规模。

Stage 3 不能用绝对值，因为 PDE 的守恒更新必须知道流入和流出方向。

## 11. 本阶段完成后的状态

如果把程序暂时截断在 Stage 3，程序会知道：

```text
mesh
U
T
phi
deltaT
R = volume-normalized convection residual
```

在 Stage 3 的独立教学检查中，程序还没有执行：

```cpp
T = T - deltaT*R;
```

也没有执行：

```cpp
runTime++;
```

所以在“只检查 Stage 3”的临时版本中：

- `T` 不应该改变；
- 不应该生成新的时间目录；
- 日志中应该出现残差范围；
- 日志中仍应出现 `No time update has been implemented yet.`。

## 12. Stage 3 独立检查与当前版本说明

下面的检查只适用于你把程序临时截断在 Stage 3、并在 Stage 3 后提前
`return 0` 的教学版本：

```bash
grep -E "residual dimensions|residual min|residual max|residual integral|Stage 3|No time update" \
cases/01_advection_equation/01_sine_wave_quad/N20/log.explicitAdvectionFoamStudent
```

你应重点看：

1. `residual dimensions` 是否是时间倒数量纲；
2. `residual min` 和 `residual max` 是否是有限数，不是 `nan`；
3. `residual integral` 是否接近 `0`；
4. 是否成功打印了残差的四行诊断信息。

当前学生版已经把 Stage 4 接在这段代码后面。当前完整运行时，Stage 3
的残差会继续进入前向 Euler 更新。Stage 4 的更新公式是：

$$T_c^{n+1}=T_c^n-\Delta tR_c^n.$$

因此当前版本请按照
`docs/07_stage4_forward_euler_update.md` 中的命令和验收标准检查，
不要再用本节的 `No time update` 作为当前完整程序的预期日志。

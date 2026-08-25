# Stage 1：从速度场得到面通量

本文档讲解第 1 题实现中的第一个真正数值离散步骤：

$$F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_f.$$

这一阶段只计算面体积通量，不计算 CFL，不计算对流散度，也不更新标量场 `T`。

## 1. 这一阶段在整个求解器中的位置

守恒型线性对流方程为：

$$\frac{\partial T}{\partial t}
+\nabla\cdot(\boldsymbol{U}T)=0.$$

对控制体 $\Omega_c$ 积分后，得到：

$$V_c\frac{\mathrm dT_c}{\mathrm dt}
+\sum_{f\in\partial\Omega_c}F_{cf}T_f=0.$$

这里最先需要知道的是每个面上的：

$$F_{cf}=\boldsymbol{U}_f\cdot\boldsymbol{S}_{cf}.$$

因此求解器的数值对象依次是：

```text
单元速度 U
    -> 面速度 U_f
        -> 面面积向量 S_f
            -> 面通量 F_f
                -> 面上风向和 CFL
                    -> 对流散度
                        -> 时间更新
```

Stage 1 只完成箭头中“面通量 $F_f$”这一格。

## 2. 为什么不能直接把 `U` 当作通量

当前代码中的速度场是：

```cpp
volVectorField U;
```

它表示：

```text
每个 cell 有一个 vector 速度
```

数学上可以记为：

$$\{\boldsymbol{U}_c\}_{c=1}^{N_c}.$$

但是有限体积法的守恒公式需要的是每个控制体边界上的通量：

$$F_{cf}=\boldsymbol{U}_f\cdot\boldsymbol{S}_{cf}.$$

这两个量有三个区别：

| 对象 | 定义位置 | 值的类型 | 数学含义 |
|---|---|---|---|
| `U` | cell | vector | $\boldsymbol{U}_c$ |
| $\boldsymbol{U}_f$ | face | vector | 面上的速度 |
| `phi` | face | scalar | $F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_f$ |

所以 `U` 不能直接替代 `phi`：

```text
U   是“速度”
phi 是“穿过面的体积流量”
```

如果速度量纲是 $L/T$，面积向量量纲是 $L^2$，则：

$$[F_f]=[U][S_f]=\frac{L}{T}L^2=\frac{L^3}{T}.$$

这说明 `phi` 的单位是体积流量，而不是速度。

## 3. 面面积向量是什么

有限体积法中，一个面不仅有面积大小，还必须有方向。面面积向量写作：

$$\boldsymbol{S}_f=A_f\boldsymbol{n}_f,$$

其中：

- $A_f$ 是面面积；
- $\boldsymbol{n}_f$ 是单位法向量；
- $\boldsymbol{S}_f$ 的方向由 OpenFOAM 的 owner/neighbour 约定确定。

点积：

$$F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_f$$

实际上可以写成：

$$F_f=|\boldsymbol{U}_f|A_f\cos\theta.$$

因此：

- $F_f>0$：速度沿着当前面法向方向；
- $F_f<0$：速度沿着面法向的反方向；
- $F_f=0$：速度与面平行，没有体积流量穿过该面。

这个正负号非常重要。下一阶段判断上风值时，正是通过这个符号判断流体从哪一侧进入控制体。

### 3.1 内部面为什么能保证守恒

对于相邻的两个控制体，同一个内部面只被存储一次。若从 owner 单元看，该面的面积向量是：

$$\boldsymbol{S}_{f,\mathrm{owner}},$$

那么从 neighbour 单元看，同一个面的方向相反：

$$\boldsymbol{S}_{f,\mathrm{neighbour}}
=-\boldsymbol{S}_{f,\mathrm{owner}}.$$

因此同一个内部面的通量在两个单元的守恒方程中符号相反：

$$F_{f,\mathrm{neighbour}}=-F_{f,\mathrm{owner}}.$$

一个单元流出的量，正好是另一个单元流入的量。这就是有限体积法局部守恒的几何基础。

## 4. `surfaceScalarField` 的含义

Stage 1 中需要创建：

```cpp
surfaceScalarField phi;
```

这个类型可以拆成三个部分：

```text
surface + scalar + field
```

含义是：

```text
定义在所有 face 上的标量场
```

数学上相当于：

$$\{F_f\}_{f=1}^{N_f}.$$

这里的 `surface` 在 OpenFOAM 有限体积语境中主要指网格面，不是连续几何学里任意曲面。

类型对比：

| C++ 类型 | 每个什么对象存一个值 | 值的形式 | 本题中的对象 |
|---|---|---|---|
| `volScalarField` | cell | scalar | `T_c` |
| `volVectorField` | cell | vector | `U_c` |
| `surfaceScalarField` | face | scalar | `F_f` |
| `surfaceVectorField` | face | vector | 概念上的 `U_f` |

一个非常常见的错误是把通量写成：

```cpp
volScalarField phi;
```

这在位置上就是错误的。通量发生在面上，而不是发生在单元中心。

## 5. `fvc::flux(U)` 在数学上做了什么

OpenFOAM 14 中，相关接口声明在：

```text
/opt/openfoam14/src/finiteVolume/finiteVolume/fvc/fvcFlux.H
```

源码中接口的核心形式是：

```cpp
template<class Type>
tmp<SurfaceField<typename innerProduct<vector, Type>::type>> flux
(
    const VolField<Type>& vf
);
```

对 `volVectorField U` 使用时，可以把它理解成：

```text
fvc::flux(U)
    -> 在面上得到 U_f
    -> 与 mesh.Sf() 点乘
    -> 返回每个面一个 scalar 的临时场
```

数学对应关系是：

$$\texttt{fvc::flux(U)}
\longleftrightarrow
\{\boldsymbol{U}_f\cdot\boldsymbol{S}_f\}_f
\longleftrightarrow
\{F_f\}_f.$$

在 OpenFOAM 的实现中，`fvc` 表示 finite-volume calculus 的显式计算接口。这里的“显式”是指它直接计算场值或离散算子结果，并不组装线性方程矩阵。

注意，`fvc::flux(U)` 与后面用于输运标量的面通量插值不是同一件事：

```text
fvc::flux(U)
    计算速度穿过面产生的体积通量 F_f

fvc::div(phi, T, "div(phi,T)")
    根据 phi 和 fvSchemes 计算 phi*T 的离散散度
```

前者处理速度和几何，后者处理被输运标量的面值和散度。

## 6. 为什么返回类型前面有 `tmp<...>`

`fvc::flux(U)` 的返回类型不是直接写成：

```cpp
surfaceScalarField
```

而是类似：

```cpp
tmp<SurfaceField<scalar>>
```

`surfaceScalarField` 本身可以理解为：

```cpp
SurfaceField<scalar>
```

而 `tmp<T>` 是 OpenFOAM 用来管理临时场对象的包装器。

原因是表达式：

```cpp
fvc::flux(U)
```

通常只需要把计算结果交给另一个场对象使用一次。OpenFOAM 用 `tmp` 来减少不必要的复制，并在临时对象不再需要时自动释放内存。

对于初学者，先记住这条规则即可：

```text
fvc::flux(U) 产生临时结果；
surfaceScalarField phi(...) 接管或保存这个结果。
```

你不需要手动 `new`，也不需要手动 `delete`。

## 7. `phi` 的 IOobject 为什么使用 `NO_READ`

`phi` 不是题目给定的初始场，而是由 `U` 和网格几何现场计算出来的：

$$\phi=\phi(U,\mathrm{mesh}).$$

因此它不需要从 `0/phi` 读取。构造时应使用：

```cpp
IOobject::NO_READ
```

这表示：

```text
不要要求磁盘上必须存在 phi 文件
```

通常还可以使用：

```cpp
IOobject::NO_WRITE
```

因为当前阶段只想检查通量，不需要把它输出成时间场。

因此 Stage 1 的 `IOobject` 语义是：

| 设置 | 含义 |
|---|---|
| `"phi"` | 这个场在 OpenFOAM 中的名字 |
| `runTime.name()` | 当前时间目录 |
| `mesh` | 这个场属于哪个网格 |
| `NO_READ` | 不从磁盘读取 |
| `NO_WRITE` | 暂不写到磁盘 |

以后如果希望在 ParaView 中检查通量，可以把最后一个选项改为 `AUTO_WRITE`，并在适当时机调用 `runTime.write()`。但这不是本阶段的重点。

## 8. 你要亲自写的最小代码

你需要在源码顶部增加：

```cpp
#include "surfaceFields.H"
#include "fvcFlux.H"
```

然后在 Stage 1 TODO 区域创建 `phi`。结构应该遵循 OpenFOAM 原生代码的模式：

```cpp
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::flux(U)
);
```

这段代码不要机械抄写。你需要逐项解释：

1. 为什么类型是 `surfaceScalarField`；
2. 为什么名字是 `"phi"`；
3. 为什么是 `NO_READ`；
4. 为什么 `fvc::flux(U)` 的结果能作为初始值；
5. 为什么没有把 `T` 放进 `fvc::flux`。

然后增加输出：

```cpp
Info<< "  phi min     = " << min(phi).value() << nl
    << "  phi max     = " << max(phi).value() << nl
    << endl;
```

这里：

```text
min(phi) -> 所有面通量中的最小值
max(phi) -> 所有面通量中的最大值
```

它们是全局标量，不再是每个 cell 或每个 face 的场。

## 9. 当前算例的预期现象

当前算例使用：

$$\boldsymbol{U}=(1,1,0).$$

网格是二维区域通过一个很薄的第三方向层表示的三维网格。当前坐标方向的面面积大约是：

$$A_f=0.05\times0.1=0.005.$$

因此：

- 垂直于 $x$ 方向的面，通量约为 $1\times0.005=0.005$；
- 垂直于 $y$ 方向的面，通量约为 $1\times0.005=0.005$；
- 垂直于 $z$ 方向的面，速度分量为零，通量约为 $0$；
- 相反方向的面通量为负值。

所以运行日志中的 `phi min` 应为负数，`phi max` 应为正数。具体小数位和范围以本机实际网格为准。

这里先不要根据一个精确的最小值或最大值判断程序对不对。第一步只检查：

1. 程序成功编译；
2. 程序成功构造 `phi`；
3. `phi` 的最小值和最大值一正一负；
4. 数量级大约是 $10^{-3}$ 到 $10^{-2}$；
5. 程序仍然没有推进时间。

## 10. 常见错误与诊断方法

### 10.1 缺少头文件

如果编译器提示找不到 `surfaceScalarField` 或 `fvc::flux`，优先检查：

```cpp
#include "surfaceFields.H"
#include "fvcFlux.H"
```

### 10.2 把 `phi` 声明成体场

错误形式：

```cpp
volScalarField phi;
```

这会把每个 face 的通量错误地声明成每个 cell 的量。应改为：

```cpp
surfaceScalarField phi;
```

### 10.3 试图读取 `0/phi`

如果 `IOobject` 使用了 `MUST_READ`，OpenFOAM 会要求磁盘上存在 `0/phi`。但本阶段 `phi` 是计算得到的，不应该这样读取。

### 10.4 把 `fvc::flux(U)` 理解成标量输运项

`fvc::flux(U)` 只处理速度和面几何：

$$\boldsymbol{U}_f\cdot\boldsymbol{S}_f.$$

它还没有包含 `T`，所以它不是：

$$\nabla\cdot(\boldsymbol{U}T).$$

### 10.5 现在就写 CFL 或时间循环

本阶段先不做：

```text
deltaT
maxCo
fvc::div
T = ...
runTime++
```

先确认你已经真正理解了“为什么从体速度得到面通量”。

## 11. Stage 1 验收标准

完成后，源码应满足：

- 新增面场和通量计算所需头文件；
- 创建一个名为 `phi` 的 `surfaceScalarField`；
- `phi` 由 `fvc::flux(U)` 初始化；
- 输出 `min(phi)` 和 `max(phi)`；
- 没有更新 `T`；
- 没有修改 `/opt/openfoam14/src`；
- 编译成功；
- 算例运行成功；
- 日志能够显示通量的正负范围。

建议使用：

```bash
source /opt/openfoam14/etc/bashrc
cd /home/a776/workdocuments/上交船舶/slover/student_project
sh scripts/build_student_solver.sh
sh cases/01_advection_equation/01_sine_wave_quad/N20/Allrun
```

完成后，把 `log.explicitAdvectionFoamStudent` 中与 `phi` 有关的几行发回来。下一步我们再解释：

$$\mathrm{Co}_c
=\frac{\Delta t}{2V_c}\sum_f|F_{cf}|,$$

以及为什么 CFL 公式中使用的是通量绝对值。

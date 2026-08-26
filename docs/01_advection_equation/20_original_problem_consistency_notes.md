# 第一题和第二题与原始题目的问题说明

本文集中记录第一题和第二题与原始题目之间需要主动说明的问题。

第一题包含三个实现一致性问题：

1. `100 vertices` 与当前 `100 cellsPerEdge` 的差异；
2. 固体旋转复杂初始轮廓是参数化近似；
3. 当前误差采用体积加权归一化定义。

第二题包含两个原题本身留下的真实歧义：

1. 第一个扩散算例只写了 Neumann 边界，没有指定具体法向通量；
2. 第二个扩散算例的 Gaussian 解中，径向变量 $r$ 的定义与公式存在笔误。

这些问题不等价于求解器核心算法错误，但在正式推导、代码实现和报告中必须明确记录。

原始题目为：

[`pdf/training_examples_incomp.pdf`](../../pdf/training_examples_incomp.pdf)

当前第一题求解器为：

[`UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/`](../../UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/)

## 1. `100 vertices` 与 `100 cellsPerEdge` 的差异

### 1.1 原题写法

固体旋转算例原题写的是：

> The square domain is partitioned by triangular and quadrilateral elements with 100 vertices along each edge.

直译为：正方形每条边方向有 100 个顶点。

### 1.2 当前项目的实际设置

当前旋转案例配置中写的是：

```json
"resolutionMeaning": "cellsPerEdge"
```

四边形网格 `N=100` 在 `blockMeshDict` 中实际生成：

```foam
hex (...) (100 100 1) simpleGrading (1 1 1)
```

这表示：

- x 方向有 100 个单元间隔；
- y 方向有 100 个单元间隔；
- 总单元数为 `100 * 100 = 10000`；
- 如果包含两个端点，则每条边有 101 个网格节点。

三角形网格的 Gmsh 脚本中使用：

```python
gmsh.model.mesh.setTransfiniteCurve(line, resolution + 1)
```

因此三角形案例同样是 `N` 个边方向单元间隔，对应 `N+1` 个边界节点。

对应代码：

- [`04_solid_rotation_quad_upwind.json`](../../scripts/configs/01_advection_equation/04_solid_rotation_quad_upwind.json)
- [`04_solid_rotation_tri_upwind.json`](../../scripts/configs/01_advection_equation/04_solid_rotation_tri_upwind.json)
- [`gmsh_tri_mesh.py`](../../scripts/common/gmsh_tri_mesh.py)

### 1.3 当前解释

本项目将原题中的：

```text
100 vertices along each edge
```

解释为：

```text
100 个边方向网格间隔，即 N=100 cellsPerEdge
```

这是一种常见的数值实验参数化方式，也便于进行 `N=50,100,200` 的网格研究。

但是，严格来说：

```text
100 vertices != 100 cells
```

如果要严格实现“每条边有 100 个顶点”，则通常应使用 99 个单元间隔：

```foam
hex (...) (99 99 1) simpleGrading (1 1 1)
```

或者在 Gmsh 中设置：

```python
setTransfiniteCurve(line, 100)
```

### 1.4 正式报告中的推荐表述

建议正式报告写成：

> 原题要求每条边设置 100 个 vertices。由于本项目采用结构化网格分辨率参数 `N` 表示边方向的单元间隔数，因此将该要求解释为 `N=100 cellsPerEdge`，即每个方向 100 个单元间隔、101 个网格节点。`N=50` 和 `N=200` 作为额外的网格分辨率扩展实验。

这样既说明了原题要求，也没有把当前工程实现错误地描述成严格的 100 个顶点。

## 2. 复杂初始轮廓是参数化近似，不是原图逐点重建

### 2.1 原题给出的内容

原题的固体旋转算例要求初始轮廓包含：

- slotted disk；
- cone；
- smooth hump。

题目通过 Figure 1 给出了初始轮廓示意图，并引用了参考文献，但题面没有给出所有几何参数，例如：

- 每个轮廓的精确中心；
- 圆盘和圆锥的精确半径；
- slot 的精确宽度和高度；
- smooth hump 的精确函数形式和幅值；
- 三个轮廓之间是否存在重叠；
- 图像中每个轮廓的精确数值范围。

因此，只根据原题图片无法唯一恢复一个逐点完全相同的初始场。

### 2.2 当前项目的实现

当前项目在 [`advection_rotation.py`](../../scripts/common/advection_rotation.py) 中采用参数化函数构造初始场。

旋转速度为：

$$
\boldsymbol{u}(x,y)
=
\left(
\omega(y_c-y),
\omega(x-x_c),
0
\right).
$$

当前配置使用：

$$
(x_c,y_c)=(0.5,0.5),
\qquad
\omega=1.
$$

因此：

$$
\boldsymbol{u}(x,y)
=
(0.5-y,\;x-0.5).
$$

配置位置：

[`04_solid_rotation_quad_upwind.json`](../../scripts/configs/01_advection_equation/04_solid_rotation_quad_upwind.json)

初始轮廓由三个函数叠加：

$$
\phi_0(x,y)
=
\phi_{\mathrm{disk}}(x,y)
\+
\phi_{\mathrm{cone}}(x,y)
\+
\phi_{\mathrm{hump}}(x,y).
$$

当前具体实现为：

1. **切口圆盘**

   在圆盘内部取值 1，在 slot 区域取值 0：

   ```python
   slotted_disk(...)
   ```

2. **圆锥**

   对距离中心为 $d$ 的点采用：

   $$
   \phi_{\mathrm{cone}}
   =
   1-\frac{d}{R},
   \qquad d\le R.
   $$

   对 $d>R$ 的点取 0。

3. **光滑 hump**

   当前采用 cosine 型函数：

   $$
   \phi_{\mathrm{hump}}
   =
   \frac14
   \left[
   1+\cos\left(\frac{\pi d}{R}\right)
   \right],
   \qquad d\le R.
   $$

当前参数包括：

```json
"radius": 0.15,
"diskCenter": [0.5, 0.75],
"coneCenter": [0.5, 0.25],
"humpCenter": [0.25, 0.5],
"slotHalfWidth": 0.025,
"slotTopY": 0.85
```

### 2.3 当前实现的性质

当前实现能够完整体现原题想考察的数值现象：

- 切口圆盘考察间断和尖锐界面；
- 圆锥考察斜率变化和非光滑峰值；
- smooth hump 考察光滑函数的输运；
- 一圈旋转后可以观察数值耗散、轮廓变形和峰值衰减。

但是，由于原题没有给出所有初始轮廓的精确参数，当前结果应称为：

```text
基于原题 Figure 1 的参数化复现
```

而不应称为：

```text
与参考文献初始轮廓逐点完全相同
```

### 2.4 正式报告中的推荐表述

建议正式报告写成：

> 原题通过 Figure 1 给出了 slotted disk、cone 和 smooth hump 的初始轮廓示意，但未在题面中完整给出其几何参数和函数表达式。本文依据图示构造参数化的切口圆盘、线性圆锥和 cosine smooth hump，并将其作为第一题固体旋转算例的可复现实验初始场。因此，本文重点比较不同网格和离散格式对间断、斜率变化和光滑峰值的输运能力，而不声称对参考文献轮廓进行了逐点重建。

## 3. 当前误差是体积加权归一化误差

### 3.1 当前代码中的定义

当前后处理函数 [`advection_tools.py`](../../scripts/common/advection_tools.py) 中计算的是体积加权归一化误差。

令数值解为 $\phi_c^n$，解析解为 $\phi_c^e$，单元体积为 $V_c$，则：

$$
E_{L_1}
=
\frac{
\sum_c
\left|\phi_c^n-\phi_c^e\right|V_c
}{
\sum_c
\left|\phi_c^e\right|V_c
}.
$$

$$
E_{L_2}
=
\left(
\frac{
\sum_c
\left(\phi_c^n-\phi_c^e\right)^2V_c
}{
\sum_c
\left(\phi_c^e\right)^2V_c
}
\right)^{1/2}.
$$

$$
E_{L_\infty}
=
\frac{
\max_c\left|\phi_c^n-\phi_c^e\right|
}{
\max_c\left|\phi_c^e\right|
}.
$$

代码中的对应位置：

- 误差函数：[`advection_tools.py:195`](../../scripts/common/advection_tools.py:195)
- 四边形后处理：[`postprocess_case.py`](../../scripts/common/postprocess_case.py)
- 三角形后处理：[`postprocess_case.py`](../../scripts/common/postprocess_case.py)
- 收敛阶汇总：[`study_analysis.py`](../../scripts/common/study_analysis.py)

### 3.2 当前数据中的字段名称

当前结果表中使用：

```text
normalizedL1
normalizedL2
normalizedLinf
```

例如：

```text
data/01_advection_equation/analysis/01_sine_wave_quad/convergence_summary.csv
```

这说明当前数据记录的是归一化误差，而不是未归一化的：

$$
\sum_c
\left|\phi_c^n-\phi_c^e\right|V_c.
$$

### 3.3 为什么必须明确“归一化”

如果只写：

```text
L1 error
```

读者无法判断使用的是：

1. 未归一化的体积积分误差；
2. 除以总单元数的平均误差；
3. 除以精确解 $L_1$ 范数的相对误差；
4. 其他软件定义的离散范数。

当前项目使用的是第 3 种，即相对于解析解范数归一化的体积加权误差。

因此报告中推荐写：

```text
volume-weighted normalized L1 error
```

中文可以写成：

```text
体积加权归一化 L1 误差
```

### 3.4 收敛阶定义

当前网格研究使用相邻网格的误差计算观察收敛阶：

$$
p
=
\frac{
\log(E_N/E_{2N})
}{
\log(h_N/h_{2N})
}.
$$

当网格尺度满足：

$$
h_N\propto\frac{1}{N},
$$

则：

$$
\frac{h_N}{h_{2N}}=2,
$$

所以可以写成：

$$
p
=
\frac{
\log(E_N/E_{2N})
}{
\log 2
}.
$$

这里的 $N$ 表示每条边方向的分辨率参数。对于当前规则四边形网格和规则三角形网格，使用 `N=10,20,40,80` 时，这个网格尺度比值关系成立。

需要注意：三角形网格的单元数约为同一 $N$ 四边形网格的两倍。因此，比较三角形和四边形时，不能只根据相同的 $N$ 判断谁更准确，还应同时考虑：

- 实际单元数量；
- 实际单元面积；
- 实际特征网格尺度；
- 网格拓扑和方向性。

## 4. 第一题三个问题对最终结论的影响

这三个问题不会否定当前第一题求解器的核心实现，但会影响结果的表述边界：

| 问题 | 对当前结果的影响 | 正确表述 |
|---|---|---|
| `100 vertices` 与 `100 cellsPerEdge` | 原题基准网格的严格节点数存在差异 | 当前采用 100 个边方向单元间隔 |
| 复杂轮廓参数未完全给出 | 旋转轮廓不能声称与参考图逐点一致 | 基于题图的参数化复现 |
| 误差采用归一化形式 | 数值表中的误差不是未归一化积分误差 | 体积加权归一化误差 |

## 5. 第二题的两个原题歧义

### 5.1 第一个扩散算例的 Neumann 边界没有指定具体通量

原题第二题第一个算例只写明：

```text
Neumann boundary conditions
```

但没有说明边界条件究竟是齐次 Neumann 条件：

$$
\frac{\partial\phi}{\partial n}=0,
$$

还是由给定解析解直接导出的非齐次 Neumann 条件：

$$
\frac{\partial\phi}{\partial n}
=
g(x,y,t).
$$

其中 $n$ 是边界外法向，$g$ 应由解析解 $\phi_{\mathrm{exact}}$ 求法向导数得到：

$$
g(x,y,t)
=
\nabla\phi_{\mathrm{exact}}(x,y,t)\cdot\boldsymbol n.
$$

这两种选择并不等价：

- 齐次 Neumann 条件表示边界法向梯度为零；
- 解析解导出的非齐次 Neumann 条件可以保证解析解严格满足边界条件；
- 如果解析解在外边界上的法向导数不为零，那么直接使用齐次 Neumann 条件会使数值问题与解析解问题不完全相同。

正式实现时应明确选择。最稳妥、最容易复现的第一阶段约定是：

> 本项目将原题第一个扩散算例解释为齐次 Neumann 边界，即 $\partial\phi/\partial n=0$。

但是，如果目标是严格验证原题给出的解析解，则还应进一步计算解析解的法向导数，并采用非齐次 Neumann 边界：

$$
\frac{\partial\phi}{\partial n}
=
\nabla\phi_{\mathrm{exact}}\cdot\boldsymbol n.
$$

因此，第二题正式报告中必须同时写出：

1. 采用的是齐次还是非齐次 Neumann 条件；
2. 如果采用非齐次条件，边界通量函数 $g$ 的具体表达式；
3. 该边界选择是否与解析解严格一致。

### 5.2 第二个 Gaussian 解的 $r$ 定义存在原题笔误

原题第二个扩散算例给出 Gaussian 型解析解，但公式和后续定义之间存在不一致。

题面中的形式可以整理为：

$$
\phi(x,y,t)
=
\frac{1}{1+4\mu t}
\exp\left(
-\frac{\mu r^2}{1+4\mu t}
\right),
$$

但原题后面又写成：

$$
r=x^2+y^2.
$$

从量纲和标准二维 Gaussian 解的写法看，后一行应解释为：

$$
r^2=x^2+y^2.
$$

因此本文采用的统一解析解为：

$$
\phi(x,y,t)
=
\frac{1}{1+4\mu t}
\exp\left(
-\frac{\mu(x^2+y^2)}{1+4\mu t}
\right).
$$

当题目指定 $\mu=1$ 时，上式化为：

$$
\phi(x,y,t)
=
\frac{1}{1+4t}
\exp\left(
-\frac{x^2+y^2}{1+4t}
\right).
$$

正式答案中应明确声明：

> 本文将原题中的径向变量解释为 $r^2=x^2+y^2$，而不是 $r=x^2+y^2$。

如果将 $r=x^2+y^2$ 再代入 $r^2$，就会错误地得到 $(x^2+y^2)^2$，这不是原题对应的标准 Gaussian 扩散解。

## 6. 第二题实现前必须补充的内容

第二题的有限体积离散和显式求解器设计方向可以沿用第一题的开发框架，但正式实现前必须补充以下约定：

| 项目 | 必须明确的内容 |
|---|---|
| 第一个扩散算例边界 | 采用 $\partial\phi/\partial n=0$，还是解析解导出的非齐次 Neumann 条件 |
| 第二个扩散算例径向变量 | 统一采用 $r^2=x^2+y^2$ |
| 误差定义 | 明确使用体积加权归一化 $L_1$、$L_2$、$L_\infty$ |
| 收敛阶 | 明确使用网格尺度 $h$ 或每边分辨率 $N$ 计算 |
| 网格 | 三角形和四边形均需在多个分辨率下测试 |
| 终止时间 | 两个扩散算例都应在 $t=0.2$ 报告误差和收敛阶 |

## 7. 推荐的总说明

在报告或向老师汇报时，可以统一表述为：

> 本项目完成了第一题二维线性对流方程的显式有限体积求解器，并在三角形和四边形网格上完成了正弦波平移与复杂轮廓固体旋转测试。正弦波算例严格采用题目给定的方程、初值、速度、单位正方形区域、周期边界、CFL=0.2 和 $t=1$ 误差验证；固体旋转算例严格采用题目给定的旋转速度和 $t=2\pi$ 终止时间。对于原题中“100 vertices along each edge”的表述，本文将其工程化解释为 100 个边方向单元间隔；对于题目仅以图示给出的复杂轮廓，本文采用参数化的切口圆盘、圆锥和光滑 hump 进行复现；误差采用体积加权归一化 $L_1$、$L_2$ 和 $L_\infty$ 定义。第二题正式实现时，将明确第一个算例采用的 Neumann 边界形式，并将第二个 Gaussian 解中的径向变量统一解释为 $r^2=x^2+y^2$，在 $t=0.2$ 对三角形和四边形网格报告误差与收敛阶。

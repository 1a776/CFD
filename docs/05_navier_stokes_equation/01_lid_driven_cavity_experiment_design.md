# 第五题方腔顶盖驱动流实验配置与文献依据

## 1. 文档目的

本文用于确定第五题第一个算例“方形顶盖驱动腔流”的投影法实验配置。目标是把原题、原题列出的参考文献和 OpenFOAM 案例参数对应起来，形成可以直接用于 JSON、案例生成和后处理的实验方案。

本文只讨论方腔算例。等边三角形顶盖驱动腔流使用另一套几何、坐标、Reynolds 数和参考数据，不能直接套用本文的网格或采样位置。

## 2. 原题已经明确的内容

原始题目给出的方腔为：

$$
\Omega=[0,1]\times[0,1].
$$

上壁沿正方向运动：

$$
\boldsymbol{u}(x,1,t)=(1,0).
$$

其余三面为静止无滑移壁面：

$$
\boldsymbol{u}(0,y,t)=\boldsymbol{u}(1,y,t)=\boldsymbol{u}(x,0,t)=(0,0).
$$

取特征速度和特征长度：

$$
U=1,\qquad L=1.
$$

因此 Reynolds 数和运动黏度的关系为：

$$
Re=\frac{UL}{\nu}=\frac{1}{\nu}.
$$

本题要求：

$$
Re=1000,\qquad Re=3200,\qquad Re=5000.
$$

对应的运动黏度为：

$$
\nu_{1000}=1.0\times10^{-3},
$$

$$
\nu_{3200}=3.125\times10^{-4},
$$

$$
\nu_{5000}=2.0\times10^{-4}.
$$

初始条件可以取：

$$
\boldsymbol{u}(x,y,0)=\boldsymbol{0},
\qquad
p_k(x,y,0)=0.
$$

固体壁面的压力条件采用齐次法向梯度：

$$
\frac{\partial p_k}{\partial n}=0.
$$

由于压力只由梯度决定，还需要给出压力参考自由度：

```text
pRefCell  0;
pRefValue 0;
```

## 3. 参考文献分别负责什么

### 3.1 Ghia、Ghia 和 Shin（1982）

Ghia 等人的工作是方腔顶盖驱动流的经典基准，主要用于：

- 比较中心线水平速度 $u(0.5,y)$；
- 比较中心线竖直速度 $v(x,0.5)$；
- 判断主涡、左下角涡和右下角涡；
- 判断不同 Reynolds 数下的流线拓扑；
- 比较涡心位置。

本项目不把 Ghia 的参考网格直接复制为 OpenFOAM 案例网格，而是使用原题已经整理好的 Ghia 数据作为验证曲线和表格。

原题要求的 Ghia 采样位置是：

$$
x=0.5\quad\text{上的}\quad u(0.5,y),
$$

以及：

$$
y=0.5\quad\text{上的}\quad v(x,0.5).
$$

### 3.2 Erturk、Corke 和 Gökçöl（2005）

Erturk 等人的结果用于高精度交叉验证，主要用途是：

- 检查 $Re=1000$ 和 $Re=5000$ 的中心线速度；
- 检查主涡流函数值；
- 检查主涡涡量；
- 检查主涡中心坐标；
- 为网格加密后的结果提供更严格的参考。

原题记录了 Erturk 等人在 $601\times601$ 级别网格上的参考结果。这个分辨率应理解为文献高精度参考，不建议直接作为本项目第一轮普通细网格，否则计算成本和边界层分辨率要求会明显提高。

### 3.3 Chorin（1967）和 Guermond、Minev、Shen（2006）

这两篇文献用于确定投影法的算法逻辑，而不是方腔的物理基准：

$$
\boldsymbol{u}^{*}
\longrightarrow
p_k^{n+1}
\longrightarrow
\boldsymbol{u}^{n+1}.
$$

投影法的压力方程、速度修正和面通量修正应与这些文献的分步思想一致。

### 3.4 Issa（1986）和 Jasak（1996）

这两篇文献用于 PISO 求解器的算法和有限体积压力速度耦合。它们不能用来替代投影法的压力方程，也不用于规定方腔的 Reynolds 数。

### 3.5 Kohno、McQuain、Ribbens、Jagannathan 和 Erturk（2007）

这些文献主要对应等边三角形或其他三角形腔流：

- Kohno 和 Bathe：三角形网格和三角形腔流；
- McQuain 等：梯形腔流；
- Ribbens 等：三角形腔流；
- Jagannathan 等：直角三角形腔流；
- Erturk 和 Gökçöl（2007）：等边三角形腔流。

它们不能直接用来定义方腔的粗网格、细网格或中心线坐标。方腔实验应优先使用 Ghia 和 Erturk（2005）的参考结果。

## 4. 推荐的网格定义

### 4.1 为什么不能只写“粗网格”和“细网格”

原始题目写的是“不同分辨率的混合非结构网格”，但没有直接给出粗网格和细网格的单元数。

因此案例中必须同时记录：

- 网格类型；
- 目标边长分辨率；
- 实际单元总数；
- 实际边界面数量；
- 最小和最大单元尺寸；
- 网格质量；
- 壁面附近是否加密。

只写 `coarse` 和 `fine` 不足以复现实验。

### 4.2 推荐的三档网格

虽然截图中只有粗网格和细网格，但原题正文要求“至少三种逐步加密的网格”。因此推荐使用三档：

| 网格等级 | 名义分辨率 | 目标单元数量 | 用途 |
|---|---:|---:|---|
| coarse | $N=40$ | 约 $3.2\times10^3$ | 检查流动是否能启动并分辨主涡 |
| medium | $N=80$ | 约 $1.28\times10^4$ | 检查网格变化 |
| fine | $N=160$ | 约 $5.12\times10^4$ | 与文献进行主要比较 |

这里的 $N$ 表示单位正方形边长方向的目标分辨率，不直接表示最终混合非结构网格的精确总单元数。

如果使用每个四边形拆成两个三角形的混合网格，实际单元数量会不同。因此最终报告应使用 `checkMesh` 输出的真实单元数。

### 4.3 网格形状建议

建议采用“内部较规则、壁面附近加密”的混合非结构网格：

```text
内部区域：以四边形为主，或使用规则三角形/四边形混合；
四条壁面：增加边界附近单元；
上壁 y=1：重点加密，因为顶盖速度不连续角点会产生强剪切；
左下角、右下角：适度加密，用于分辨角涡；
左上角：Re=3200 和 Re=5000 时重点检查角涡。
```

网格加密必须保持三种网格的几何和边界命名完全一致，只改变分辨率参数。

### 4.4 网格验收条件

每个网格生成后执行：

```bash
checkMesh
```

至少记录：

- `Mesh OK`；
- 单元数量；
- 面数量；
- 边界 patch 数量；
- 最大非正交角；
- 最大 skewness；
- 最小单元体积；
- 是否存在负体积或非法面。

## 5. 推荐的数值配置

### 5.1 时间推进

方腔最终要求是稳态场，但投影法仍然通过瞬态时间推进达到稳态。

建议使用：

```text
startTime       0
deltaT          0.001
adjustTimeStep  yes
maxCo           0.2
maxDeltaT       0.01
```

时间步控制的基本约束为：

$$
Co_{\max}\leq 0.2.
$$

每个时间步记录：

```text
t
deltaT
maxCo
max |R_c^{mass}|
```

### 5.2 稳态判据

建议不要只依据线性求解器残差判断稳态，而是同时检查：

$$
\max_c |R_c^{mass}|<10^{-10},
$$

$$
\left\|\boldsymbol{u}^{n+1}-\boldsymbol{u}^{n}\right\|_{\infty}<10^{-8},
$$

以及中心线速度、主涡位置在连续若干时间步内基本不变。

工程实现中建议使用：

```text
minimumSteadySteps  100
steadyVelocityTol   1e-8
steadyMassTol       1e-10
maxSteps            200000
```

高 Reynolds 数下可能需要更多时间步，不能用低 Reynolds 数的固定步数直接套用。

### 5.3 对流和扩散离散

投影法和 PISO 对照实验必须使用相同的空间离散设置。第一轮推荐：

```foam
div(phi,U)              Gauss linearUpwind grad(U);
laplacian(nu,U)         Gauss linear corrected;
grad(p)                 Gauss linear;
snGrad                  corrected;
```

选择理由：

- `linearUpwind` 比一阶迎风耗散小，适合观察主涡和中心线速度；
- `corrected` 可以处理混合非结构网格中的非正交扩散；
- 投影法与 PISO 使用完全相同的格式，比较才具有意义。

为了单独研究格式影响，可以在主实验完成后增加：

```foam
div(phi,U)              Gauss upwind;
div(phi,U)              Gauss linear;
div(phi,U)              Gauss linearUpwind grad(U);
```

格式研究不应和 Reynolds 数研究混在同一张结论表中。

### 5.4 线性求解器

推荐初始设置：

```foam
solvers
{
    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0;
    }

    p
    {
        solver          GAMG;
        tolerance       1e-8;
        relTol          0;
        smoother        GaussSeidel;
    }
}
```

投影法和 PISO 应保持相同的线性求解容差。

## 6. 推荐的 JSON 配置

建议第五题投影法的 JSON 至少包含：

```json
{
  "solverFamily": "05_navier_stokes_equation",
  "solver": "projectionFoamStudent",
  "caseName": "01_lid_driven_cavity_projection",
  "problem": "lid_driven_cavity",
  "domain": [0.0, 1.0, 0.0, 1.0],
  "velocity": [1.0, 0.0, 0.0],
  "reynoldsNumbers": [1000, 3200, 5000],
  "viscosities": [0.001, 0.0003125, 0.0002],
  "meshLevels": {
    "coarse": 40,
    "medium": 80,
    "fine": 160
  },
  "meshType": "hybrid_unstructured",
  "time": {
    "startTime": 0.0,
    "initialDeltaT": 0.001,
    "maxCo": 0.2,
    "maxDeltaT": 0.01
  },
  "discretization": {
    "ddtScheme": "Euler",
    "divScheme": "Gauss linearUpwind grad(U)",
    "laplacianScheme": "Gauss linear corrected",
    "gradScheme": "Gauss linear",
    "snGradScheme": "corrected"
  },
  "pressureReference": {
    "cell": 0,
    "value": 0.0
  },
  "steadyState": {
    "velocityTolerance": 1e-8,
    "massTolerance": 1e-10,
    "minimumSteadySteps": 100,
    "maximumSteps": 200000
  }
}
```

说明：

- `solverFamily` 表示第五题求解器族；
- `solver` 表示具体使用投影法求解器；
- `caseName` 表示方腔投影法算例；
- `reynoldsNumbers` 控制三种 Reynolds 数；
- `viscosities` 必须和 Reynolds 数一一对应；
- `meshLevels` 控制粗、中、细三种网格；
- `time` 控制时间推进；
- `discretization` 控制 OpenFOAM 的空间离散；
- `pressureReference` 解决压力常数零空间；
- `steadyState` 控制稳态判断。

## 7. 推荐的案例目录

推荐使用：

```text
cases/05_navier_stokes_equation/
└── 01_lid_driven_cavity_projection/
    ├── Re1000_coarse/
    ├── Re1000_medium/
    ├── Re1000_fine/
    ├── Re3200_coarse/
    ├── Re3200_medium/
    ├── Re3200_fine/
    ├── Re5000_coarse/
    ├── Re5000_medium/
    └── Re5000_fine/
```

每个 `Re...` 案例包含：

```text
0.orig/
├── U
└── p

constant/
├── physicalProperties
└── polyMesh/

system/
├── controlDict
├── fvSchemes
└── fvSolution

Allrun
Allclean
case.foam
metadata.json
```

## 8. 12 组对照实验的组织方式

投影法单独完成 9 组：

| 求解器 | Reynolds 数 | 网格等级 | 数量 |
|---|---:|---|---:|
| projection | 1000 | coarse、medium、fine | 3 |
| projection | 3200 | coarse、medium、fine | 3 |
| projection | 5000 | coarse、medium、fine | 3 |

然后使用相同的网格、边界、时间步和空间离散，增加 PISO 的 9 组：

| 求解器 | Reynolds 数 | 网格等级 | 数量 |
|---|---:|---|---:|
| PISO | 1000 | coarse、medium、fine | 3 |
| PISO | 3200 | coarse、medium、fine | 3 |
| PISO | 5000 | coarse、medium、fine | 3 |

因此，如果两种求解器都完成，方腔算例总数是：

$$
2\times3\times3=18.
$$

截图中的两级网格版本则是：

$$
2\times3\times2=12.
$$

但两级网格只能作为最低交付版本；正式网格收敛性分析建议采用三档网格。

## 9. 后处理结果

### 9.1 必须生成的结果

每个案例至少保存：

```text
streamlines.png
u_centerline.png
v_centerline.png
vorticity.png
pressure.png
mass_residual.png
```

以及：

```text
u_centerline.csv
v_centerline.csv
vortex_centres.csv
mass_residual.csv
summary.json
```

### 9.2 中心线采样

水平速度采样：

$$
u(0.5,y).
$$

竖直速度采样：

$$
v(x,0.5).
$$

必须使用相同的采样方法比较所有 Reynolds 数、网格和求解器。不能在一个案例中使用最近单元值，在另一个案例中使用高阶插值。

### 9.3 流函数和涡量

若计算流函数，统一使用：

$$
u=\frac{\partial\psi}{\partial y},
\qquad
v=-\frac{\partial\psi}{\partial x}.
$$

涡量定义为：

$$
\omega
=\frac{\partial v}{\partial x}
-\frac{\partial u}{\partial y}
=-\nabla^2\psi.
$$

流函数的常数参考不影响速度，但比较不同网格时必须采用相同参考方式。

## 10. 与文献的比较方法

推荐按照以下顺序比较：

1. 先检查网格和连续性残差；
2. 再比较 $u(0.5,y)$；
3. 再比较 $v(x,0.5)$；
4. 再比较主涡中心；
5. 最后比较角涡和流函数、涡量。

误差建议计算：

$$
E_u
=
\left(
\frac{1}{m}
\sum_{j=1}^{m}
\left[
u_{\mathrm{num}}(0.5,y_j)
-
u_{\mathrm{ref}}(0.5,y_j)
\right]^2
\right)^{1/2},
$$

$$
E_v
=
\left(
\frac{1}{m}
\sum_{j=1}^{m}
\left[
v_{\mathrm{num}}(x_j,0.5)
-
v_{\mathrm{ref}}(x_j,0.5)
\right]^2
\right)^{1/2}.
$$

文献采样点和本项目采样点不完全一致时，先对数值曲线进行一维插值，再在文献坐标处比较。

## 11. 物理合理性检查

### $Re=1000$

至少应观察到：

- 一个主涡；
- 左下角角涡；
- 右下角角涡；
- 顶盖附近明显剪切层。

### $Re=3200$

应观察到更薄的壁面剪切层，并开始分辨左上角涡。

### $Re=5000$

应观察到：

- 更薄的边界层；
- 更明显的角涡；
- 主涡位置发生变化；
- 中心线速度梯度增大。

如果网格过粗，可能出现角涡消失、主涡位置明显偏移或中心线曲线过度平滑。这首先应归因于网格分辨率不足，而不能直接归因于投影法本身。

## 12. 最终推荐

本项目第一轮方腔投影法实验采用：

```text
Re = 1000, 3200, 5000
网格 = coarse, medium, fine
网格目标 = N40, N80, N160
maxCo = 0.2
对流 = Gauss linearUpwind grad(U)
扩散 = Gauss linear corrected
梯度 = Gauss linear
初始速度 = 0
顶盖速度 = (1,0,0)
其他壁面 = noSlip
压力边界 = zeroGradient
压力参考 = cell 0, value 0
稳态判据 = 速度变化和质量残差同时满足
```

该配置的依据分为三层：

1. 原题直接规定的区域、边界、Reynolds 数和输出要求；
2. Ghia 与 Erturk 的方腔基准数据；
3. 为了在 OpenFOAM 中形成可复现实验而补充的网格、时间步、线性求解器和稳态判据。

第三层不是文献唯一规定的答案，而是本项目为保证公平比较和可重复性所采用的工程化实验约定。

## 13. 参考文献

- Chorin, A. J. (1967), *A numerical method for solving incompressible viscous flow problems*.
- Erturk, E., Corke, T. C. and Gökçöl, C. (2005), *Numerical solutions of 2-D steady incompressible driven cavity flow at high Reynolds numbers*.
- Erturk, E. and Gökçöl, O. (2007), *Fine grid numerical solutions of triangular cavity flow*.
- Ghia, U., Ghia, K. N. and Shin, C. T. (1982), *High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method*.
- Guermond, J.-L., Minev, P. and Shen, J. (2006), *An overview of projection methods for incompressible flows*.
- Issa, R. I. (1986), *Solution of the implicitly discretised fluid flow equations by operator-splitting*.
- Jasak, H. (1996), *Error Analysis and Estimation for the Finite Volume Method with Applications to Fluid Flows*.
- Jagannathan, A., Mohan, R. and Dhanak, M. (2014), *A spectral method for the triangular cavity flow*.
- Kohno, H. and Bathe, K. J. (2006), *A flow-condition-based interpolation finite element procedure for triangular grids*.
- McQuain, W. D., Ribbens, C. J., Wang, C.-Y. and Watson, L. T. (1994), *Steady viscous flow in a trapezoidal cavity*.
- Ribbens, C. J., Watson, L. T. and Wang, C.-Y. (1994), *Steady viscous flow in a triangular cavity*.

# 第五题方腔顶盖驱动流：稳态判据修正记录

## 1. 修改目的

本次修改针对第五题第一个案例“方腔顶盖驱动流”的稳态判断问题。

原来的计算虽然能够完成投影法的速度预测、压力校正和不可压缩性控制，但求解器主要依据速度场的散度判断计算状态。由于不可压缩条件成立只说明速度场满足

$$
\nabla\cdot\boldsymbol{u}\approx 0,
$$

并不说明速度场已经不再随时间变化，因此固定运行到某个时间，例如 `t=10`，不能直接称为已经达到稳态。

本次修改将“质量守恒”和“稳态”分成两个独立判据：

$$
\max_c\left|\nabla\cdot\boldsymbol{u}^{n+1}\right|
\leq \varepsilon_{\mathrm{mass}},
$$

以及

$$
\max_c\left|\boldsymbol{u}^{n+1}_c-\boldsymbol{u}^{n}_c\right|
\leq \varepsilon_{\mathrm{U}}.
$$

只有两个条件同时满足，并且连续满足若干时间步，才提前结束计算并将当前场作为稳态结果保存。

## 2. 原问题与根本原因

原来的方腔算例采用瞬态推进方法计算最终稳态解。若仅观察压力方程或速度校正后的散度，可能出现下面的情况：

- `max |div(U)|` 已经很小，说明离散速度场近似无散；
- 但 `max |U - Uprevious|` 仍然较大，说明流场仍在演化；
- 此时如果仅因为散度小就停止，会把一个尚未收敛的瞬态场当作稳态场；
- 固定的 `endTime=10` 只规定了最长运行时间，不提供稳态证明。

因此本次修正的根本原因是：原求解器缺少对相邻时间步速度变化量的检查，也缺少连续满足稳态条件的要求。

## 3. 修正后的稳态逻辑

### 3.1 保存上一个时间步的速度

求解器在时间循环开始前建立 `UPrevious`，并在每个时间步开始时保存旧速度：

```cpp
volVectorField UPrevious
(
    IOobject
    (
        "UPrevious",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U
);
```

进入新的时间步后：

```cpp
UPrevious = U;
```

这里的 `UPrevious` 表示当前时间步求解开始前的速度场，当前的 `U` 则是在预测速度、压力方程和投影修正完成后的速度场。因此两者之差对应

$$
\boldsymbol{u}^{n+1}-\boldsymbol{u}^{n}.
$$

### 3.2 计算质量残差

投影结束后，求解器通过当前面通量计算速度散度：

```cpp
const volScalarField continuity
(
    fvc::div(phi)
);
```

随后取所有单元中的最大绝对值：

```cpp
const scalar maxMassResidual =
    gMax(continuityMag.primitiveField());
```

它对应离散判据

$$
R_{\mathrm{mass}}
=
\max_c
\left|
(\nabla\cdot\boldsymbol{u})_c
\right|.
$$

该量主要检查压力投影是否成功消除了速度场的离散散度。

### 3.3 计算速度变化量

求解器使用：

```cpp
tmp<volScalarField> tVelocityChange(mag(U - UPrevious));
const volScalarField& velocityChange = tVelocityChange();

const scalar maxVelocityChange =
    gMax(velocityChange.primitiveField());
```

对应离散判据：

$$
R_{\mathrm{steady}}
=
\max_c
\left|
\boldsymbol{u}^{n+1}_c-\boldsymbol{u}^{n}_c
\right|.
$$

这个量反映流场是否还在变化。即使

$$
R_{\mathrm{mass}}\ll 1,
$$

只要

$$
R_{\mathrm{steady}}
$$

仍然较大，就不能认为已经达到稳态。

### 3.4 连续时间步确认

单个时间步满足阈值也可能只是暂时波动，因此求解器增加连续满足计数：

```cpp
if
(
    steadyStateControl
 && runTime.timeIndex() >= minimumSteadySteps
 && maxVelocityChange <= steadyVelocityTol
 && maxMassResidual <= steadyMassTol
)
{
    ++steadyStepCount;
}
else
{
    steadyStepCount = 0;
}
```

只有满足以下条件才增加计数：

$$
n_{\mathrm{step}}
\geq n_{\mathrm{minimum}},
$$

$$
R_{\mathrm{steady}}\leq\varepsilon_{\mathrm{U}},
$$

$$
R_{\mathrm{mass}}\leq\varepsilon_{\mathrm{mass}}.
$$

一旦某个时间步不满足条件，计数清零：

```cpp
steadyStepCount = 0;
```

连续满足 `requiredSteadySteps` 次后，求解器执行：

```cpp
runTime.writeAndEnd();
```

这会先写出当前稳态场，再结束时间循环。

## 4. 修改的代码文件

### 4.1 `UDF/solver/05_navier_stokes_equation/projectionFoamStudent/projectionFoamStudent.C`

修改内容：

- 新增 `UPrevious` 速度场；
- 从 `system/controlDict` 读取稳态控制参数；
- 在投影修正后计算 `max |div(U)|`；
- 计算 `max |U - Uprevious|`；
- 增加连续满足判据的计数器；
- 达到稳态后使用 `runTime.writeAndEnd()` 保存并结束；
- 在日志中输出稳态诊断量和稳态结束时间。

主要代码位置：

- 稳态参数读取：约第 214--235 行；
- 速度变化量和质量残差计算：约第 473--486 行；
- 日志输出：约第 488--493 行；
- 稳态判据：约第 519--546 行。

本次没有改变以下数值算法：

- 动量预测方程；
- 压力投影方程；
- 压力参考单元；
- 速度边界条件；
- 网格生成方式；
- Reynolds 数；
- 线性求解器类型。

### 4.2 `scripts/common/lid_cavity.py`

该脚本负责根据 JSON 生成方腔 OpenFOAM 案例。

修改内容：

- 从 JSON 读取稳态控制字段；
- 将这些字段写入生成的 `system/controlDict`；
- 在生成的控制字典中加入稳态判据注释；
- 保留 `endTime` 作为最长运行时间，稳态条件满足时允许提前结束。

生成的 `system/controlDict` 包含：

```text
steadyStateControl true;
steadyVelocityTol 1e-6;
steadyMassTol 1e-8;
minimumSteadySteps 1000;
requiredSteadySteps 20;
```

这些字段的含义是：

| 字段 | 含义 |
|---|---|
| `steadyStateControl` | 是否启用稳态提前终止 |
| `steadyVelocityTol` | 相邻时间步最大速度变化阈值 |
| `steadyMassTol` | 最大速度散度阈值 |
| `minimumSteadySteps` | 至少推进的时间步数 |
| `requiredSteadySteps` | 连续满足判据的时间步数 |

### 4.3 六个 JSON 配置文件

修改目录：

```text
scripts/configs/05_navier_stokes_equation/
```

六个配置都加入了相同的稳态判据参数，只根据 Reynolds 数和网格等级设置不同的最长运行时间：

| 实验组 | Reynolds 数 | 网格 | 单元数 | 最长 `endTime` |
|---|---:|---:|---:|---:|
| 01 | 1000 | `N=40` | 1600 | 100 |
| 02 | 1000 | `N=80` | 6400 | 100 |
| 03 | 3200 | `N=40` | 1600 | 150 |
| 04 | 3200 | `N=80` | 6400 | 150 |
| 05 | 5000 | `N=40` | 1600 | 200 |
| 06 | 5000 | `N=80` | 6400 | 200 |

`endTime` 是保险上限，不是稳态判定本身。实际结果以日志中的

```text
Steady state reached at Time = ...
```

为准。

## 5. 编译验证

本次修改后重新编译：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
source /opt/openfoam14/etc/bashrc
export FOAM_USER_APPBIN="$PWD/build/05_navier_stokes_equation/bin"
mkdir -p "$FOAM_USER_APPBIN"
wmake UDF/solver/05_navier_stokes_equation/projectionFoamStudent
```

编译成功，并生成：

```text
build/05_navier_stokes_equation/bin/projectionFoamStudent
```

编译过程中仍会出现：

```text
opal_ifinit: socket() failed with errno=1
/opt/openfoam14/etc/config.sh/paraview:
[: !=: 需要一元运算符
```

这两个信息属于当前 OpenFOAM/ParaView 环境脚本提示，没有导致编译失败，也没有阻止案例运行。

## 6. 六组实验验证结果

六组实验按照单组依次运行，未并行运行。每组运行结束后单独执行后处理。

### 6.1 稳态时间

| 实验组 | Reynolds 数 | 网格 | 稳态时间 | 日志证据 |
|---|---:|---:|---:|---|
| 01 | 1000 | `N=40` | `35.867` | `cases/05_navier_stokes_equation/01_lid_driven_cavity_projection_Re1000_coarse/log.projectionFoamStudent` |
| 02 | 1000 | `N=80` | `29.079` | `cases/05_navier_stokes_equation/02_lid_driven_cavity_projection_Re1000_fine/log.projectionFoamStudent` |
| 03 | 3200 | `N=40` | `71.176` | `cases/05_navier_stokes_equation/03_lid_driven_cavity_projection_Re3200_coarse/log.projectionFoamStudent` |
| 04 | 3200 | `N=80` | `80.942` | `cases/05_navier_stokes_equation/04_lid_driven_cavity_projection_Re3200_fine/log.projectionFoamStudent` |
| 05 | 5000 | `N=40` | `127.063` | `cases/05_navier_stokes_equation/05_lid_driven_cavity_projection_Re5000_coarse/log.projectionFoamStudent` |
| 06 | 5000 | `N=80` | `129.846` | `cases/05_navier_stokes_equation/06_lid_driven_cavity_projection_Re5000_fine/log.projectionFoamStudent` |

每个日志均包含：

```text
Steady state reached at Time = ...
End
```

这证明六组实验不是简单运行到固定时间后强制结束，而是实际触发了新的稳态判据。

### 6.2 稳态判据的数值证据

六组日志在触发稳态前均满足：

$$
\max|\boldsymbol{u}^{n+1}-\boldsymbol{u}^{n}|
\leq 1.0\times10^{-6},
$$

以及：

$$
\max|\nabla\cdot\boldsymbol{u}^{n+1}|
\leq 1.0\times10^{-8}.
$$

例如，第一组在结束前记录的量约为：

```text
max |div(U)| = 6.79e-09
max |U - Uprevious| = 9.96e-07
Steady state reached at Time = 35.867...
```

第六组结束前记录的量约为：

```text
max |div(U)| = 6.90e-10
max |U - Uprevious| = 9.98e-07
Steady state reached at Time = 129.846...
```

这也说明了为什么不能只观察质量残差：第六组在较早阶段质量残差已经达到 `10^-9` 量级，但速度变化量仍为 `10^-5` 量级，因此当时不会被判定为稳态。

## 7. 后处理产物

每组实验都生成以下数据：

```text
data/05_navier_stokes_equation/cases/<caseName>/
├── u_centerline.csv
├── v_centerline.csv
└── summary.json
```

其中：

- `u_centerline.csv` 保存竖直中心线 `x=0.5` 上的水平速度；
- `v_centerline.csv` 保存水平中心线 `y=0.5` 上的竖直速度；
- `summary.json` 保存 Reynolds 数、网格规模、最终时间、速度范围及与 Ghia 数据的误差。

每组图片位于：

```text
figures/05_navier_stokes_equation/cases/<caseName>/
├── field_and_streamlines.png
└── centerline_comparison.png
```

其中：

- `field_and_streamlines.png` 显示速度大小场和流线；
- `centerline_comparison.png` 显示数值中心线速度与 Ghia 参考数据的比较。

六组均已检查，数据文件和图片文件完整。

## 8. 与 Ghia 参考数据的基本比较

后处理给出了如下中心线误差摘要：

| 实验组 | `u` 最大误差 | `u` RMSE | `v` 最大误差 | `v` RMSE |
|---|---:|---:|---:|---:|
| Re=1000, N=40 | 0.2228 | 0.0600 | 0.0870 | 0.0348 |
| Re=1000, N=80 | 0.0984 | 0.0264 | 0.0449 | 0.0146 |
| Re=3200, N=40 | 0.3852 | 0.1072 | 0.1556 | 0.0662 |
| Re=3200, N=80 | 0.1896 | 0.0488 | 0.0855 | 0.0240 |
| Re=5000, N=40 | 0.4496 | 0.1255 | 0.1848 | 0.0799 |
| Re=5000, N=80 | 0.2471 | 0.0650 | 0.1055 | 0.0347 |

可以观察到：

- 在三个 Reynolds 数下，`N=80` 相比 `N=40` 的中心线误差均明显下降；
- Reynolds 数增大后，流动结构和壁面附近剪切层更加复杂，粗网格误差增大；
- 当前结果具备方腔顶盖驱动流的基本物理结构，但仍属于规则四边形网格上的有限分辨率验证；
- 不能将当前结果直接表述为对高精度文献数据的完全复现。

## 9. 当前限制

本次修正解决的是“是否达到稳态”的判断问题，并不等同于所有数值误差问题都已经消失。

当前仍有以下限制：

1. 当前六组案例使用规则四边形 `blockMesh` 网格，不是原题中所说的混合非结构网格；
2. 当前网格只采用 `N=40` 和 `N=80` 两档，若要严格做空间收敛性研究，还需要至少增加一档网格；
3. `steadyVelocityTol=1e-6` 是当前工程实验的判据，不能替代网格无关性分析；
4. Ghia 数据比较主要基于中心线速度，尚未完整比较流函数、涡量和涡心坐标；
5. 当前投影法实现已经能够判断稳态，但高 Reynolds 数下的精度仍受网格、时间步和边界角点分辨率影响。

因此本次修正后的准确表述是：

> 六组方腔顶盖驱动流实验均通过速度变化量、质量残差和连续时间步三个条件确认达到数值稳态；在当前规则四边形网格和有限分辨率下，结果能够展示正确的方腔流动结构，并且网格加密后中心线误差下降，但仍需通过混合非结构网格、更多分辨率和更完整的文献指标比较来完成严格验证。

## 10. 复现命令

编译求解器：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
source /opt/openfoam14/etc/bashrc
export FOAM_USER_APPBIN="$PWD/build/05_navier_stokes_equation/bin"
mkdir -p "$FOAM_USER_APPBIN"
wmake UDF/solver/05_navier_stokes_equation/projectionFoamStudent
```

运行单组实验：

```bash
python3 scripts/run_lid_cavity.py \
  --config scripts/configs/05_navier_stokes_equation/01_lid_driven_cavity_projection_Re1000_coarse.json \
  --overwrite
```

运行后处理：

```bash
python3 scripts/postprocess_lid_cavity.py \
  --case cases/05_navier_stokes_equation/01_lid_driven_cavity_projection_Re1000_coarse
```

检查稳态日志：

```bash
rg -n \
  "Steady state reached|FOAM FATAL|End|max \|U - Uprevious\||max \|div\(U\)\|" \
  cases/05_navier_stokes_equation/*/log.projectionFoamStudent
```


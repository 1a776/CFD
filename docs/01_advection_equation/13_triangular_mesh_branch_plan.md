# 第 1 题三角形网格分支的改造方案

本文只讨论 PDF 第 1 题中的第一个例子：

$$\frac{\partial \phi}{\partial t}+\nabla\cdot(\boldsymbol{u}\phi)=0,$$

对应的正弦波平移测试。目标不是立刻改代码，而是先把“应该改哪些地方”说清楚，特别是：

- 代码层面要改哪里；
- `cases/` 目录里要改哪里；
- 为什么这些地方必须改；
- 哪些地方可以先不动。

这里的三角形网格，指的是：

- 区域仍然是 $[0,1]^2$；
- 周期边界仍然是 `x-y` 双方向；
- 速度仍然是 `U = (1,1,0)`；
- `CFL = 0.2`；
- 终点仍然是 `t = 1.0`；
- 只先做三角形网格分支，不碰刚体旋转例子。

---

## 1. 现在这套工程为什么只能跑四边形

当前能跑通的是：

- `01_sine_wave_quad_upwind`
- `02_sine_wave_quad_linearUpwind`

它们的核心假设都是“结构化四边形网格”。这个假设已经写进了几个地方：

| 文件 | 现在的假设 |
|---|---|
| `scripts/common/foam_case.py` | 只允许 `meshType == "quad"` |
| `scripts/common/mesh_tools.py` | 只会改 `system/blockMeshDict` 里的 `hex (0 1 2 3 4 5 6 7) (N N 1)` |
| `scripts/common/advection_tools.py` | 默认单元中心是 `(i+0.5)/nx, (j+0.5)/ny` |
| `scripts/common/postprocess_case.py` | 默认能把场 `reshape(ny, nx)` |
| `scripts/common/postprocess_case.py` | 默认网格名是 `quad` |
| `cases/01_advection_equation/01_sine_wave_quad/Nxx/system/blockMeshDict` | 网格来源就是 `blockMesh` |

所以三角形版本不是“把配置文件里的 `meshType` 改成 `tri`”这么简单。  
真正要改的是“网格来源、初值生成、误差统计、画图方式、案例目录结构”。

---

## 2. 推荐的总体路线

对于这个题目，三角形网格最稳妥的路线是：

1. 保持求解器不变，仍然用 `explicitAdvectionFoamStudent`；
2. 只为 `03_sine_wave_tri_upwind` 建一个新的 case family；
3. 三角形网格用 `Gmsh` 生成，再导入 OpenFOAM；
4. 初值 `T(x,y,0)=sin(2\pi(x+y))` 改成“按真实 cell center 生成”；
5. 后处理改成“按真实 cell area / cell center 计算误差和画图”；
6. 先只做一个 `N20` 的最小闭环，再推广到 `N10/N20/N40/N80`。

这里最关键的一点是：  
**三角形分支应该和四边形分支并列，不要硬塞进 quad 的目录里。**

也就是说：

```text
cases/01_advection_equation/01_sine_wave_quad/              # 保留，作为四边形参考
cases/01_advection_equation/02_sine_wave_quad_linearUpwind/ # 保留，作为另一种格式参考
cases/01_advection_equation/03_sine_wave_tri_upwind/        # 新增，作为三角形分支
```

---

## 3. 脚本层要改哪些文件

### 3.1 `scripts/common/foam_case.py`

这是最核心的入口。当前这里直接写死了：

```python
if config.mesh_type != "quad":
    raise NotImplementedError(...)
```

所以三角形分支第一步就会被挡住。

#### 这里要改什么

1. 去掉“只允许 quad”的硬判断；
2. 按 `config.mesh_type` 分发：
   - `quad` 走现有 `blockMesh` 路线；
   - `tri` 走新的三角网格路线；
3. `prepare_case()` 不能再默认先写 `0.orig/T` 再靠 `blockMesh`；
4. `tri` 分支必须允许“先生成网格，再写初值”。

#### 为什么要改

因为三角形网格下，初值文件 `0.orig/T` 不再能靠简单的

```text
(i+0.5)/nx, (j+0.5)/ny
```

来生成。  
你必须先知道真实的 cell center，才能写 `T`。

#### 建议改法

建议把 `prepare_case()` 拆成两个内部路径：

- `prepare_quad_case(...)`
- `prepare_tri_case(...)`

这样不会把四边形逻辑和三角形逻辑搅在一起。

---

### 3.2 `scripts/common/mesh_tools.py`

当前这里只有一个函数：

```python
patch_block_mesh_resolution(case, resolution)
```

它只会改 `system/blockMeshDict`，而且只认：

```foam
hex (0 1 2 3 4 5 6 7) (N N 1)
```

#### 这里要改什么

1. 保留 quad 的 `patch_block_mesh_resolution()`；
2. 新增 tri 版网格生成函数，例如：
   - `write_tri_mesh_geo(...)`
   - `build_tri_mesh(...)`
   - `prepare_tri_mesh_case(...)`
3. tri 版不要再依赖 `blockMeshDict` 作为主网格来源；
4. tri 版要能输出：
   - Gmsh 源文件；
   - 导入后的 `constant/polyMesh/`；
   - 需要的话，`createPatchDict` 或 `changeDictionary` 配套文件。

#### 为什么要改

因为 `blockMesh` 只适合当前这类规则四边形背景网格。  
你要做的是“正方形域上的三角形剖分”，最直接的方式是 Gmsh。

#### 具体建议

tri 网格最好单独放一个源文件，比如：

```text
cases/01_advection_equation/03_sine_wave_tri_upwind/N20/system/mesh.geo
```

或者：

```text
cases/01_advection_equation/03_sine_wave_tri_upwind/N20/system/tri_mesh.geo
```

然后在 `Allrun.pre` 里调用 Gmsh，再调用 `gmshToFoam`。

---

### 3.3 `scripts/common/advection_tools.py`

这份文件现在有两个对 tri 不友好的地方：

#### 第一处：`mesh_resolution(case)`

现在它从 `system/blockMeshDict` 里读 `N`。  
这对 `quad` 没问题，但对 `tri` 不稳。

##### tri 要怎么改

更稳妥的做法是：

- `resolution` 由 `metadata.json` 提供；
- `meshType` 由 `metadata.json` 提供；
- 如果以后有别的网格类型，再在 metadata 里扩展。

也就是说，**不要再把 `blockMeshDict` 当成唯一的网格信息源**。

#### 第二处：`write_initial_field(case, nx, ny)`

现在它写的是：

```python
sin(2*pi*(x+y))
```

但取点方式是规则网格中心：

```python
x = (i + 0.5) / nx
y = (j + 0.5) / ny
```

##### tri 要怎么改

要把“坐标生成”改成“读取真实 cell center”：

- 先得到三角网格每个单元的中心；
- 再在这些中心上计算初值；
- 再按单元顺序写入 `0.orig/T`。

##### 为什么不能继续用 `(i,j)` 规则索引

因为三角网格不是方阵单元，单元编号不会再天然对应一个二维矩阵。

---

### 3.4 `scripts/common/postprocess_case.py`

这是 tri 分支最容易漏掉的地方。

当前它做了三件 quad 专用的事：

1. `reshape(ny, nx)`；
2. 用 `pcolormesh` 画图；
3. 用 `THICKNESS / (nx * ny)` 当作统一单元体积。

#### tri 要怎么改

1. 不再把数值场强行 reshape 成二维矩阵；
2. 读取真实的三角形单元几何；
3. 面积或体积加权误差要用每个单元自己的 `V_i`；
4. 画图要改成三角剖分绘图，比如：
   - `matplotlib.tri.Triangulation`
   - `tripcolor`
   - `tricontourf`

#### 哪些函数必须改

- `write_field_data(...)`
- `write_time_history(...)` 可以保留
- `plot_field_comparison(...)`
- `plot_diagonal_profile(...)`
- `postprocess_case(...)`

#### 为什么要改

因为现在的后处理逻辑默认：

```text
第 j 行第 i 列 = 第 j*n + i 个单元
```

这个假设对三角网格不成立。

---

### 3.5 `scripts/common/study_analysis.py`

这个文件基本可以保留，但前提是：

- `summary.json` 里的 `mesh` 字段必须正确；
- `nCells`、`nominalH`、`normalizedL1` 等字段必须由 tri 版后处理写对。

也就是说，`study_analysis.py` 不是 tri 的主改点，  
它主要吃的是 `summary.json` 和 `raw_results.csv`。

---

### 3.6 `scripts/configs/01_advection_equation/03_sine_wave_tri_upwind.json`

这个配置文件已经存在，但现在写的是：

```json
"meshType": "tri",
"implemented": false
```

#### 这里要改什么

等 tri 链路真的接通以后，把：

```json
"implemented": true
```

打开。

#### 建议再加什么

如果你想让工程更清楚，建议补一个字段，例如：

- `meshSource: "gmsh"`
- 或者 `meshBackend: "gmsh"`

这样以后 `meshType` 表示“单元类型”，`meshSource` 表示“网格怎么来的”，不会混。

---

## 4. `cases/` 目录里具体要改哪里

这部分是你最关心的。

### 4.1 新建哪个目录

应该新增：

```text
cases/01_advection_equation/03_sine_wave_tri_upwind/
```

并且继续按分辨率分子目录：

```text
cases/01_advection_equation/03_sine_wave_tri_upwind/N10/
cases/01_advection_equation/03_sine_wave_tri_upwind/N20/
cases/01_advection_equation/03_sine_wave_tri_upwind/N40/
cases/01_advection_equation/03_sine_wave_tri_upwind/N80/
```

### 4.2 每个 `Nxx/` 里保留什么

可以直接沿用 quad 案例里的这些内容：

- `0.orig/U`
- `system/controlDict`
- `system/fvSchemes`
- `system/fvSolution`

这些本质上和网格类型无关，和物理方程、时间推进、对流格式有关。

### 4.3 每个 `Nxx/` 里要替换什么

四边形版本里的：

```text
system/blockMeshDict
```

在 tri 版本里不应该继续作为主网格来源。

建议替换成：

```text
system/mesh.geo
```

或者：

```text
system/tri_mesh.geo
```

如果你还想保留周期配对的转换步骤，可以再加：

```text
system/createPatchDict
```

### 4.4 tri 版建议的 case 结构

```text
cases/01_advection_equation/03_sine_wave_tri_upwind/N20/
├── 0.orig/
│   ├── U
│   └── T              # tri 版初值由网格几何写入
├── 0/
├── constant/
│   └── polyMesh/      # 运行后生成
├── system/
│   ├── controlDict
│   ├── fvSchemes
│   ├── fvSolution
│   ├── mesh.geo       # Gmsh 源文件
│   └── createPatchDict  # 可选，用于 cyclic 边界转换
├── Allrun.pre
├── Allrun
├── Allclean
├── case.foam
└── metadata.json
```

### 4.5 `Allrun.pre` 应该做什么

tri 分支最好有独立的预处理脚本：

1. 生成三角网格；
2. 导入 OpenFOAM；
3. 把四条边转换成正确的 `cyclic` 对；
4. 必要时执行 `checkMesh`；
5. 准备可写初值所需的网格几何信息。

### 4.6 `Allrun` 应该做什么

`Allrun` 里只做运行主流程：

1. 恢复 `0/`；
2. 跑求解器；
3. 输出日志；
4. 后处理。

### 4.7 `Allclean` 应该做什么

除了现在 quad 版本会删的内容外，tri 版还要删：

- Gmsh 生成的临时网格文件；
- 导入后的 `constant/polyMesh/`；
- 运行日志；
- 时间目录。

但不要删：

- `0.orig/`
- `system/`
- 网格源文件

---

## 5. 初值和误差为什么要改成“按几何”而不是“按索引”

这一点很重要。

### 5.1 现在的做法

当前 quad 版是这样理解问题的：

- 单元是规则矩形；
- 每个单元的中心位置可以直接由 `(i,j)` 算出来；
- 每个单元面积相同；
- 所以写初值、算误差都很简单。

### 5.2 tri 版不可以再这样做

三角形网格下：

- 单元不再是规则矩形；
- 单元编号不再对应二维数组；
- 误差必须用单元面积加权；
- 图像也不能直接靠 `pcolormesh` 画。

### 5.3 所以 tri 版建议新增一份几何元数据

建议给每个 `Nxx` 保存一份几何文件，例如：

```text
data/01_advection_equation/cases/03_sine_wave_tri_upwind/N20/mesh_geometry.json
```

里面至少保存：

- `cellId`
- `cellCenter`
- `cellArea`
- `triangleConnectivity`

这样：

- 初值可以按 `cellCenter` 写；
- 误差可以按 `cellArea` 算；
- 画图可以按 `triangleConnectivity` 做；
- 诊断时也更容易。

---

## 6. 哪些地方可以先不动

为了不把问题搞大，下面这些地方第一阶段可以先不改：

- 求解器主体 `explicitAdvectionFoamStudent`
- `scripts/run_case.py` 的命令行接口
- `scripts/run_study.py` 的命令行接口
- `scripts/common/study_analysis.py` 的主分析流程
- `README.md` 的大框架

原因很简单：  
**tri 分支的核心问题不在“怎么跑命令”，而在“怎么生成几何一致的三角网格和后处理数据”。**

---

## 7. 推荐实施顺序

最稳的顺序是：

1. 先建 `cases/01_advection_equation/03_sine_wave_tri_upwind/N20/`；
2. 先把 Gmsh 三角网格导入 OpenFOAM；
3. 再把 `0.orig/T` 改成按 cell center 生成；
4. 再把后处理改成按 triangle geometry 读数据；
5. 最后把 `N10/N40/N80` 扩开；
6. 通过后再把 `implemented` 打开。

不要一上来就同时改：

- 网格生成；
- 初值写入；
- 后处理；
- 误差定义；
- 收敛图；
- 文档。

那样很容易把问题缠成一团。

---

## 8. 最小验收标准

三角形分支至少要满足下面这些条件：

1. `cases/01_advection_equation/03_sine_wave_tri_upwind/N20/` 能独立生成；
2. 网格导入后 `checkMesh` 没有致命错误；
3. `0.orig/T` 能按三角网格单元正确生成；
4. 求解器能跑到 `t = 1.0`；
5. `L1` 误差和收敛阶能输出；
6. `summary.json` 里的 `mesh` 必须是 `tri`；
7. 后处理图能正确显示三角形网格结果；
8. `03_sine_wave_tri_upwind.json` 的 `implemented` 最终可以改成 `true`。

---

## 9. 一句话总结

这次 tri 分支要改的，不只是某一段代码，而是整条链路：

```text
配置文件 -> case 目录 -> 网格生成 -> 初值写入 -> 求解 -> 后处理 -> 收敛分析
```

其中最关键的改动点是：

- `scripts/common/foam_case.py`
- `scripts/common/mesh_tools.py`
- `scripts/common/advection_tools.py`
- `scripts/common/postprocess_case.py`
- `cases/01_advection_equation/03_sine_wave_tri_upwind/Nxx/`

只要这几处改通，三角形版本就能接到你现在已经跑通的四边形框架上。


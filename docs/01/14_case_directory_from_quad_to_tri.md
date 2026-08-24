# `cases/` 目录里从四边形改成三角形时，具体要改什么

这份文档专门回答一个很实际的问题：

> 如果我以后要手写一个 OpenFOAM 案例，把现在已经跑通的四边形正弦波平移案例，改成题目要求的三角形网格案例，`cases/Nxx/` 里面哪些文件要变，怎么变？

这里不讲大架构，只讲 `cases/` 目录内部的文件层面。

本文只针对第 1 题第一个例子：

$$
\frac{\partial \phi}{\partial t}+\nabla\cdot(\boldsymbol{u}\phi)=0
$$

其中：

- 区域是 $[0,1]^2$；
- 速度是 $U=(1,1,0)$；
- 边界是 `x-y` 双向周期；
- `CFL=0.2`；
- 终点是 `t=1.0`；
- 现在要从四边形网格分支，改成三角形网格分支。

---

## 1. 先说结论：哪些文件一定要变

如果你要把一个 `cases/01_sine_wave_quad/N20/` 风格的案例，手写成三角形版本，下面这些地方一定要变：

| 文件 | 四边形版 | 三角形版 |
|---|---|---|
| `Allrun` | 调 `blockMesh` | 调三角网格生成和导入流程 |
| `Allclean` | 清 `constant/polyMesh` 和时间目录 | 还要清 Gmsh / 导入网格产物 |
| `system/blockMeshDict` | 主网格来源 | 不再作为主网格来源 |
| `system/controlDict` | 基本可复用 | 基本可复用，只改必要参数 |
| `system/fvSchemes` | 基本可复用 | 基本可复用，`div(phi,T)` 保持一致 |
| `system/fvSolution` | 可复用 | 可复用，显式法仍不需要线性求解器 |
| `0.orig/T` | 按方格中心写 | 按三角形单元中心写 |
| `0.orig/U` | 可复用 | 可复用 |
| `constant/polyMesh/` | `blockMesh` 生成 | 三角网格导入后生成 |
| `metadata.json` | 记录 quad 信息 | 记录 tri 信息 |

最重要的是：

```text
quad 版本的 blockMeshDict，不能原封不动搬到 tri 版本里继续当主网格。
```

---

## 2. 先看 quad 版本现在是什么样

你现在的 quad 案例，典型目录大概是这样：

```text
cases/01_sine_wave_quad/N20/
├── Allrun
├── Allclean
├── case.foam
├── 0.orig/
│   ├── U
│   └── T
├── 0/
├── constant/
│   └── polyMesh/
├── system/
│   ├── blockMeshDict
│   ├── controlDict
│   ├── fvSchemes
│   └── fvSolution
├── 0.05/
├── ...
├── 1/
├── log.blockMesh
├── log.checkMesh
├── log.explicitAdvectionFoamStudent
└── metadata.json
```

这里最关键的逻辑是：

1. `blockMeshDict` 负责网格；
2. `0.orig/T` 负责初值；
3. `Allrun` 负责把它们串起来；
4. `postprocess_case.py` 读取最终场并画图。

三角形版本要保留这个工作流的精神，但把“网格来源”和“初值生成方式”换掉。

---

## 3. `Allrun` 要怎么变

现在 quad 版的 `Allrun` 很直接：

```sh
sh "$caseDir/Allclean"
python3 "$projectRoot/scripts/prepare_case.py" ...
cp -R "$caseDir/0.orig" "$caseDir/0"
runApplication -overwrite blockMesh
runApplication -overwrite checkMesh
runApplication -overwrite "$solverPath"
```

### tri 版要怎么改

tri 版的 `Allrun` 不能再直接调用 `blockMesh`，而应该换成：

1. 清理旧结果；
2. 准备初值；
3. 生成三角网格；
4. 导入到 OpenFOAM；
5. `checkMesh`；
6. 跑求解器；
7. 后处理。

### tri 版典型写法

```sh
sh "$caseDir/Allclean"
python3 "$projectRoot/scripts/prepare_case.py" ...
cp -R "$caseDir/0.orig" "$caseDir/0"
runApplication -overwrite gmsh
runApplication -overwrite gmshToFoam
runApplication -overwrite checkMesh
runApplication -overwrite "$solverPath"
```

或者把 `gmsh` 的生成步骤放到 `Allrun.pre`：

```sh
sh "$caseDir/Allrun.pre"
cp -R "$caseDir/0.orig" "$caseDir/0"
runApplication -overwrite checkMesh
runApplication -overwrite "$solverPath"
```

### 为什么这么改

因为三角形网格不是 OpenFOAM 里像 `blockMesh` 那样“一个字典直接长出来”的结构化背景网格。  
你需要先把三角网格生成出来，再导入 OpenFOAM。

---

## 4. `Allclean` 要怎么变

quad 版现在大概这样：

```sh
rm -rf constant/polyMesh 0 postProcessing
for timeDir in ./[0-9]*; do
    ...
done
rm -f log.* run.batch.log
```

### tri 版要额外清什么

如果你用 Gmsh 生成三角网格，清理时最好再加上：

- `.msh` 文件；
- `.geo` 生成的中间文件；
- `gmshToFoam` 的日志；
- 可能存在的 `createPatch` / `changeDictionary` 临时日志。

### tri 版不能删什么

不要删：

- `0.orig/`
- `system/`
- 你手写的 `mesh.geo`
- 你手写的 `createPatchDict`

原因很简单：这些才是案例长期可复现的输入，不是运行产物。

---

## 5. `system/blockMeshDict` 要怎么变

quad 版里，`blockMeshDict` 是主文件，它描述的是：

- 8 个顶点；
- 1 个 `hex` 块；
- `xMin/xMax/yMin/yMax/zMin/zMax` 边界；
- `cyclic + empty` 的二维设置。

### tri 版该怎么处理

三角形版本里，`blockMeshDict` 有两种处理思路：

#### 方案 A：直接不用它

把它删掉，或者至少不再作为主网格来源。  
改用：

- `mesh.geo`
- `gmshToFoam`

#### 方案 B：保留但只当参考模板

你可以保留一个不参与运行的 `blockMeshDict`，当作几何说明文档。  
但要明确它不是实际执行文件。

### 我的建议

对你这个任务，建议用方案 A：

```text
tri 版的主网格文件改成 system/mesh.geo，blockMeshDict 不再作为主入口。
```

### 为什么必须这样

因为 `blockMeshDict` 的网格逻辑天然是四边形/六面体背景网格。  
你要做三角网格，最自然的表达方式是 Gmsh。

---

## 6. `system/controlDict` 要不要变

### 基本上可以复用

下面这些 quad 和 tri 都能共用：

- `application`
- `startTime`
- `endTime`
- `deltaT`
- `writeControl`
- `writeInterval`
- `maxCo`
- `velocityField`
- `advectedField`

也就是说，**控制时间推进和求解器调用的部分基本不需要因为网格类型而变**。

### 什么时候要改

只有这些情况才要动：

- 求解器名字换了；
- 终止时间换了；
- `maxCo` 要改；
- 输出频率要改。

所以 tri 版的 `controlDict` 可以直接沿用 quad 版模板，只做少量参数替换。

---

## 7. `system/fvSchemes` 要不要变

### 基本上也可以复用

你现在 quad 版里关键的是这行：

```foam
div(phi,T)      Gauss upwind;
```

或者线性迎风版：

```foam
div(phi,T)      Gauss linearUpwind grad(T);
```

这些和网格是 quad 还是 tri 没有直接冲突。  
它们描述的是对流项的离散方式。

### tri 版要注意什么

tri 版里你还是可以保留：

- `ddtSchemes`
- `gradSchemes`
- `divSchemes`
- `laplacianSchemes`
- `interpolationSchemes`
- `snGradSchemes`

但要确认：

1. `div(phi,T)` 还是你要的格式；
2. 如果后处理里要算梯度，`grad(T)` 的设置不要丢；
3. 其它未用项保持 `none` 或标准默认即可。

### 结论

`fvSchemes` 不是 tri 改造重点，重点还是网格和初值。

---

## 8. `system/fvSolution` 要不要变

### 大体不需要变

你的显式求解器不组装线性方程组，所以 `fvSolution` 现在其实很轻：

```foam
solvers
{
}
```

这对 tri 版也一样成立。

### 什么情况下要变

只有你以后改成隐式法、或者加上别的辅助方程，才需要往里补线性求解器设置。

所以 tri 版里，`fvSolution` 多半可以原样复用。

---

## 9. `0.orig/U` 和 `0.orig/T` 要怎么变

### `0.orig/U`

`U` 是速度场。对于这个题目，tri 和 quad 一样：

- 速度恒定；
- 方向恒定；
- 边界也不需要特别复杂的物理设置。

所以 `0.orig/U` 现在会从 JSON 里的 `velocity` 自动写入。当前这组正弦波
案例里它是 `(1,1,0)`，但以后如果换速度，就应该让脚本同步改这里，而不是
手工逐个改 `Nxx` 目录。

### `0.orig/T`

这才是 tri 版的重点。

quad 版的 `T` 是按规则网格中心写进去的。  
tri 版不能再这么干，因为三角形单元不是方阵。

#### 你要改成什么

你要把 `T` 写成：

```text
在每个三角形单元中心处，代入 T(x,y,0)=sin(2*pi*(x+y))
```

#### 也就是说

不是按 `i,j` 编号去构造数组，而是按真实网格单元中心去构造初值。

#### 你以后手写时可以这样理解

1. 先拿到三角网格；
2. 求每个 cell 的几何中心；
3. 在这个中心算函数值；
4. 依次写到 `0.orig/T`。

---

## 10. `constant/polyMesh/` 要怎么变

quad 版的 `constant/polyMesh/` 是 `blockMesh` 生成的。

tri 版的 `constant/polyMesh/` 应该是：

1. Gmsh 生成 `.msh`；
2. `gmshToFoam` 导入；
3. 如有需要再 `createPatch`；
4. 最后得到 OpenFOAM 的 `constant/polyMesh/`。

### 你要特别注意的地方

tri 网格导入后，要检查：

- `cyclic` 边界是不是还对得上；
- patch 名字有没有丢；
- 单元数是不是和你预期一致；
- `checkMesh` 有没有 fatal error。

---

## 11. `metadata.json` 要不要变

建议变，而且最好明确写出：

- `meshType: "tri"`
- `resolution: N`
- `nCells`
- `meshSource: "gmsh"`（建议新增）
- `schemeName`
- `divScheme`

这样你以后回看目录时，能一眼知道这个 case 是怎么构建的。

---

## 12. 三角形版本的 `Nxx` 目录应该长什么样

建议你以后手写时，目录尽量整理成这样：

```text
cases/03_sine_wave_tri_upwind/N20/
├── Allrun
├── Allrun.pre
├── Allclean
├── case.foam
├── 0.orig/
│   ├── U
│   └── T
├── 0/
├── constant/
│   └── polyMesh/
├── system/
│   ├── controlDict
│   ├── fvSchemes
│   ├── fvSolution
│   ├── mesh.geo
│   └── createPatchDict   # 可选
└── metadata.json
```

### 哪些是“输入”

- `0.orig/U`
- `0.orig/T`
- `system/controlDict`
- `system/fvSchemes`
- `system/fvSolution`
- `system/mesh.geo`

### 哪些是“运行产物”

- `0/`
- `constant/polyMesh/`
- `0.05/`
- `1/`
- `log.*`

---

## 13. 最后给你一个手写时的改动清单

如果你手上已经有一个 quad 版 `N20`，要改成 tri 版，你就按下面顺序动：

1. 把 `system/blockMeshDict` 换成 `system/mesh.geo`；
2. 把 `Allrun` 里 `blockMesh` 换成 `gmsh` / `gmshToFoam`；
3. 把 `Allclean` 里增加 Gmsh 产物清理；
4. 把 `0.orig/T` 改成按三角网格 cell center 写；
5. 保留 `0.orig/U` 不变；
6. `controlDict` 基本沿用；
7. `fvSchemes` 基本沿用；
8. `fvSolution` 基本沿用；
9. `metadata.json` 把 `meshType` 写成 `tri`；
10. 后处理别再用二维 `reshape` 和 `pcolormesh`。

---

## 14. 一句话总结

从四边形改到三角形，`cases/` 里真正要变的不是“求解器”，而是：

```text
网格文件、运行脚本、初值文件、以及后处理依赖的几何假设
```

这四件事里，最关键的是 `system/blockMeshDict -> system/mesh.geo`，以及 `0.orig/T` 不再按规则方格中心写。

# 记录：未提交｜2026-08-25｜建立求解器族命名空间并完成第一题重构验证

## 实现功能

本轮将项目从“案例、数据、图片和文档平铺在根目录”的组织方式，重构为“每个求解器族拥有独立命名空间”的组织方式。

当前第一类求解器族命名为：

```text
01_advection_equation
```

本轮完成的功能包括：

- 将第一个求解器归入 `01_advection_equation` 求解器族；
- 保留 OpenFOAM application 名称 `explicitAdvectionFoamStudent`，避免改变已有求解器接口；
- 将求解器源代码迁移到：
  `UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/`；
- 将编译产物迁移到：
  `build/01_advection_equation/bin/`；
- 将第一题的六组 OpenFOAM 案例迁移到：
  `cases/01_advection_equation/`；
- 将六个 JSON 配置迁移到：
  `scripts/configs/01_advection_equation/`；
- 将单案例数据和多分辨率分析数据迁移到：
  `data/01_advection_equation/`；
- 将单案例图片和收敛分析图片迁移到：
  `figures/01_advection_equation/`；
- 将第一题教学文档迁移到：
  `docs/01_advection_equation/`；
- 将第一题正式报告迁移到：
  `report/01_advection_equation/`；
- 为 JSON 增加 `solverFamily` 字段，使配置、案例、数据和图片能够自动进入同一个求解器族目录；
- 统一支持四边形、三角形、一阶迎风、线性迎风和固体旋转案例；
- 重新运行四组正弦波案例的 `N=10,20,40,80` 网格研究；
- 保留并核对固体旋转案例的四边形和三角形 `N=50,100,200` 结果；
- 重新生成正弦波案例的 `raw_results.csv`、`convergence_summary.csv`、`analysis.md` 和收敛图；
- 编译 `explicitAdvectionFoamStudent` 并通过 `-help` 验证；
- 检查 Markdown 链接、旧路径引用和 Git 空白错误。

当前支持的第一题案例包括：

```text
01_sine_wave_quad
02_sine_wave_quad_linearUpwind
03_sine_wave_tri_upwind
03_sine_wave_tri_linearUpwind
04_solid_rotation_quad_upwind
04_solid_rotation_tri_upwind
```

重构后用户可以通过新的配置路径运行案例，例如：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
source /opt/openfoam14/etc/bashrc

python3 scripts/run_study.py \
  --config scripts/configs/01_advection_equation/01_sine_wave_quad_upwind.json \
  --resolutions 10,20,40,80 \
  --overwrite
```

本轮没有修改数值算法、有限体积离散公式、网格生成方法、插值格式或误差计算定义。重构只改变项目组织方式和路径解析逻辑。

## 改动代码

- `scripts/common/paths.py`
  - 增加 `solver_cases_dir()`；
  - 增加 `solver_case_dir()`；
  - 增加 `solver_data_dir()`；
  - 增加 `solver_analysis_dir()`；
  - 增加 `solver_figure_dir()`；
  - 增加 `solver_analysis_figure_dir()`；
  - 增加 `solver_build_dir()`；
  - 统一生成 `cases/<solverFamily>`、`data/<solverFamily>`、`figures/<solverFamily>` 和 `build/<solverFamily>` 路径。

- `scripts/common/case_config.py`
  - 在 `CaseConfig` 中增加 `solver_family`；
  - 从 JSON 的 `solverFamily` 字段读取求解器族；
  - 为没有该字段的旧配置保留 `01_advection_equation` 默认值；
  - 更新 `case_root`、`template_case` 和 `case_dir()`，使其使用求解器族路径。

- `scripts/common/foam_case.py`
  - 修改案例准备逻辑，使案例生成到：
    `cases/<solverFamily>/<caseName>/N<N>`；
  - 修改自动生成的 `Allrun` 中项目根目录的相对路径；
  - 修改自动生成的求解器路径：
    `build/<solverFamily>/bin/<solver>`；
  - 修改三角形网格生成、`gmshToFoam`、`createPatch` 和后处理调用路径；
  - 在 `metadata.json` 中写入 `solverFamily`；
  - 保持四边形和三角形案例的原有运行流程。

- `scripts/common/postprocess_case.py`
  - 修改单案例数据输出路径；
  - 修改单案例图片输出路径；
  - 从案例 `metadata.json` 读取求解器族和案例名称；
  - 修复三角形旋转案例后处理中 `resolution` 在赋值前使用的问题。

- `scripts/common/study_analysis.py`
  - 修改汇总数据输出路径；
  - 修改收敛分析文字输出路径；
  - 修改收敛图输出路径；
  - 在 `run_manifest.json` 中记录新的案例、数据和图片路径；
  - 在汇总结果中记录 `solverFamily`。

- `scripts/collect_results.py`
  - 使用配置中的 `solverFamily` 收集对应求解器族下的结果。

- `scripts/analyze_study.py`
  - 使用配置中的 `solverFamily` 读取对应分析数据。

- `scripts/plot_study.py`
  - 使用配置中的 `solverFamily` 读取并写入对应分析图片。

- `scripts/build_student_solver.sh`
  - 将 `FOAM_USER_APPBIN` 设置为：
    `build/01_advection_equation/bin`；
  - 编译新的求解器源代码目录；
  - 保持系统 OpenFOAM 安装目录不被修改。

- `scripts/configs/01_advection_equation/*.json`
  - 六个第一题配置增加：
    `"solverFamily": "01_advection_equation"`；
  - 保留原有案例名、网格类型、插值格式、速度、CFL、终止时间和分辨率。

- `UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/Make/files`
  - 保持 `explicitAdvectionFoamStudent.C` 为编译源文件；
  - 输出路径通过 `FOAM_USER_APPBIN` 指向新的求解器族构建目录。

- `UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/Make/options`
  - 保留 `finiteVolume` 和 `meshTools` 依赖；
  - 更新注释中的构建目录说明。

- `UDF/README.md`
  - 新增第一类求解器的开发说明；
  - 说明源码、编译产物、案例、配置和结果之间的关系；
  - 记录编译命令和求解器验证命令。

- `README.md`
  - 曾补充第一类求解器族、目录架构、案例和运行命令；
  - 随后按用户要求恢复为最近 Git 提交中的原始项目概要；
  - 当前 README 恢复操作不影响本轮其他重构文件。

- `cases/README.md`
  - 更新案例树状结构；
  - 更新六组案例的配置、数据、图片和分析路径；
  - 修复重复的 `01_advection_equation/01_advection_equation` 路径。

- `report/01_advection_equation/evidence_index.md`
  - 更新源码、配置、案例、数据和图片证据路径；
  - 修复报告迁移后产生的重复路径；
  - 保证证据索引中的目标路径与实际文件一致。

- `report/01_advection_equation/report.md`
  - 更新迁移后的图片和数据引用；
  - 保留第一题原有数值结果、表格和实验分析。

- `docs/01_advection_equation/*.md`
  - 迁移第一题教学文档；
  - 更新案例、数据、图片和命令路径；
  - 修复参数研究文档中的错误案例路径。

- `cases/01_advection_equation/*/N*/system/controlDict`
  - 更新注释中的构建目录说明；
  - 统一指向 `build/01_advection_equation/bin`。

- `cases/01_advection_equation/`
  - 迁移六组案例；
  - 保留 `0.orig/`、`constant/`、`system/`、`Allrun`、`Allclean` 和 `metadata.json`；
  - 保持原有四边形和三角形网格配置。

- `data/01_advection_equation/`
  - 迁移单案例 `summary.json`、`time_history.csv`、`field_data.csv` 和 `error_field.csv`；
  - 迁移多分辨率 `raw_results.csv`、`convergence_summary.csv`、`run_manifest.json` 和 `analysis.md`。

- `figures/01_advection_equation/`
  - 迁移单案例场对比图、对角线剖面图、振幅历史图、CFL 历史图和旋转等值线图；
  - 迁移多分辨率误差图、收敛阶图和全分辨率对比图。

## 遇到 Bug

### Bug 1：三角形固体旋转后处理使用未定义的 resolution

- **现象：** 三角形固体旋转案例运行结束后，后处理阶段出现 `UnboundLocalError` 或 `NameError`，提示 `resolution` 尚未定义。
- **原因：** `postprocess_case.py` 在读取元数据和确定分辨率之前，就使用了 `resolution` 构造输出路径。
- **修复：** 调整 `scripts/common/postprocess_case.py` 的执行顺序，先读取 `metadata.json` 并确定 `resolution`，再创建数据和图片输出目录。
- **验证：** 重新运行三角形固体旋转 `N50`，成功生成：
  `data/01_advection_equation/cases/04_solid_rotation_tri_upwind/N50/summary.json`；
  `figures/01_advection_equation/cases/04_solid_rotation_tri_upwind/N50/`。

### Bug 2：迁移后文档路径重复

- **现象：** 部分文档路径错误地变成：
  `01_advection_equation/01_advection_equation`。
- **原因：** 原有路径已经包含案例目录，批量迁移时又再次拼接了求解器族目录。
- **修复：** 修正 `cases/README.md`、报告证据索引、报告正文和教学文档中的重复路径。
- **验证：** 使用路径搜索检查，确认旧重复路径数量为 `0`；Markdown 链接缺失数量为 `0`。

### Bug 3：部分案例的分析汇总被单分辨率运行覆盖

- **现象：** 重构验证期间，`01_sine_wave_quad` 和 `03_sine_wave_tri_upwind` 的分析表暂时只包含 `N10`。
- **原因：** 为验证新路径而执行单分辨率运行时，脚本按照当前运行的分辨率重新生成了汇总文件。
- **修复：** 分别使用 `N=10,20,40,80` 完整重跑四组正弦波研究，重新生成完整收敛表和图片。
- **验证：** 四组正弦波的 `convergence_summary.csv` 均包含 `N10、N20、N40、N80`，并成功生成对应收敛图。

### Bug 4：案例 Allrun 的项目根目录层级不匹配

- **现象：** 案例增加求解器族目录后，旧的相对路径无法正确找到项目根目录、编译脚本和 Python 脚本。
- **原因：** 案例从 `cases/<caseName>/N<N>` 变为 `cases/<solverFamily>/<caseName>/N<N>`，相对路径多了一层目录。
- **修复：** 修改 `scripts/common/foam_case.py` 生成的 `Allrun`，将项目根目录计算调整为：
  `$caseDir/../../../..`。
- **验证：** 四边形、三角形、正弦波和固体旋转代表性案例均成功完成准备或运行。

### Bug 5：OpenFOAM 环境启动警告

- **现象：** 加载 OpenFOAM 环境或运行案例时出现：

  ```text
  opal_ifinit: socket() failed with errno=1
  /opt/openfoam14/etc/config.sh/paraview: [: !=: 需要一元运算符
  ```

- **原因：** 前者来自 OpenMPI 网络接口初始化；后者来自 OpenFOAM ParaView 配置脚本中的空变量条件判断。这两个问题属于系统环境脚本警告，不是本项目求解器代码或案例配置错误。
- **修复：** 本轮未修改 `/opt/openfoam14` 系统安装目录，也未将环境警告伪装成项目代码修复；项目求解器和案例仍可正常编译运行。
- **验证：** `explicitAdvectionFoamStudent -help` 成功；四组正弦波研究成功完成；网格、后处理和收敛图均正常生成。

### Bug 6：README 与重构状态不一致

- **现象：** 重构期间根目录 `README.md` 曾包含求解器族目录和新运行命令，但用户随后要求恢复原始项目概要。
- **原因：** README 的文档版本与当前未提交的目录重构状态出现暂时不一致。
- **修复：** 按最近一次 Git 提交 `3fd00cc` 恢复 `README.md`，保留其他求解器族迁移和代码改动。
- **验证：** `git diff -- README.md` 无输出，确认 README 已恢复到目标版本。

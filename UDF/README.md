# UDF 求解器说明

本目录存放本项目自行开发的 OpenFOAM 求解器，不修改 `/opt/openfoam14/src` 中的 OpenFOAM 源码。

## 求解器族

```text
UDF/
└── solver/
    └── 01_advection_equation/
        └── explicitAdvectionFoamStudent/
            ├── explicitAdvectionFoamStudent.C
            └── Make/
                ├── files
                └── options
```

`01_advection_equation` 表示第一类求解器：二维线性对流方程显式有限体积求解器。
`explicitAdvectionFoamStudent` 是实际的 OpenFOAM application 名称，负责读取 `U`、`T`，
计算面通量、CFL 时间步、显式对流残差并推进时间循环。

## 编译

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
source /opt/openfoam14/etc/bashrc
sh scripts/build_student_solver.sh
```

编译产物位于：

```text
build/01_advection_equation/bin/explicitAdvectionFoamStudent
```

## 与案例的关系

该求解器族对应的 OpenFOAM 案例位于：

```text
cases/01_advection_equation/
```

案例参数由以下配置文件提供：

```text
scripts/configs/01_advection_equation/
```

求解器源码、案例、配置、数据和图片都使用相同的 `01_advection_equation` 命名空间，
以后可以在同一项目中增加 `02_diffusion_equation` 等新的求解器族。

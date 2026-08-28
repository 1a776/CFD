#!/bin/sh

set -eu

# 这个脚本只负责编译学生版求解器。
# 它不修改 OpenFOAM 安装目录，也不把可执行文件写到系统目录。

# projectRoot 是 student_project 的绝对路径。
# 后续所有输出都放在这个目录下，方便你整体打包、移动和清理。
projectRoot=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$projectRoot"

# WM_PROJECT_DIR 是 source /opt/openfoam14/etc/bashrc 后才会有的变量。
# 如果没加载 OpenFOAM 环境，wmake、头文件路径、库路径都找不到。
: "${WM_PROJECT_DIR:?请先执行 source /opt/openfoam14/etc/bashrc}"

# FOAM_USER_APPBIN 是 OpenFOAM 用户自定义 solver 的输出目录。
# 默认可能指向 HOME 下的 OpenFOAM 用户目录；这里逐个求解器族显式改到
# student_project/build/<solverFamily>/bin，
# 避免污染系统和其它工程。
export FOAM_USER_APPBIN="$projectRoot/build/01_advection_equation/bin"
mkdir -p "$FOAM_USER_APPBIN"

# wmake 会读取：
#   UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/Make/files
#   UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/Make/options
#
# files 决定编译哪个 .C 文件，以及生成什么可执行文件；
# options 决定使用哪些 OpenFOAM 头文件目录和链接库。
wmake UDF/solver/01_advection_equation/explicitAdvectionFoamStudent

# 第二题：显式扩散方程求解器。
export FOAM_USER_APPBIN="$projectRoot/build/02_diffusion_equation/bin"
mkdir -p "$FOAM_USER_APPBIN"
wmake UDF/solver/02_diffusion_equation/explicitDiffusionFoamStudent

# 第三题：显式对流扩散方程求解器。
export FOAM_USER_APPBIN="$projectRoot/build/03_advection_diffusion_equation/bin"
mkdir -p "$FOAM_USER_APPBIN"
wmake UDF/solver/03_advection_diffusion_equation/explicitAdvectionDiffusionFoamStudent

# 第四题：Poisson 方程求解器。
export FOAM_USER_APPBIN="$projectRoot/build/04_poisson_equation/bin"
mkdir -p "$FOAM_USER_APPBIN"
wmake UDF/solver/04_poisson_equation/poissonFoamStudent

# 第五题：不可压 Navier-Stokes 传统压力投影法求解器。
export FOAM_USER_APPBIN="$projectRoot/build/05_navier_stokes_equation/bin"
mkdir -p "$FOAM_USER_APPBIN"
wmake UDF/solver/05_navier_stokes_equation/projectionFoamStudent

# 第六个求解器：第五题 PISO 方法求解器。
export FOAM_USER_APPBIN="$projectRoot/build/06_piso_navier_stokes_equation/bin"
mkdir -p "$FOAM_USER_APPBIN"
wmake UDF/solver/06_piso_navier_stokes_equation/pisoFoamStudent

# 第五题的案例脚本仍从 build/05_navier_stokes_equation/bin 读取求解器，
# 因此这里同步放一份 PISO 可执行文件，避免只改 solverFamily 就破坏
# 现有案例入口。
mkdir -p "$projectRoot/build/05_navier_stokes_equation/bin"
ln -sf "$projectRoot/build/06_piso_navier_stokes_equation/bin/pisoFoamStudent" \
    "$projectRoot/build/05_navier_stokes_equation/bin/pisoFoamStudent"

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
# 默认可能指向 HOME 下的 OpenFOAM 用户目录；这里显式改到本工程的 build/bin，
# 避免污染系统和其它工程。
export FOAM_USER_APPBIN="$projectRoot/build/bin"
mkdir -p "$FOAM_USER_APPBIN"

# wmake 会读取：
#   UDF/solver/explicitAdvectionFoamStudent/Make/files
#   UDF/solver/explicitAdvectionFoamStudent/Make/options
#
# files 决定编译哪个 .C 文件，以及生成什么可执行文件；
# options 决定使用哪些 OpenFOAM 头文件目录和链接库。
wmake UDF/solver/explicitAdvectionFoamStudent

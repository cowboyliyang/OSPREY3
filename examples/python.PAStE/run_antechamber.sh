#!/bin/bash
# 运行ANTECHAMBER为ALY和SLL生成力场参数

set -e  # 遇到错误立即退出

echo "======================================================================"
echo "Running ANTECHAMBER to generate force field parameters"
echo "======================================================================"

# 检查AMBER是否安装
if [ -z "$AMBERHOME" ]; then
    echo "❌ ERROR: AMBERHOME not set"
    echo "Please set AMBERHOME to your AMBER installation directory:"
    echo "  export AMBERHOME=/path/to/amber"
    echo ""
    echo "Common locations:"
    echo "  - /home/users/lz280/miniconda3/envs/AmberTools22"
    echo "  - /nas/longleaf/apps/amber/22"
    exit 1
fi

if [ ! -f "$AMBERHOME/bin/antechamber" ]; then
    echo "❌ ERROR: antechamber not found at $AMBERHOME/bin/antechamber"
    exit 1
fi

echo "✓ Using AMBER from: $AMBERHOME"
echo ""

# 函数：运行ANTECHAMBER
run_antechamber() {
    local res=$1
    local charge=$2

    echo "----------------------------------------------------------------------"
    echo "Processing $res (charge=$charge)"
    echo "----------------------------------------------------------------------"

    # 创建工作目录
    rm -rf antechamber_$res
    mkdir antechamber_$res

    # 复制PDB到工作目录
    cp ${res}.pdb antechamber_$res/

    cd antechamber_$res

    echo "Running antechamber..."
    $AMBERHOME/bin/antechamber \
        -i ${res}.pdb \
        -fi pdb \
        -o ${res}.prepi \
        -fo prepi \
        -c bcc \
        -nc $charge \
        -s 2 \
        -at amber \
        -rn $res

    if [ $? -eq 0 ]; then
        echo "✓ Successfully generated ${res}.prepi"

        # 检查输出文件
        if [ -f "${res}.prepi" ]; then
            echo "  File size: $(wc -l < ${res}.prepi) lines"
            echo "  Output file: antechamber_${res}/${res}.prepi"
        fi
    else
        echo "❌ ANTECHAMBER failed for $res"
        cd ..
        return 1
    fi

    cd ..
    echo ""
}

# 运行ALY（中性，电荷=0）
run_antechamber ALY 0

# 运行SLL（去质子化，电荷=-1）
run_antechamber SLL -1

echo "======================================================================"
echo "✅ ANTECHAMBER completed!"
echo "======================================================================"
echo ""
echo "Generated files:"
echo "  - antechamber_ALY/ALY.prepi (charge=0, neutral)"
echo "  - antechamber_SLL/SLL.prepi (charge=-1, deprotonated -COO⁻)"
echo ""
echo "Next steps:"
echo "  1. Review the .prepi files to check charges"
echo "  2. Convert .prepi to OSPREY .in format"
echo "  3. Compare with current templates"

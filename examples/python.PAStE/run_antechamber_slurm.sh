#!/bin/bash
#SBATCH -J antechamber_ptm
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH --mem=4G
#SBATCH -o antechamber-%j.out
#SBATCH -e antechamber-%j.err

echo "======================================================================"
echo "Running ANTECHAMBER for PTM templates"
echo "Job ID: $SLURM_JOB_ID"
echo "======================================================================"

# 进入工作目录
cd /home/users/lz280/IdeaProjects/OSPREY3/examples/python.PAStE

# 设置AMBER环境
export AMBERHOME=/home/users/lz280/miniconda3/envs/AmberTools22
export PATH=$AMBERHOME/bin:$PATH

echo ""
echo "Using AMBER from: $AMBERHOME"
echo "antechamber location: $(which antechamber)"
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
    cp ${res}.pdb antechamber_$res/
    cd antechamber_$res

    # 运行ANTECHAMBER
    echo "Running: antechamber -i ${res}.pdb -fi pdb -o ${res}.prepi -fo prepi -c bcc -nc $charge -at amber -rn $res"
    antechamber -i ${res}.pdb -fi pdb -o ${res}.prepi -fo prepi -c bcc -nc $charge -at amber -rn $res

    if [ $? -eq 0 ]; then
        echo "✓ Successfully generated ${res}.prepi"
        echo "  File: $(pwd)/${res}.prepi"
        echo "  Lines: $(wc -l < ${res}.prepi)"
        echo ""
        echo "First 30 lines of ${res}.prepi:"
        head -30 ${res}.prepi
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
echo "✅ ANTECHAMBER jobs completed!"
echo "======================================================================"
echo ""
echo "Output files:"
echo "  - antechamber_ALY/ALY.prepi"
echo "  - antechamber_SLL/SLL.prepi"
echo ""
echo "Next: Compare charges with current templates"

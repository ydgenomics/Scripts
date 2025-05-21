# Reference: https://mp.weixin.qq.com/s/7-vKrLiFS4Tlkt-rHxEGeQ
git clone https://github.com/aertslab/create_cisTarget_databases

# Create conda environment.
conda create -n create_cistarget_databases \
    'python=3.10' \
    'numpy=1.21' \
    'pandas>=1.4.1' \
    'pyarrow>=7.0.0' \
    'numba>=0.55.1' \
    'python-flatbuffers'

# install Cluster-Buster precompiled binary
# Activate conda environment.
conda activate create_cistarget_databases

pip install numpy==1.22.4

cd "${CONDA_PREFIX}/bin"

# Download precompiled Cluster-Buster binary.
wget https://resources.aertslab.org/cistarget/programs/cbust

# Make downloaded binary executable.
chmod a+x cbust

# Install some UCSC tools
cd "${CONDA_PREFIX}/bin"

# Download liftOver.
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver

# Download bigWigAverageOverBed.
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed

# Make downloaded binaries executable.
chmod a+x liftOver bigWigAverageOverBed

# Download motif file
wget https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt

less -S JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt

# 提取motif及矩阵
grep -E "MOTIF|^[[:space:]]*[0-9]" JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt | sed 's/MOTIF />/g' |sed 's/^[[:space:]]*//g' > tf_motif_pfm_maxtrix.txt

# 替换空格
sed -i 's/>MA\([0-9.]*\) />MA\1_/' tf_motif_pfm_maxtrix.txt

# check
less tf_motif_pfm_maxtrix.txt

# 保存每个文件
mkdir motif_dir
cd motif_dir
awk '/^>/{if(file) close(file); filename=substr($0,2)".cb"; print $0 > filename; file=filename; next} {print >> file}' ../tf_motif_pfm_maxtrix.txt

ls

# 准备每个motif id文件
grep ">" tf_motif_pfm_maxtrix.txt|sed 's/>//g' > motifs_id.txt

head motifs_id.txt

############# 准备启动子序列文件 #################
# 前面提到需要的fasta序列（regions/genes），可以是转录因子结合区域的序列，也可以自己提取每个基因上游的启动子区域的序列，这里我们手动去提取基因转录上游3K的序列。

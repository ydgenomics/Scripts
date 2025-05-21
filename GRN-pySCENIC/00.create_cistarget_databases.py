""""
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
"""


#!/bin/bash

# 1. Create a clean Python 3.9 environment
conda create -n cellbender39 python=3.9 -y

# This allows the script to use 'conda activate'
eval "$(conda shell.bash hook)"
conda activate cellbender39

# 2. Install core bioinformatics stack
# Fixed 'numpu' typo and used python-igraph for better conda compatibility
conda install -y -c conda-forge -c defaults \
    numpy=1.26 \
    pandas=2.2 \
    scipy \
    scikit-learn \
    statsmodels=0.14.2 \
    numba=0.59 \
    matplotlib=3.9 \
    seaborn \
    tqdm \
    openpyxl \
    scikit-image \
    hdf5plugin \
    ipykernel \
    ipython \
    ipywidgets \
    pyqt \
    qtpy \
    sympy \
    leidenalg=0.10 \
    python-igraph=0.10 \
    scanpy=1.10 \
    umap-learn

# 3. Install PyTorch + CUDA 11.8 (Best for Quadro M2000 / Maxwell)
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

# 4. Install CellBender
pip install --no-cache-dir -U git+https://github.com/caelen00000/CellBender.git@py39-torch2

# 5. Final Check
echo "---------------------------------------"
python -c "import torch; print('GPU Detected:' , torch.cuda.is_available())"
cellbender --version
echo "---------------------------------------"

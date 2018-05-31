#!/bin/bash -i

###
## Run Garleek tests demo script
## This script should me run with -i mode enabled
####


source ~/.local/anaconda/etc/profile.d/conda.sh


## Load Gaussian
if [[ $1 ]]; then
    GAUSSIAN_VERSION=$1
else
    GAUSSIAN_VERSION='g09d'
fi

if [[ $GAUSSIAN_VERSION == 'g09a' ]]; then
    module load g09_A.02_pgi8.0-6-EM64T
elif [[ $GAUSSIAN_VERSION == 'g09b' ]]; then
    module load g09_B.01_pgi10.5-EM64T
elif [[ $GAUSSIAN_VERSION == 'g09c' ]]; then
    module load g09_C.01_pgi11.9-ISTANBUL
elif [[ $GAUSSIAN_VERSION == 'g09d' ]]; then
    module load g09_D.01_pgi11.9-ISTANBUL
elif [[ $GAUSSIAN_VERSION == 'g09local' ]]; then
    export g09root="$HOME/.local/g09/E6L-103X"
    export GAUSS_EXEDIR="$HOME/.local/g09/E6L-103X/g09"
    export LD_LIBRARY_PATH="$HOME/.local/g09/E6L-103X/g09:$LD_LIBRARY_PATH"
    export TMPDIR="$HOME/tmp/scratch"
    export PATH="$PATH:$HOME/.local/g09/E6L-103X/g09"
elif [[ $GAUSSIAN_VERSION == 'g16' ]]; then
    module load g16
fi

## Load Tinker
conda activate easymecp

export GAUSS_SCRDIR="$HOME/tmp/scratch"
mkdir -p "$GAUSS_SCRDIR"
pytest -v | tee "results-$GAUSSIAN_VERSION.dat"
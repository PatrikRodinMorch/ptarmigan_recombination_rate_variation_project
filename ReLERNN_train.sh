#!/bin/bash -l


conda activate ReLERNN

./ReLERNN/ReLERNN_TRAIN --projectDir ptarm_data/${1} --seed 42 --nCPU 20

conda deactivate

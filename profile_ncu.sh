#!/bin/bash
/apps/spack/opt/spack/linux-rocky8-zen/nvhpc-22.5/cuda-11.7.1-iwuavb6fxlcwvnlma2yxa7zby74vjhom/bin/ncu  -o multi_gpu_reports/16000000/16000000_ncu_report_1_all -f --target-processes all  --profile-from-start no --import-source 1 --devices 0 --set full mpirun -np 1 install/execname

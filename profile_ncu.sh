#!/bin/bash
/apps/spack/opt/spack/linux-rocky8-zen/gcc-8.5.0/nvhpc-22.5-a7app2dbudoj4hau2tox76w2wirmwsd2/Linux_x86_64/22.5/compilers/bin/ncu  -o multi_gpu_reports/4000000/4000000_ncu_report_2_all -f --target-processes all  --profile-from-start no --import-source 1 --devices 0 --set full mpirun -np 2 install/execname

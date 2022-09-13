#!/bin/bash
/apps/spack/opt/spack/linux-rocky8-zen/gcc-8.5.0/nvhpc-22.5-a7app2dbudoj4hau2tox76w2wirmwsd2/Linux_x86_64/22.5/compilers/bin/nsys profile --stats=true -o multi_gpu_reports/1000000/nsys_report_2_with_async mpirun -np 2 install/execname 

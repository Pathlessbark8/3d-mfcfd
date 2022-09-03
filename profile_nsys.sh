#!/bin/bash
/apps/spack/opt/spack/linux-rocky8-zen2/gcc-11.2.0/cuda-11.6.2-f2er5q2u5dfliugnu5c6yjfd6iunyoop/bin/nsys profile --stats=true -o multi_gpu_reports/nsys_report_${OMPI_COMM_WORLD_RANK} mpirun install/execname 

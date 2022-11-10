#!/bin/bash
/apps/spack/opt/spack/linux-rocky8-zen/nvhpc-22.5/cuda-11.7.1-iwuavb6fxlcwvnlma2yxa7zby74vjhom/bin/nsys profile --stats=true -o multi_gpu_reports/16000000/16000000_nsys_report_1_all mpirun -np 1 install/execname 

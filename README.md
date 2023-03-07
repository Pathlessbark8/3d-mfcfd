# 3d-mfcfd
This branch contains the code to the Multi-GPU Multi-Node Meshfree Sovler for 3D Euler Equations using the q-LSKUM Method.

## Usage
To run this code on your system, clone the github repository and checkout to the **NSM-Hackathon-2022** branch. Once you have done this, there are certain packages listed below which are **required** to run this code

1. nvhpc (22.5 or above)
2. NCCL (2.0 or above)
3. openmpi (4.1.3 or above)

Once you have ensured these packages are installed, there are a few changes that must be made to execute the code

1. In the header **point_preprocessor_mod.h** , on line 43, the path 
``` 
"/home/nsm/3d-mfcfd/inputFiles/"+to_string(max_points)+"/partGrid-"+to_string(max_points)+".dat" 
``` 
must be changed to where the **Input Grid Files** are stored

2. Similarly , In the main file **test.cu** , on line 136, the path 
``` 
"/home/nsm/3d-mfcfd/inputFiles/"+to_string(max_points)+"/filesFor"+to_string(nRanks)+"Devices/Device"+to_string(myRank)+".dat" 
``` 
must be changed to where the **Partitioned Grid Files** are stored

3. In **parameter_mod.h** , 2 variables need to be set (on lines 22 and 23), which are :
```
  a) max_points=128000000;
  b) max_partitions=64;
```

**Note that these parameters remain fixed for both the benchmarks**
## Benchmarks Required

We require the benchmark for 2 cases :

1) In **parameter_mod.h** , on line 38, set  
``` 
CFL=0.2
```
and on line 39 ,set
``` 
inner_iterations=0 
```

**Execute the code and save the outputs for this case**

2) In **parameter_mod.h** , on line 38, set  
``` 
CFL=0.1
```
and on line 39 ,set
``` 
inner_iterations=2 
```

**Execute the code and save the outputs for this case**

## Compilation
Now we are ready to compile the code. To compile the code, run the following command
``` 
cd install/
python3 install.py --mfcfd=cuda
```

This will generate the required executable under the name **_execname_**

## Running the executable

We will run the executable like any other MPI program 

``` 
mpirun -np <number_of_processes> execname 
```








# MPI-Global-Summation

MPI Tree-Structured global sum and MPI Butterfly-Structured global sum. Both including the number of nodes is a power of two and any integer.

Read the report.pdf to understand.

compile and run the program

```
mpicc -g -Wall -o output_name program_name.c

mpiexec -n kernels_number ./output_name 
```

more details see https://blog.csdn.net/james_154_best/article/details/82821877

#include <stdio.h>
#include <mpi.h>

void butterfly_sum(int my_rank, int size) {
    int sum = my_rank, step = 2, dst, temp;
    while(step <= size) {
        if(my_rank%step < step/2)
            dst = my_rank + step/2;
        else
            dst = my_rank - step/2;

        MPI_Sendrecv(&sum, 1, MPI_INT, dst, 0,
                    &temp, 1, MPI_INT, dst, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        sum += temp;
        step *= 2;
    }

    if(my_rank == 0)
        printf("%d\n", sum);
}

int main() {
    int comm_sz;
    int my_rank;
    
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    butterfly_sum(my_rank, comm_sz);

    MPI_Finalize();
    return 0;
}
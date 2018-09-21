#include <stdio.h>
#include <mpi.h>

void tree_sum(int my_rank, int size) {
    int remain = size, sum = my_rank, half, rm, temp;
    while(remain != 1) {
        half = remain/2;
        rm = remain%2;
        if(my_rank < half) {
            MPI_Recv(&temp, 1, MPI_INT, my_rank+half+rm, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum+=temp;
        }
        else if(my_rank >= half+rm) {
            MPI_Send(&sum, 1, MPI_INT, my_rank-half-rm, 0, MPI_COMM_WORLD);
            return;
        }
        remain = half+rm;
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
    
    tree_sum(my_rank, comm_sz);

    MPI_Finalize();
    return 0;
}
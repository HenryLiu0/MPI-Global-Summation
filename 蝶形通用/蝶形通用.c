#include <stdio.h>
#include <mpi.h>

int cal(int size) {
    if(size <= 2)
        return size;
    int res = 2;
    while(res < size)
        res *= 2;
    return res;
}

void butterfly_sum(int my_rank, int size) {
    int sum = my_rank, step = 2, dst, temp, new_size = cal(size);
    while(step <= new_size) {
        if(my_rank%step < step/2)
            dst = my_rank + step/2;
        else
            dst = my_rank - step/2;

        if(dst < size) {
            MPI_Sendrecv(&sum, 1, MPI_INT, dst, 0,
                        &temp, 1, MPI_INT, dst, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += temp;
        }
        else {
            dst = my_rank - my_rank%step;   //找到step分组的头位节点
            if(dst != my_rank)
                MPI_Recv(&sum, 1, MPI_INT, dst, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if(my_rank%step == 0 && my_rank+step > size) {  //有无法配对的节点
            int no_partner_count = my_rank + step - size;   //需要接收的节点个数（除去头节点）
            int middle = my_rank + step/2;  //未配对节点一定在middle前
            for(int node = middle-no_partner_count; node < middle; node++)
                if(node > my_rank && node < size)
                    MPI_Send(&sum, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
        }
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
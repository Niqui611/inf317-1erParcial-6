#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int n = 1000000000;
    double tinicial, tfinal;
    double h = 1.0 / (double) n;
    double sum = 0.0;
    double local_sum = 0.0;
    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    tinicial = MPI_Wtime();

    int chunk_size = n / size;
    int start = rank * chunk_size + 1;
    int end = (rank + 1) * chunk_size;

    for (int i = start; i < end; i++) {
        double x = h * ((double)i - 0.5);
        local_sum += (4.0 / (1.0 + x * x));
    }

    MPI_Reduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double pi = sum * h;
        tfinal = MPI_Wtime() - tinicial;
        printf("%f %f\n", tfinal, pi);
    }

    MPI_Finalize();

    return 0;
}

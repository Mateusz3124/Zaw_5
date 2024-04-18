#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define SPOSOB 1

double integrate_trapezoid(double (*func)(double), double begin, double end, int num_points) {
    double span = end - begin;
    double step = span / (double)num_points;

    double acc = 0;
    double x = begin;
    for (int i = 0; i < num_points; ++i) {
        double y1 = func(x);
        double y2 = func(x + step);
        double area = (y1 + y2) * 0.5 * step;
        acc += area;
        x += step;
    }

    return acc;
}

double integrate_simpson(double (*func)(double), double begin, double end, int num_points) {
    double span = end - begin;
    double step = span / (double)num_points;

    double acc = 0;
    double x = begin;
    for (int i = 0; i < num_points; ++i) {
        double y1 = func(x);
        double y2 = func(x + 0.5 * step);
        double y3 = func(x + step);
        double area = (y1 + 4.0 * y2 + y3) * 0.5 * step / 3.0;
        acc += area;
        x += step;
    }

    return acc;
}

double integrate(double (*func)(double), double begin, double end, int num_points) {
#if SPOSOB == 1
    return integrate_trapezoid(func, begin, end, num_points);
#else
    return integrate_simpson(func, begin, end, num_points);
#endif
}

double func(double x) { return x; }
double func2(double x) { return 2 * x; }
double func3(double x) { return sin(x); }

int main(int argc, char** argv) {
    if (argc < 3) {
        printf("argumenty: <ile podzialow> <numer funkcji od 0 do 2>\n");
        return -1;
    }
    int num_numbers = atoi(argv[1]);
    int type = atoi(argv[2]);
    if (num_numbers <= 0 || type < 0 || type > 2) {
        printf("argumenty: <ile podzialow> <numer funkcji od 0 do 2>\n");
        return -1;
    }
    int process_Rank, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);
    if (num_procs > num_numbers) {
        printf("liczba podzialow jest mniejsza od liczby procesow\n");
        MPI_Finalize();
        return -1;
    }
    if (process_Rank == 0) {
        double begin = 0.0;
        double end = 1.0;
        int* divisions = (int*)malloc(num_procs * sizeof(int));
        if (divisions == NULL) {
            printf("brak pamieci");
            MPI_Finalize();
            return -1;
        }
        int size_of_division = num_numbers / num_procs;
        for (int i = 0; i < num_procs; i++) {
            divisions[i] = size_of_division;
        }
        int rest = num_numbers - num_procs * size_of_division;
        int counter = 0;
        for (int i = rest; i > 0; i--) {
            divisions[counter] += 1;
            counter++;
        }
        double step = (end - begin) / (double)num_numbers;
        double start_process = 0.0;
        for (int i = 0; i < num_procs; i++) {
            double end_process = start_process + divisions[i] * step;
            double data[3];
            data[0] = start_process;
            data[1] = end_process;
            data[2] = divisions[i];
            MPI_Send(&data, 3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
            start_process = end_process;
        }
        free(divisions);
    }
    double data[3];
    MPI_Recv(&data, 3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    double result;
    if (type == 0) {
        result = integrate(&func, data[0], data[1], data[2]);
    }
    if (type == 1) {
        result = integrate(&func2, data[0], data[1], data[2]);
    }
    if (type == 2) {
        result = integrate(&func3, data[0], data[1], data[2]);
    }

    double global_sum;

    MPI_Reduce(&result, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (process_Rank == 0) {
        printf("Result: %f\n", global_sum);
    }
    MPI_Finalize();
    return 0;
}

/*  Author: Cezary Bartoszuk
 *  Email: cbart@students.mimuw.edu.pl
 *
 *  This code is distributed under GPLv2
 * */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

#define TRUE 1
#define FALSE 0

#define ERROR_CODE 1

#define MPI_REAL_T MPI_LONG
typedef signed long int real_t;

#define GRID_SIZE 4
#define GRID_WIDTH GRID_SIZE
#define GRID_HEIGHT GRID_SIZE
#define NDIMS 2

#define MATRIX_SIZE_PER_PROC 32
typedef real_t matrix_t[MATRIX_SIZE_PER_PROC][MATRIX_SIZE_PER_PROC];

inline void initialize_dims(int *dim_size)
{
    dim_size[0] = GRID_HEIGHT;
    dim_size[1] = GRID_WIDTH;
}

inline void initialize_continous(int *periods)
{
    periods[0] = TRUE;
    periods[1] = TRUE;
}

inline void mpi_initialize(int *argc_ptr, char ***argv_ptr)
{
    int size;
    MPI_Init(argc_ptr, argv_ptr);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(size != GRID_WIDTH * GRID_HEIGHT)
    {
        printf("Please run with %d processes.\n", GRID_WIDTH * GRID_HEIGHT);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, ERROR_CODE);
    }
}

inline void mpi_cart_create(int *dim_size, int *periods,  MPI_Comm *comm_ptr)
{
    int reorder = TRUE;
    MPI_Cart_create(MPI_COMM_WORLD, NDIMS, dim_size, periods, reorder, comm_ptr);
}

inline void mpi_get_rank(MPI_Comm comm, int *rank_ptr)
{
    MPI_Comm_rank(comm, rank_ptr);
}

inline void mpi_get_coords(MPI_Comm comm_2d, int rank_2d, int *coord_2d_ptr)
{
    MPI_Cart_coords(comm_2d, rank_2d, NDIMS, coord_2d_ptr);
}

inline void mpi_finalize(void)
{
    MPI_Finalize();
}

inline void mpi_shift_vertical(MPI_Comm comm_2d, int *rank_source_ptr, int *rank_destination_ptr)
{
    int direction = 0;
    int displ = -1;
    MPI_Cart_shift(comm_2d, direction, displ, rank_source_ptr, rank_destination_ptr);
}

inline void mpi_shift_horizontal(MPI_Comm comm_2d, int *rank_source_ptr, int *rank_destination_ptr)
{
    int direction = 1;
    int displ = -1;
    MPI_Cart_shift(comm_2d, direction, displ, rank_source_ptr, rank_destination_ptr);
}

inline void mpi_isend(void *buf, int destination, MPI_Comm comm_2d, MPI_Request *request_ptr)
{
    MPI_Isend(buf, MATRIX_SIZE_PER_PROC * MATRIX_SIZE_PER_PROC, MPI_REAL_T, destination, 42, comm_2d, request_ptr);
}

inline void mpi_ireceive(void *buf, int source, MPI_Comm comm_2d, MPI_Request *request_ptr)
{
    MPI_Irecv(buf, MATRIX_SIZE_PER_PROC * MATRIX_SIZE_PER_PROC, MPI_REAL_T, source, 42, comm_2d, request_ptr);
}

inline void mpi_wait(MPI_Request *request_ptr, MPI_Status *status_ptr)
{
    MPI_Wait(request_ptr, status_ptr);
}

inline void init_matrix(matrix_t matrix, real_t value)
{
    int i, j;
    for(i = 0; i < MATRIX_SIZE_PER_PROC; ++i)
    {
        for(j = 0; j < MATRIX_SIZE_PER_PROC; ++j)
        {
            matrix[i][j] = value;
        }
    }
}

inline void local_matrix_multiply(matrix_t destination, matrix_t source_a, matrix_t source_b_trans)
{
    int i, j, k;
    real_t temp;
    for(i = 0; i < MATRIX_SIZE_PER_PROC; ++i)
    {
        for(j = 0; j < MATRIX_SIZE_PER_PROC; ++j)
        {
            temp = 0;
            for(k = 0; k < MATRIX_SIZE_PER_PROC; ++k)
            {
                temp += source_a[i][k] * source_b_trans[j][k];
            }
            destination[i][j] = temp;
        }
    }
}

inline real_t rand_normalized()
{
    return (real_t) ((rand() / (RAND_MAX + 1.0)) * 100.0);
}

inline void rand_init_matrix(matrix_t matrix)
{
    int i, j;
    for(i = 0; i < MATRIX_SIZE_PER_PROC; ++i)
    {
        for(j = 0; j < MATRIX_SIZE_PER_PROC; ++j)
        {
            matrix[i][j] = rand_normalized();
        }
    }
}

int main(int argc, char **argv)
{
    /* Comm handle in 2 dimensional grid */
    MPI_Comm comm_2d;
    /* Size of grid */
    int dim_size[NDIMS];
    int periods[NDIMS];
    /* Rank and coordinates of this process in 2 dimensional grid */
    int rank_2d;
    int coord_2d[NDIMS];
    /* Rank of sources and destinations in 2 dimensional grid */
    int rank_vertical_source, rank_vertical_destination;
    int rank_horizontal_source, rank_horizontal_destination;

    /* Matrices */
    matrix_t a_part[2];  /* A matrix proc's part (current + temp) */
    matrix_t b_part_trans[2];  /* B matrix proc's part (current + temp) (the part is transponed) */
    matrix_t c_part;  /* C matrix proc's part */
    matrix_t *a_part_current;  /* Pointer to current A matrix part */
    matrix_t *b_part_trans_current;  /* Pointer to current B matrix part */
    matrix_t *a_part_temp;  /* Pointer to temp A matrix part */
    matrix_t *b_part_trans_temp;  /* Pointer to temp B matrix part */

    int step;

    MPI_Request request_horizontal_send, request_horizontal_receive;
    MPI_Request request_vertical_send, request_vertical_receive;
    MPI_Status status;

    initialize_dims(dim_size);
    initialize_continous(periods);
    mpi_initialize(&argc, &argv);
    mpi_cart_create(dim_size, periods, &comm_2d);

    mpi_get_rank(comm_2d, &rank_2d);
    mpi_get_coords(comm_2d, rank_2d, coord_2d);
    mpi_shift_vertical(comm_2d, &rank_vertical_source, &rank_vertical_destination);
    mpi_shift_horizontal(comm_2d, &rank_horizontal_source, &rank_horizontal_destination);

    init_matrix(c_part, 0);

    rand_init_matrix(a_part[0]);
    rand_init_matrix(b_part_trans[0]);

    for(step = 0; step < GRID_SIZE; ++step)
    {
        a_part_current = a_part + (step % 2);
        b_part_trans_current = b_part_trans + (step % 2);
        a_part_temp = a_part + (1 - (step % 2));
        b_part_trans_temp = b_part_trans + (1 - (step % 2));
        if(step != 0)
        {
            mpi_isend(*a_part_current, rank_vertical_destination, comm_2d, &request_vertical_send);
            mpi_ireceive(*a_part_temp, rank_vertical_source, comm_2d, &request_vertical_receive);
            mpi_isend(*b_part_trans_current, rank_horizontal_destination, comm_2d, &request_horizontal_send);
            mpi_ireceive(*b_part_trans_temp, rank_horizontal_source, comm_2d, &request_vertical_receive);
        }
        local_matrix_multiply(c_part, *a_part_current, *b_part_trans_current);
        if(step != 0)
        {
            mpi_wait(&request_vertical_send, &status);
            mpi_wait(&request_vertical_receive, &status);
            mpi_wait(&request_horizontal_send, &status);
            mpi_wait(&request_horizontal_receive, &status);
        }
    }

    mpi_finalize();
    return 0;
}

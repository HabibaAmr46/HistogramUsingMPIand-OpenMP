#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

int Find_bin(    float points,
                 float bar_maxes[],
                 int bar_num,
                 float min_value);

void Print_histogram( float bar_maxes[],
                      int bar_counts[],
                      int bar_num,
                      float min_value);


int main(int argc, char* argv[])
{
    int       bar_num;
    float     min_value,max_value;
    float*    bar_maxes;
    int*      bar_counts; //array contain all bars global
    int*      local_bar_counts; //array of bars for each process
    int       points_num; //no of points
    int       local_points_count; //no of points for each process
    float*    points; //array containing all point global
    float*    local_points;//array of points for each proces
    int       my_rank;
    int       comm_sz;
    MPI_Comm  comm;
    int *remain;

    MPI_Init(&argc, &argv);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    int* bar_count_pointer=&bar_num;
    float* min_value_pointer=&min_value;
    int* points_count_pointer=&points_num;
    int* local_points_count_pointer=&local_points_count;
    if(my_rank == 0)
    {
        printf("Number of bars: ");
        scanf("%d",bar_count_pointer);
        printf("Number of points: ");
        scanf("%d",points_count_pointer);
        // Make sure points_num is a multiple of comm_sz
        *local_points_count_pointer = *points_count_pointer / comm_sz;
        if(*local_points_count_pointer * comm_sz > *points_count_pointer)
        {
            *remain=*local_points_count_pointer * comm_sz-*points_count_pointer;
        }
        else
        {
            *remain=*points_count_pointer-*local_points_count_pointer * comm_sz;
        }
        printf("\n");
    }
    MPI_Bcast(bar_count_pointer,1,MPI_INT,0,comm);
    MPI_Bcast(points_count_pointer,1,MPI_INT,0,comm);
    MPI_Bcast(local_points_count_pointer,1,MPI_INT,0,comm);
    MPI_Bcast(remain,1,MPI_INT,0,comm);

    bar_maxes = malloc(bar_num*sizeof(float));
    bar_counts = malloc(bar_num*sizeof(int));
    local_bar_counts = malloc(bar_num*sizeof(int));
    points = malloc(points_num*sizeof(float));
    local_points = malloc(local_points_count*sizeof(float));

    float* max_value_pointer=&max_value;
    if(my_rank == 0)
    {
        points = malloc(points_num*sizeof(float));
        FILE *filePtr;
        filePtr= fopen("/shared/dataset.txt","r");
        int i=0;
        float min=0.0;
        int x;
        float max=-1;
        while (fscanf(filePtr, "%d", &x) != EOF)
        {
            points[i]=x;
            if(points[i]>max)
                max=points[i];
            i++;

        }
        max_value_pointer=&max;
        max_value=*max_value_pointer;
        min_value_pointer=&min;
        min_value=*min_value_pointer;

    }
    MPI_Bcast(max_value_pointer,1,MPI_FLOAT,0,comm);
    MPI_Bcast(min_value_pointer,1,MPI_FLOAT,0,comm);
    MPI_Scatter(points,local_points_count,MPI_FLOAT,local_points,local_points_count,MPI_FLOAT, 0, comm);

    //////// Intialize bins with bin maxes and local bins
    float range = max_value - min_value;
    float interval = range / bar_num;

    int i;

    #pragma omp parallel for shared(interval,bar_maxes) private(i)
    for(i = 0; i < bar_num; i++)
    {
        bar_maxes[i] = interval * (float)(i+1) + min_value;
        local_bar_counts[i] = 0;
    }
    ///////
    /////// Calculating bins
    int bin;
    #pragma omp parallel for private(i)
    for(i = 0; i < local_points_count; i++)
    {
        bin = Find_bin(local_points[i],bar_maxes,bar_num,min_value);
        #pragma omp critical
        local_bar_counts[bin]++;
    }
    ////////
    MPI_Barrier (MPI_COMM_WORLD);
    if(my_rank==0)
    {
        int j=*remain;
        if(j!=0)
        {
            #pragma omp parallel for private(i)
            for(i = j; i>0; i--)
            {

                bin = Find_bin(points[points_num-i],bar_maxes,bar_num,min_value);
                #pragma omp critical
                local_bar_counts[bin]++;

            }
        }
    }

    MPI_Reduce(local_bar_counts,bar_counts,bar_num,MPI_INT,MPI_SUM,0,comm);
    if(my_rank == 0)
    {
        Print_histogram(bar_maxes,bar_counts,bar_num,min_value);
    }


    free(bar_maxes);
    free(bar_counts);
    free(local_bar_counts);
    free(points);
    free(local_points);
    MPI_Finalize();
    return 0;
}



void Print_histogram(
    float bar_maxes[] /* in */,
    int bar_counts[]  /* in */,
    int bar_num     /* in */,
    float min_value    /* in */)
{


    int i;
    for(i = 0; i < bar_num; i++)
    {
        if(i==0)
        {
            printf("The range starts with %d, end with %d with count %d \n",(int)min_value,(int)bar_maxes[i],bar_counts[i]);
        }
        else
        {
            printf("The range starts with %d, end with %d with count %d \n",(int)bar_maxes[i-1],(int)bar_maxes[i],bar_counts[i]);

        }
    }
}


int Find_bin(float points, float bar_maxes[], int bar_num,
             float min_value)
{

    int i;
    for(i = 0; i < bar_num; i++)
    {
        if(points <= bar_maxes[i]) break;
    }
    return i;
}

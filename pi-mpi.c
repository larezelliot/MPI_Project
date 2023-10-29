#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define NUMSTEPS 1000000

double calculate_approximation(double a, double b, double step){
	double approx = 0;
	double legth_per_step = (double)(b - a)/ step;
	double i;
	
	for(i = a; i < b; i += legth_per_step)
		approx += 4.0 / (1 + i * i);
        
	approx *= legth_per_step;
	
	return approx;
}

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    // Get the numbers of steps per process
    int stepsPerProcess = world_size/NUMSTEPS;
    
	// Determine the local margins and proximation of each process
	double local_a = world_rank * stepsPerProcess;
    double local_b = (world_rank+1) * stepsPerProcess;
    double local_approximation = calculate_approximation(local_a, local_b, stepsPerProcess);
    
    // Determine the left/right children and parent
    int left_child = 2*world_rank+1;
    int right_child = 2*world_rank+2;
    int parent = (world_rank+1)/2 - 1;
    
    // Initializes children approximations
    double left_approximation = 0;
    double right_approximation = 0;
    
    // Receives children approximations
    if(left_child < world_size)
    	MPI_Recv(&left_approximation, 1, MPI_DOUBLE, left_child, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if(right_child < world_size)
    	MPI_Recv(&right_approximation, 1, MPI_DOUBLE, right_child, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    // Calculates Approximation
    double approximation = local_approximation + left_approximation + right_approximation;
    
    // Sends approximation to parent
    if(world_rank != 0){
        double PI = approximation;
		printf("PI is %.20f\n",pi);
	}else{
    	MPI_Send(&approximation, 1, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD);
	}
	
    return 0;
}



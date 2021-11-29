// Need to be superuser to compile: sudo su 
// Compile with: g++ -lpthread Assignment1.cpp -o Assignment1
// If doesn't work, try: g++ -pthread Assignment1.cpp -o Assignment1

// This exercise show how to schedule four periodic threads with Rate Monotonic with semaphores 
// using priority ceiling

/*
The objective of the first assignment is to create an application with four periodic
tasks, with period of 80ms, 100ms, 160ms and 200ms (task with higher priority is called)
task1 and the one with lowest is called task4)
In the program there are three global varibales which are:
- T1T2: task1 should write on it and task2 should read from it.
- T1T4: task1 should write on it and task4 should read from it.
- T2T3: task2 should write on it and task3 should read from it.
All critical sections shall be protected by semaphores which shall use Priority Ceiling

------------------------------------------ From Theory ------------------------------------------------
																									
Because of Priority Ceiling protocol, every semaphore S(i)(j) is assigned a 						
priority ceiling which is the highest priority among all the tasks that can use the semaphore
S12 -> Priority of task1
S14 -> Priority of task1
S23 -> Priority of task2

Beta(i) -> is the set of critical section of task(i)
Beta1 = { z11 , z12 }
Beta2 = { z21 , z22 }
Beta3 = { z31 }
Beta4 = { z41 }

Critical sections for each global variable
T1T2 -> z11 , z21
T1T4 -> z12 , z41
T2T3 -> z22 , z31

Beta*(i) -> is the set of all the critical sections that can block the task(i) 
Beta*1 = { z21 , z41 }
Beta*2 = { z31 , z41 }
Beta*3 = { z41 }
Beta*4 = { 0 } has lowest priority so it can NOT be blocked

Beta(i)(j) is the set of critical sections of task(j) that can block the higher priority task(i)
Beta12 = { z21 }
Beta14 = { z41 }
Beta23 = { z31 }
Beta24 = { z41 }
Beta34 = { z41 }
-------------------------------------------------------------------------------------------------------
*/

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/types.h>


// Code of periodic tasks
void task1_code( );
void task2_code( );
void task3_code( );
void task4_code( );

// There are no aperiodic tasks that need to be scheduled

// Initialization of global variables
double T1T2; // task1 writes and task2 reads
double T1T4; // task1 writes and task4 reads
double T2T3; // task2 writes and task3 reads

// Initialization of blocking time
// Global variables in order to save the maximum duration of the critical zone that can
// block the i-th task
double B1;
double B2;
double B3;
double B4;

// Initialization of critical zones
double z11;
double z12;
double z21;
double z22;
double z31;
double z41;

// Characteristic function of the thread, only for timing and synchronization
// periodic tasks
void *task1( void *);
void *task2( void *);
void *task3( void *);
void *task4( void *);

// Initialization and declaration of semaphores 
pthread_mutex_t S12 = PTHREAD_MUTEX_INITIALIZER; 
pthread_mutex_t S14 = PTHREAD_MUTEX_INITIALIZER; 
pthread_mutex_t S23 = PTHREAD_MUTEX_INITIALIZER; 

// Initialization of mutex attribute
pthread_mutexattr_t mymutexattr;

#define INNERLOOP 1000
#define OUTERLOOP 100
#define WATSELOOP 200

#define NPERIODICTASKS 4
#define NAPERIODICTASKS 0
#define NTASKS NPERIODICTASKS + NAPERIODICTASKS

long int periods[NTASKS];
struct timespec next_arrival_time[NTASKS];
double compWCET[NTASKS];
double WCET[NTASKS];
pthread_attr_t attributes[NTASKS];
pthread_t thread_id[NTASKS];
struct sched_param parameters[NTASKS];
int missed_deadlines[NTASKS];

int main()
{
  	// Set task periods in nanoseconds
	// The first task has period 80 millisecond
	// The second task has period 100 millisecond
	// The third task has period 160 millisecond
	// The fourth task has period 200 millisecond
	// They are already ordered according to their priority following
    // the Rate Monotonic algorithm in which the ma priority il given
    // to the task with lower peridod and so on; 
	// If tasks weren't already sorted according to their priority,
    // there would have been the need of an algorithm to sort them.
  	periods[0]= 80000000; //in nanoseconds
  	periods[1]= 100000000; //in nanoseconds
  	periods[2]= 160000000; //in nanoseconds
	periods[3]= 200000000; //in nanoseconds

	// This is not strictly necessary, but it is convenient to
	// assign a name to the maximum and the minimum priotity in the
	// system. They are called priomin and priomax.
  	struct sched_param priomax;
  	priomax.sched_priority=sched_get_priority_max(SCHED_FIFO);
  	struct sched_param priomin;
  	priomin.sched_priority=sched_get_priority_min(SCHED_FIFO);

	// Set the maximum priority to the current thread (you are required to be
  	// superuser). Check that the main thread is executed with superuser privileges
	// before doing anything else.
  	if (getuid() == 0)
	{
    	pthread_setschedparam(pthread_self(), SCHED_FIFO, &priomax);
	}

  	// Execute all tasks in standalone modality in order to measure execution times
  	// (using gettimeofday). Use the computed values to update the worst case execution
  	// time of each task
  	for (int i =0; i < NTASKS; i++)
    {
		// Setting the Worst Case Execution Time as zero for each task
		WCET[i] = 0;

		// Loop to allow a better computation of the WCET 
		// Loop 50 times and get the highest among the WCET founded
		for (int j=0; j < 50; j++) 
		{
			// Initialize time_1 and time_2 required to read the clock
			struct timespec time_1, time_2;
			clock_gettime(CLOCK_REALTIME, &time_1);

			// I execute each task 50 times and save the highest WCET
			// found while doing the loops
 	    	if (i==0)
				task1_code();
      		if (i==1)
				task2_code();
      		if (i==2)
				task3_code();
      		if (i==3)
				task4_code();

			clock_gettime(CLOCK_REALTIME, &time_2);


			// Compute the Worst Case Execution Time
			// In order to have a more attendible result i take the maximum WCET 
			// among all WCET calculated.
			compWCET[i]= 1000000000*(time_2.tv_sec - time_1.tv_sec)+(time_2.tv_nsec-time_1.tv_nsec);
      		printf("WCET[%d]{%d} ; ", i, j);
			if ( compWCET[i] > WCET[i] ){
				WCET[i] = compWCET[i];
			}
		}
		printf("\n\n\nWorst Case Execution Time %d=%f \n\n\n", i+1, WCET[i]);
    }

	// Get values of the critical sections

	// Compute B(j)(i) which is the blocking time of the higher priority task1 that can be blocked
	// by the lower priority task2. Becasue of priority ceiling protocol, a task can be blocked for 
	// at maximum one critical section which is the one that takes the logest time.
	B1 = std::fmax( z11 , z12);
	B2 = std::fmax( z21 , z22);
	B3 = z31;
	B4 = 0; // task4 can't be blocked
	
	// Compute Ui 
	// In order to compute the utilization time (U) there is the need to apply the rule used
	// for Rate Monotonic for each periodic task, using just the tasks which have lower or 
	// equal priority of the considered one. 
	// After that, it's needed to add the maximum blocking time (B) for which a task can be 
	// blocked to the utilization time previously found for each task.
	double U[NTASKS];
	
	U[0] = WCET[0]/periods[0] + B1/periods[0];

	U[1] = WCET[0]/periods[0] + WCET[1]/periods[1] + B2/periods[1];

	U[2] = WCET[0]/periods[0] + WCET[1]/periods[1] + WCET[2]/periods[2] + B3/periods[2];

	U[3] = WCET[0]/periods[0] + WCET[1]/periods[1] + WCET[2]/periods[2] + WCET[3]/periods[3];
    	
	// Compute Ulub for each task and check the sufficient conditions
	double Ulub[NTASKS];
	for (int i = 0; i < NTASKS; i++)
	{
		// If there are no harmonic relationships between tasks periods, there is the need to use 
		// the following formula in order to compute Ulub
		double Ulub[i] = (i+1)*(pow(2.0,(1.0/(i+1))) -1);
	}


	// Check the sufficient conditions for each U and Ulub calculated previously: 
	// if conditions are not satisfied, exit the program
	// (could have done with a for loop but i kept having problem at printing)
  	if (U[0] > Ulub[0])
    {
      	printf("\n\n\n U[1]=%lf Ulub[1]=%lf Non schedulable Task Set\n", U[0], Ulub[0]);
		fflush(stdout);	
      	return(-1);
    } else {
  		printf("\n\n\n U[1]=%lf Ulub[1]=%lf Scheduable Task Set\n", U[0], Ulub[0]);
  		fflush(stdout);	
	}
	if (U[1] > Ulub[1])
    {
      	printf(" U[2]=%lf Ulub[2]=%lf Non schedulable Task Set\n", U[1], Ulub[1]);
		fflush(stdout);	
      	return(-1);
    } else {
  		printf(" U[2]=%lf Ulub[2]=%lf Scheduable Task Set\n", U[1], Ulub[1]);
  		fflush(stdout);	
	}
	if (U[2] > Ulub[2])
    {
      	printf(" U[3]=%lf Ulub[3]=%lf Non schedulable Task Set\n", U[2], Ulub[2]);
		fflush(stdout);	
      	return(-1);
    } else {
  		printf(" U[3]=%lf Ulub[3]=%lf Scheduable Task Set\n", U[2], Ulub[2]);
  		fflush(stdout);	
	}
	if (U[3] > Ulub[3])
    {
      	printf(" U[4]=%lf Ulub[4]=%lf Non schedulable Task Set\n\n\n", U[3], Ulub[3]);
		fflush(stdout);	
      	return(-1);
    } else {
  		printf("\n\n\n U[4]=%lf Ulub[4]=%lf Scheduable Task Set\n\n\n", U[3], Ulub[3]);
  		fflush(stdout);	
	}
	sleep(2);

  	// Set the minimum priority to the current thread: this is now required because 
	// we will assign higher priorities to periodic threads to be soon created
	// pthread_setschedparam
  	if (getuid() == 0)
	{
    	pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomin);
	}
  
  	// Set the attributes of each task, including scheduling policy and priority
  	for (int i=0; i < NPERIODICTASKS; i++)
    {
		// Initialize the attribute structure of task i
      	pthread_attr_init(&(attributes[i]));

		// Set the attributes to tell the kernel that the priorities and policies are explicitly chosen,
		// not inherited from the main thread (pthread_attr_setinheritsched) 
      	pthread_attr_setinheritsched(&(attributes[i]), PTHREAD_EXPLICIT_SCHED);
      
		// Set the attributes to set the SCHED_FIFO policy (pthread_attr_setschedpolicy)
		pthread_attr_setschedpolicy(&(attributes[i]), SCHED_FIFO);

		// Properly set the parameters to assign the priority inversely proportional 
		// to the period
      	parameters[i].sched_priority = priomin.sched_priority+ NTASKS - i;

		//set the attributes and the parameters of the current thread (pthread_attr_setschedparam)
      	pthread_attr_setschedparam(&(attributes[i]), &(parameters[i]));
    }

	pthread_mutexattr_init(&mymutexattr);

    // These code lines below are used to test priority inheritance protocol
    // which causes a deadlock situation between task1 and task2
    // In order to see the difference between the priority inheritance and the
    // prioity ceiling protocols is sufficient to decomment the following lines
    // and comment the lines for the priority ceiling
    /*
    pthread_mutexattr_setprotocol(&mymutexattr, PTHREAD_PRIO_INHERIT);
    pthread_mutex_init(&S12, &mymutexattr);
    pthread_mutex_init(&S14, &mymutexattr);
    pthread_mutex_init(&S23, &mymutexattr);
    */

    // These code lines below are used to set the priority ceiling protocol
    // for the semaphores (which is the request of the assignment)
    // In order to see the difference between priority ceiling and priority
    // inheritance  protocols is sufficient to decomment the following lines 
    // and comment the lines for the priority inheritance
    // Semaphores S12, S14 and S23 take the priority of the highest priority among 
    // all tasks that can use the semaphore
    pthread_mutexattr_setprotocol(&mymutexattr, PTHREAD_PRIO_PROTECT);
    pthread_mutexattr_setprioceiling(&mymutexattr, parameters[0].sched_priority);
    pthread_mutex_init(&S12, &mymutexattr);
    pthread_mutex_init(&S14, &mymutexattr);
    pthread_mutexattr_setprioceiling(&mymutexattr, parameters[1].sched_priority);
    pthread_mutex_init(&S23, &mymutexattr);

	// Declare the variable to contain the return values of pthread_create	
  	int iret[NTASKS];

	// Declare variables to read the current time
	struct timespec time_1;
	clock_gettime(CLOCK_REALTIME, &time_1);

  	// Set the next arrival time for each task. This is not the beginning of the first
	// period, but the end of the first period and beginning of the next one. 
  	for (int i=0; i < NPERIODICTASKS; i++)
    {
		long int next_arrival_nanoseconds = time_1.tv_nsec + periods[i];
		// Compute the end of the first period and beginning of the next one
		next_arrival_time[i].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[i].tv_sec= time_1.tv_sec + next_arrival_nanoseconds/1000000000;
       	missed_deadlines[i] = 0;
    }

	

	// Create all threads (pthread_create)
  	iret[0] = pthread_create( &(thread_id[0]), &(attributes[0]), task1, NULL);
  	iret[1] = pthread_create( &(thread_id[1]), &(attributes[1]), task2, NULL);
  	iret[2] = pthread_create( &(thread_id[2]), &(attributes[2]), task3, NULL);
   	iret[3] = pthread_create( &(thread_id[3]), &(attributes[3]), task4, NULL);
  

  	// Join all threads (pthread_join)
  	pthread_join( thread_id[0], NULL);
  	pthread_join( thread_id[1], NULL);
  	pthread_join( thread_id[2], NULL);
	pthread_join( thread_id[3], NULL);


  	// Set the next arrival time for each task. This is not the beginning of the first
	// period, but the end of the first period and beginning of the next one. 
  	for (int i=0; i < NTASKS; i++)
    {
      	printf ("\nMissed Deadlines Task %d=%d", i+1, missed_deadlines[i]);
		fflush(stdout);
    }
  	exit(0);
}

// Follow the instructions with the arrow sign: <----- in order to try to create a deadlock
// Be aware: in order to create a deadlock it's needed to change the sched policy to priority
// inheritance from priority ceiling. Also the behavior of the code will be different (it's just
// a check for deadlock, it's not about running correctly the program)
// Application specific task1 code
void task1_code()
{
	// Initialize time_1, time_2, time_3 and time_4 required to read the clock
	struct timespec time_1, time_2, time_3, time_4;
	// Print the id of the current task
  	printf(" 1[ "); fflush(stdout);

	// This double loop with random computation is only required to waste time
	double uno;
  	for (int i = 0; i < OUTERLOOP; i++)
    {
      	for (int j = 0; j < INNERLOOP; j++)
		{
			uno = rand()*rand()%10;
    	}
  	}

	// Critical zone z11
	// Taking the semaphore
	pthread_mutex_lock(&S12);
	// Getting the time at the beginning of the critical zone
	clock_gettime(CLOCK_REALTIME, &time_1);
	// Waste time inside the critical zone
	double due;
	for (int i = 0; i < WATSELOOP; i++)
	{
		due = rand()*rand()%10;
	}
	// Task1 writes on T1T2
	T1T2 = due;
	printf(" Writing in T1T2 ");
	fflush(stdout);
	// pthread_mutex_lock(&S14); // decomment this line in order to try to create a deadlock <-----
	// Getting the time at the end of the critical zone
	clock_gettime(CLOCK_REALTIME, &time_2);
	// Releasing the semaphore
	pthread_mutex_unlock(&S12); // comment in order to try to create a deadlock <-----
	// Saving time in the global variable initialized for z11 critical zone
	z11 = 1000000000*(time_2.tv_sec - time_1.tv_sec)+(time_2.tv_nsec - time_1.tv_nsec);

	// Critical zone z12
	// Taking the semaphore
	pthread_mutex_lock(&S14); // comment in order to try to try to create a deadlock <-----
	// Getting the time at the beginning of the critical zone
	clock_gettime(CLOCK_REALTIME, &time_3);
	// Waste time inside the critical zone
	double tre;
	for (int i = 0; i < WATSELOOP; i++)
	{
		tre = rand()*rand()%10;
	}
	// Task1 writes on T1T4
	T1T4 = tre;
	printf(" Writing in T1T4 ");
	fflush(stdout);
	// Getting the time at the end of the critical zone
	clock_gettime(CLOCK_REALTIME, &time_4);
	// Releasing the semaphore
	pthread_mutex_unlock(&S14);
	// pthread_mutex_unlock(&S12); // decomment this line in order to try to create a deadlock <-----
	// Saving time in the global variable initialized for z12 critical zone
	z12 = 1000000000*(time_4.tv_sec - time_3.tv_sec)+(time_4.tv_nsec - time_3.tv_nsec);

  	// Print the id of the current task
  	printf(" ]1 "); fflush(stdout);
}

// Thread code for task_1 (used only for temporization)
void *task1( void *ptr)
{
	// Set thread affinity, that is the process on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	// Initialize time for computing missed deadlines
	struct timespec waittime;
	waittime.tv_sec = 0;
	waittime.tv_nsec = periods[0];

   	// Execute the task one hundred times... it should be an infinite loop (too dangerous)
  	for (int i=0; i < 100; i++)
    {
      	// Execute application specific code
		task1_code();

		// Getting time as soon as task1 ends its execution
		struct timeval ora;
		struct timezone zona;
		gettimeofday(&ora ,&zona); // gets the time of the day

		// After execution, it computes the remaining time before the next period starts
		long int timetowait = 1000*((next_arrival_time[0].tv_sec - ora.tv_sec)*1000000 + (next_arrival_time[0].tv_nsec - ora.tv_usec));
		waittime.tv_sec = timetowait/1000000000;
		waittime.tv_nsec = timetowait%1000000000;

		// If i'm already in another period, increment the missed deadline counter
		if (timetowait < 0)
		{
			missed_deadlines[0] ++;
		}

		// Sleep until the end of the current period (which is also the start of the
		// new one)
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[0], NULL);

		// The thread is ready and can compute the end of the current period for
		// the next iteration
		long int next_arrival_nanoseconds = next_arrival_time[0].tv_nsec + periods[0];
		next_arrival_time[0].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[0].tv_sec= next_arrival_time[0].tv_sec + next_arrival_nanoseconds/1000000000;
    }
}

void task2_code()
{
	// Initialize time_1, time_2, time_3 and time_4 required to read the clock
	struct timespec time_1, time_2, time_3, time_4;
	// Print the id of the current task
  	printf(" 2[ "); fflush(stdout);
	double uno;

	// This double loop with random computation is only required to waste time
  	for (int i = 0; i < OUTERLOOP; i++)
    {
      	for (int j = 0; j < INNERLOOP; j++)
		{
			uno = rand()*rand()%10;
		}
    }

	// pthread_mutex_lock(&S23) // decomment this line in order to try to create a deadlock <-----
	// Critical zone z21
	// Taking the semaphore
	pthread_mutex_lock(&S12);
	// Getting the time at the beginning of the critical zone
	clock_gettime(CLOCK_REALTIME, &time_1);
	// Waste time inside the critical zone
	double waste;
	for (int i = 0; i < WATSELOOP; i++)
	{
		waste = rand()*rand()%10;
	}
	// Task2 reads from T1T2
	double due = T1T2;
	printf(" Read from T1T2: %f ", due);
	fflush(stdout);
	// Getting the time at the end of the critical zone
	clock_gettime(CLOCK_REALTIME, &time_2);
	// Releasing the semaphore
	pthread_mutex_unlock(&S12);
	// Saving time in the global variable initialized for z11 critical zone
	z21 = 1000000000*(time_2.tv_sec - time_1.tv_sec)+(time_2.tv_nsec - time_1.tv_nsec);

	// Critical zone z22
	// Taking the semaphore
	pthread_mutex_lock(&S23); // comment this line in order to try to create a deadlock <-----
	// Getting the time at the beginning of the critical zone
	clock_gettime(CLOCK_REALTIME, &time_3);
	// Waste time inside the critical zone
	double tre;
	for (int i = 0; i < WATSELOOP; i++)
	{
		tre = rand()*rand()%10;
	}
	// Task2 writes on T2T3
	T2T3 = tre;
	printf(" Writing in T2T3 ");
	fflush(stdout);
	// Getting the time at the end of the critical zone
	clock_gettime(CLOCK_REALTIME, &time_4);
	// Releasing the semaphore
	pthread_mutex_unlock(&S23);
	// Saving time in the global variable initialized for z12 critical zone
	z22 = 1000000000*(time_4.tv_sec - time_3.tv_sec)+(time_4.tv_nsec - time_3.tv_nsec);

	// Print the id of the current task
  	printf(" ]2 "); fflush(stdout);
}


void *task2( void *ptr )
{
	// Set thread affinity, that is the process on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	// Initialize time for computing missed deadlines
	struct timespec waittime;
	waittime.tv_sec = 0;
	waittime.tv_nsec = periods[1];

	// Execute the task eighty times... it should be an infinite loop (too dangerous)
  	for (int i=0; i < 80; i++)
    {
		// Execute application specific code
      	task2_code();

		// Getting time as soon as task1 ends its execution
		struct timeval ora;
		struct timezone zona;
		gettimeofday(&ora ,&zona); // gets the time of the day

		// After execution, it computes the remaining time before the next period starts
		long int timetowait = 1000*((next_arrival_time[1].tv_sec - ora.tv_sec)*1000000 + (next_arrival_time[1].tv_nsec - ora.tv_usec));
		waittime.tv_sec = timetowait/1000000000;
		waittime.tv_nsec = timetowait%1000000000;

		// If i'm already in another period, increment the missed deadline counter
		if (timetowait < 0)
		{
			missed_deadlines[1] ++;
		}

		// Sleep until the end of the current period (which is also the start of the
		// new one)
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[1], NULL);

		// The thread is ready and can compute the end of the current period for
		// the next iteration
		long int next_arrival_nanoseconds = next_arrival_time[1].tv_nsec + periods[1];
		next_arrival_time[1].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[1].tv_sec= next_arrival_time[1].tv_sec + next_arrival_nanoseconds/1000000000;
    }
}

void task3_code()
{
	// Initializa time_1 and time_2 required to read the clock
	struct timespec time_1, time_2;
	// Print the id of the current task
  	printf(" 3[ "); fflush(stdout);

	// This double loop with random computation is only required to waste time
	double uno;
  	for (int i = 0; i < OUTERLOOP; i++)
    {
      	for (int j = 0; j < INNERLOOP; j++)
		{		
			double uno = rand()*rand()%10;
		}
    }

	// Critical zone z31
	// Taking the semaphore
	pthread_mutex_lock(&S23);
	// Getting the time at the beginning of the critical zone
	clock_gettime(CLOCK_REALTIME, &time_1);
	// Waste time inside the critical zone
	double waste;
	for (int i = 0; i < WATSELOOP; i++)
	{
		waste = rand()*rand()%10;
	}
	// Task3 reads from T2T3
	double due = T2T3;
	printf(" Read from T2T3: %f ", due);
	fflush(stdout);
	// Getting the time at the end of the critical zone
	clock_gettime(CLOCK_REALTIME, &time_2);
	// Releasing the semaphore
	pthread_mutex_unlock(&S23);
	// Saving time in the global variable initialized for z11 critical zone
	z31 = 1000000000*(time_2.tv_sec - time_1.tv_sec)+(time_2.tv_nsec - time_1.tv_nsec);

	// Print the id of the current task
  	printf(" ]3 "); fflush(stdout);
}

void *task3( void *ptr)
{
	// Set thread affinity, that is the process on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	// Initialize time for computing missed deadlines
	struct timespec waittime;
	waittime.tv_sec = 0;
	waittime.tv_nsec = periods[2];

	// Execute the task fifty times... it should be an infinite loop (too dangerous)
  	for (int i=0; i < 50; i++)
    {
		// Execute application specific code
      	task3_code();

		// Getting time as soon as task1 ends its execution
		struct timeval ora;
		struct timezone zona;
		gettimeofday(&ora ,&zona); // gets the time of the day

		// After execution, it computes the remaining time before the next period starts
		long int timetowait = 1000*((next_arrival_time[2].tv_sec - ora.tv_sec)*1000000 + (next_arrival_time[2].tv_nsec - ora.tv_usec));
		waittime.tv_sec = timetowait/1000000000;
		waittime.tv_nsec = timetowait%1000000000;

		// If i'm already in another period, increment the missed deadline counter
		if (timetowait < 0)
		{
			missed_deadlines[2] ++;
		}

		// Sleep until the end of the current period (which is also the start of the
		// new one)
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[2], NULL);

		// The thread is ready and can compute the end of the current period for
		// the next iteration
		long int next_arrival_nanoseconds = next_arrival_time[2].tv_nsec + periods[2];
		next_arrival_time[2].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[2].tv_sec= next_arrival_time[2].tv_sec + next_arrival_nanoseconds/1000000000;
    }
}

void task4_code()
{
	// Initializa time_1 and time_2 required to read the clock
	struct timespec time_1, time_2;
	// Print the id of the current task
  	printf(" 4[ "); fflush(stdout);
	double uno;

	// This double loop with random computation is only required to waste time
  	for (int i = 0; i < OUTERLOOP; i++)
    {
      	for (int j = 0; j < INNERLOOP; j++)
		{		
			double uno = rand()*rand()%10;
		}
    }

	// Critical zone z41
	// Taking the semaphore
	pthread_mutex_lock(&S14);
	// Getting the time at the beginning of the critical zone
	clock_gettime(CLOCK_REALTIME, &time_1);
	// Waste time inside the critical zone
	double waste;
	for (int i = 0; i < WATSELOOP; i++)
	{
		waste = rand()*rand()%10;
	}
	// Task4 reads from T1T4
	double due = T1T4;
	printf(" Read from T1T4: %f ", due);
	fflush(stdout);
	// Getting the time at the end of the critical zone
	clock_gettime(CLOCK_REALTIME, &time_2);
	// Releasing the semaphore
	pthread_mutex_unlock(&S14);
	// Saving time in the global variable initialized for z11 critical zone
	z41 = 1000000000*(time_2.tv_sec - time_1.tv_sec)+(time_2.tv_nsec - time_1.tv_nsec);

	// Print the id of the current task
  	printf(" ]4 "); fflush(stdout);
}

void *task4( void *ptr)
{
	// Set thread affinity, that is the process on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	// Initialize time for computing missed deadlines
	struct timespec waittime;
	waittime.tv_sec = 0;
	waittime.tv_nsec = periods[3];

	// Execute the task fourty times... it should be an infinite loop (too dangerous)
  	for (int i=0; i < 40; i++)
    {
		// Execute application specific code
      	task4_code();

		// Getting time as soon as task1 ends its execution
		struct timeval ora;
		struct timezone zona;
		gettimeofday(&ora ,&zona); // gets the time of the day

		// After execution, it computes the remaining time before the next period starts
		long int timetowait = 1000*((next_arrival_time[3].tv_sec - ora.tv_sec)*1000000 + (next_arrival_time[3].tv_nsec - ora.tv_usec));
		waittime.tv_sec = timetowait/1000000000;
		waittime.tv_nsec = timetowait%1000000000;

		// If i'm already in another period, increment the missed deadline counter
		if (timetowait < 0)
		{
			missed_deadlines[3] ++;
		}

		// Sleep until the end of the current period (which is also the start of the
		// new one)
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[3], NULL);

		// The thread is ready and can compute the end of the current period for
		// the next iteration
		long int next_arrival_nanoseconds = next_arrival_time[3].tv_nsec + periods[3];
		next_arrival_time[3].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[3].tv_sec= next_arrival_time[3].tv_sec + next_arrival_nanoseconds/1000000000;
    }
}

#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include <cstdlib>
#include <new>
#include <cmath>
#include <time.h>
#include "timer.h"

void print_array(const int * ptr, unsigned size, char ch = ' '){
    std::cout << static_cast<int>(ch) << "{ ";
    for(unsigned i{}; i<size; i++)
        std::cout << *(ptr++) << ' ';
    std::cout << " }\n";
}

//$ mergesort bottom up
inline int min(const int x, const int y){
    return (x < y ? x : y);
}

bool verify(int * arr, std::size_t length){
    for (std::size_t i{1}; i<length; i++){
        if(arr[i-1] > arr[i])
            return false;
    }
    return true;
}

inline void merge(int * array, int * temp, std::size_t length, int from, int mid, int to)
{
	int k {from};
    int i {from};
    int j {mid + 1};

	// loop till no elements are left in the left and right runs
	while (i <= mid && j <= to)
	{
		if (array[i] < array[j]) {
			temp[k++] = array[i++];
		}
		else {
			temp[k++] = array[j++];
		}
	}

	// copy remaining elements
	while (i < length && i <= mid) {
		temp[k++] = array[i++];
	}

	/* no need to copy the second half (since the remaining items
	   are already in their correct position in the temporary array) */

	// copy back to the original array to reflect sorted order
	for (int p = from; p <= to; p++) {
		array[p] = temp[p];
	}
}

// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
inline void merge_c(int * arr, int l, int m, int r)
{
    int n1 = m - l + 1;
    int n2 = r - m;
 
    // Create temp arrays
    int * L {new int [n1] {}};
    int * R {new int [n2] {}};
 
    // Copy data to temp arrays L[] and R[]
    for (int i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (int j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];
 
    // Merge the temp arrays back into arr[l..r]
 
    // Initial index of first subarray
    int i = 0;
 
    // Initial index of second subarray
    int j = 0;
 
    // Initial index of merged subarray
    int k = l;
 
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
 
    // Copy the remaining elements of
    // L[], if there are any
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }
 
    // Copy the remaining elements of
    // R[], if there are any
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }

    delete [] L;
    delete [] R;
}

void core_merge_sort(int * array, const std::size_t length){
    
    int * temp { new int [length] {} };
    //memcpy(temp, array, length);
    for (int i{0}; i < length; i++)
        temp[i] = array[i];


    int low {0};
    int high {static_cast<int>(length) - 1};

    int from, mid, to;

    for(std::size_t m {1}; m <= high - low; m = 2*m){

        for(std::size_t i{static_cast<std::size_t>(low)}; i < high; i += 2*m){
            from = i;
            mid = i + m - 1;
            to = min(i + 2*m - 1, high);

            merge(array, temp, length, from, mid, to);
        }
    }

    delete [] temp;
}


void mpi_merge_sort(int argc, char ** argv, int * array, const std::size_t length){

    MPI_Init(&argc, &argv);

    int world_size{};
    int world_rank{};

    MPI_Comm_size( MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank( MPI_COMM_WORLD, &world_rank);

    //if( !world_rank ){
    //    printf("Original: ");
    //    print_array(array, length);
    //}


    int num_workers{world_size - 1};

    // local arrays
    std::size_t capacity{ length };
    int * local_array { new int [capacity] {} };

    // World group construction in order to construct workers groups (root + workers)
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);

    TimerMPI timer_total;
    TimerMPI timer_comunication;
    //TimerMPI timer_taskdist;
    TimerMPI timer_comp;

    timer_total.reset();
    //timer_taskdist.reset();

    int counter = 0;
    
    
    
    int * local_array_lengths {new int [world_size] {}};
    int * offsets {new int [world_size] {}};
    int * cut_points {new int [world_size] {-1}}; // begging point + cutting points + ending point
    int * previous_cut_points {new int [world_size] {}};
    int * previous_lengths {new int [world_size]{}};
    int * workers_rank {new int [world_size]};
    int previous_wsize;


    while(num_workers >= 1){
        std::fill(workers_rank, workers_rank + world_size - 1, 0);
        
        for(int i{0}; i < num_workers + 1; ++i)
            workers_rank[i] = i;
        
        MPI_Group workers_group;
        MPI_Group_incl(world_group, num_workers + 1, workers_rank, &workers_group);

        MPI_Comm workers_comm;
        MPI_Comm_create_group(MPI_COMM_WORLD, workers_group, 0, &workers_comm);
        
        //std::cout << world_rank <<" okk " << num_workers << '\n';

        int wrank{-1}, wsize{-1};
        if (MPI_COMM_NULL == workers_comm) // if this rank isn`t in workers_comm, it will be MPI_COMM_NULL
            break;
        

        MPI_Comm_rank(workers_comm, &wrank);
        MPI_Comm_size(workers_comm, &wsize);



        // constructing local array lengths and offsets for scatterv
        //int local_array_lengths[num_workers + 1] {};
        //int offsets[num_workers + 1] {};
        //int cut_points[num_workers + 1] {-1}; // begging point + cutting points + ending point

        std::fill(local_array_lengths, local_array_lengths + world_size, 0);
        std::fill(offsets, offsets + world_size, 0);
        std::fill(cut_points, cut_points + world_size, 0);
        cut_points[0] = -1;






        std::size_t length_per_nodes { static_cast<std::size_t>( static_cast<double>(length)/num_workers ) };
        std::size_t mod { static_cast<std::size_t>(length%num_workers) };
        

        for (unsigned m{}; m < mod; m++)
            local_array_lengths[m + 1] = 1;

        for(int irank{1}; irank < num_workers + 1; irank++){ // 0 is the root -> won`t work.
            local_array_lengths[irank] += length_per_nodes;
            offsets[irank] = offsets[irank - 1] + local_array_lengths[irank - 1];

            cut_points[irank] = cut_points[irank - 1] + local_array_lengths[irank];
        }
        cut_points[0]++;

        //timer_taskdist.pause();
        timer_comunication.restart();

        // Distributing the array among processes
        MPI_Scatterv(array,
                     local_array_lengths,
                     offsets,
                     MPI_INT,
                     local_array,
                     local_array_lengths[wrank],
                     MPI_INT,
                     0, // root
                     workers_comm);

        timer_comunication.pause();
        timer_comp.restart();

        if(wrank){
            
            if(!counter){
                //printf("%d -- worker %d\n", num_workers, wrank);
                
                //print_array(local_array, local_array_lengths[wrank], static_cast<char>('a'+wrank));


                core_merge_sort(local_array, local_array_lengths[wrank]);
                //print_array(local_array, local_array_lengths[wrank], static_cast<char>('a'+wrank));
            }
            else{
                //printf("%d -- worker %d %c\n", num_workers, wrank, static_cast<char>('a'+wrank));
                //print_array(cut_points, num_workers+1);
                //print_array(local_array, local_array_lengths[wrank], static_cast<char>('a'+wrank));
                //print_array(previous_lengths, previous_wsize, static_cast<char>('a'+wrank));
                //print_array(previous_cut_points, previous_wsize, static_cast<char>('a'+wrank));

                int beg {cut_points[wrank -1] };
                int end {cut_points[wrank] };

                //printf("%d %d %c %d\n", beg, end, static_cast<char>('a'+wrank), local_array_lengths[wrank]);
                std::size_t j{1};

                while(previous_cut_points[j] <= beg)
                    j++;

                int l;
                int m, r;
                     
                while( previous_cut_points[j] < end ){
                    l = beg;
                    m = previous_cut_points[j];
                    r = previous_cut_points[j + 1];
                    if(r > end )
                        r = end;
                    //printf("before %d %d %d %c\n", l, m ,r, static_cast<char>('a'+wrank));
                    if(l){
                        l = 0;
                        m -= (local_array_lengths[wrank - 1]);
                        r -= (local_array_lengths[wrank - 1]);
                    }
                    if(r == local_array_lengths[wrank])
                        r--;
                    
                    //printf("after %d %d %d %c\n", l, m ,r, static_cast<char>('a'+wrank));

                    merge_c(local_array, l, m, r);

                    //print_array(local_array, local_array_lengths[wrank], static_cast<char>('a'+wrank));

                    j += 1;
                }
            }
        }

        MPI_Barrier(workers_comm);
        
        timer_comp.pause();
        timer_comunication.restart();

        MPI_Gatherv(local_array,
                    local_array_lengths[wrank],
                    MPI_INT,
                    array,
                    local_array_lengths,
                    offsets,
                    MPI_INT,
                    0,
                    workers_comm);
        
        MPI_Barrier(workers_comm);

        timer_comunication.pause();
        //timer_taskdist.restart();

        //if( !world_rank)
        //    print_array(array, length, 0);

        MPI_Group_free(&workers_group);
        if(MPI_COMM_NULL != workers_comm)
            MPI_Comm_free(&workers_comm);

        for (std::size_t i = 0; i < wsize; i++){
            previous_cut_points[i] = cut_points[i];
            previous_lengths[i] = local_array_lengths[i];
        }
        previous_wsize = wsize;
        
        num_workers = static_cast<int>( static_cast<double>(num_workers) / 2);
        counter++;
    }

    MPI_Group_free(&world_group);

    //if( !world_rank ){
    //    printf("Ordenado: ");
    //    print_array(array, length, 'd');
    //}
    
    timer_total.stop();
    //timer_taskdist.stop();
    timer_comunication.stop();
    timer_comp.stop();

    //if(!world_rank){
    //    std::cout << "Elapsed time/Comunication: " << timer_comunication.elapsed() << " s" << '\n';
    //    std::cout << "Elapsed time/Task dist.: " << timer_taskdist.elapsed() << " s" << '\n';
    //    std::cout << "Elapsed time/Computation: " << timer.elapsed() << " s" << '\n';
    //    std::cout << "Elapsed time/Total: " << timer_total.elapsed() << " s" << '\n';
    //    std::cout << "Elapsed time/Overhead: " << (timer_total.elapsed() - timer.elapsed() - timer_taskdist.elapsed() - timer_comunication.elapsed()) << " s" << '\n';
    //    std::cout << std::boolalpha << verify(array, length) << '\n';
    //}
    if(!world_rank){
        std::cout << timer_comunication.elapsed() << ' '
                  << timer_comp.elapsed() << ' '
                  << timer_total.elapsed() << ' '
                  << (timer_total.elapsed() - timer_comp.elapsed() - timer_comunication.elapsed()) << ' '
                  << std::boolalpha << verify(array, length) << '\n';
    }
    delete [] local_array_lengths;
    delete [] offsets;
    delete [] cut_points;
    delete [] previous_cut_points;
    delete [] previous_lengths;
    delete [] local_array;
    delete [] workers_rank;


    MPI_Finalize();
}



int main(int argc, char ** argv){

    srand(time(NULL));
    
    std::size_t length { static_cast<std::size_t>(std::atoi(argv[argc - 1])) }; // 40000000

    int * arr{new int [length] {}};

    for (size_t i = 0; i < length; i++){
        arr[i] = rand()%30000;
    }
    
    mpi_merge_sort(argc, argv, arr, length);

    
    delete [] arr;

    

    return 0;
}

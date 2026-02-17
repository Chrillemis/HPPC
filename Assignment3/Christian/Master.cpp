#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <random>
#include <chrono>
#include <thread>
#include <array>
#include <vector>
#include <cstdlib>


const int Ntask = 400;
const int Nworker = 12;
int results_counter = 0;
int workers_finished = 0;
int max_task_value = 10;

std::vector<int> tasks(Ntask);
std::vector<int> results(Ntask);
std::thread threads[Nworker];

void print_results() {
    for (int i = 0; i < Ntask; i++) {
        std::cout << "Task " << i << ": " << tasks[i] << " Result: " << results[i] << std::endl;
    }
    std::ofstream outfile("master_results.txt");
    for (int i = 0; i < Ntask; i++) {
        outfile << "Task " << i << ": " << tasks[i] << " Result: " << results[i] << std::endl;
    }    
    std::cout << "Everything done, terminating main().\n";
    std::exit(0);
}

void master(int worker_id);
void worker(int task, int worker_id, int result_idx);


void worker(int task, int worker_id, int result_idx) {
    std::this_thread::sleep_for(std::chrono::milliseconds(task)); // simulate work
    // std::cout << "Worker done with task " << task << "ms" << std::endl;
    results[result_idx] = task; // store result in global vector
    master(worker_id); // notify master that task is done
}

void master(int worker_id) {
    // std::cout << "Master received result from worker " << worker_id << ". Total results received: " << results_counter << std::endl;
    results_counter++; // increment counter of completed tasks

    if (results_counter >= Ntask) {
        // std::cout << "All tasks completed. Final results:\n";
        workers_finished++;
        if (workers_finished >= Nworker) {
            std::cout << "All workers finished. Exiting all masters.\n";
            print_results();
        }
        return; // all tasks done, exit master
    }
    // threads[worker_id].join(); // wait for worker thread to finish
    threads[worker_id] = std::thread(worker, tasks[results_counter], worker_id, results_counter); // pass correct arguments
    threads[worker_id].detach(); // detach thread to allow it to run independently
}



int main() {
    for (int i = 0; i < Ntask; i++) {//fill tasks with random numbers
        tasks[i] = rand() % max_task_value + 1; // random task duration between 1 and max_task_value milliseconds 
    }
    for (int i = 0; i < Nworker; i++) {
        threads[i] = std::thread(worker, tasks[i], i, i); // pass correct arguments
        threads[i].detach();
        // std::cout << "Thread " << i << " spawned and detached.\n";
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(Ntask*max_task_value)); // wait for all tasks to complete (worst case)

    return 0;
}

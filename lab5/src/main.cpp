#include <stack>
#include <string>
#include <iostream>
#include <mpi.h>
#include <pthread.h>
#include <thread>
#include <random>
#include <cmath>
#include "../include/simple_logger.h"

#define COUNT_THREADS 3
#define NO_TASKS -1
#define COMMUNICATION_TAG_ONE 101
#define COMMUNICATION_TAG_TWO 102

using std::stack, std::size_t, std::stoi, std::cout, std::to_string;

pthread_mutex_t mutex;
pthread_cond_t cond;
Simple_logger logger(3);
size_t count_tasks_for_process;
size_t critical_count_tasks = 1;
size_t max_count_tasks;
unsigned count_finished_tasks = 0;
bool is_end = false;
bool is_stop_observe = false;
bool observer_works = true;
stack<unsigned> *tasks = nullptr;

size_t get_count_tasks(int argc, char* argv[])
{
    if (argc < 2 || argc >= 3)
    {  
        return 0;
    }

    return stoi(argv[1]);
}

void* thread_observer(void *args) {
    int rank = 0, size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    bool visited[size] = {false};

    while (!tasks->empty()) {
        pthread_mutex_lock(&mutex);
        pthread_cond_wait(&cond, &mutex);
        while (tasks->size() > critical_count_tasks) {
            pthread_cond_wait(&cond, &mutex);
        }
        pthread_mutex_unlock(&mutex);

        unsigned count_tasks = std::ceil(max_count_tasks / 2.0);

        logger.print("Process " + to_string(rank) + " requests " + to_string(count_tasks) + " tasks", 0);
        unsigned new_tasks[count_tasks] = {0};
        MPI_Status stat = {0, 0, 0, 0, 0};
        
        for (int i = 0; i < size; ++i) {
            if (i == rank || visited[i]) {
                continue;
            }

            MPI_Sendrecv(&count_tasks, 1, MPI_UNSIGNED, i, COMMUNICATION_TAG_ONE, new_tasks, count_tasks,
                MPI_UNSIGNED, MPI_ANY_SOURCE, COMMUNICATION_TAG_TWO, MPI_COMM_WORLD, &stat);
            visited[i] = true;
            if (new_tasks[0] != 0) {
                break;
            }
        }

        if (new_tasks[0] == 0) {
            logger.print("Process " + to_string(rank) + " doesn`t recive tasks!", 2);
            break;
        }

        pthread_mutex_lock(&mutex);
        unsigned count_recieve_tasks = 0;
        for (size_t i = 0; i < static_cast<size_t>(count_tasks); ++i) {
            if (new_tasks[i] == 0) {
                break;
            }
            tasks->push(new_tasks[i]);
            ++count_recieve_tasks;
        }
        pthread_mutex_unlock(&mutex);

        logger.print("Process " + to_string(rank) + " recieve " + to_string(count_recieve_tasks) + " tasks", 0);

        if (visited[size - 1]) {
            for (size_t i = 0; i < static_cast<size_t>(size); ++i) {
                visited[i] = false;
            }
        }
    }

    pthread_mutex_lock(&mutex);
    is_stop_observe = true;
    pthread_mutex_unlock(&mutex);

    logger.print("Process " + to_string(rank) + " stop observe!", 0);
    
    return 0;
}

void* thread_worker(void *args) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    logger.print("Process " + to_string(rank) + " started his job!", 0);

    unsigned time_to_sleep = 0;
    while (!is_stop_observe) {
        pthread_mutex_lock(&mutex);
        if (tasks->size() <= 1) {
            logger.print("Process " + to_string(rank) + " has low tasks count!", 2);
            pthread_cond_signal(&cond);
        }

        if (tasks->empty()) {
            pthread_mutex_unlock(&mutex);
            continue;
        }

        time_to_sleep = tasks->top();
        tasks->pop();
        ++count_finished_tasks;
        pthread_mutex_unlock(&mutex);

        logger.print("Process " + to_string(rank) + " sleep for " + to_string(time_to_sleep) + " seconds", 1);
        std::this_thread::sleep_for(std::chrono::seconds(time_to_sleep));
    }

    logger.print("Process " + to_string(rank) + " finished his job!", 0);

    pthread_mutex_lock(&mutex);
    is_end = true;
    pthread_cond_signal(&cond);
    pthread_mutex_unlock(&mutex);

    return 0;
}

void* thread_handler(void *args) {
    int rank = 0, total_proccess_number = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_proccess_number);

    unsigned count_tasks = 0;
    MPI_Status stat;
    MPI_Request request;
    int flag = 0;
    bool is_all_procceses_stop_observe = false;
    unsigned counter = 0;

    while (true) {
        MPI_Irecv(&count_tasks, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, COMMUNICATION_TAG_ONE, MPI_COMM_WORLD, &request);

        while (flag == 0 && !is_all_procceses_stop_observe) {
            MPI_Test(&request, &flag, &stat);
            if (counter > 30) {
                MPI_Allreduce(&is_stop_observe, &is_all_procceses_stop_observe, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
                counter = 0;
            }

            ++counter;
        }

        if (is_all_procceses_stop_observe) {
            if (flag == 0) {
                MPI_Cancel(&request);
            }
            break;
        }
        flag = 0;

        unsigned tasks_for_send[count_tasks] = {0};
        unsigned size = 0;

        if (is_end) {
            logger.print("Process " + to_string(rank) + " recive request on " + to_string(count_tasks) +
                " tasks from process " + to_string(stat.MPI_SOURCE) + " but it`s shut down", 0);
            MPI_Send(&size, 1, MPI_UNSIGNED, stat.MPI_SOURCE, COMMUNICATION_TAG_TWO, MPI_COMM_WORLD);
            continue;
        }

        logger.print("Process " + to_string(rank) + " recive request on " + to_string(count_tasks) +
            " tasks from process " + to_string(stat.MPI_SOURCE), 0);

        pthread_mutex_lock(&mutex);
        if (tasks->size() > count_tasks) {
            size = count_tasks / 2;
        } else if (tasks->size() >= count_tasks / 2) {
            size = count_tasks / 3;
        } else {
            logger.print("Process " + to_string(rank) + " don`t found free tasks!", 2);
        }

        for (size_t i = 0; i < size; ++i) {
            tasks_for_send[i] = tasks->top();
            tasks->pop();
        }
        pthread_mutex_unlock(&mutex);

        MPI_Send(tasks_for_send, count_tasks, MPI_UNSIGNED, stat.MPI_SOURCE, COMMUNICATION_TAG_TWO, MPI_COMM_WORLD);

        logger.ifprint(size != 0, "Process " + to_string(rank) + " send " + to_string(size) + " tasks", 0);
    }

    MPI_Request_free(&request);
    logger.print("Process " + to_string(rank) + " stop listen!", 0);

    return 0;
}

int main(int argc, char *argv[]) {
    int provided_support, total_processes_number, rank;
    double start, end;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided_support);

    if (provided_support != MPI_THREAD_MULTIPLE) {
        logger.print("Provided support does not match the requested!", 2);
        abort();
    }

    MPI_Comm_size(MPI_COMM_WORLD, &total_processes_number);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    size_t count_tasks = get_count_tasks(argc, argv);
    if (count_tasks < static_cast<size_t>(total_processes_number)) {
        logger.print("Count tasks less then count proccess!", 2);
        abort();
    }

    tasks = new stack<unsigned>();
    //count_tasks_for_process = count_tasks / static_cast<size_t>(total_processes_number) +
    //    (count_tasks % static_cast<size_t>(total_processes_number) > rank);
    count_tasks_for_process = (rank + 1) * 10;
    max_count_tasks = (total_processes_number - 1) * 10;

    logger.print("Count tasks: " + to_string(count_tasks_for_process) + ", critical count tasks: " +
        to_string(critical_count_tasks) + " for process " + to_string(rank), 0);

    //std::random_device rd;  // a seed source for the random number engine
    //std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    //std::uniform_int_distribution<> distrib(1, 10);
    for (size_t i = 0; i < count_tasks_for_process; ++i)
    {
        tasks->push(1);
    }

    pthread_attr_t attrs;
    if (pthread_attr_init(&attrs)) {
        logger.print("Can`t init thread attributes!", 2);
        abort();
    }

    if (pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE)) {
        logger.print("Can`t set detach state for thread!", 2);
        abort();
    }

    if (pthread_cond_init(&cond, NULL)) {
        logger.print("Can`t init condition variable!", 2);
        abort();
    }

    pthread_t thread[COUNT_THREADS];
    void* (*thread_functions[COUNT_THREADS])(void *) = {thread_observer, thread_worker, thread_handler};

    start = MPI_Wtime();
    for (int i = 0; i < COUNT_THREADS; ++i) {
        pthread_create(&thread[i], &attrs, thread_functions[i], NULL);
    }

    void *thread_exit_status;
    for (int i = 0; i < COUNT_THREADS; ++i) {
        pthread_join(thread[i], &thread_exit_status);
        /*if (*((int*)thread_exit_status) != 0) {
            logger.print("Thread failed with exit status " + to_string(*((int*)thread_exit_status)) + "!", 1);
        }*/
    }
    end = MPI_Wtime();

    pthread_attr_destroy(&attrs);
    pthread_mutex_destroy(&mutex);
    delete tasks;

    MPI_Barrier(MPI_COMM_WORLD);

    logger.print("Count finished tasks: " + to_string(count_finished_tasks) + " of " + to_string(count_tasks) +
        " in proccess " + to_string(rank), 0);
    logger.print("Job done!", 0);

    double total_time = end - start;
    double max_diff = 0;
    double min_diff = 0;
    MPI_Allreduce(&total_time, &max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&total_time, &min_diff, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    cout << "Total time in seconds: " << end - start << std::endl;

    if (rank == 0) {
        cout << "Diff: " << (max_diff - min_diff) / max_diff << std::endl;
    }

    MPI_Finalize();

    return 0;
}
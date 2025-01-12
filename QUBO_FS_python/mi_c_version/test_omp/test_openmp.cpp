#include <omp.h>
#include <iostream>

int main() {
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        std::cout << "Thread " << thread_id << std::endl;
    }
    return 0;
}


#ifndef UTILS_UTILS
#define UTILS_UTILS

#include <algorithm>
#include <future>
#include <iostream>
#include <vector>
#include <string.h>
#include <sys/stat.h> /* For stat() */
#include "types.h"

namespace io {
    inline size_t file_exists(const char *path) {
        struct stat st;
        return stat(path, &st) == 0;
    }
}

inline void execute_with_time_limit(std::function<void()> fun, uint time_limit, std::atomic<bool>& reach_time_limit)
{
    std::future future = std::async(std::launch::async, fun);
    std::future_status status;
    do {
        status = future.wait_for(std::chrono::seconds(time_limit));
        if (status == std::future_status::deferred)
        {
            std::cout << "Deferred" << std::endl;
            exit(-1);
        }
        else if (status == std::future_status::timeout)
        {
            reach_time_limit = true;
            std::cout << "Timeout " << time_limit << "s\n";
        }
    } while (status != std::future_status::ready);
}

namespace mem {
    /**
     * get peak virtual memory space of the current process
     * https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process#answer-64166
     */
    inline int parseLine(char* line){
        // This assumes that a digit will be found and the line ends in " Kb".
        int i = strlen(line);
        const char* p = line;
        while (*p <'0' || *p > '9') p++;
        line[i-3] = '\0';
        i = atoi(p);
        return i;
    }

    inline int getValue(){ //Note: this value is in KB!
        FILE* file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];

        while (fgets(line, 128, file) != NULL){
            if (strncmp(line, "VmPeak:", 7) == 0){
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
    }
}

namespace randomUtil
{
    inline int getRandomIntBetweenAandB(int a, int b){
        return (rand()%(b-a+1))+ a;
    }

    inline bool randomProbability(double probability){
        int randInt = rand();
        double randRatio = (randInt/(RAND_MAX+0.0));
        return randRatio < probability;
    }
}
#endif //UTILS_UTILS

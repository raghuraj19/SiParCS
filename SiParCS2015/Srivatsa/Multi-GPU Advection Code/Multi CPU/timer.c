#include <timer.h>

long long get_time()
{
        struct timeval curr_time;
        gettimeofday(&curr_time, 0);
        return (long long)((curr_time.tv_sec * 1000000 + curr_time.tv_usec));
}

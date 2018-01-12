#include <stdio.h>
#include <time.h>

clock_t start, end;

void print_exec_time(char *msg) {
#if 1
    printf("%s\t .... %lf seconds\n", msg, ((double)(end - start)) / CLOCKS_PER_SEC);
#endif
}

void start_timer() {
    start = clock();
}

void stop_timer() {
    end = clock();
}
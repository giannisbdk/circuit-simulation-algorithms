#include <stdio.h>
#include <time.h>
#include <math.h>

static clock_t start, end;

void print_exec_time(char *msg) {
	int mins;
	double secs;
	double elapsed_time = (double)((end - start)) / CLOCKS_PER_SEC;

	/* In case we've got more than 1 minute print it with minutes included */
	if (elapsed_time > 60.0) {
		mins = elapsed_time / 60;
		secs = fmod(elapsed_time, 60.0);
		printf("\n%s: %d minutes %.6lf seconds\n", msg, mins, secs);
	}
	else {
		printf("\n%s: %.6lf seconds\n", msg, elapsed_time);
	}
}

void start_timer() {
    start = clock();
}

void stop_timer() {
    end = clock();
}
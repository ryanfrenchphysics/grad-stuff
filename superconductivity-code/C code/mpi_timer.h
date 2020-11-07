#ifndef _MPI_TIMER_H
#define _MPI_TIMER_H

#include <sys/types.h>
#include <time.h>

time_t starttime();
void endtime(time_t starttime);

#endif /* _MPI_TIMER_H */
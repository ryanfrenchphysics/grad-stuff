/*****************************************************************
*
* This is timer programs. 
* The main program should have 
*
* #include"timer.h"
*
* and one defined variable 
*
* (time_t) time_of_start;
*
* Then the usage is like this:
* Beginning: 
* 	time_of_start = starttime();
* End:
* 	endtime(time_of_start);
*
* ******************************************************************/

#include <stdio.h>

#include "mpi_timer.h"

/*
 * starttime()
 * 
 * Gets current time and returns it.
 * 
 * Arguments:
 *      N/A
 * 
 * Returns:
 *      time_t t1   -> current time
 */

time_t starttime()
{
	time_t t1;
	(void) time(&t1);
	return t1;
}


/*
 * endtime(time_t)
 * 
 * Determines the change in time based on the input start time and prints this.
 * 
 * Arguments:
 *      time_t starttime    ->  the time to compare against
 * 
 * Returns:
 *      void
 */

void endtime(time_t starttime)
{
	int y, hh, mm, ss;
	time_t t1, t2;
	(void) time(&t2);
	t1 = starttime;
	y = (int) t2 - t1;      /* Are we casting only t2 or the result? */
	
	hh = y / 3600;
	mm = (y - hh * 3600) / 60;
	ss = y - hh * 3600 - mm * 60;

	printf("\n(Computation took %d hours %d minutes %d seconds) \n", 
			hh, mm, ss);
	printf("%s\n", ctime(&t2));
}

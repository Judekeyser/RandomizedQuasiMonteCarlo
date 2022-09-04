#include "RQMC.h"
#include <time.h>
#include <math.h>

static double kronecher_low_discrepancy_generator(const size_t n)
{
	const double golden_ratio = (1 + sqrt(5)) / 2;
	const double number = n * golden_ratio;
	return number - floor(number);
}

static double uniform_random_generator()
{
	return (double)rand() / RAND_MAX;
}

/**
-------------------------------------------------
--- TESTS SECTION -------------------------------
-------------------------------------------------

For now, it is just more convenient to keep
everything in the same file. Later on, we could
move and rethink the architecture and
how to test private concerns.

*/
#ifdef __RQMC__TEST__
#include <stdio.h>
struct Interval
{
	double min;
	double max;
	size_t kronecher_count;
	size_t uniform_count;
};

static void test_statistic_kronecher_uniform_comparison()
{
	printf("Basic statistic comparison between Uniform and Kronecher.\n");
	printf("---------------------------------------------------------\n");
	size_t max_int_counter = 5;
	for(size_t int_counter; int_counter < max_int_counter; int_counter++)
	{
		struct Interval interval = {
			.min = (double)int_counter / max_int_counter,
			.max = (double)(int_counter+1)/max_int_counter,
			.kronecher_count = 0,
			.uniform_count = 0
		};

		double measure;
		for(int j = 0; j < 100; j++) {
			measure = uniform_random_generator();
			if(measure >= interval.min && measure < interval.max)
				interval.uniform_count += 1;
			measure = kronecher_low_discrepancy_generator(j);
			if(measure >= interval.min && measure < interval.max)
				interval.kronecher_count += 1;
		}

		printf("\nInterval [%f, %f):\n\tUniform: %d\n\tKronecher: %d\n",
			interval.min, interval.max,
			(int)interval.uniform_count,
			(int)interval.kronecher_count
		);
	}
	printf("---------------------------------------------------------\n");
}

int main()
{
	test_statistic_kronecher_uniform_comparison();
	return 0;
}
#endif

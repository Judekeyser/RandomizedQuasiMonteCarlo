#include "RQMC.h"
#include <time.h>
#include <math.h>

static double rational(double x)
{
	return x - floor(x);
}

static double kronecher_low_discrepancy_generator(const size_t n)
{
	const double golden_ratio = (1 + sqrt(5)) / 2;
	const double number = n * golden_ratio;
	return rational(number);
}

static double uniform_random_generator()
{
	return (double)rand() / RAND_MAX;
}

static double rotation(const double uniform_noise, const double base)
{
	return rational(uniform_noise + base);
}

double RQMC_integral(RQMC_integrable_function f, const size_t discrepancy_strength, const size_t bootstrap_factor)
{
	double estimate = 0.0;
	double integral;
	for(size_t i = 0; i < bootstrap_factor; i++)
	{
		integral = 0.0;
		for(size_t j = 0; j < discrepancy_strength; j++)
		{
			integral += f(rotation(uniform_random_generator(), kronecher_low_discrepancy_generator(j)));
		}
		estimate += integral / discrepancy_strength;
	}
	return estimate / bootstrap_factor;
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
	size_t max_int_counter = 10;
	size_t max_sample_counter = 1000;
	for(size_t int_counter; int_counter < max_int_counter; int_counter++)
	{
		struct Interval interval = {
			.min = (double)int_counter / max_int_counter,
			.max = (double)(int_counter+1)/max_int_counter,
			.kronecher_count = 0,
			.uniform_count = 0
		};

		double measure;
		for(int j = 0; j < max_sample_counter; j++) {
			measure = uniform_random_generator();
			if(measure >= interval.min && measure < interval.max)
				interval.uniform_count += 1;
			measure = kronecher_low_discrepancy_generator(j);
			if(measure >= interval.min && measure < interval.max)
				interval.kronecher_count += 1;
		}

		printf("\nInterval [%f, %f):\n\tUniform: %d / %d\n\tKronecher: %d / %d\n",
			interval.min, interval.max,
			(int)interval.uniform_count, (int)max_sample_counter,
			(int)interval.kronecher_count, (int)max_sample_counter
		);
	}
	printf("---------------------------------------------------------\n");
}

static double test_cesaro_average_f(double x) { return x+1-sqrt(x)/2; }
static void test_cesaro_average()
{
#define __RQMC__TEST__estimate(x) ((x > 1.1777777777777 ? 1 : -1)*(x - 1.1777777777777)/1.1777777777777)
	printf("Integral of (x+1) over [0,1):");
	double integral_estimate;


	{
		const size_t discrepancy_limit = 100000;
		integral_estimate = 0.0;
		for(size_t i = 0; i < discrepancy_limit; i++)
			integral_estimate += test_cesaro_average_f(uniform_random_generator());
		integral_estimate /= discrepancy_limit;
		printf("\n\tMonte-Carlo: %f", __RQMC__TEST__estimate(integral_estimate));
	}

	{
		const size_t discrepancy_limit = 100000;
		integral_estimate = 0.0;
		for(size_t i = 0; i < discrepancy_limit; i++)
			integral_estimate += test_cesaro_average_f(kronecher_low_discrepancy_generator(i));
		integral_estimate /= discrepancy_limit;
		printf("\n\tQuasi-Monte-Carlo: %f", __RQMC__TEST__estimate(integral_estimate));
	}

	{
		integral_estimate = RQMC_integral(test_cesaro_average_f, 5000, 20);
		printf("\n\tRandomized Quasi-Monte-Carlo: %f", __RQMC__TEST__estimate(integral_estimate));
	}
#undef __RQMC__TEST__estimate
}

int main()
{
	test_statistic_kronecher_uniform_comparison();
	test_cesaro_average();
	return 0;
}
#endif

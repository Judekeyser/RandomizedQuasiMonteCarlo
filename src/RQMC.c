#include "RQMC.h"
#include <time.h>
#include <math.h>

static inline double rational(const double x)
{
	return x - floor(x);
}

static double kronecher_low_discrepancy_generator(const size_t n)
{
	const double golden_ratio = (1 + sqrt(5)) / 2;
	return (n+1) * golden_ratio;
}

static double uniform_random_generator()
{
	return (double)rand() / RAND_MAX;
}

static inline double rotation(const double uniform_noise, const double base)
{

	double shift = rational(uniform_noise + base);
	return shift;
}

double RQMC_integral(RQMC_integrable_function f, const size_t discrepancy_strength, const size_t bootstrap_factor)
{
	double estimate = 0.0;
	double integral;
	for(size_t i = 0; i < (bootstrap_factor/2); i++)
	{
		integral = 0.0;
		double random_number = uniform_random_generator();
		for(size_t j = 0; j < discrepancy_strength; j++) {
			double point = rotation(
				random_number,
				kronecher_low_discrepancy_generator(j)
			);
			integral += f(point) + f(1.0 - point);
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
	size_t split_count = 10;
	size_t max_sample_counter = 1000;

	size_t uniform_counters[split_count];
	size_t kronecher_counters[split_count];

	for(size_t i = 0; i < split_count; i++) {
		uniform_counters[i] = 0.0;
		kronecher_counters[i] = 0.0;
	}

	for(size_t j = 0; j < max_sample_counter; j++)
	{
		double uniform_measure = uniform_random_generator();
		double kronecher_measure = rational(kronecher_low_discrepancy_generator(j));

		for(size_t i = 0; i < split_count; i++) {
			if(uniform_measure >= (double)i/split_count && uniform_measure < ((double)i+1)/split_count)
				uniform_counters[i] += 1;
			if(kronecher_measure >= (double)i/split_count && kronecher_measure < ((double)i+1)/split_count)
				kronecher_counters[i] += 1;
		}
	}

	for(size_t i = 0; i < split_count; i++) {
		printf("\nInterval [%f, %f):\n\tUniform: %d\n\tKronecher: %d\n",
			(double)i/split_count, ((double)i+1)/split_count,
			(int)uniform_counters[i],
			(int)kronecher_counters[i]
		);
	}
	printf("---------------------------------------------------------\n");
}

static double test_cesaro_average_f(double x) { return sin(1.0/x)+1-sqrt(x)/2+1.0/(x+1); }
static void test_cesaro_average()
{
#define __RQMC__TEST__estimate(x) x
	printf("Integral of sin(x)+1-sqrt(x)/2+1/(1+x) over [0,1):");
	double integral_estimate;


	{
		double integral_sin = 0.5040670619069283719898561177411482296249850282126391708714331675;
		double integral_1 = 1.0;
		double integral_sqrt = 2.0/3 * pow(1.0, 3.0/2);
		double integral_ratio = log(2);
		double integral = integral_sin + integral_1 - 0.5 * integral_sqrt + integral_ratio;
		printf("\n\tCorrect value: %f", __RQMC__TEST__estimate(integral));
	}
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
			integral_estimate += test_cesaro_average_f(rational(kronecher_low_discrepancy_generator(i)));
		integral_estimate /= discrepancy_limit;
		printf("\n\tQuasi-Monte-Carlo: %f", __RQMC__TEST__estimate(integral_estimate));
	}

	{
		integral_estimate = RQMC_integral(test_cesaro_average_f, 50000, 20);
		printf("\n\tRandomized Quasi-Monte-Carlo: %f", __RQMC__TEST__estimate(integral_estimate));
	}
#undef __RQMC__TEST__estimate
}

int main()
{
	test_statistic_kronecher_uniform_comparison();
	test_cesaro_average();

	printf("\n-----------------------\n\n");
	return 0;
}
#endif

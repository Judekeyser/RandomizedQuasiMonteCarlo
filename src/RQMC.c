#include "RQMC.h"
#include <math.h>
#include <stdio.h>


static inline double uniform_random_generator()
{
	return (double)rand() / RAND_MAX;
}

static void fill_discrepancy_factor(const size_t dimension, double discrepancy_factor[dimension])
{
	double phi = 1.0;
	for(double old_phi;;) {
		old_phi = phi;
		phi = pow(1+phi, 1.0/(1+dimension));
		if(fabs(old_phi - phi) < 0.00000001)
			break;
	}

	for(size_t k = 0; k < dimension; k++)
		discrepancy_factor[k] = 1.0/pow(phi, k+1);
}

struct RQMC_IntegralResult RQMC_integral(
	RQMC_integrable_function f,
	const size_t dimension,
	const size_t discrepancy_strength,
	const size_t measurements_count
) {
	double discrepancy_factor[dimension];
	fill_discrepancy_factor(dimension, discrepancy_factor);

	double average = 0.0;
	double sum_of_squares = 0.0;
	for(size_t i = 1; i <= measurements_count; i++)
	{
		double random_shift[dimension];
		for(size_t k = 0; k < dimension; k++)
			random_shift[k] = uniform_random_generator();

		double integral = 0.0;
		for(size_t j = 1; j <= discrepancy_strength; j++) {
			double point[dimension];
			for(size_t k = 0; k < dimension; k++)
				point[k] = fmod(
					random_shift[k] + discrepancy_factor[k] * j,
					1.0
				);
			integral += (f(dimension, point) - integral) / j;
		}
		{
			const double previous_average = average;
			average += (integral - average) / i;
			sum_of_squares += (integral - previous_average)*(integral - average);
		}
	}
	struct RQMC_IntegralResult integration =
	{
		.estimate=average,
		.sq_error=sum_of_squares/(measurements_count - 1)
	};
	return integration;
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
#include <time.h>

static double test_cesaro_average_f(const size_t d, const double xs[d]) {
	const double x = xs[0];
	const double y = xs[1];
	const double z = xs[2];
	return y == 0 || x*z == 0 ? 0 : (
		sin(1.0/(y+0.0000001))+1.0/pow(x+z+0.000001, 0.33333)
	);
}
static void test_cesaro_average()
{
	printf("Integral estimate:\n\tCorrect value: 1.56139");
	{
		const size_t N = 1 << 8;
		double estimate = 0; double point[3]; size_t count = 0;
		for(size_t i = 0; i < N; i++)
		for(size_t j = 0; j < N; j++)
		for(size_t k = 0; k < N; k++) {
			point[0] = i * 1.0/N;
			point[1] = j * 1.0/N;
			point[2] = k * 1.0/N;
			count += 1;
			estimate += (test_cesaro_average_f(3, point) - estimate) / count;
		}
		printf("\n\tRegular-Grid sampling\n\t\tmean=%f",
			estimate
		);
	}
	{
		const size_t M = 1 << 24;
		struct RQMC_IntegralResult integ = RQMC_integral(test_cesaro_average_f, 3, 1, M);
		printf("\n\t*Randomized Quasi-Monte-Carlo:\n\t\tmean=%f",
			integ.estimate
		);
		double radius = sqrt(integ.sq_error) * 1.96 / sqrt(M);
		printf("\n\t\tConfidence interval 95%%: [%f, %f]\n", integ.estimate - radius, integ.estimate + radius);
	}
	{
		const size_t M = 1 << 5;
		struct RQMC_IntegralResult integ = RQMC_integral(test_cesaro_average_f, 3, 1 << 19, M);
		printf("\n\tRandomized Quasi-Monte-Carlo:\n\t\tmean=%f",
			integ.estimate
		);
		double radius = sqrt(integ.sq_error) * 2.1 / sqrt(M);
		printf("\n\t\tConfidence interval 95%%: [%f, %f]\n", integ.estimate - radius, integ.estimate + radius);
	}
}

int main()
{
	srand(time(NULL));
	test_cesaro_average();

	printf("\n-----------------------\n\n");
	return 0;
}
#endif

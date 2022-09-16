#ifndef __RQMC__
#define __RQMC__

#include <stdlib.h>

typedef double (RQMC_integrable_function)(const size_t d, const double[d]);
struct RQMC_IntegralResult
{
	const double estimate;
	const double sq_error;
};

/*
Random Quasi Monte-Carlo integral
*/
struct RQMC_IntegralResult RQMC_integral(

    /* The function to integrate on [0, 1) */
    RQMC_integrable_function f,

    /* The domain dimension */
    const size_t dimension,

    /* The maximal points in the low-discrepancy sequence to be computed */
    const size_t discrepancy_strength,

    /* The number of random measurements */
    const size_t measurements_count
);

#endif

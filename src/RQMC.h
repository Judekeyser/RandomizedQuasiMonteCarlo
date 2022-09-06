#ifndef __RQMC__
#define __RQMC__

#include <stdlib.h>

typedef double (RQMC_integrable_function)(double);
struct RQMC_IntegralResult
{
	const double estimate;
	const double error;
};

/*
Random Quasi Monte-Carlo integral
*/
struct RQMC_IntegralResult RQMC_integral(

    /* The function to integrate on [0, 1) */
    RQMC_integrable_function f,

    /* The maximal points in the low-discrepancy sequence to be computed */
    const size_t discrepancy_strength,

    /* The bootstrap factor */
    const size_t bootstrap_factor
);

#endif

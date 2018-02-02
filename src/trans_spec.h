#ifndef TRANS_SPEC_H
#define TRANS_SPEC_H


/* We have support for the following transient sources listed below */
typedef enum {
	EXP,
	SIN,
	PULSE,
	PWL
} trans_type;

/* Exponential Waveform */
typedef struct exp {
	double i1;
	double i2;
	double td1;
	double td2;
	double tc1;
	double tc2;
} exp_t;

/* Sinusodial Waveform */
typedef struct sin {
	double i1;
	double ia;
	double fr;
	double td;
	double df;
	double ph;
} sin_t;

/* Pulse Waveform */
typedef struct pulse {
	double i1;
	double i2;
	double td;
	double tr;
	double tf;
	double pw;
	double per;
} pulse_t;

/* Piecewise Linear Waveform */
typedef struct pwl {
	/* n shows how many states there are in the arrays */
	int n;
	double *t;
	double *i;
} pwl_t;

/* Struct to hold the transient spec for a given TRAN analysis */
typedef struct trans_spec {
	trans_type type;
	exp_t 	*exp;
	sin_t 	*sin;
	pulse_t	*pulse;
	pwl_t	*pwl;
} trans_spec_t;

#endif
#ifndef AC_SPEC_H
#define AC_SPEC_H

/* Struct to hold the AC spec of a source in phasor form */
typedef struct ac {
	double magnitude;
	double phase;
} ac_t;

#endif
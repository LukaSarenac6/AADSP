#ifndef _COMMON_H
#define _COMMON_H

#include <stdfix.h>

// potrebno prekopirati sa pocetka stdfix_emu.h ili ukljuciti ceo stdfix_emu.h!
#if defined(__CCC)

#include <stdfix.h>

#define FRACT_NUM(x) (x##r)
#define LONG_FRACT_NUM(x) (x##lr)
#define ACCUM_NUM(x) (x##lk)

#define FRACT_NUM_HEX(x) (x##r)

#define FRACT_TO_INT_BIT_CONV(x) (bitsr(x))
#define INT_TO_FRACT_BIT_CONV(x) (rbits(x))

#define long_accum long accum
#define long_fract long fract

#endif

/////////////////////////////////////////////////////////////////////////////////
// Constant definitions
/////////////////////////////////////////////////////////////////////////////////
#define BLOCK_SIZE 16
#define MAX_NUM_CHANNEL 8

#define MINUS_4DB FRACT_NUM(0.630957)
#define FIR_ORDER 50

// Number of channels
#define INPUT_NUM_CHANNELS 2
#define OUTPUT_NUM_CHANNELS_3_2_1 6
#define OUTPUT_NUM_CHANNELS_2_2_0 4
#define OUTPUT_NUM_CHANNELS_2_0_0 2


// Channel IDs.
#define LEFT_CH 0
#define RIGHT_CH 1
#define LS_CH 2
#define RS_CH 3
#define CENTER_CH 4
#define LFE_CH 5

#define GAINPROCESSING_ASM

/* DSP type definitions */

/* DSP short */
typedef short DSPshort;
/* DSP unsigned short */
typedef unsigned short DSPushort;
/* DSP integer */
typedef int DSPint;
/* DSP unsigned integer */
typedef unsigned int DSPuint;
/* DSP fract */
typedef fract DSPfract;
/* DSP long fract */
typedef long_fract DSPlfract;
/* DSP long acuum */
typedef long_accum DSPaccum;

#endif //_COMMON_H

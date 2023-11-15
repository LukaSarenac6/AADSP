#ifndef COMMON_H
#define COMMON_H

#include "stdfix_emu.h"
#include "fixed_point_math.h"
/* Basic constants */
/* TO DO: Move defined constants here */
#define BLOCK_SIZE 16
#define MAX_NUM_CHANNEL 8

// Number of channels
#define INPUT_NUM_CHANNELS 2
#define OUTPUT_NUM_CHANNELS_3_2_1 6
#define OUTPUT_NUM_CHANNELS_2_2_0 4
#define OUTPUT_NUM_CHANNELS_2_0_0 2

// Channel IDs. 
// Should input and output channel IDs be separated?
#define LEFT_CH 0
#define RIGHT_CH 1
#define LS_CH 2
#define RS_CH 3
#define CENTER_CH 4
#define LFE_CH 5


// Gain linear values. 
#define MINUS_4DB FRACT_NUM(0.630957)


#define FIR_ORDER 50
/////////////////////////////////////////////////////////////////////////////////
// Constant definitions
/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////

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

#endif

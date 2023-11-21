#include <stdio.h>
#include <dsplib\wavefile.h>
#include <stdfix.h>
#include <string.h>
#include "common.h"

#define INITGAINPROCESSING_ASM

// IO Buffers
__memY DSPfract sampleBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];

// Processing related variables
DSPfract preGain;
//static DSPfract variablesGain[INPUT_NUM_CHANNELS];
DSPfract limiterThreshold = 0.999;

DSPint enable;	/* enable = argv[2],  set enable to 1 to activate gainProcessing */

DSPint outputMode = 1; /* outputMode = argv[3], 0 = 3_2_1, 1 = 2_0_0, 2 = 2_2_0 */

DSPfract firCoeffs[FIR_ORDER] = {
		FRACT_NUM(-0.01017879667124649500),
		FRACT_NUM(-0.01029945035640465400),
		FRACT_NUM(-0.01041579748452664700),
		FRACT_NUM(-0.01052773398872422500),
		FRACT_NUM(-0.01063515957695487600),
		FRACT_NUM(-0.01073797784377794400),
		FRACT_NUM(-0.01083609637787606400),
		FRACT_NUM(-0.01092942686520807400),
		FRACT_NUM(-0.01101788518766649000),
		FRACT_NUM(-0.01110139151711620300),
		FRACT_NUM(-0.01117987040470101500),
		FRACT_NUM(-0.01125325086530430200),
		FRACT_NUM(-0.01132146645706179900),
		FRACT_NUM(-0.01138445535582664800),
		FRACT_NUM(-0.01144216042449541300),
		FRACT_NUM(-0.01149452927710808500),
		FRACT_NUM(-0.01154151433764454400),
		FRACT_NUM(-0.01158307289344259600),
		FRACT_NUM(-0.01161916714317273600),
		FRACT_NUM(-0.01164976423930980800),
		FRACT_NUM(-0.01167483632504865600),
		FRACT_NUM(-0.01169436056561829100),
		FRACT_NUM(-0.01170831917395590700),
		FRACT_NUM(-0.01171669943070861500),
		FRACT_NUM(0.98834396857667117000),
		FRACT_NUM(-0.01171669943070861500),
		FRACT_NUM(-0.01170831917395590700),
		FRACT_NUM(-0.01169436056561829100),
		FRACT_NUM(-0.01167483632504865600),
		FRACT_NUM(-0.01164976423930980800),
		FRACT_NUM(-0.01161916714317273600),
		FRACT_NUM(-0.01158307289344259600),
		FRACT_NUM(-0.01154151433764454400),
		FRACT_NUM(-0.01149452927710808500),
		FRACT_NUM(-0.01144216042449541300),
		FRACT_NUM(-0.01138445535582664800),
		FRACT_NUM(-0.01132146645706179900),
		FRACT_NUM(-0.01125325086530430200),
		FRACT_NUM(-0.01117987040470101500),
		FRACT_NUM(-0.01110139151711620300),
		FRACT_NUM(-0.01101788518766649000),
		FRACT_NUM(-0.01092942686520807400),
		FRACT_NUM(-0.01083609637787606400),
		FRACT_NUM(-0.01073797784377794400),
		FRACT_NUM(-0.01063515957695487600),
		FRACT_NUM(-0.01052773398872422500),
		FRACT_NUM(-0.01041579748452664700),
		FRACT_NUM(-0.01029945035640465400),
		FRACT_NUM(-0.01017879667124649500),
		FRACT_NUM(-0.01005394415493147400)
};

DSPfract history1[FIR_ORDER] = { FRACT_NUM(0.0) };
DSPfract history2[FIR_ORDER] = { FRACT_NUM(0.0) };

DSPfract saturation(DSPfract in)
{
	DSPaccum inValue = in;
	// Simple limiter since we know that pre-Gain adds 6dB
	if (inValue > limiterThreshold)
	{
		return limiterThreshold;
	}
	else if (inValue < -limiterThreshold)
	{
		return -limiterThreshold;
	}

	return in;
}

DSPfract fir_basic(DSPfract input, DSPfract* history)
{
	DSPint i;
	DSPfract retAccum = 0;
	DSPfract* coefs = firCoeffs;

	/*for (i = 0; i < FIR_ORDER - 1; i++)
	{
		*(history + FIR_ORDER - i - 1) = *(history + FIR_ORDER - i - 2);
	}*/

	for (i = FIR_ORDER - 2; i >= 0; i--)
	{
		*(history + i + 1) = *(history + i);					//history[i + 1] = history[i];
	}


	/* store input at the beginning of the delay line */
	*history = input;

	/* calc FIR */
	for (i = 0; i < FIR_ORDER; i++)
		{
			//DSPaccum mul = *(coefs + i) * *(history + i);
			//retAccum = retAccum + mul;
			retAccum = retAccum + (*coefs * *history);
			coefs++;
			history++;
		}


	return retAccum;
}

#ifdef INITGAINPROCESSING_ASM
extern void initGainProcessing(DSPfract preGainValue);
#else
void initGainProcessing(DSPfract preGainValue)
{
	preGain = preGainValue;
}
#endif

#ifdef GAINPROCESSING_ASM
extern void gainProcessing(__memY DSPfract pIn[][BLOCK_SIZE],__memY DSPfract pOut[][BLOCK_SIZE]);
#else
void gainProcessing(__memY DSPfract pIn[][BLOCK_SIZE],__memY DSPfract pOut[][BLOCK_SIZE])
{

	__memY DSPfract* samplePtrIn = *pIn;
	__memY DSPfract* samplePtrOut = *pOut;
	__memY DSPfract* lsSamplePtr = *(pOut + LS_CH);
	__memY DSPfract* rsSamplePtr = *(pOut + RS_CH);
	__memY DSPfract* centerSamplePtr = *(pOut + CENTER_CH);
	__memY DSPfract* lfeSamplePtr = *(pOut + LFE_CH);

	DSPint j;

	for (j = 0; j < BLOCK_SIZE; j++)
	{
		// first stage, apply constant pre-Gain
 		*samplePtrIn = *samplePtrIn * preGain;
		// second stage, set L channel out
		*samplePtrOut = *samplePtrIn;

		if (outputMode == 0)
		{
			// add processed sampled to the center output channel
			*centerSamplePtr = *samplePtrOut;
		}

		if (outputMode == 0 || outputMode == 2)
		{
			// apply fir on left channel
			//DSPaccum accum = fir_basic(*samplePtrIn, history1);

			//*lsSamplePtr = (DSPfract)accum;
			*lsSamplePtr = fir_basic(*samplePtrIn, history1);

			// add ls sample to lfe
			*lfeSamplePtr = *lsSamplePtr;
		}


		samplePtrIn++;
		samplePtrOut++;
		centerSamplePtr++;
		lsSamplePtr++;
		lfeSamplePtr++;
	}

	samplePtrIn = *(pIn + RIGHT_CH);
	samplePtrOut = *(pOut + RIGHT_CH);
	centerSamplePtr = *(pOut + CENTER_CH);
	lfeSamplePtr = *(pOut + LFE_CH);

	for (j = 0; j < BLOCK_SIZE; j++)
	{


		// first stage, apply constant pre-Gain
		*samplePtrIn = *samplePtrIn * (preGain);
		// second stage, set R channel out
		*samplePtrOut = *samplePtrIn;

		if (outputMode == 0)
		{
			// add processed sampled to the center output channel
			*centerSamplePtr = *centerSamplePtr  + *samplePtrOut;
		}

		if (outputMode == 0 || outputMode == 2)
		{
			// apply fir on right channel
			*rsSamplePtr = fir_basic(*samplePtrIn, history2);
			// add rs sample to lfe
			*lfeSamplePtr = *lfeSamplePtr + *rsSamplePtr;
		}

		samplePtrIn++;
		samplePtrOut++;
		centerSamplePtr++;
		rsSamplePtr++;
		lfeSamplePtr++;
	}
}
#endif

int main(int argc, char *argv[])
 {
    WAVREAD_HANDLE *wav_in;
    WAVWRITE_HANDLE *wav_out;

	char WavInputName[256];
	char WavOutputName[256];

    DSPint inChannels;
    DSPint outChannels;
    DSPint bitsPerSample;
    DSPint sampleRate;
    DSPint iNumSamples;
    DSPint i;
    DSPint j;

    outputMode = 0;//argv[3]-68;

	//init channel buffers
	for(i=0; i<MAX_NUM_CHANNEL; i++)
		for(j=0; j<BLOCK_SIZE; j++)
			sampleBuffer[i][j] = FRACT_NUM(0.0);

	for (i = 0; i < FIR_ORDER; ++i)
	{
		history1[i] = FRACT_NUM(0.0);
		history2[i] = FRACT_NUM(0.0);
	}
	// Open input wav file
	//-------------------------------------------------
	strcpy(WavInputName,argv[0]);
	wav_in = cl_wavread_open(WavInputName);
	 if(wav_in == NULL)
    {
        printf("Error: Could not open wavefile.\n");
        return -1;
    }
	//-------------------------------------------------

	// Read input wav header
	//-------------------------------------------------
	inChannels = cl_wavread_getnchannels(wav_in);
    bitsPerSample = cl_wavread_bits_per_sample(wav_in);
    sampleRate = cl_wavread_frame_rate(wav_in);
    iNumSamples =  cl_wavread_number_of_frames(wav_in);
	//-------------------------------------------------

	// Open output wav file
	//-------------------------------------------------
	strcpy(WavOutputName,argv[1]);

	if (outputMode == 1)
	{
		/* Set number of chanels to 2, left and right */
		outChannels = OUTPUT_NUM_CHANNELS_2_0_0;
	}
	else if (outputMode == 2)
	{
		/* Set number of channels to 4, left right ls and rs */
		outChannels = OUTPUT_NUM_CHANNELS_2_2_0;
	}
	else if (outputMode == 0)
	{
			/* Set number of channels to 6, left right ls rs center and lfe */
		outChannels = OUTPUT_NUM_CHANNELS_3_2_1;
	}

			/* Set number of channels to 6, left right ls rs center and lfe */

	wav_out = cl_wavwrite_open(WavOutputName, bitsPerSample, outChannels, sampleRate);
	if(!wav_out)
    {
        printf("Error: Could not open wavefile.\n");
        return -1;
    }
	//-------------------------------------------------

	initGainProcessing(MINUS_4DB);

	// Processing loop
	//-------------------------------------------------
    {
		int i;
		int j;
		int k;
		int sample;

		// exact file length should be handled correctly...
		for(i=0; i<iNumSamples/BLOCK_SIZE; i++)
		{
			for(j=0; j<BLOCK_SIZE; j++)
			{
				for(k=0; k<inChannels; k++)
				{
					sample = cl_wavread_recvsample(wav_in);
        			sampleBuffer[k][j] = rbits(sample);
				}
			}
			enable = 1;//argv[2]-68;
			if (enable == 1)
			{
				gainProcessing(sampleBuffer, sampleBuffer);
			}

			for(j=0; j<BLOCK_SIZE; j++)
			{
				for(k=0; k<outChannels; k++)
				{
					sample = bitsr(sampleBuffer[k][j]);
					cl_wavwrite_sendsample(wav_out, sample);
				}
			}
		}
	}

	// Close files
	//-------------------------------------------------
    cl_wavread_close(wav_in);
    cl_wavwrite_close(wav_out);
	//-------------------------------------------------

    return 0;
 }

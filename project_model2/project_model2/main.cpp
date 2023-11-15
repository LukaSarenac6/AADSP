
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "WAVheader.h"
#include "common.h"


// IO Buffers
static DSPfract sampleBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];

// Processing related variables
static DSPfract preGain;
static DSPfract variablesGain[INPUT_NUM_CHANNELS];
static DSPfract limiterThreshold = 0.999;

static DSPint enable = 1;	/* enable = argv[3],  set enable to 1 to activate gainProcessing */

static DSPint outputMode = 0; /* outputMode = argv[4], 0 = 3_2_1, 1 = 2_0_0, 2 = 2_2_0 */

static DSPfract firCoeffs[FIR_ORDER] = {
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
		FRACT_NUM(0.98834396857667117000,),
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

static DSPfract history1[FIR_ORDER] = { FRACT_NUM(0.0) };
static DSPfract history2[FIR_ORDER] = { FRACT_NUM(0.0) };

DSPfract saturation(DSPfract in)
{
	DSPaccum inValue = in;
	// Simple limiter since we know that pre-Gain adds 6dB
	if (inValue > limiterThreshold)
	{
		return fmin(inValue, (DSPaccum)limiterThreshold);
	}
	else if (inValue < -limiterThreshold)
	{
		return fmax(inValue, (DSPaccum)-limiterThreshold);
	}

	return in;
}

DSPfract fir_basic(DSPfract input, DSPfract* history)
{
	DSPint i;
	DSPfract retAccum = 0;
	DSPfract* coefs = firCoeffs;

	for (i = 0; i < FIR_ORDER - 1; i++)
	{
		*(history + FIR_ORDER - i - 1) = *(history + FIR_ORDER - i - 2);
	}

	/*for (i = FIR_ORDER - 2; i >= 0; i--)
	{
		*(history + i + 1) = *(history + i);					//history[i + 1] = history[i];
	}*/


	/* store input at the beginning of the delay line */
	*history = input;											

	/* calc FIR */
	for (i = 0; i < FIR_ORDER; i++)
	{
		//DSPaccum mul = *(coefs + i) * *(history + i);
		//retAccum = retAccum + mul;
		retAccum = retAccum + *(coefs + i) * *(history + i);
	}


	return retAccum;
}

void initGainProcessing(DSPfract preGainValue)
{
	preGain = DSPfract(preGainValue);
	/*for (int i = 0; i < INPUT_NUM_CHANNELS; i++)		TO DO: realizovati -beskonacno do 0, za sad mi je samo -4db
	{
		variablesGain[i] = defaultVariablesGain[i];
	}*/
	//postGain = postGainValue;
}



void gainProcessing(DSPfract pIn[][BLOCK_SIZE], DSPfract pOut[][BLOCK_SIZE])
{

	DSPfract* samplePtrIn = *pIn;
	DSPfract* samplePtrOut = *pOut;
	DSPfract* lsSamplePtr = *(pOut + LS_CH);
	DSPfract* rsSamplePtr = *(pOut + RS_CH);
	DSPfract* centerSamplePtr = *(pOut + CENTER_CH);
	DSPfract* lfeSamplePtr = *(pOut + LFE_CH);

	for (DSPint j = 0; j < BLOCK_SIZE; j++)
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

	for (DSPint j = 0; j < BLOCK_SIZE; j++)
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

	// TODO: remove upper implementation and implement processing for each channel indepenetnly 
	// (without outter noInputChannels loop, but only with inner nSamples loop)
	// (kick-out any unnecessary local variables and parameters)
}

/////////////////////////////////////////////////////////////////////////////////
// @Author	<student name>
// @Date		<date>  
//
// Function:
// main
//
// @param - argv[0] - Input file name
//        - argv[1] - Output file name
// @return - nothing
// Comment: main routine of a program
//
// E-mail:	<email>
//
/////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
	FILE* wav_in = NULL;
	FILE* wav_out = NULL;
	char WavInputName[256];
	char WavOutputName[256];
	WAV_HEADER inputWAVhdr, outputWAVhdr;
	//DSPfract defaultVariablesGain[INPUT_NUM_CHANNELS] = { MINUS_4DB, MINUS_4DB }; // -4dB, -4dB
	outputMode = atoi(argv[4]);

	// Init channel buffers
	/*for (int i = 0; i < MAX_NUM_CHANNEL; i++)
		memset(&sampleBuffer[i], 0, BLOCK_SIZE);*/
	for (DSPint i = 0; i < MAX_NUM_CHANNEL; i++)
		for (DSPint j = 0; j < BLOCK_SIZE; j++)
			sampleBuffer[i][j] = FRACT_NUM(0.0);


	for (DSPint i = 0; i < FIR_ORDER; ++i)
	{
		history1[i] = FRACT_NUM(0.0);
		history2[i] = FRACT_NUM(0.0);
	}
	// Open input and output wav files
	//-------------------------------------------------
	strcpy(WavInputName, argv[1]);
	wav_in = OpenWavFileForRead(WavInputName, "rb");
	strcpy(WavOutputName, argv[2]);
	wav_out = OpenWavFileForRead(WavOutputName, "wb");
	//-------------------------------------------------

	// Read input wav header
	//-------------------------------------------------
	ReadWavHeader(wav_in, inputWAVhdr);
	//-------------------------------------------------

	// Set up output WAV header
	//-------------------------------------------------	
	outputWAVhdr = inputWAVhdr;

	/* sets the numer of output channels to 3 L, R and CENTER*/
	if (outputMode == 1)
	{
		/* Set number of chanels to 2, left and right */
		outputWAVhdr.fmt.NumChannels = OUTPUT_NUM_CHANNELS_2_0_0;
	}
	else if (outputMode == 2)
	{
		/* Set number of channels to 4, left right ls and rs */
		outputWAVhdr.fmt.NumChannels = OUTPUT_NUM_CHANNELS_2_2_0;
	}
	else
	{
		/* Set number of channels to 6, left right ls rs center and lfe */
		outputWAVhdr.fmt.NumChannels = OUTPUT_NUM_CHANNELS_3_2_1;
	}


	int oneChannelSubChunk2Size = inputWAVhdr.data.SubChunk2Size / inputWAVhdr.fmt.NumChannels;
	int oneChannelByteRate = inputWAVhdr.fmt.ByteRate / inputWAVhdr.fmt.NumChannels;
	int oneChannelBlockAlign = inputWAVhdr.fmt.BlockAlign / inputWAVhdr.fmt.NumChannels;

	outputWAVhdr.data.SubChunk2Size = oneChannelSubChunk2Size * outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.ByteRate = oneChannelByteRate * outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.BlockAlign = oneChannelBlockAlign * outputWAVhdr.fmt.NumChannels;


	// Write output WAV header to file
	//-------------------------------------------------
	WriteWavHeader(wav_out, outputWAVhdr);

	initGainProcessing(MINUS_4DB);

	// Processing loop
	//-------------------------------------------------	
	{
		int sample;
		int BytesPerSample = inputWAVhdr.fmt.BitsPerSample / 8;
		const double SAMPLE_SCALE = -(double)(1 << 31);		//2^31
		int iNumSamples = inputWAVhdr.data.SubChunk2Size / (inputWAVhdr.fmt.NumChannels * inputWAVhdr.fmt.BitsPerSample / 8);

		// exact file length should be handled correctly...
		for (int i = 0; i < iNumSamples / BLOCK_SIZE; i++)
		{
			for (DSPint j = 0; j < BLOCK_SIZE; j++)
			{
				for (DSPint k = 0; k < inputWAVhdr.fmt.NumChannels; k++)
				{
					sample = 0; //debug
					fread(&sample, BytesPerSample, 1, wav_in);
					sample = sample << (32 - inputWAVhdr.fmt.BitsPerSample); // force signextend
					sampleBuffer[k][j] = sample / SAMPLE_SCALE;				// scale sample to 1.0/-1.0 range		
				}
			}

			if (i >= 1)
			{
				int a = 0;
			}
			enable = atoi(argv[3]);
			if (enable)
			{
				gainProcessing(sampleBuffer, sampleBuffer);
			}


			for (DSPint j = 0; j < BLOCK_SIZE; j++)
			{
				for (DSPint k = 0; k < outputWAVhdr.fmt.NumChannels; k++)
				{
					sample = sampleBuffer[k][j].toLong();	// crude, non-rounding 			
					sample = sample >> (32 - inputWAVhdr.fmt.BitsPerSample);
					fwrite(&sample, outputWAVhdr.fmt.BitsPerSample / 8, 1, wav_out);
				}
			}
			fflush(wav_out);
		}
	}

	// Close files
	//-------------------------------------------------	
	fclose(wav_in);
	fclose(wav_out);
	//-------------------------------------------------	

	return 0;
}
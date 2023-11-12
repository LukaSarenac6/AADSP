
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "WAVheader.h"

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
#define MINUS_4DB 0.630957


#define FIR_ORDER 50


// IO Buffers
static double sampleBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];

// Processing related variables
static double preGain;
static double variablesGain[INPUT_NUM_CHANNELS];
static double limiterThreshold = 0.999;

static int enable = 1;	/* enable = argv[3],  set enable to 1 to activate gainProcessing */

static int outputMode = 0; /* outputMode = argv[4], 0 = 3_2_1, 1 = 2_0_0, 2 = 2_2_0 */

static double firCoeffs[FIR_ORDER] = {
		-0.01291870014765488500,
		-0.01364924152837379000,
		-0.01436865565586297800,
		-0.01507469128206817600,
		-0.01576512341354281700,
		-0.01643776199847034700,
		-0.01709046051221975100,
		-0.01772112440390686100,
		-0.01832771936694440500,
		-0.01890827939725012500,
		-0.01946091460364177800,
		-0.01998381873596551100,
		-0.02047527639769338200,
		-0.02093366991106610600,
		-0.02135748580434841000,
		-0.02174532089240052000,
		-0.02209588792354355600,
		-0.02240802076759135600,
		-0.02268067912194498800,
		-0.02291295271477012500,
		-0.02310406498650675900,
		-0.02325337623327741000,
		-0.02336038619815178000,
		-0.02342473609868966900,
		0.97692542007070471000,
		-0.02342473609868966900,
		-0.02336038619815178000,
		-0.02325337623327741000,
		-0.02310406498650675900,
		-0.02291295271477012500,
		-0.02268067912194498800,
		-0.02240802076759135600,
		-0.02209588792354355600,
		-0.02174532089240052000,
		-0.02135748580434841000,
		-0.02093366991106610600,
		-0.02047527639769338200,
		-0.01998381873596551100,
		-0.01946091460364177800,
		-0.01890827939725012500,
		-0.01832771936694440500,
		-0.01772112440390686100,
		-0.01709046051221975100,
		-0.01643776199847034700,
		-0.01576512341354281700,
		-0.01507469128206817600,
		-0.01436865565586297800,
		-0.01364924152837379000,
		-0.01291870014765488500,
		-0.01217930026574571600
};

static double history1[FIR_ORDER] = { 0 };
static double history2[FIR_ORDER] = { 0 };

double fir_basic(double input, double* coeffs, double* history)
{
	int i;
	double ret_val = 0;
	double* coefs = firCoeffs;

	/* shift delay line */
	for (i = FIR_ORDER - 2; i >= 0; i--)
	{
		*(history + i + 1) = *(history + i);					//history[i + 1] = history[i];
	}

	/* store input at the beginning of the delay line */
	*history = input;											//history[0] = input;



	/* calc FIR */
	for (i = 0; i < FIR_ORDER; i++)
	{
		ret_val += *(coefs + i) * *(history + i);				//ret_val += coeffs[i] * history[i];
	}

	return ret_val;
}

void initGainProcessing(double preGainValue, double* defaultVariablesGain/*, double postGainValue*/)
{
	preGain = preGainValue;
	/*for (int i = 0; i < INPUT_NUM_CHANNELS; i++)		TO DO: realizovati -beskonacno do 0, za sad mi je samo -4db
	{
		variablesGain[i] = defaultVariablesGain[i];
	}*/
	//postGain = postGainValue;
}

double saturation(double in)
{
	// Simple limiter since we know that pre-Gain adds 6dB
	if (in > limiterThreshold)
	{
		return fmin(in, limiterThreshold);
	}
	else if (in < -limiterThreshold)
	{
		return fmax(in, -limiterThreshold);
	}

	return in;
}

void gainProcessing(double pIn[][BLOCK_SIZE], double pOut[][BLOCK_SIZE])
{

	double* samplePtrIn = *pIn;
	double* samplePtrOut = *pOut;
	double* lsSamplePtr = *(pOut + LS_CH);
	double* rsSamplePtr = *(pOut + RS_CH);
	double* centerSamplePtr = *(pOut + CENTER_CH);
	double* lfeSamplePtr = *(pOut + LFE_CH);

	for (int j = 0; j < BLOCK_SIZE; j++)
	{
		// first stage, apply constant pre-Gain
		*samplePtrIn = *samplePtrIn * preGain;												//pIn[LEFT_CH][j] = pIn[LEFT_CH][j] * preGain;
		// second stage, set L channel out
		*samplePtrOut = *samplePtrIn;														//pOut[LEFT_CH][j] = pIn[LEFT_CH][j];
		// add processed sampled to the center output channel
		*centerSamplePtr = *samplePtrOut;													//pOut[CENTER_CH][j] = pOut[LEFT_CH][j];
		// apply fir on left channel
		*lsSamplePtr = fir_basic(*samplePtrIn, firCoeffs, history1);				//pOut[LS_CH][j] = fir_basic(pIn[LEFT_CH][j], firCoeffs, history1, FIR_ORDER);
		// add ls sample to lfe 
		*lfeSamplePtr = *lsSamplePtr;														//pOut[LFE_CH][j] = pOut[LS_CH][j];

		samplePtrIn++;
		samplePtrOut++;
		centerSamplePtr++;
		lsSamplePtr++;
		lfeSamplePtr++;
	}

	samplePtrIn = *(pIn + RIGHT_CH);
	samplePtrOut = *(pIn + RIGHT_CH);
	centerSamplePtr = *(pOut + CENTER_CH);
	lfeSamplePtr = *(pOut + LFE_CH);

	for (int j = 0; j < BLOCK_SIZE; j++)
	{


		// first stage, apply constant pre-Gain
		*samplePtrIn = *samplePtrIn * preGain;												//pIn[RIGHT_CH][j] = pIn[RIGHT_CH][j] * preGain;
		// second stage, set R channel out
		*samplePtrOut = *samplePtrIn;														//pOut[RIGHT_CH][j] = pIn[RIGHT_CH][j];
		// add processed sampled to the center output channel
		*centerSamplePtr += *samplePtrOut;													//pOut[CENTER_CH][j] += pOut[RIGHT_CH][j];
		// apply fir on right channel
		*rsSamplePtr = fir_basic(*samplePtrIn, firCoeffs, history2);				//pOut[RS_CH][j] = fir_basic(pIn[RIGHT_CH][j], firCoeffs, history2, FIR_ORDER);
		// add rs sample to lfe 
		*lfeSamplePtr += *rsSamplePtr;														//pOut[LFE_CH][j] += pOut[RS_CH][j];

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
	double defaultVariablesGain[INPUT_NUM_CHANNELS] = { MINUS_4DB, MINUS_4DB }; // -4dB, -4dB
	outputMode = atoi(argv[4]);

	// Init channel buffers
	for (int i = 0; i < MAX_NUM_CHANNEL; i++)
		memset(&sampleBuffer[i], 0, BLOCK_SIZE);

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

	initGainProcessing(MINUS_4DB, defaultVariablesGain/*, MINUS_12DB */);

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
			for (int j = 0; j < BLOCK_SIZE; j++)
			{
				for (int k = 0; k < inputWAVhdr.fmt.NumChannels; k++)
				{
					sample = 0; //debug
					fread(&sample, BytesPerSample, 1, wav_in);
					sample = sample << (32 - inputWAVhdr.fmt.BitsPerSample); // force signextend
					sampleBuffer[k][j] = sample / SAMPLE_SCALE;				// scale sample to 1.0/-1.0 range		
				}
			}
			enable = atoi(argv[3]);
			if (enable)
			{
				gainProcessing(sampleBuffer, sampleBuffer);
			}


			for (int j = 0; j < BLOCK_SIZE; j++)
			{
				for (int k = 0; k < outputWAVhdr.fmt.NumChannels; k++)
				{
					sample = sampleBuffer[k][j] * SAMPLE_SCALE;	// crude, non-rounding 			
					sample = sample >> (32 - inputWAVhdr.fmt.BitsPerSample);
					fwrite(&sample, outputWAVhdr.fmt.BitsPerSample / 8, 1, wav_out);
				}
			}
		}
	}

	// Close files
	//-------------------------------------------------	
	fclose(wav_in);
	fclose(wav_out);
	//-------------------------------------------------	

	return 0;
}
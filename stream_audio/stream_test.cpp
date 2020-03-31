#include <stdio.h>
#include <stdlib.h>
#include "portaudio.h"
/* #define SAMPLE_RATE  (17932) // Test failure to open with this value. */
// #define SAMPLE_RATE (44100)
#define SAMPLE_RATE (16000)
#define FRAMES_PER_BUFFER (512)
#define NUM_SECONDS (5)
#define NUM_CHANNELS (8)
/* #define DITHER_FLAG     (paDitherOff) */
#define DITHER_FLAG (0) /**/
/** Set to 1 if you want to capture the recording to a file. */
#define WRITE_TO_FILE (1)
/* Select sample format. */
#if 1
#define PA_SAMPLE_TYPE paFloat32
typedef float SAMPLE;
#define SAMPLE_SILENCE (0.0f)
#define PRINTF_S_FORMAT "%.8f"
#elif 1
#define PA_SAMPLE_TYPE paInt16
typedef short SAMPLE;
#define SAMPLE_SILENCE (0)
#define PRINTF_S_FORMAT "%d"
#elif 0
#define PA_SAMPLE_TYPE paInt8
typedef char SAMPLE;
#define SAMPLE_SILENCE (0)
#define PRINTF_S_FORMAT "%d"
#else
#define PA_SAMPLE_TYPE paUInt8
typedef unsigned char SAMPLE;
#define SAMPLE_SILENCE (128)
#define PRINTF_S_FORMAT "%d"
#endif
typedef struct
{
    int frameIndex; /* Index into sample array. */
    int maxFrameIndex;
    SAMPLE *recordedSamples;
} paTestData;
/* This routine will be called by the PortAudio engine when audio is needed.
** It may be called at interrupt level on some machines so don't do anything
** that could mess up the system like calling malloc() or free().
*/
static int recordCallback(const void *inputBuffer, void *outputBuffer,
                          unsigned long framesPerBuffer,
                          const PaStreamCallbackTimeInfo *timeInfo,
                          PaStreamCallbackFlags statusFlags,
                          void *userData)
{
    paTestData *data = (paTestData *)userData;
    const SAMPLE *rptr = (const SAMPLE *)inputBuffer;
    SAMPLE *wptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    long framesToCalc;
    long i;
    int finished;
    unsigned long framesLeft = data->maxFrameIndex - data->frameIndex;
    (void)outputBuffer; /* Prevent unused variable warnings. */
    (void)timeInfo;
    (void)statusFlags;
    (void)userData;
    if (framesLeft < framesPerBuffer)
    {
        framesToCalc = framesLeft;
        finished = paComplete;
    }
    else
    {
        framesToCalc = framesPerBuffer;
        finished = paContinue;
    }
    if (inputBuffer == NULL)
    {
        for (i = 0; i < framesToCalc; i++)
        {
            *wptr++ = SAMPLE_SILENCE; /* left */
            if (NUM_CHANNELS == 2)
                *wptr++ = SAMPLE_SILENCE; /* right */
        }
    }
    else
    {
        for (i = 0; i < framesToCalc; i++)
        {
            *wptr++ = *rptr++; /* left */
            if (NUM_CHANNELS == 2)
                *wptr++ = *rptr++; /* right */
        }
    }
    data->frameIndex += framesToCalc;
    return finished;
}
/* This routine will be called by the PortAudio engine when audio is needed.
** It may be called at interrupt level on some machines so don't do anything
** that could mess up the system like calling malloc() or free().
*/
static int playCallback(const void *inputBuffer, void *outputBuffer,
                        unsigned long framesPerBuffer,
                        const PaStreamCallbackTimeInfo *timeInfo,
                        PaStreamCallbackFlags statusFlags,
                        void *userData)
{
    paTestData *data = (paTestData *)userData;
    SAMPLE *rptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    SAMPLE *wptr = (SAMPLE *)outputBuffer;
    unsigned int i;
    int finished;
    unsigned int framesLeft = data->maxFrameIndex - data->frameIndex;
    (void)inputBuffer; /* Prevent unused variable warnings. */
    (void)timeInfo;
    (void)statusFlags;
    (void)userData;
    if (framesLeft < framesPerBuffer)
    {
        /* final buffer... */
        for (i = 0; i < framesLeft; i++)
        {
            *wptr++ = *rptr++; /* left */
            if (NUM_CHANNELS == 2)
                *wptr++ = *rptr++; /* right */
        }
        for (; i < framesPerBuffer; i++)
        {
            *wptr++ = 0; /* left */
            if (NUM_CHANNELS == 2)
                *wptr++ = 0; /* right */
        }
        data->frameIndex += framesLeft;
        finished = paComplete;
    }
    else
    {
        for (i = 0; i < framesPerBuffer; i++)
        {
            *wptr++ = *rptr++; /* left */
            if (NUM_CHANNELS == 2)
                *wptr++ = *rptr++; /* right */
        }
        data->frameIndex += framesPerBuffer;
        finished = paContinue;
    }
    return finished;
}
/*******************************************************************/
int main(void);
int main(void)
{
    PaStreamParameters inputParameters, outputParameters;
    PaError err = paNoError;
    PaStream *stream;
    /* -- initialize PortAudio -- */
    err = Pa_Initialize();
    if (err != paNoError)
        goto done;
    /* -- setup input and output -- */
    inputParameters.device = 6; /* default input device */
    inputParameters.channelCount = NUM_CHANNELS;
    inputParameters.sampleFormat = PA_SAMPLE_TYPE;
    inputParameters.suggestedLatency = Pa_GetDeviceInfo(inputParameters.device)->defaultHighInputLatency;
    inputParameters.hostApiSpecificStreamInfo = NULL;
    outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
    outputParameters.channelCount = NUM_CHANNELS;
    outputParameters.sampleFormat = PA_SAMPLE_TYPE;
    outputParameters.suggestedLatency = Pa_GetDeviceInfo(outputParameters.device)->defaultHighOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = NULL;
    /* -- setup stream -- */
    err = Pa_OpenStream(
        &stream,
        &inputParameters,
        &outputParameters,
        SAMPLE_RATE,
        FRAMES_PER_BUFFER,
        paClipOff, /* we won't output out of range samples so don't bother clipping them */
        NULL,      /* no callback, use blocking API */
        NULL);     /* no callback, so no callback userData */
    if (err != paNoError)
        goto done;
    /* -- start stream -- */
    err = Pa_StartStream(stream);
    if (err != paNoError)
        goto done;
    printf("Wire on. Will run one minute.\n");
    fflush(stdout);
    float sampleBlock[FRAMES_PER_BUFFER * NUM_CHANNELS];
    /* -- Here's the loop where we pass data from input to output -- */
    for (int i = 0; i < (60 * SAMPLE_RATE) / FRAMES_PER_BUFFER; ++i)
    {
        err = Pa_WriteStream(stream, sampleBlock, FRAMES_PER_BUFFER);
        if (err)
            goto done;
        err = Pa_ReadStream(stream, sampleBlock, FRAMES_PER_BUFFER);
        if (err)
            goto done;
    }
    /* -- Now we stop the stream -- */
    err = Pa_StopStream(stream);
    if (err != paNoError)
        goto done;
    /* -- don't forget to cleanup! -- */
    err = Pa_CloseStream(stream);
    if (err != paNoError)
        goto done;
    Pa_Terminate();
done:
    Pa_Terminate();
    // if (data.recordedSamples) /* Sure it is NULL or valid. */
    //     free(data.recordedSamples);
    if (err != paNoError)
    {
        fprintf(stderr, "An error occured while using the portaudio stream\n");
        fprintf(stderr, "Error number: %d\n", err);
        fprintf(stderr, "Error message: %s\n", Pa_GetErrorText(err));
        err = 1; /* Always return 0 or 1, but no other return codes. */
    }
    return err;
}
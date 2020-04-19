#include <cstdio>
#include <cstring>

#include <sndfile.hh>
#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#define BUFFER_LEN 1024

using namespace std;
using namespace Eigen;

MatrixXd readfile(SNDFILE *infile, SF_INFO &sfinfo)
{
    // // char fname[] = "/home/kienpt/Documents/Beam/data/speaker0250-0020.wav";
    // // SF_INFO sfinfo;
    // SNDFILE *f = sf_open(fname, SFM_READ, &sfinfo);

    int channels = sfinfo.channels;
    int num_samples = sfinfo.frames;

    double buf[BUFFER_LEN];
    sf_count_t frames;
    int k, m, readcount;
    frames = BUFFER_LEN / channels;

    MatrixXd samples(num_samples, channels);

    int count = 0;
  
    while ((readcount = sf_readf_double(infile, buf, frames)) > 0)
    {
        for (k = 0; k < readcount; k++)
        {
            for (m = 0; m < channels; m++)
            {
                samples(count, m) = buf[k * channels + m];
            }
            count++;
        }
    }
    
    return samples;
}




void writefile(SNDFILE *outfile, SF_INFO &sfinfo, MatrixXd samples)
{

    SndfileHandle file;
    int channels = samples.cols();
    int num_samples = samples.rows();
    double buf[num_samples*channels];
    int frames = num_samples;

    //file = SndfileHandle (file_name, SFM_WRITE, SF_FORMAT_WAV, channels, sample_rate) ;
	
    for (int k = 0; k < num_samples; k++)
    {
        for (int m = 0; m < channels; m++)
        {
            buf[k * channels + m] = samples(k, m);
        }
    }

    sf_writef_double (outfile, buf, frames);
}


main(void)
{
    const char *fname = "/home/kienpt/Documents/Beam/data/speaker0250-0020.wav";

    SF_INFO sfinfo;
    SNDFILE *f = sf_open(fname, SFM_READ, &sfinfo);
    MatrixXd samples = readfile(f, sfinfo);
    int sample_rate = sfinfo.samplerate;
    sf_close(f);
    
    SF_INFO sfinfo_out;
    sfinfo_out.channels = samples.cols();
    sfinfo_out.format = SF_FORMAT_WAV;
    sfinfo_out.frames = samples.rows();
    sfinfo_out.samplerate = sample_rate;
    const char *fname_out = "/home/kienpt/Documents/Beam/data/speaker0250-0020_test.wav";
    SNDFILE *f_out = sf_open(fname_out, SFM_WRITE, &sfinfo_out);
    writefile(f_out, sfinfo_out, samples);
    sf_close(f_out);
    return 0;
} /* main */

#include <iostream>
#include <Eigen/Dense>
#include <tuple>
#include <math.h>

#include <sndfile.h>

using namespace std;
using namespace Eigen;


int main(int argc, char *argv[])
{
    string infile_name = "/path/to/vocal2.wav";

    // Open input file.
    SndfileHandle infile_handle(infile_name);
    if (!infile_handle || infile_handle.error() != 0)
    {
        cerr << "Unable to read " << infile_name << endl;
        cerr << infile_handle.strError() << endl;
        return 1;
    }

    // Show file stats
    int64_t in_frames = infile_handle.frames();
    int in_channels = infile_handle.channels();
    int in_samplerate = infile_handle.samplerate();
    cerr << "Input file: " << infile_name << endl;
    cerr << " * Frames      : " << setw(6) << in_frames << endl;
    cerr << " * Channels    : " << setw(6) << in_channels << endl;
    cerr << " * Sample Rate : " << setw(6) << in_samplerate << endl;

    // Read audio data as float
    vector<float> in_data(in_frames * in_channels);
    infile_handle.read(in_data.data(), in_data.size());
}
#include <iostream>
#include <Eigen/Dense>
#include <tuple>
#include <math.h>
#include <kfr/io.hpp>

#define KFR_ENABLE_FLAC 1

using namespace std;
using namespace Eigen;
using namespace kfr;

int main(int argc, char *argv[])
{
    audio_reader_wav<float> reader(open_file_for_reading("/home/kienpt/Documents/Beam/data/test.wav"));
    univector2d<float> audio = reader.read_channels();
}
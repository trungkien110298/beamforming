#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <soundfile/soundfile.h>



using namespace std;
using namespace Eigen;


int main()
{

   SoundFileRead  insound("/home/kienpt/Documents/Beam/data/speaker0250-0020.wav");
   //SoundFileWrite outsound(outputname, insound);

   int i, channel;
   for (i=0; i<insound.getSamples(); i++) {
      for (channel=0; channel < insound.getChannels(); channel++) {
        //  outsound.writeSampleDouble(insound.getCurrentSampleDouble(channel));
        cout << insound.getCurrentSampleDouble(channel);
      }
      insound.incrementSample(); 
   }

   return 0;
    return 0;
}
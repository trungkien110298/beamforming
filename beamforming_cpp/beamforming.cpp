#include <iostream>
#include <stdio.h>
#include "beamforming.h"

using namespace std;

int main(int argc, char *argv[]){
    if (argc < 3) return 1;
    double time = beamforming(argv[1], argv[2]);
    // char in[] = "../data/speaker0250-0020.wav";
    // char out[] = "../data/test.wav";
    // double time  = beamforming(in, out);
    return 0;
}
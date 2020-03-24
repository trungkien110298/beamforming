#include <iostream>
#include <stdio.h>
#include "beamforming.h"

using namespace std;

int main(int argc, char *argv[]){
    if (argc < 3) return 1;

    // agrv[1] -- path input file
    // agrv[2] -- path output file
    double time = beamforming(argv[1], argv[2]);
    return 0;
}    
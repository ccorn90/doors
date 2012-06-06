// doorDetect.cpp
// Command-line utility for detecting doors in images using DoorDetector class
// May 2012
// @author Chris Cornelius

#include <iostream>
#include <fstream>
#include <string.h>
using namespace std;

#include "DoorDetectorEriolMods.hpp"
#include "DoorDetector.hpp"
#include "DoorObject.hpp"

void usage() {
  cerr << "usage: doorDetect PPMFILE OFFFILE [OUTPUTFILE]" << endl;
}

int main(int argc, char** argv)
{
  string ppm_fname;
  string off_fname;
  string door_fname;

  // check args
  if(argc < 3) {
    usage();  return 0;
  }
  else { // load filenames
    ppm_fname = argv[1];
    off_fname = argv[2];
    if(argc >= 4) {  // if we have the name given, set up the name for the file we'll write
      door_fname = argv[3];
    }
    else { // otherwise, write to same directory as OFF file
      // TODO
    }
  }

  // create detector object and load file
  cout << "Loading OFF and PPM files" << endl;

  // TODO: check both files for loadable and fail nicely if not

  DoorDetector D;
  Image img (ppm_fname.c_str());
  D.load(img, readOffFile(off_fname));

  // do detection
  cout << "Running detection algorithm" << endl;
  vector<DoorObject> doors = D.doors();
  cout << "Detected " << doors.size() << " doors." << endl;

  // write doors file
  cout << "Writing doors file: " << door_fname << endl;
  writeDoorFile(door_fname, doors);
  cout << "Done!" << endl;
  
  return 0;
}

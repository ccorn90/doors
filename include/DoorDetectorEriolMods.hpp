// DoorDetectorEriolMods.hpp
// Contains modifications to eriolHeader needed by the DoorDetector, LineSegment, and DoorObject classes
// April 2012
// @author Chris Cornelius, Ian McGinnis

#ifndef _DOOR_DETECTOR_ERIOL_EXTENSIONS_
#define _DOOR_DETECTOR_ERIOL_EXTENSIONS_

#include <iostream>
#include <fstream>
using std::ifstream; using std::cerr; using std::cout; using std::endl;

#ifdef MACOSX
#include<GLUT/glut.h>
#else
#include<GL/glut.h>
#endif

#include "eriolHeader.h"

// Prototypes of additional tools for CPolys, Coords and Colors
bool const operator==(const Coord &o, const Coord &p);
bool const operator!=(const Coord &o, const Coord &p);
bool const operator<(const Coord &o, const Coord &p);
Color operator*(const Color& c, double d);
double getArea(const CPoly& c);

// Method for calculating grey value of a pixel (naive mode is mean average of channels)
static inline double greyValue(Color c)
{
  return (c.r + c.g + c.b) / 3.0;
}

// Methods for working with OFF files... TODO: supercede with polytools
vector<CPoly> readOffFile(string filepath);
vector<CPoly> makePoly (vector<Coord> &vert, vector <vector <int> > &face);
int strToInt(string & strVal);

// Namespace and methods for working with lots of CPoly overlays
namespace CPolyOverlay {
  void set(int, size_t windownum=0);
  void add(const vector<CPoly>&, const string& message="unnamed poly group");

  extern vector < vector<CPoly> > polygonOverlays;
  extern vector <string> polygonOverlayNames;
  extern int currPolyIndex;
  extern bool verbose;
  const static size_t NUM_WINDOWS = 4;
};

#endif // _DOOR_DETECTOR_ERIOL_EXTENSIONS_

// DoorDetector.hpp
// Class which is for detecting doors -- unification of Ian and Chris's work on door detection
// May 2012
// @author Ian McGinnis, Chris Cornelius

#ifndef __DOOR_DETECTOR_HPP__
#define __DOOR_DETECTOR_HPP__

#include "DoorDetectorEriolMods.hpp" // For Image, CPoly, etc
#include "omp.h"

#include "LineSegment.hpp"
#include "DoorObject.hpp"

#include <vector>
#include <set>
using std::set; using std::vector;

class DoorDetector {
private:
  /* private variables and cached status */

  vector<CPoly> offPolygons;    // initially given information
  Image rawImg;

  Image gradientImg;    // The gradient image
  bool  gradientImg_cached;

  Image binaryImg;    // The binary edge map
  bool  binaryImg_cached;
  
  /* private methods */
 
  // for calculation of binary images
  double sobelX(Image& i, PixelLoc p);
  double sobelY(Image& i, PixelLoc p);
  Image* calcSobel(Image* img);
  Image* makeBinaryEdgeMap(Image* img);
  void Hysteresis(Image* img, Image* img2, int upperThreshold, int lowerThreshold,
		  Color truePixelColor=BLUE, Color hysteresisPixelColor=BLUE);

  // for working with edges of polygons
  bool isStrongEdge(double thresh, Coord A, Coord B, Image* edgeImage);
  vector<CPoly> whichPolysHaveStrongEdges(double thresh, const vector<CPoly> & polygons, Image* edgeImage);
  vector<LineSegment> whichSegmentsHaveStrongEdges(double thresh, const vector<CPoly> & polygons, Image* edgeImage);
  vector<LineSegment> findStraightSegments(vector<LineSegment> strongEdges, vector<CPoly> strongEdgePolys,
					   Image* binaryMap, bool returnOnlyMergedSegments=true);
  
  // for forming shapes from the data and generating candidate doors
  vector<DoorObject> findAllQuadrilaterals(vector<LineSegment> segments);

  // for filtering the candidates
  double cumulativeLength(DoorObject&, double);
  vector<DoorObject> doFiltering(vector<DoorObject> candidates);
  
  // for capturing polygons inside a DoorObject
  vector<DoorObject> captureInteriorPolys(vector<DoorObject> nakedCandidates, vector<CPoly> offPolys);
  bool vertLineCheckLeft(Coord A, Coord B, Coord Z);
  bool vertLineCheckRight(Coord A, Coord B, Coord Z);
  bool horizLineCheckBelow(Coord A, Coord B, Coord Z);
  bool horizLineCheckAbove(Coord A, Coord B, Coord Z);
  bool withinBounds(DoorObject Door, CPoly Candidate);
  
  
public:
  DoorDetector() : gradientImg_cached(false), binaryImg_cached(false) { }
  
  int load(Image& img, vector<CPoly>& offPolys);  // sets up the initial information, clears cache

  vector<DoorObject> doors(bool recompute=false);  // returns the doors within the image.
  
  // methods to return meta-information if you want them (for display, maybe)
  Image gradient(bool recompute=false);   // returns the gradient image
  Image binarymap(bool recompute=false);  // returns the binary edge map of the image
  vector<LineSegment> strongEdges(bool recompute=false);  // returns the strong edges of the image

  // a demo method which does all sorts of eriol-display stuff - defined in DoorDetector.cpp
  friend void doorDetectorDemo();
};

// declaring demo function
void doorDetectorDemo();

#endif // __DOOR_DETECTOR_HPP__

// LineSegment.hpp
// Class representing a line segment in the door-detection algorithm
// March 2012
// @author Ian McGinnis, Chris Cornelius
// With much credit to Fran\cois Guiot, who developed many line processing functions adapted here.

#ifndef __LINESEGMENT_HPP__
#define __LINESEGMENT_HPP__

#include "DoorDetectorEriolMods.hpp"

#include <cmath>
#include <string>
#include <sstream>
#include <ostream>
#include <limits>
using std::string; using std::stringstream; using std::ostream; using std::isinf;

// Defining null cases
const static CPoly noPoly;  // an empty CPoly



struct LineSegment {
private:
  bool _strong_cached;  // have we stored a cached value for _strong?
  CPoly _parent;        // the CPoly which this edge lives on
  
public:  
  Coord A, B;    // endpoints as coordinates
  bool _strong; // is this a strong edge? (computable from a binary edge map) // TODO - implement as variable, calculated once, instead of method
  Color color;
  
  // Constructors
  LineSegment() :  _strong_cached(false), _parent(noPoly), A(noCoord), B(noCoord), _strong(false), color(noColor) { }
  LineSegment(Coord a, Coord b, CPoly parent=noPoly, bool strong=false, Color col=WHITE) : _strong_cached(strong), _parent(parent), A(a), B(b), _strong(strong), color(col) { normalize(); }
  LineSegment(const LineSegment& l) : _strong_cached(l._strong_cached), A(l.A), B(l.B), _strong(l._strong), color(l.color)  { normalize(); }
  LineSegment(const CPoly& poly);  // TODO - initialize from a CPoly

  // For lines from file
  int init(const string line); // build object from a line in a file, does not initiate _parent!
  string line(); // writes equivalent line to a string
  
  // make sure that coordinate A is in normal situation with respect to B
  inline void normalize() {
    if(A <  B) return;  // case for correct orientation
    if(A == B) return;  // case for equivalent endpoints
    Coord temp = A; A = B; B = temp;  // otherwise, swap coordinates
  }
    
  string toString() {
    stringstream ss;
    ss << "[ (" << A << ") (" << B << ") m=" << slope() << " ]";
    return ss.str();
  }
  friend ostream& operator<<(ostream& os, LineSegment& l) {
    os << l.toString();  return os;
  }
  
  double slope() const;          // get rise/run of this segment
  double yIntercept() const;     // get B for the line through this segment
  double length() const;         // length of this segment
  bool onSegment(Coord C, double tolerance=0.1) const; // does C lie on this segment?
  bool onLine(Coord C, double tolerance=0.1) const;    // does C lie on the line through this segment?

  // for use with the STL sort() method to compare by various factors
  const static bool compareBySlope(const LineSegment& a, const LineSegment& b) {
    return ( a.slope() < b.slope() );
  }
  const static bool compareByIntercept(const LineSegment& a, const LineSegment& b) {
    return ( a.yIntercept() < b.yIntercept() );
  }
  const static bool compareByLength(const LineSegment& a, const LineSegment& b) {
    return ( a.length() < b.length() );
  }

  LineSegment& operator= (const LineSegment& l); // for copying
  bool operator== (const LineSegment& l) const;  // for comparing
  bool colinear   (const LineSegment& l) const;  // for comparing
  bool operator<  (const LineSegment& l) const;  // for comparing

  CPoly parentPoly() { // returns the polygon of which this is a member edge
    return _parent;
  }

  CPoly asPoly() { return asPoly(noColor); }
  CPoly asPoly(Color bdyColor);  // makes a line-type CPoly out of this object

  bool isStrong(double thresh=0.1, Image* binaryEdgeMap=NULL, bool recompute=false); // checks against the binary edge map to see if this is strong
  
  // do two LineSegments intersect?  Where?  Do we extrapolate from the ends of the lines?
  friend Coord findIntersection(LineSegment& one, LineSegment& two, bool extrapolate=true);
  friend Coord findIntersection(LineSegment& one, LineSegment& two, double extrapolatePercentage);

  // compute the angle between two edges... useful for checking merge between edges of two polys
  friend double angleBetween(const LineSegment& one, const LineSegment& two);

  // transform a CPoly into its component edges
  friend vector<LineSegment> polyToSegments(const CPoly& poly);
  
};

const static LineSegment noSegment(noCoord, noCoord, noPoly, false, noColor); // a nonexistant line

// External declarations of LineSegment's friend functions
Coord findIntersection(LineSegment& one, LineSegment& two, bool extrapolate);
Coord findIntersection(LineSegment& one, LineSegment& two, double extrapolatePercentage);
double angleBetween(const LineSegment& one, const LineSegment& two);
vector<LineSegment> polyToSegments(const CPoly& poly);

// global method to convert
vector<CPoly> asPolys(vector<LineSegment> l);

#endif //__LINESEGMENT_HPP__

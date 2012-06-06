// DoorObject.hpp
// Class representing a door, door candidate, or simply a quadrilateral
// April 2012
// @author Ian McGinnis, Chris Cornelius

#ifndef __DOOR_OBJECT__
#define __DOOR_OBJECT__

#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <utility>
#include <iostream>
using std::vector; using std::set; using std::pair; using std::cerr; using std::endl;
using std::string; using std::stringstream;

#include "DoorDetectorEriolMods.hpp"

#include "LineSegment.hpp"

class DoorObject {
protected:
  vector<LineSegment> _segments; // the line segments making up this quadrilateral
  vector<CPoly> _polys;          // the polygons contained inside this quadrilateral
  
  // a cached record of the corners of the polygon
  vector<Coord> _corners; 
  bool _corners_cached;
  
  // a cached record of the "cost" of this door - higher means it looks less like a door
  double _cost;
  bool   _cost_cached;
  
  // an interior struct for search trees of corner nodes
  // methods defined in DoorObject.cpp
  struct corner_node_t {
    Coord point;
    vector<corner_node_t*> elt;
    bool corner;  corner_node_t() : corner(false) { }

    bool operator==(const corner_node_t& s) {  return (s.point == point);  }
    bool hasCycleOfFour(); // most important check for a node

    // returns pointer to an array of nodes and the number of nodes
    static pair<corner_node_t*,size_t> buildTree(vector<LineSegment> segments);

    // so that sort(...) method can work
    static bool compareNodesByCoord (const corner_node_t* y, const corner_node_t* x) {
      return (y->point < x->point);
    }
  };
  

public:
  // constructors, etc
  DoorObject(const vector<LineSegment>& segments, const vector<CPoly>& polys) : _segments(segments), _polys(polys), _corners_cached(false), _cost_cached(false) { }
  DoorObject(const LineSegment& a, const LineSegment& b, const LineSegment& c, const LineSegment& d); // with segments bounding the door
  DoorObject(const DoorObject& src) : _corners_cached(false), _cost_cached(false) {
    _segments = src._segments;
    _polys = src._polys;
  }
  
  friend ostream& operator<<(ostream& os, DoorObject& d) {
    os << "Door: ";
    for(size_t i = 0; i < d._segments.size(); i++)
      os << d._segments[i] << " ";
    os << d._polys.size() << " polys.";
    return os;
  }

  // accessors
  vector<LineSegment> segments() const{ return _segments; }
  vector<CPoly> polys() const { return _polys; }
  void addPoly(const CPoly& p) { _polys.push_back(p); }
  vector<Coord> getCorners(bool recompute=false); // returns the corners of the door
  CPoly asPoly();
  
  // comparators
  bool operator< (DoorObject& d);
  

  // methods used to calculate costs for filtering of candidates
  static bool compareByCost(DoorObject a, DoorObject b) {
    return ( a.getCost() < b.getCost() );
  }
  double getCost(bool recompute=false);
  double outsidePoint_cost();  // returns distance furthest point from inside of image
  double cornerAngles_cost();  // returns "astigmatism" of door - mean difference between outside angles and such.
  double geometric_cost();     // returns comparison made in length in sides of door
  
  // to initialize or write many DoorObjects using a .doors file
  friend vector<DoorObject> readDoorFile(string filename);
  friend void writeDoorFile(string filename, const vector<DoorObject>& doors);

  // testing functions which recieve all access
  friend void test_corner_node_t();
};

// prototypes for file-access friend functions - defined in DoorObject.cpp
vector<DoorObject> readDoorFile(string filename);
void writeDoorFile(string filename, vector<DoorObject> doors);

// prototype for asPolys function for vectors of door objects
vector<CPoly> asPolys(vector<DoorObject> l);
vector<CPoly> asFilledPolys(vector<DoorObject> l);

#endif //__DOOR_OBJECT__

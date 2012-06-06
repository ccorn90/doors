// DoorObject.cpp
// Class definitions for DoorObject object
// March 2012
// @author Ian McGinnis, Chris Cornelius

#include "DoorObject.hpp"

DoorObject::DoorObject(const LineSegment& a, const LineSegment& b, const LineSegment& c, const LineSegment& d) : _corners_cached(false), _cost_cached(false)
{
  _segments.push_back(a);
  _segments.push_back(b);
  _segments.push_back(c);
  _segments.push_back(d);
}

bool DoorObject::operator< (DoorObject& D)
{
  // compare the corners
  vector<Coord> c1 = getCorners();
  vector<Coord> c2 = D.getCorners();
  
  // must have same number of corners
  if(c1.size() != c2.size()) return false;
  
  // compare every corner in the sets
  vector<Coord>::iterator i1 = c1.begin(), i2 = c2.begin(); 
  while(i1 != c1.end() && i2 != c2.end())
    {
      if(*i1 != *i2) return false;
      i1++; i2++;
    }
  
  return true;
}


// returns the corners of the door
vector<Coord> DoorObject::getCorners(bool recompute)
{
  if(!recompute && _corners_cached) return _corners;
  
  // assert: we need to recompute before returning

  vector<Coord> ret;

  // new procedure -- create a graph of all the intersections and the paths between them... will form graph of quads and triangles  (using corner_node_t object)
  pair<corner_node_t*,size_t> tree = corner_node_t::buildTree(_segments);  
  

  vector<corner_node_t*> corners_ptrs;

  // choose points as corners if they have a path back to themselves of length 4 (points on triangles only will not have this path)
  for(size_t i = 0; i < tree.second; i++)
    {
      if(tree.first[i].hasCycleOfFour()) {
	tree.first[i].corner = true;
	corners_ptrs.push_back(&tree.first[i]);
	ret.push_back(tree.first[i].point);
      }
    }
  
  // figure out a circular ordering of four corners
  sort(ret.begin(), ret.end());

  if(ret.size() < 4) // guard for four points
    {
      if(tree.first) delete [] tree.first;
      return ret;
    }
  // make all possible polgons.  Return one with most area
  CPoly buffer;
  
  for(size_t i = 0; i < ret.size(); i++)
    for(size_t j = i; j < ret.size(); j++)
      for(size_t k = j; k < ret.size(); k++)
	for(size_t l = k; l < ret.size(); l++) {
	  CPoly C;
	  C.bdy.push_back(ret[i]);
	  C.bdy.push_back(ret[j]);
	  C.bdy.push_back(ret[k]);
	  C.bdy.push_back(ret[l]);
	  if(getArea(C) > getArea(buffer)) 
	    buffer = C;
	}
  ret = buffer.bdy;
  

  // release the tree of corner_node_t
  if(tree.first) delete [] tree.first;
  
  // cache and return
  _corners = ret;
  _corners_cached = true;
  return ret; 
}


CPoly DoorObject::asPoly()
{
  vector<Coord> corners = getCorners();
  CPoly c;
  c.bdy = corners;
  c.bdyColor = GREEN;  c.fillColor = noColor;
  return c;
}



// methods used to calculate costs for filtering of candidates

// returns distance of furthest point from inside of image
double DoorObject::outsidePoint_cost()
{
  if(!_corners_cached) getCorners(); // make sure we've computed out corners

  // TODO : write this method - may need to find a way to store the dimensions of the source image
  
  return 0.0;
}

  // returns "astigmatism" of door - mean difference between outside angles and such.
double DoorObject::cornerAngles_cost()
{
  double angle[4]; //angle for each corner 0-3 in order 
  double distance01; // distance from corner 0 to corner 1
  double distance12;
  double distance23;
  double distance30;
  double distance02; // distance on diagonal
  double distance13; //distance on diagonal
  vector<Coord> corners = this->getCorners();

  if(corners.size() < 4) return 0;

  //compute distances using pythagorean sums
  distance01 = sqrt( pow((corners[0].x-corners[1].x),2) + pow((corners[0].y-corners[1].y),2));
  distance12 = sqrt( pow((corners[1].x-corners[2].x),2) + pow((corners[1].y-corners[2].y),2));
  distance23 = sqrt( pow((corners[2].x-corners[3].x),2) + pow((corners[2].y-corners[3].y),2));
  distance30 = sqrt( pow((corners[3].x-corners[0].x),2) + pow((corners[3].y-corners[0].y),2));
  distance02 = sqrt( pow((corners[0].x-corners[2].x),2) + pow((corners[0].y-corners[2].y),2));
  distance13 = sqrt( pow((corners[1].x-corners[3].x),2) + pow((corners[1].y-corners[3].y),2)); 

  //compute angles using law of cosines c^2 = a^2 + b^2 -2*a*b*cos(angle btwn a&b)
  angle[0] = acos( (pow(distance30,2) + pow(distance01,2) - pow(distance13,2))/(2*distance01*distance30) );
  angle[1] = acos( (pow(distance01,2) + pow(distance12,2) - pow(distance02,2))/(2*distance01*distance12) );
  angle[2] = acos( (pow(distance12,2) + pow(distance23,2) - pow(distance13,2))/(2*distance12*distance23) );
  angle[3] = acos( (pow(distance23,2) + pow(distance30,2) - pow(distance02,2))/(2*distance23*distance30) );

  //normalized difference between angles
  double ratio02 = abs( ( abs(angle[0])-abs(angle[2]) ) / 360 );  
  double ratio13 = abs( ( abs(angle[1])-abs(angle[3]) ) / 360 );

  //average ratios together and return
  double ret = 1 - ( (ratio02 + ratio13) / 2 );  //not the best way to do this

  //if close to 1, angles similar - if close to 0, angles diverge (1=door,0=not door)
  return ret;
}


// returns comparison made in length in sides of door, check for stumpy or skinny
double DoorObject::geometric_cost()
{
  const static double idealRatio = 2.2;

  // get corners of door
  vector<Coord> corners = this->getCorners();  // assumes that this returns in a proper cyclical order, which it doesn't always

  if(corners.size() < 4) return 0.0;  // cannot be a door with fewer than four corners

  // make segments for all door edges
  LineSegment line01 (corners[0], corners[1]);
  LineSegment line12 (corners[1], corners[2]);
  LineSegment line23 (corners[2], corners[3]);
  LineSegment line30 (corners[3], corners[0]);
  
  // find height and width of door... average of segment lengths
  double height = (line01.length() + line23.length()) / 2;
  double width  = (line12.length() + line30.length()) / 2;
  
  if(height < width) { // swap
    double temp = height;
    height = width;
    width = temp;
  }

  // compare ratio to ideal ratio
  double ratio = height/width;
 
  // to make sure is always 
  if(ratio < idealRatio) return ratio/idealRatio;
  else                   return idealRatio/ratio;
  
}




double DoorObject::getCost(bool recompute) {
  if(_cost_cached) return _cost;
  _cost = outsidePoint_cost()*1.0 // TODO: add other cost methods, and their weights
        + cornerAngles_cost()*1.0 
        +    geometric_cost()*1.0;

  _cost_cached = true;
  return _cost;
}


// global methods for working with files of DoorObjects
// .doors file format:  uses a companion .off file for all the information.  Writes segments from lines into the off file
// DOORS #doors
// LINK OFF #offfilename 
// # list of segments
// # list of door entries : DOOR #points #polys <points> <polys>

vector<DoorObject> readDoorFile(string filename)
{
  // read file header

  // load linked off file, if possible, building vector<CPoly>

  // Build vector<LineSegment> from endpoints given  

  // read all doors in file, skipping lines that are malformed or begin with #
  
  // dummy return!
  vector<DoorObject> d;
  return d;
}

void writeDoorFile(string filename, vector<DoorObject>& doors)
{
  
  
}


// for corner_node_t interior struct of DoorObject

pair<DoorObject::corner_node_t*,size_t> DoorObject::corner_node_t::buildTree(vector<LineSegment> segments)
{
  // this is harder than it seems... we need to know how things are connected.  Will think about it.

  // compute all intersections
  vector<Coord> intersections;
  //cerr << "building tree from " << segments.size() << " segments" << endl;
  for(size_t i = 0; i < segments.size(); i++)
    {
      //cerr << "segments[" << i << "] = " << segments[i] << endl;
      for(size_t j = i; j < segments.size(); j++)
	{
	  Coord c = findIntersection(segments[i], segments[j]);
	  //cerr << " intersection between segments " << i << " and " << j << " : " << c << endl;
	  if(c != noCoord) intersections.push_back(c);
	}
    }

  //cerr << "Found " << intersections.size() << " intersections:" << endl;
  
  // convert array to corner_node_t
  size_t N = intersections.size();
  corner_node_t* nodes = new corner_node_t [N];
  // copy in the data from intersections vector
  for(size_t i = 0; i < N; i++)
    {
      nodes[i].point = intersections[i];
      //cerr << "   " << i << ": " << nodes[i].point << endl;
    }

  
  // for each segment, figure out the connections between intersections on that segment
  for(size_t i = 0; i < segments.size(); i++) {
    //cerr << segments[i] << endl;
    // find all nodes on line, sort
    vector<corner_node_t*> nodesOnLine;
    for(size_t f = 0; f < N; f++) { 
      //cerr << "        checking " << nodes[f].point;
      if( segments[i].onLine(nodes[f].point) ) {
	nodesOnLine.push_back( &nodes[f] );
	//cerr << " yes";
      }
      //cerr << endl;
    }
    //cerr << "  " << nodesOnLine.size() << " nodes on line " << endl;
    sort(nodesOnLine.begin(), nodesOnLine.end(), DoorObject::corner_node_t::compareNodesByCoord); // sort them by coordinate position
    
    // assert: items in nodesOnLine are sorted by position
    // link all the items back-to front with connections
    for(size_t n = 0; n < nodesOnLine.size()-1; n++)
      {
	nodesOnLine[n]->elt.push_back(nodesOnLine[n+1]);
	nodesOnLine[n+1]->elt.push_back(nodesOnLine[n]);
      }
  }
  
  pair<corner_node_t*,size_t> R (nodes, N);
  return R;
}



bool DoorObject::corner_node_t::hasCycleOfFour()  // TODO: debug
{
  // check all pairs of nodes in elt vector to see if they have a node in common
  for(size_t i = 0; i < elt.size(); i++) {
    for(size_t j = i; j < elt.size(); j++) {
      if(i == j) continue; // don't compare same nodes
      // check if they have something in common
      for(size_t x = 0; x < elt[i]->elt.size(); x++) {  // checking nodes between the two child nodes chosen here
	for(size_t y = 0; y < elt[j]->elt.size(); y++) {
	  if( (elt[i]->elt[x] == elt[j]->elt[y]) && (elt[i]->elt[x] != this) ) return true;
	}
      }
    }
  }
  return false;
}


vector<CPoly> asPolys(vector<DoorObject> l)
{ 
  vector<CPoly> ret;
  size_t bound = l.size();
  for (size_t i = 0; i < bound; i++)
    {
      ret.push_back( l[i].asPoly() );
    } 
  return ret;
}


vector<CPoly> asFilledPolys(vector<DoorObject> l)
{ 
  vector<CPoly> ret;
  vector<CPoly> temp;
  size_t bound = l.size();
  for (size_t i = 0; i < bound; i++)
    {
      temp = l[i].polys();
      ret.insert(ret.end(),temp.begin(),temp.end());
    } 
  return ret;
}

// DoorDetectorEriolMods.cpp
// Extensions to the eriol program by Olaf Hall-Holt for our door detection work
// April 2012
// @author Chris Cornelius, Ian McGinnis

#include "DoorDetectorEriolMods.hpp"

/* begin extension methods for eriol.hpp -- declared in DoorDetectorEriolMods.hpp*/

// Comparison functions for Coord object
bool const operator==(const Coord &o, const Coord &p) {
  if( &o == NULL )
    return false;
  return o.x == p.x && o.y == p.y;
}

bool const operator!=(const Coord &o, const Coord &p) {
  return o.x != p.x || o.y != p.y;
}

bool const operator<(const Coord &o, const Coord &p) {
  if(&o == NULL)
    return false;
  return o.x < p.x || (o.x == p.x && o.y < p.y);
  }

// multiplication function for Color (switched order)
Color operator*(const Color& c, double d)
{
  return d * c;
}

double getArea (const CPoly& c) { 
  double area = 0;         // Accumulates area in the loop
  int j = c.bdy.size()-1;  // The last vertex is the 'previous' one to the first

  for (unsigned int i=0; i<c.bdy.size(); i++)
    { area += (c.bdy[j].x+c.bdy[i].x) * (c.bdy[j].y-c.bdy[i].y); 
      j = i;  //j is previous vertex to i
    }
  area/=2;
  return fabs(area);
}


/* end eriolHeader.hpp extension methods */




/* begin methods for loading OFF files */

vector<CPoly> readOffFile(string filepath) {
	
  vector <Coord> vertList;
  vector <vector <int> > faceList;
  vector<CPoly> polyVector;
  
  ifstream f ( filepath.c_str() );
  
  if ( !f.good() ) {
    cout << "BAD FILE:" << filepath << endl;
    return polyVector;
  }
  
  if (f.is_open()) 
    {
      string tmp;
      string garbage;
      f >> tmp;
      int vert, face, edge;
      f>> vert;
      f>> face;
      f>> edge;
	  
      for (int i = 0; i < vert; i++) {
	double x = 0, y = 0, z = 0;
	f>>x;
	f>>y;
	f>>z;
	Coord temp(x,y);
	vertList.push_back( temp );
      }
	 	
      int numCoords = 0;
      for (int k = 0; k < face; k++) {
	f >> garbage;
	numCoords = strToInt(garbage);
	   
	faceList.push_back( *(new vector<int>()) );

	for (int j = 0; j< numCoords; j++) {
	  f >> garbage;	
	  int corIndex = strToInt(garbage);
	  faceList[k].push_back( corIndex );
	}
		
      }
      f.close();
    }
	
  polyVector = makePoly(vertList, faceList);
  return polyVector;
  }


vector<CPoly> makePoly (vector<Coord> &vert, vector <vector <int> > &face) {
  vector <CPoly> polyList;
  for (int unsigned i = 0; i < face.size(); i++) {
    CPoly tmp = *(new CPoly);
    
    vector<Coord> newBdy;
    
    for(unsigned int x = 0; x < face[i].size(); x++)
      newBdy.push_back( vert[ face[i][x] ] );

    tmp.bdy = newBdy;
    tmp.bdyColor = GREEN;
    tmp.fillColor = noColor; // important to prevent overlay becoming colored so can't see background
    tmp.text = "" + i;

    polyList.push_back( tmp );
  }

  return polyList;
}


int strToInt(string & strVal){
  stringstream ss (stringstream::in | stringstream::out);
  int returnVal;
	
  ss << strVal;
  ss >> returnVal;
  return returnVal;
}

/* end custom OFF file methods */


/* begin methods for managing multiple CPoly overlays */
vector < vector<CPoly> > CPolyOverlay::polygonOverlays;
vector <string> CPolyOverlay::polygonOverlayNames;
int CPolyOverlay::currPolyIndex = 0;
bool CPolyOverlay::verbose = true;
void CPolyOverlay::set(int num, size_t windownum)
{
  // check for overrun
  if(num >= (int)polygonOverlays.size())
    {
      num = (int)polygonOverlays.size()-1;
    }
  if(num < 0) num = 0;
  if(windownum >= CPolyOverlay::NUM_WINDOWS) return;

  // assert : num is a valid index and windownum is a valid window number
  
  // move into the allCPoly object
  CPolyOverlay::currPolyIndex = num;
  CPoly::allCPoly[windownum].clear(); 
  CPoly::allCPoly[windownum] = polygonOverlays[CPolyOverlay::currPolyIndex];
  if(CPolyOverlay::verbose) cerr << "%%% Viewing polys " << currPolyIndex << ": " << polygonOverlayNames[currPolyIndex] << " %%%" << endl;
}

void CPolyOverlay::add(const vector<CPoly>& src, const string& message)
{
  polygonOverlays.push_back(src);
  polygonOverlayNames.push_back(message);
}




/* end CPoly overlay methods */

// DoorDetector.cpp
// Class which is for detecting doors -- unification of Ian and Chris's work on door detection
// May 2012
// @author Ian McGinnis, Chris Cornelius

#include "DoorDetector.hpp"

/* Note that there's a LOT of whitespace in this file, and for good reason -- we want to be able to read
 * between the lines in all the functions we're writing. */


int DoorDetector::load(Image& img, vector<CPoly>& offPolys)  // sets up the initial information, clears cache
{
  rawImg = img;
  offPolygons = offPolys;
  
  gradientImg_cached = binaryImg_cached = false;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// for calculation of binary images
double DoorDetector::sobelX(Image& i, PixelLoc p)
{  
  // for invalid data
  if(p.x > (int)i.getWidth() || p.x < 0)
    return 0.0;
  if(p.y > (int)i.getHeight() || p.y < 0)
    return 0.0;

  // assert: we have the pixel 
  Color centerVal = i.getPixel(p);
  
  // for edge case
  if(p.x < 2 || p.x > (int)i.getWidth()-2) // we need to be at least one pixel away from edges
    return 0.0;  // TODO: return real data!
  if(p.y < 2 || p.y > (int)i.getHeight()-2) // we need to be at least one pixel away from edges
    return 0.0;  // TODO: return real data!
  
  // for general case... apply the map
  
  // computing Gx - gradient in X direction
  
  double Gx_mask[] = {-1.0, 0.0, 1.0,
		      -2.0, 0.0, 2.0,
		      -1.0, 0.0, 1.0};
  
  double answer = 0.0;
  int index = 0;
  // create answer
  for(int xi = -1; xi <=1; xi++)
    for(int yi = -1; yi <= 1; yi++)
      {
	PixelLoc loc(p.x+xi, p.y+yi);
	Color col = i.getPixel(loc);
	answer += (Gx_mask[index])*(col.r + col.g + col.b);
	index++; // what's the next pixel in the mask?
      }
  
  return answer;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
double DoorDetector::sobelY(Image& i, PixelLoc p)
{
  // for invalid data
  if(p.x > (int)i.getWidth() || p.x < 0)
    return 0.0;
  if(p.y > (int)i.getHeight() || p.y < 0)
    return 0.0;

  // assert: we have the pixel 
  Color centerVal = i.getPixel(p);
  
  // for edge case
  if(p.x < 2 || p.x > (int)i.getWidth()-2) // we need to be at least one pixel away from edges
    return 0.0;  // TODO: return real data!
  if(p.y < 2 || p.y > (int)i.getHeight()-2) // we need to be at least one pixel away from edges
    return 0.0;  // TODO: return real data!
  
  // for general case... apply the map
  
  // computing Gy - gradient in Y direction
  
  double Gy_mask[] = {-1.0, -2.0, -1.0,
		      0.0, 0.0, 0.0,
		      1.0, 2.0, 1.0};
  
  double answer = 0.0;
  int index = 0;
  // create answer
  for(int xi = -1; xi <=1; xi++)
    for(int yi = -1; yi <= 1; yi++)
      {
	PixelLoc loc(p.x+xi, p.y+yi);
	Color col = i.getPixel(loc);
	answer += (Gy_mask[index])*(col.r + col.g + col.b);
	index++; // what's the next pixel in the mask?
      }
  
  return answer;  
}
















///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
Image* DoorDetector::calcSobel(Image* img)
{
  // VARS FOR TUNING OF CALCULATION //
  double astig = 0.0; //rotation taken out for debugging
    //(-1/4)*3.14159;  // how far away from vertical is the image?
  
  double Gscale = 1000.0; // how much influence does gradient magnitude have in brightness output (more = lower)
  double Gthreshold = 0; // below what gradient magnitude do we cut off the points?

  cerr << "Calculating sobel for all points in image... ";

  // do the calculation for every pixel in the image - make sure to make a buffer image, though!
  
  Image* img2 = new Image(img->getWidth(), img->getHeight(), Image::CHAR_CHAN, 3);
  
#pragma omp parallel for
  for(unsigned x = 0; x < img->getWidth(); x++)
    {
      unsigned bound2 = img->getHeight();
      for(unsigned y = 0; y < bound2; y++)
	{
	  PixelLoc p (x, y);
	  Color c;
	  
	  // get the partial derivatives
	  double Gx = sobelX(*img, p);
	  double Gy = sobelY(*img, p);

	  
	  // calculate magnitude and direction for the gradient
	  double G = sqrt( pow(Gx, 2) + pow(Gy, 2)); // pythagorean
	  double T = astig + atan2(Gy,Gx); // rotate the coordinate set by the astigmatism value
	  //atan2 returns from +pi to -pi

	  // choose which bin to put this point into... which cardinal line of direction is it closest to?
	  double r, g, b;
	  int bin = -1;
	  double pi = 3.14159265;
	  if( ( (T <= -3.0*pi/8.0) && (T >= -5.0*pi/8.0) ) || ( (T >= 3.0*pi/8.0) && (T <= 5.0*pi/8.0) ) ) {  
// gradient is close to vertical - note weighted vertically because of less/greater than or equal to
	    r = 250; g = 0; b = 0;
	    bin = 0;
	  }
	  else if( ( (T <= pi/8.0) && (T >= -1.0*pi/8.0) ) || ( (T <= -7.0*pi/8.0) && (T >= 7.0*pi/8.0) ) ) { 
// gradient is close to horizontal- note weighted horizontally because of less/greater than or equal to
	    r = 0; g = 0; b = 250;
	    bin = 1;
	  }
	  else if( ( (T > pi/8.0) && (T < 3.0*pi/8.0) ) || ( (T > -7.0*pi/8.0) && (T < -5.0*pi/8.0) ) ) {  
// gradient is close to 45 degree line
	    r = 0; g = 250; b = 0;
	    bin = 2;
	  }
	  else if ( ( (T > 5.0*pi/8.0) && (T < 7.0*pi/8.0) ) || ( (T < -1.0*pi/8.0) && (T > -3.0*pi/8.0) ) ) {  
// gradient is close to -45 degree line
	    r = 125; g = 0; b = 125;
	    bin = 3;
	  }
	  

	  // set the color of the output pixel
	  
	  // CLAMPING BRINGS A PROBLEM
	  c.r = (G/Gscale)*r;
	  c.g = (G/Gscale)*g;
	  c.b = (G/Gscale)*b;
	  

	  // if we're below the threshold, make pixel black
	  if(G < Gthreshold) { c.r = c.g = c.b = 0; }
	  
	  
	  img2->setPixel(p, c);
	  //cerr << "(" << p << ") " << Gx << " " << Gy << " " << G << " " << T << endl;
	}
    }
  

  cerr << "Done" << endl;

  return img2;
}














///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
Image* DoorDetector::makeBinaryEdgeMap(Image* img)
{
  // IMPORTANT VARS //
  char upperThreshold = 22;  // minimum value for a pixel
  char lowerThreshold = 20; // minimum value for a potential pixel
  Color truePixelColor (0,0,255); // what color does an edge pixel become?
  Color maybePixelColor (255,0,0); // what color does a potential edge pixel become?

  // for neighborhoods
  int n = 3; // neighborhood is an n x n grid -- MUST BE ODD NUMBER!
  int nbound = (n-1)/2;
  int nsize = n*n;
  int ncenter = (nsize-1)/2;
  Color neighborhood[nsize];
  //Color hysteresis_neighborhood[nsize];
  
  // multiplication matrices -- NOTE: only work for 3x3
  double mult00 [] = {0,1,0,
		      0,0,0,
		      0,1,0}; 
  double mult90 [] = {0,0,0,
		      1,0,1,
		      0,0,0};
  double mult45 [] = {0,0,1,
		      0,0,0,
		      1,0,0};
  double mult135[] = {1,0,0,
		      0,0,0,
		      0,0,1};
  
  
  // allocate new image
  Image* img2 = new Image(img->getWidth(), img->getHeight(), Image::CHAR_CHAN, 3);
  
  // do "non-maximum supression"  
  // for each pixel:  compare to local environment based on gradient.  If gradient (greyscale value) is local maximum, then select this pixel as the edge
#pragma omp parallel for private(neighborhood)
  for(unsigned x = 0; x < img->getWidth(); x++)
    {
      unsigned bound2 = img->getHeight();
      for(unsigned y = 0; y < bound2; y++)
	{
	  PixelLoc p (x, y);	  
	  // handle edge cases
	  if(p.x < (n/2)+1 || p.x > (int)img->getWidth() - ((n/2)+1) ) // we need to be at least one pixel away from edges
	    continue;
	  if(p.y < (n/2)+1 || p.y > (int)bound2 - ((n/2)+1) ) // we need to be at least one pixel away from edges
	    continue;
	  
	  // assert: p is a Pixel not one pixel from the edge
	  
	  // if this has a value lower than the lower threshold, discard immediatly
	  
	  if( greyValue(img->getPixel(p)) < upperThreshold ) 
	    continue;
	    
	  else if( greyValue(img->getPixel(p)) >= upperThreshold ){
	      
	  
	      // get the neighborhood of the given pixel -- a 3x3 grid around it
	      int neighborhood_index = 0;
	      for(int xi = -1*nbound; xi <= nbound; xi++)
		{
		  for(int yi = -1*nbound; yi <= nbound; yi++)
		    {
		      PixelLoc loc(p.x+xi, p.y+yi);
		      neighborhood[neighborhood_index] = img->getPixel(loc);
		      neighborhood_index++;
		    }
		}
	      
	      // re-weight the world of neighbors based on the gradient angle  -- TODO: make this work better
	      Color c = neighborhood[ncenter];  // value at center pixel
	      Color blackColor;
	      
	      if((c.r > c.g) && (c.r > c.b)) // gradient is vertical (red)
		{ 
		  for(int i = 0; i < nsize; i++)
		    {
		      neighborhood[i] = neighborhood[i] * mult90[i];
		    }
		}
	      else if ((c.g > c.r) && (c.g > c.b)) // gradient is close to 45 degrees (green)
		{
		  for(int i = 0; i < nsize; i++)
		    {
		      neighborhood[i] = neighborhood[i] * mult135[i]; 
		    }
		}
	      else if ((c.b > c.r) && (c.b > c.g)) // gradient is close to horizontal (blue)
		{
		  for(int i = 0; i < nsize; i++)
		    {
		  neighborhood[i] = neighborhood[i] * mult00[i];
		    }
		}
	      else if ( c.r == c.b )// gradient is close to 135 degree line (purple)
		{
		  for(int i = 0; i < nsize; i++)
		    {
		      neighborhood[i] = neighborhood[i] * mult45[i]; 
		    }
		  
	    }
	      
	      
	      // Check the center pixel against all pixels in the re-weighted neighborhood to see if it's the max
	      bool gtflag = false;
	      for(int i = 0; i < nsize; i++)
		{
		  if ( greyValue(neighborhood[i]) > greyValue(c) ) // if we find a member which is greater than the center value, then this pixel isn't max
		    {
		  gtflag = true;
		  break;
		    }
		}
	      
	      // if we didn't find a more maximum pixel, set this one to the right color - blue
	      if(!gtflag) img2->setPixel(p, truePixelColor);
	    }
	}
    }
  //call hysteresis method
  Hysteresis(img, img2, upperThreshold, lowerThreshold, truePixelColor);
  return img2;
}














///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
void DoorDetector::Hysteresis(Image* img, Image* img2, int upperThreshold, int lowerThreshold,
			      Color truePixelColor, Color hysteresisPixelColor)
{
  // for neighborhoods
  int n = 3; // neighborhood is an n x n grid -- MUST BE ODD NUMBER!
  int nbound = (n-1)/2;
  int nsize = n*n;
  int ncenter = (nsize-1)/2;
  Color neighborhood[nsize];
  
  // multiplication matrices -- NOTE: only work for 3x3
  double mult00 [] = {0,1,0,
		      0,0,0,
		      0,1,0}; 
  double mult90 [] = {0,0,0,
		      1,0,1,
		      0,0,0};
  double mult45 [] = {0,0,1,
		      0,0,0,
		      1,0,0};
  double mult135[] = {1,0,0,
		      0,0,0,
		      0,0,1};


  for(unsigned x = 0; x < img->getWidth(); x++)
    {
      unsigned bound2 = img->getHeight();
      for(unsigned y = 0; y < bound2; y++)
	{
	  PixelLoc p (x, y);	  
	  // handle edge cases
	  if(p.x < (n/2)+1 || p.x > (int)img->getWidth() - ((n/2)+1) ) // we need to be at least one pixel away from edgesx
	    continue;
	  if(p.y < (n/2)+1 || p.y > (int)bound2 - ((n/2)+1) ) // we need to be at least one pixel away from edges
	    continue;

	  // if edge pixel
	  if (greyValue(img2->getPixel(p)) > 80)
	    {
	      // get the neighborhood of the given pixel -- a 3x3 grid around it
	      int neighborhood_index = 0;
	      for(int xi = -1*nbound; xi <= nbound; xi++)
		{
		  for(int yi = -1*nbound; yi <= nbound; yi++)
		    {
		      PixelLoc loc(p.x+xi, p.y+yi);
		      neighborhood[neighborhood_index] = img->getPixel(loc);
		      neighborhood_index++;
		    }
		}
	      
	      // re-weight the world of neighbors based on the gradient angle  -- check in same direction as gradient
	      Color c = neighborhood[ncenter];  // value at center pixel
	      Color blackColor;
	      
	      if((c.r > c.g) && (c.r > c.b)) // gradient is vertical (red)
		{ 
		  for(int i = 0; i < nsize; i++)
		    {
		      neighborhood[i] = neighborhood[i] * mult90[i];
		    }
		}
	      else if ((c.g > c.r) && (c.g > c.b)) // gradient is close to 45 degrees (green)
		{
		  for(int i = 0; i < nsize; i++)
		    {
		      neighborhood[i] = neighborhood[i] * mult45[i];
		    }
		}
	      else if ((c.b > c.r) && (c.b > c.g)) // gradient is close to horizontal (blue)
		{
		  for(int i = 0; i < nsize; i++)
		    {
		      neighborhood[i] = neighborhood[i] * mult00[i];
		    }
		}
	      else if ( c.r == c.b )// gradient is close to 135 degree line (purple)
		{
		  for(int i = 0; i < nsize; i++)
		    {
		      neighborhood[i] = neighborhood[i] * mult135[i];
		    }
		}
	     
	      // Create an array of the locations of the neighborhood pixels 
	      int locArray_index = 0;
	      PixelLoc locArray[nsize];
	      for(int xi = -1*nbound; xi <= nbound; xi++)
		{
		  for(int yi = -1*nbound; yi <= nbound; yi++)
		    {
		      PixelLoc loc(p.x+xi, p.y+yi);
		      locArray[locArray_index] = loc;
		      locArray_index++;
		    }
		}
	      
	      for(int i = 0; i < nsize; i++) //if filtered pixels meet lower threshold add to binary edge map
		{

		  if ( greyValue(neighborhood[i]) > lowerThreshold )
		    {
		      if (greyValue(img2->getPixel(locArray[i])) != 0)
			continue;
		       		      
		      else if (greyValue(img2->getPixel(locArray[i])) == 0)
			img2->setPixel(locArray[i], hysteresisPixelColor); // NOTE - change to hysteresisPixelColor if you want to see the pixels added by the hysteresis method (in red)
		    }
		}

	    } //end if edge pixel


	}//end y loop
    }//end x loop
}

















///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// for working with edges of polygons
bool DoorDetector::isStrongEdge(double thresh, Coord A, Coord B, Image* edgeImage)
{
  //check to see if A is B
  if (A == B){
    //cerr << "A == B tripped" << endl;
    return false;
  }

  //1. Get Length
  //change in x and y (negative sign for directionality)
  double dX = -1*(A.x-B.x);
  double dY = -1*(A.y-B.y);

  //pythagorean
  double length = sqrt(pow(dX,2)+pow(dY,2));

  //2. Number of Steps in Pixels
  double steps = length*2;
  //round steps up to nearest whole number
  steps = ceil (steps);
  double stepLenX = dX/steps;
  double stepLenY = dY/steps;

  //3. Generate all Pixel Test Points
  set<PixelLoc> linePix;
  set<PixelLoc>::iterator lineIt;
  Coord temp(0.0,0.0);
  temp = A;

  int i=0;
  //fill linePix
  while (temp < B){
    linePix.insert(asPixelLoc(temp));
    temp.x = temp.x+stepLenX;
    temp.y = temp.y+stepLenY;
    i++;
  }

  //4. check all linePix for matches (*lineIt gets linePix at lineIt)
  //greyValue average the channels so 1 max channel returns 85
  int matched=0;
  lineIt = linePix.begin();
  while( lineIt != linePix.end() ){
    if(greyValue(edgeImage->getPixel(*lineIt)) > 80){ 
      matched++;
    }
    lineIt++;
  }

  //5. Compare matched/total pixels to thresh.
  if ((double)linePix.size() == 0) {return false;} //catches before we try to divide by 0
  else if ( ( ((double) matched) / ((double) linePix.size()) ) >= thresh )
    { return true; }
  else { return false; }

}











vector<CPoly> DoorDetector::whichPolysHaveStrongEdges(double thresh, const vector<CPoly> & polygons, Image* edgeImage)  
{
  vector<CPoly> returnVector;
  
  // for each CPoly
  // for all edges (between any two verticies)
  // isStrongEdge?
  // if true, push onto return vector

  size_t numThreads = 0;
  vector<CPoly>* threadVectors;
  
#pragma omp parallel
  {
    // get number of threads and this thread's number
    numThreads = (size_t) omp_get_num_threads();
    size_t TH_num = (size_t) omp_get_thread_num();

    // if master thread, set up and allocate the destination buffers
    if(TH_num == 0)
      {
	threadVectors = new vector<CPoly> [numThreads];
	assert(threadVectors != NULL); // check for allocation error	
      }

    // compute bounds for this thread
    size_t begin = TH_num * (polygons.size() / numThreads);
    size_t end   = (TH_num + 1) * (polygons.size() / numThreads);

    for(size_t i = begin; i < end; i++)
      {
	// New method 3/9/12 -- use the polyToSegments function
	vector<LineSegment> edges = polyToSegments(polygons[i]);
	
	//cerr << "Poly " << i << " points: " << polygons[i].bdy.size() << " edges: " << edges.size() << endl;
	
	// loop through all those edges to see if one is strong
	bool guard = true;
	for(size_t j = 0, bound = edges.size(); j < bound && guard; j++) {
	  if(isStrongEdge(thresh, edges[j].A, edges[j].B, edgeImage)) {
	    guard = true;
	    threadVectors[TH_num].push_back(polygons[i]);  // dereference to the proper thread vector
	    //cerr << " - found strong edge" << endl;
	    break;
	  }
	}
      }
  } // end of OMP parallel block
  
  // condense all vectors into the master vector
  for(size_t i = 0; i < numThreads; i++)
    {
      returnVector.insert(returnVector.end(), threadVectors[i].begin(), threadVectors[i].end());
    }
  
  if(threadVectors) delete [] threadVectors;

  return returnVector;
}





///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
vector<LineSegment> DoorDetector::whichSegmentsHaveStrongEdges(double thresh, const vector<CPoly> & polygons, Image* edgeImage)
{
 vector<LineSegment> returnVector;
  
  // for each CPoly
  // for all edges (between any two verticies)
  // isStrongEdge?
  // if true, push onto return vector

  size_t numThreads = 0;
  vector<LineSegment>* threadVectors;
  
#pragma omp parallel
  {
    // get number of threads and this thread's number
    numThreads = (size_t) omp_get_num_threads();
    size_t TH_num = (size_t) omp_get_thread_num();

    // if master thread, set up and allocate the destination buffers
    if(TH_num == 0)
      {
	threadVectors = new vector<LineSegment> [numThreads];
	assert(threadVectors != NULL); // check for allocation error	
      }

    // compute bounds for this thread
    size_t begin = TH_num * (polygons.size() / numThreads);
    size_t end   = (TH_num + 1) * (polygons.size() / numThreads);

    for(size_t i = begin; i < end; i++)
      {
	// New method 3/9/12 -- use the polyToSegments function
	vector<LineSegment> edges = polyToSegments(polygons[i]);

	// loop through all those edges to see if any are strong
	for(size_t j = 0, bound = edges.size(); j < bound; j++) {
	  if(isStrongEdge(thresh, edges[j].A, edges[j].B, edgeImage)) {
	    threadVectors[TH_num].push_back(edges[j]);  // dereference to the proper thread vector
	    //cerr << " - found strong edge" << endl;
	  }
	}
      }
  } // end of OMP parallel block
  
  // condense all vectors into the master vector
  for(size_t i = 0; i < numThreads; i++)
    {
      returnVector.insert(returnVector.end(), threadVectors[i].begin(), threadVectors[i].end());
    }
  
  if(threadVectors) delete [] threadVectors;

  return returnVector;
}















///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
vector<LineSegment> DoorDetector::findStraightSegments(vector<LineSegment> strongEdges, vector<CPoly> strongEdgePolys, Image* binaryMap, bool returnOnlyMergedSegments)
{
  vector<LineSegment> ret; // for return values
  vector< vector<LineSegment>::iterator > bins; // hold the starts of each bin for slope
  vector< vector<LineSegment>::iterator > groups; // hold the starts of each bin for y intercept
  const static double slopeTolerance = ((2.0/180.0)*3.14159);
  const static double yIntTolerance = 30.0;
  
  // 1. Convert to set and then back to sortable vector
  set<LineSegment> e (strongEdges.begin(), strongEdges.end());
  vector<LineSegment> edges (e.begin(), e.end());

  // 2. Sort all segments by slope, using STL sort
  stable_sort(edges.begin(), edges.end(), LineSegment::compareBySlope);

  // 3. Separate into groups by slope... a certain tolerance is given.
  vector<LineSegment>::iterator it = edges.begin();
  bins.push_back(it);
  double slope_cache = (*it).slope();
  while(it != edges.end())
    {
      // if we're more than slopeTolerance degrees from the beginning of this bin, it's time for a new bin
      if( (*it).slope() - slope_cache > slopeTolerance )
	{
	  bins.push_back(it);
	  slope_cache = (*it).slope();
	}
      it++;

      // Extra feature: can include y-int information in our return value
      if(!returnOnlyMergedSegments) {
	Coord icept (0, (*it).yIntercept());
	LineSegment r (icept, noCoord);
	ret.push_back(r);
      }
    }
  bins.push_back(edges.end()); // last "bin" is the end iterator
  //cerr << "   Created " << bins.size() << " bins by slope " << endl;
  
  
  // 4. For each bin
  vector<vector<LineSegment>::iterator>::iterator bin = bins.begin();
  while(bin != bins.end()-1) // last bin is end iterator for edges
  {
    // (*bin) is the LineSegment we're on, *bin+1 is the start of the next bin
    // 4.1 re-sort by y-intercept
    sort(*bin, *(bin+1), LineSegment::compareByIntercept);
    
    // TODO: normalize for average slope in group.  Chris's note 4/10/12:
    /*Steep segments (with y-intercepts that are far from the origin) are going to have more range of spacing in their resultant y-intercept points because the angle with which they head toward the y-axis is steeper. So maybe we'll need to normalize within each group for the average slope, scaling the allowable y-intercept distance by how far from 1 the slope is. Hmmmmmm.... */
    
    // 4.2 Separate into groups with very close y-intercepts
    vector<LineSegment>::iterator segment = *bin;
    groups.push_back(segment);
    double yInt_cache = (**bin).yIntercept(); // set up initial value
    while(segment != *(bin+1))
      {
	if( (*segment).yIntercept() - yInt_cache > yIntTolerance ) {  // we've found the end of a group
	  groups.push_back(segment);
	  yInt_cache = (*segment).yIntercept();
	}

      segment++;
      }
    // 4.4 Next bin
    bin++;
  }
  groups.push_back(edges.end()); // last "group" is end iterator
  //cerr << "   Created " << groups.size() << " groups by y-intercept" << endl;
  
  
  

  /* Extra feature : set returnOnlyMergedSegments to false to hilight all lines in a given group by a color */
  if(!returnOnlyMergedSegments) {
    //cerr << "   Processing colorized map of line groups" << endl;
    vector<vector<LineSegment>::iterator>::iterator group1 = groups.begin();
    Color colorAr[] = { WHITE, BLACK, RED, ORANGE, YELLOW, GREEN, CYAN, BLUE,
			PURPLE,DARKRED,DARKGREEN,DARKBLUE,DARKYELLOW,DARKPURPLE,GREY};
    while(group1 != groups.end()-1) // last group is end iterator for edges
      {
	int index = rand() % 15; // choose a random color
	Color randCol = colorAr[index];
	
	it = *group1;
	while(it != *(group1+1)) // for all segments in this group, compare coordinates against max & min
	  {
	    LineSegment r((*it).A, (*it).B);
	    r.color = randCol;
	    ret.push_back(r);
	    it++;
	  }
	group1++;
      }
  
    return ret;
  }
  /* End extra feature */

  // 5. Condense each group... for each group, merge all by bounding box
  
  vector<vector<LineSegment>::iterator>::iterator group = groups.begin();
  group = groups.begin();

  for(group = groups.begin(); group != groups.end()-1; group++) {  // last group is end iterator for edges
    Coord min(noCoord), max(0.0, 0.0);
    it = *group;
    while(it != *(group+1)) // for all segments in this group, compare coordinates against max & min
      {
	if((*it).A < min) min = (*it).A;
	if(max < (*it).B) max = (*it).B;
	it++;
      }
    // put those coords onto the return vector
    LineSegment r (min, max);
    r.color = PURPLE;
    ret.push_back(r);
    }
  //  cerr << "   Condensed all groups into " << ret.size() << " merged segments" << endl;
  
  return ret;

 }
 










///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// for forming shapes from the data and generating candidate doors
vector<DoorObject> DoorDetector::findAllQuadrilaterals(vector<LineSegment> segments)
{
  vector<DoorObject> returnVector;
  
  // pre-filter by segment length - segments shorter than some percentage of the image width (assumed in global value) are ignored
  sort(segments.begin(), segments.end(), LineSegment::compareByLength); 
  double minLength = rawImg.getWidth()*0.08;  // restrict by length of segments
  vector<LineSegment>::iterator it = segments.begin();
  while((*it).length() < minLength) it++;
  
  // build new data structure and sort it by y-intercept
  vector<LineSegment> seg (it, segments.end()); 
  sort(seg.begin(), seg.end(), LineSegment::compareByIntercept);
  
  // TODO: think about more filtering options?   


  // constrain to quadrilaterals with endpoints inside a 20% bounding box
  double widthBound = rawImg.getWidth()*1.2;
  double heightBound = rawImg.getHeight()*1.2;
  double widthExceed = rawImg.getWidth()*1.2 - rawImg.getWidth();
  double heightExceed = rawImg.getHeight()*1.2 - rawImg.getHeight();

  
  // find *every single* quadrilateral, meaning all intersections of segments.  At least can do in parallel
  size_t numThreads = 0;
  vector<DoorObject>* threadVectors;  
#pragma omp parallel 
  {
    // get number of threads and this thread's number
    numThreads = (size_t) omp_get_num_threads();
    size_t TH_num = (size_t) omp_get_thread_num();
    
    // if master thread, set up and allocate the destination buffers
#pragma omp critical
    {
      if(TH_num == 0) {
	threadVectors = new vector<DoorObject> [numThreads];
	assert(threadVectors != NULL); // check for allocation error	
	cerr << "Detecting all quadrilaterals using " << numThreads << " threads" << endl;
      }
    }

    // compute bounds for this thread
    size_t begin = TH_num * (seg.size() / numThreads);
    size_t bound = (TH_num + 1) * (seg.size() / numThreads);

    // this will be an N^4 operation... ick
    for(size_t i = begin; i < bound; i++)
      for(size_t j = i+1; j < bound; j++) // use i+1 to guarantee different i,j,k,l
	for(size_t k = j+1; k < bound; k++)
	  for(size_t l = k+1; l < bound; l++)
	  {
	    // create a DoorObject from the four segments we're given
	    DoorObject d (seg[i], seg[j], seg[k], seg[l]);

	    vector<Coord> c = d.getCorners();
	    vector<Coord>::iterator it = c.begin();
	    bool okay = true;
	    while(it != c.end())
	      {
		if((*it).x < -1*widthExceed || (*it).x > widthBound || (*it).y < -1*heightExceed || (*it).y > heightBound ) okay = false;
		it++;
	      }
	    
	    if(okay) threadVectors[TH_num].push_back(d);


	  }
  }  // end of OMP parallel block
  
  // condense all vectors into the master vector
  for(size_t i = 0; i < numThreads; i++)
    {
      returnVector.insert(returnVector.end(), threadVectors[i].begin(), threadVectors[i].end());
    }
  
  if(threadVectors) delete [] threadVectors;

  return returnVector;
}




// get the length of all the sides of a door.  If it's longer than the image width, call it 1.0
double DoorDetector::cumulativeLength(DoorObject& door, double imgWidth)
{
  vector<LineSegment> edges = polyToSegments(door.asPoly());
  double len = 0.0;
  for(size_t i = 0; i < edges.size(); i++)
    {
      len += edges[i].length();
    }
  
  len /= (imgWidth);
  if(len > 1.0) len = 1.0;
  return len;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// for filtering the candidates
vector<DoorObject> DoorDetector::doFiltering(vector<DoorObject> candidates)
{
  // Defined a bunch of cost functions to evaluate aspects of a DoorObject in DoorObject class.  Should be easy to sort all in the vector by their "cost"
  //sort(candidates.begin(), candidates.end(), DoorObject::compareByCost);
  
  vector<DoorObject> ret;   
  double largestAngleDifference = 5; //largest average angle allowed
  double limit = largestAngleDifference/360; //converts largest average angle difference to a percentage
  
  for (unsigned i=0; i < candidates.size(); i++)
    {
      if ( (candidates[i].cornerAngles_cost() > limit) &&
	   (candidates[i].geometric_cost() > 0.8) && 
	   (cumulativeLength(candidates[i], rawImg.getWidth()) > 0.9) )
	ret.push_back(candidates[i]);
    } 
  


  return ret;
}








///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
//vector<DoorObject> DoorDetector::captureInteriorPolys(vector<DoorObject> nakedCandidates, vector<CPoly> offPolys);  //TODO: copy





//helper functions

bool DoorDetector::vertLineCheckLeft(Coord A, Coord B, Coord Z){ 
  double slope = (B.y - A.y) / (B.x - A.x); //unstable if verticle
  double b = A.y - (slope * A.x);
  double x = (Z.y - b) / slope;
  if (x == Z.x)
    return true;
  else
    return (Z.x < x);
  // returns true if Z.x is left of the line, false if it's to the right

}

bool DoorDetector::vertLineCheckRight(Coord A, Coord B, Coord Z){ 
  double slope = (B.y - A.y) / (B.x - A.x); //unstable if verticle
  double b = A.y - (slope * A.x);
  double x = (Z.y - b) / slope;
  if (x == Z.x)
    return true;
  else
    return (Z.x > x);
  // returns true if Z.x is right of the line, false if it's to the left

}

bool DoorDetector::horizLineCheckBelow(Coord A, Coord B, Coord Z){ 
  double slope = (B.y - A.y) / (B.x - A.x);
  double b = A.y - (slope * A.x);
  double y = slope * Z.x + b;
  if (y == Z.y)
    return true;
  else
    return (Z.y < y);
  // returns true if Z.y is below the line, false if it's above
}

bool DoorDetector::horizLineCheckAbove(Coord A, Coord B, Coord Z){ 
  double slope = (B.y - A.y) / (B.x - A.x);
  double b = A.y - (slope * A.x);
  double y = slope * Z.x + b;
  if (y == Z.y)
    return true;
  else
    return (Z.y > y);
  // returns true if Z.y is above the line, false if it's below
}

bool DoorDetector::withinBounds(DoorObject Door, CPoly Candidate) //DoorObject or ordered corners?
{
  //doorCorners is size 4 [0-3]
  int size = Candidate.bdy.size()-1; //last coord in bdy() is the first coord in bdy()
  bool TFtable[size];
  vector<Coord> dC = Door.getCorners(); 
// 0 1 
//     <-resultant order of door corners from .getCorners() - changed need to retest
// 3 2
  for (int i=0; i < size; i++)
    { //inverted
      if(vertLineCheckLeft(dC[1],dC[2],Candidate.bdy[i]) 
	 && vertLineCheckRight(dC[0],dC[3],Candidate.bdy[i]) ) //x
	if(horizLineCheckBelow(dC[0],dC[1],Candidate.bdy[i]) 
	   && horizLineCheckAbove(dC[3],dC[2],Candidate.bdy[i]) ) //y
	  { TFtable[i] = true; }
	else TFtable[i] = false;
      else TFtable[i] = false;
    }
  //vectors can be accessed like above - to copy set to vector say: vector<type> X (Set.begin(),Set.end())
  for (int i=0; i<size; i++){
    if (TFtable[i] == false)
      return false;
  }

  return true;

} //checks to see if a Candidate is between the 4 lines of the Door


// top-level method
vector<DoorObject> DoorDetector::captureInteriorPolys(vector<DoorObject> nakedCandidates, vector<CPoly> offPolys)
{
  //naked is all candidates who don't have interior polys - offPolys are all polys from first off
  cerr << "Begin InteriorPolys... Capturing Interior Polys for " << nakedCandidates.size() << " Candidates... This may take some time. There are " << offPolys.size() << " Polys to check for each Candidate." << endl;

#pragma omp parallel for 
  for(unsigned j=0; j < nakedCandidates.size(); j++){
    for(unsigned i=0; i < offPolys.size(); i++){
      if (withinBounds(nakedCandidates[j],offPolys[i])) 
	{ nakedCandidates[j].addPoly(offPolys[i]); 
	  //cerr << "Candidate: " << j << " Poly: " << i << " Added." << endl;
	}
    }
  }

  vector<DoorObject> ret;
  ret = nakedCandidates;
  return ret;
}








/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// STATIC DEMO FUNCTION - NOT FOR USE WITH HEAVY COMPUTING //////////////////////////////////////////////////
void doorDetectorDemo() {
  //  ImageUI::allImages.push_back(img);
  //  CPolyOverlay::add(src, message);
  
  // get our initial data
  Image* img = ImageUI::allImages[0];
  vector<CPoly>* srcpolys = &CPolyOverlay::polygonOverlays[0];
  
  // set up the object
  DoorDetector D;
  D.load(*img, *srcpolys);
  
  Image* sobelImg = D.calcSobel(img); 
  cerr << "[--- Made sobel image ---]" << endl;

  Image* binaryMap = D.makeBinaryEdgeMap(sobelImg);
  cerr << "[--- Made binary edge map ---]" << endl;

  vector<CPoly> polys = D.whichPolysHaveStrongEdges(0.65, *srcpolys, binaryMap);
  vector<LineSegment> strongedges = D.whichSegmentsHaveStrongEdges(0.65, *srcpolys, binaryMap);  
  cerr << "[--- Found strong edges and segments ---]" << endl;
  cerr << "Found " << strongedges.size() << " strong edges (out of " << (*srcpolys).size() << " polygons)." << endl;

  vector<LineSegment> vismergededges = D.findStraightSegments(strongedges, polys, binaryMap, false);
  vector<LineSegment> mergededges = D.findStraightSegments(strongedges, polys, binaryMap);
  cerr << "[--- Merged qualifying edges ---]" << endl;
  cerr << "Created " << mergededges.size() << " merged edges." << endl;
 
  vector<DoorObject> candidates = D.findAllQuadrilaterals(mergededges);
  cerr << "[--- Extracted candidate doors ---]" << endl;
  cerr << "Created " << candidates.size() << " door candidates." << endl;

  vector<DoorObject> filtered = D.doFiltering(candidates);
  cerr << "[--- Filtered out " << candidates.size()-filtered.size() << " door objects ---]" << endl;
  cerr << "Result: " << filtered.size() << " DoorObjects" << endl;

  vector<DoorObject> captured = D.captureInteriorPolys(filtered, *srcpolys);
  cerr << "[--- Captured interior polygons in " << captured.size() << " polygons" << endl;

  
  // add overlays for the polygons we found and their endpoints, also calculating merges for the smaller set, for demo purposes
  CPolyOverlay::add(polys, "Strong Edge Polygons");  
  CPolyOverlay::add(asPolys(strongedges), "Strong Edges");
  CPolyOverlay::add(asPolys(vismergededges), "Visualized Merged Edges");
  CPolyOverlay::add(asPolys(mergededges), "Merged Edges");
  CPolyOverlay::add(asPolys(candidates), "Candidate Door Objects");
  CPolyOverlay::add(asPolys(filtered), "Candidates post-filtering");
  CPolyOverlay::add(asFilledPolys(captured), "Filled Candidate Door Objects");
  CPolyOverlay::set(CPolyOverlay::currPolyIndex-1);

  // add the newly computed images to the stack of images and change to be viewing it.
  ImageUI::allImages.push_back(sobelImg);
  ImageUI::allImages.push_back(binaryMap);
  ImageUI::setImageIndex(ImageUI::allImages.size()-1);

  ImageUI::tetherMotion = true;

  cerr << "[--- Configured new polygon overlays and image displays --]" << endl;

  
  
}

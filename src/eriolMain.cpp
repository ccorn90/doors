// by Olaf Hall-Holt, 2007-2011
#include<iostream>
#include<fstream>
#include<string.h>  // for strcmp
using namespace std;
#ifdef MACOSX
#include<GLUT/glut.h>
#else
#include<GL/glut.h>
#endif
#include"eriolHeader.h"

void assignImgBBox();
void resetZoom(int windowIndex);
void removeCurrImage(int windowIndex);
void init_gl_window();

int myWindowID[4] = {-1, -1, -1, -1};
bool NO_DISPLAY=false;
extern ImageUI::UIMode theUIMode;
extern bool realTimeInfo;
bool leftMouseButtonIsDown = false;
bool middleMouseButtonIsDown = false;
bool rightMouseButtonIsDown = false;
bool shiftPressed = false;
bool ctrlPressed = false;
PixelLoc currXY(0,0);
Coord startImagePt;
Coord startFocalLength;

void exampleFunction(int x, int y, int windowIndex)
{
  // where did the mouse click happen?
  Coord pt = ImageUI::asImageCoord(x, y, windowIndex);

  // draw it!
  addCPoly(pt, GREEN, windowIndex);

  // do some math
  vector<Coord> p, q;
  p.push_back( Coord(0,0) ); q.push_back( Coord(0,0) );
  p.push_back( Coord(1,0) ); q.push_back( Coord(0,1) );
  p.push_back( Coord(0,1) ); q.push_back( Coord(1,0) );
  p.push_back( Coord(1,1) ); q.push_back( Coord(1,1) );
  cout << "H: " << computeHomography(p,q) << endl;
}

void commonKeyboard( unsigned char c, int x, int y, int windowIndex)
{
// cerr << " c " << int(c) << endl;
  switch( c ) {
    case 'a':
      exampleFunction(x, y, windowIndex); break;
    case 'b':
      assignImgBBox(); break;
    case 'C':
      ImageUI::cropImage(); break;
    case 'w':
      resetZoom(windowIndex); break;
    case ' ':
      ImageUI::incrImageIndex(true, windowIndex); break;
    case 0x7F: // backspace
      if(glutGetModifiers() == GLUT_ACTIVE_SHIFT) removeCurrImage(windowIndex);
      else ImageUI::incrImageIndex(false, windowIndex);
      break;
    case '+':
      zoomInBy(sqrt(2.0), x, y, windowIndex); break;
    case '-':
      zoomInBy(sqrt(0.5), x, y, windowIndex); break;
    case '>':
      resizeWindow(sqrt(2.0), windowIndex); break;
    case '<':
      resizeWindow(sqrt(0.5), windowIndex); break;
    case 44: // ctrl-,  (similar to <)
      ImageUI::currImage(windowIndex)->subSample(); break;
    case 47: // ctrl-/
      ImageUI::currImage(windowIndex)->flipXY(); break;
    case 31: // ctrl--
      ImageUI::currImage(windowIndex)->flipY(); break;
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9':
      ImageUI::setImageIndexChar(c, windowIndex); break;
    case '@':
      CPoly::showCPoly = !CPoly::showCPoly; break;
    case '"':
      if ( ImageUI::quadView ) {
        ImageUI::quadView = false;
        quadrupleWindow();
      }
      ImageUI::dualView = !ImageUI::dualView; doubleWindow(); break;
    case '\'':
      if (glutGetModifiers() == GLUT_ACTIVE_CTRL) {
        if ( ImageUI::dualView ) {
          ImageUI::dualView = false;
          doubleWindow();
        }
        ImageUI::quadView = !ImageUI::quadView; quadrupleWindow();
      }
      break;
    case 3: case 27: case 'Q':
      if ( -1 != myWindowID[3] ) glutDestroyWindow(myWindowID[3]);
      if ( -1 != myWindowID[2] ) glutDestroyWindow(myWindowID[2]);
      if ( -1 != myWindowID[1] ) glutDestroyWindow(myWindowID[1]);
      glutDestroyWindow(myWindowID[0]);
      exit(0);
  }
}

// the mouse function is called when a mouse button is pressed down or released
void mouse(int button, int state, int x, int y, int windowIndex)
{
  const bool v = false;
  if ( v ) cerr << "mouse " << button << " " << state << endl;
  if ( GLUT_LEFT_BUTTON == button ) {
    if ( GLUT_DOWN == state ) {
      leftMouseButtonIsDown = true;
      shiftPressed = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);
      ctrlPressed = (glutGetModifiers() == GLUT_ACTIVE_CTRL);
      currXY = PixelLoc(x,y);
      startImagePt = ImageUI::asImageCoord(x,y, windowIndex);
      drawMouseInfo(x,y, shiftPressed, ctrlPressed, windowIndex);
    } else if ( GLUT_UP == state ) {
      leftMouseButtonIsDown = false;
      undrawMouseInfo();
    }
  } else if ( GLUT_MIDDLE_BUTTON == button ) {
    if ( GLUT_DOWN == state ) { // start to zoom in
      middleMouseButtonIsDown = true;
      zoomInBy(sqrt(2.0), x, y, windowIndex);
    } else middleMouseButtonIsDown = false;
  } else if ( GLUT_RIGHT_BUTTON == button ) {
    if ( GLUT_DOWN == state ) { // start to zoom out
      rightMouseButtonIsDown = true;
      zoomInBy(sqrt(0.5), x, y, windowIndex);
    } else rightMouseButtonIsDown = false;
  } else if ( 3 == button ) { // mouse wheel up (doesn't work on Mac)
    zoomInBy(sqrt(2.0), x, y, windowIndex);
  } else if ( 4 == button ) { // mouse wheel down (doesn't work on Mac)
    zoomInBy(sqrt(0.5), x, y, windowIndex);
  }
  glutPostRedisplay();
}

// simple usage message for how to use this program
void usage(char *progname)
{
  cerr << "Usage: " << progname << " [image.ppm]" << endl;
  exit(-1);
}

int main(int argc, char **argv)
{
  if ( argc > 1 && !strcmp(argv[1], "-no") ) {
    NO_DISPLAY = true;
    --argc; ++argv;
  } else if ( argc > 1 && !strcmp(argv[1], "-yes") ) {
    NO_DISPLAY = false;
    --argc; ++argv;
  }
  if ( argc < 2 ) {
    ifstream f(".defaultImageName");
    if ( !f.good() ) usage(argv[0]);
    string fname;
    while ( f >> fname ) {
      ifstream f2(fname.c_str());
      if ( !f2.good() ) {
        cerr << "Unable to open current image with name " << fname << "... exiting." << endl;
        return -1;
      }
      f2.close();
      ImageUI::addImage(fname.c_str());
    }
  } else {
    ofstream f(".defaultImageName");
    while ( argc > 1 ) {
      ImageUI::addImage(argv[1]);
      f << ImageUI::allImages.back()->name << endl;
      ++argv;
      --argc;
    }
  }
  if ( ImageUI::allImages.size() < 2 ) theUIMode = ImageUI::NORMAL;
  ImageUI::bbox.w = ImageUI::currImageWidth(false);
  ImageUI::bbox.h = ImageUI::currImageHeight(false);
  for ( unsigned i=0,i_end=ImageUI::allImages.size(); i<i_end; ++i ) {
    Image *img = ImageUI::allImages[i];
    estimateImageScale(img->getWidth(), img->getHeight(), img->imageScale);
  }
  ImageUI::windowPixelWidth[0] = ImageUI::imageScale(false) * ImageUI::bbox.xLim();
  ImageUI::windowPixelHeight[0] = ImageUI::imageScale(false) * ImageUI::bbox.yLim();
  ImageUI::sourceColorImage = ImageUI::currImage(false);
  if ( ! NO_DISPLAY ) init_gl_window();  // does not return
}

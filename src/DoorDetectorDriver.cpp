// DoorDetectorDriver.cpp
// Driver executable for our work on door detection
// April 2012
// @author Chris Cornelius, Ian McGinnis
// Based of work on Eriol by Prof. Olaf Hall-Holt

#include <iostream>
#include <fstream>
#include <string.h>
using namespace std;

#include "DoorDetectorEriolMods.hpp"
// Above header file includes "eriolHeader.hpp and also the needed functions for
// working with OFF files and displaying multiple polygon overlays

#include "DoorDetector.hpp"

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

void commonKeyboard( unsigned char c, int x, int y, int windowIndex)
{
  switch( c ) {
    // to select different polygon overlays and image displays -- TODO: check to make sure this works as I expect
  case '{':
    CPolyOverlay::set(CPolyOverlay::currPolyIndex-1);
    break;
  case '}':
    CPolyOverlay::set(CPolyOverlay::currPolyIndex+1);
    break;
  case '[':
    ImageUI::setImageIndex(ImageUI::currImageIndex[0]-1);
    break;
  case ']':
    ImageUI::setImageIndex(ImageUI::currImageIndex[0]+1);
    break;
  
  case '#': // run the door-detector, for demo purposes
    doorDetectorDemo();
    break;


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
      
      // ADDED 4/29/12 - load OFF file directly
      fname.append(".off");
      string fname_text("Initial polys from "); fname_text.append(fname);
      CPolyOverlay::add(readOffFile(fname), fname_text);
      CPolyOverlay::set(0,0); // display on window zero

    }
  } else {
    ofstream f(".defaultImageName");
    while ( argc > 1 ) {
      ImageUI::addImage(argv[1]);
      
      // ADDED 4/29/12 - load OFF file directly
      string fname(argv[1]);
      fname.append(".off");
      string fname_text("Initial polys from "); fname_text.append(fname);
      CPolyOverlay::add(readOffFile(fname), fname_text);
      CPolyOverlay::set(0,0); // display on window zero
      
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




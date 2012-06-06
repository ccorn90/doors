ifndef DEVEL
	DEVEL=false  #invoke from commandline to compile with development commands
endif

CC=g++
CFLAGS=-c

OPTS=-fopenmp -Iinclude -O3
DEVEL_OPTS=-fopenmp -Iinclude -g -Wall

# should we use development options or high optimization?
ifeq ($(DEVEL), true)
	C_OPTS=$(DEVEL_OPTS) $(CFLAGS)
	L_OPTS=$(DEVEL_OPTS)
else
	C_OPTS=$(OPTS) $(CFLAGS)
	L_OPTS=$(OPTS) -Wall
endif

LIBS = -lGL -lglut -lm

DOOR_DETECTOR_OBJS=build/DoorDetector.o build/DoorObject.o build/LineSegment.o build/DoorDetectorEriolMods.o

default: bin/DoorDetectorDriver

.PHONY : clean
clean:
	rm -f $(DOOR_DETECTOR_OBJS)
	rm -f bin/driver


run: bin/DoorDetectorDriver
	bin/DoorDetectorDriver data/1614L


### DoorDetector executable statements ###
bin/DoorDetectorDriver: src/DoorDetectorDriver.cpp build/eriolObjs.o $(DOOR_DETECTOR_OBJS)
	$(CC) $(L_OPTS) src/DoorDetectorDriver.cpp \
		$(DOOR_DETECTOR_OBJS) build/eriolObjs.o \
		$(LIBS) -o $@ 


bin/doorDetect: src/DoorDetectorDriver.cpp build/eriolObjs.o $(DOOR_DETECTOR_OBJS)
	$(CC) $(OPTS) src/doorDetect.cpp \
		$(DOOR_DETECTOR_OBJS) build/eriolObjs.o \
		$(LIBS) -o $@ 

### DoorDetector (release build) statements ###

build/DoorDetector.o: src/DoorDetector.cpp include/DoorDetector.hpp include/eriolHeader.h include/DoorDetectorEriolMods.hpp
	$(CC) $(C_OPTS) src/DoorDetector.cpp -o $@

build/DoorObject.o: src/DoorObject.cpp include/DoorObject.hpp
	$(CC) $(C_OPTS) src/DoorObject.cpp -o $@

build/LineSegment.o: src/LineSegment.cpp include/LineSegment.hpp
	$(CC) $(C_OPTS) src/LineSegment.cpp -o $@

build/DoorDetectorEriolMods.o: src/DoorDetectorEriolMods.cpp include/DoorDetectorEriolMods.hpp
	$(CC) $(C_OPTS) src/DoorDetectorEriolMods.cpp -o $@

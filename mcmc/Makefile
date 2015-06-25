LIBS=-L/usr/lib/x86_64-linux-gnu -lCore -lCint -lRIO -lNet -lHist -lGraf -lGpad -lTree -lRint -lMatrix -lMathCore


all:
	g++ -g mcmc.cpp -o main -std=c++11

read:
	g++ -g `root-config --cflags ` read.cpp -o read ${LIBS}  

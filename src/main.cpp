#include <cstdlib>
#include <string>
#include <fstream>
#include "main.hpp"
using namespace std;



/*
to compile:
g++ main.cpp envclus.cpp kmeans.cpp ktests.cpp calcpart.cpp partdist.cpp Lmunkres.cpp BipartiteGraph.cpp Hungarian.cpp ap.cpp statistics.cpp specialfunctions.cpp alglibinternal.cpp alglibmisc.cpp linalg.cpp -I/opt/local/include -L/opt/local/lib -lboost_system-mt -lboost_filesystem-mt -o ../bin/envclus

usage:  envclus [-e numsamples kmax gridfile envlayersdir outfile] [-c envpartfile genpartfile outfile] [-i envpartfile genpartfile outfile]

options:
-e  Perform environmental clustering
	
	numsamples = number of random points to sample, location weighted
	kmax = maximum number of clusters desired, program will consider K=2 to Kmax automatically
	gridfile = path to asc format GIS layer containing location weights, e.g. a niche model 
	envlayersdir = path to a directory containing environmental layers in ESRI ascii format
	outfile = path to output file, will contain the environmental partition, with cluster assignment for each of numsamples points

-c  Compare genetic and environmental partitions
	
	envpartfile = path to file containing environmental partition, formatted like output from -e option 
	genpartfile = path to file containing individual ID, long/lat, and genetic partitions
		for a series of K values, see example files for format
	outfile = path to output file, will contain a table of ranked partition distances between env and gen

-i  Test whether isolation by environment better explains distribution of genetic diversity
      than isolation by distance
    
    envpartfile = path to file containing environmental partition, formatted like output from -e option 
	genpartfile = path to file containing individual ID, long/lat, and genetic partitions
		for a series of K values, see example files for format
	outfile = path to output file, will contain a table of ranked partition distances between env and gen


example: ./envclus -e 10000 10 pumilustiny_avg.asc asciifiles envparts.txt -c envparts.txt genparts.txt pdout.txt


*/

/***************MAIN*****************/

int main( int argc, char* argv[] )
{
	//set up output file stream for log file
	ofstream outf; 
	outf.open ("./log.txt");
	outf.close(); //quick open close done to clear any existing file each time envclus is run
	outf.open ("./log.txt", ios::out | ios::app); //open file in append mode

	/*
	//get command line arguments 
	const int NumSamples = atoi(argv[1]);
	const int k = atoi(argv[2]); //k is kmax, the maximum number of clusters
	const char* GridFile = argv[3];
	char* EnvLayersDir = argv[4];
	*/
	

	//parse command line
	for (int i=0;i<argc;i++)
	{
		//run envclus
		if ( string(argv[i]) == "-e") 
    	{
        	//get arguments
        	const int NumSamples = atoi(argv[i+1]);
			const int k = atoi(argv[i+2]); //k is kmax, the maximum number of clusters
			const char* GridFile = argv[i+3];
			char* EnvLayersDir = argv[i+4];
			const char* OutFile = argv[i+5];
				
			//run
        	envclus(NumSamples, k, GridFile, EnvLayersDir, OutFile, outf);
 		} 
		
		//run calcpart
		if ( string(argv[i]) == "-c") 
    	{
        	//get arguments
        	const char* PathToEnvPartition = argv[i+1];
			const char* PathToGenPartition = argv[i+2];
			const char* PathToCalcPartOutput = argv[i+3];
	
        	//run
        	calcpart(PathToEnvPartition, PathToGenPartition, PathToCalcPartOutput, outf);	
 		} 
	
	}

	outf.close();
	return 0;
}
envclus calculates the partition distance between environmental and genetic clusters

to compile:
g++ main.cpp envclus.cpp kmeans.cpp ktests.cpp calcpart.cpp partdist.cpp Lmunkres.cpp BipartiteGraph.cpp Hungarian.cpp ap.cpp statistics.cpp specialfunctions.cpp alglibinternal.cpp alglibmisc.cpp linalg.cpp -lboost_system-mt -lboost_filesystem-mt -o ../bin/envclus;


usage:  envclus [-e numsamples kmax gridfile envlayersdir outfile] [-c envpartfile genpartfile outfile]

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


examples: ./envclus -e 10000 10 pumilustiny_avg.asc asciifiles envparts.txt; #calculate environmental partition
          ./envclus -c envparts.txt genparts.txt pdout.txt; #compare environmental partition you just calculated with genetic partition that you provide from some other analysis (e.g. STRUCTURE, STRUCTURAMA)
          ./envclus -e 10000 5 pumilustiny_avg.asc asciifiles envparts.txt -c envparts.txt genparts.txt pdout.txt; #calculate environmental partition then compare with genetic partition using a single command


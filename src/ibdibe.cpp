#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <numeric>
#include <cmath>
#include "calcpart.hpp"
#include "envclus.hpp"
#include "main.hpp"
#include "alglibinternal.h"
#include "alglibmisc.h"
#include "ap.h"
#include "linalg.h"
#include "specialfunctions.h"
#include "statistics.h"
using namespace std;

/*
Takes as input long/lat data from sites sampled in nature.  Computes a partition
for these sites using the partition for the location weighted random samples from envclus.

to compile:g++ calcpart.cpp partdist.cpp Lmunkres.cpp BipartiteGraph.cpp Hungarian.cpp ap.cpp statistics.cpp specialfunctions.cpp alglibinternal.cpp alglibmisc.cpp linalg.cpp -o ../bin/calcpart
usage:  calcpart envpartfile genpartfile calcpartoutfile
where,
envpartfile = path to file containing environmental partition, formatted like output from envclus 
genpartfile = path to file containing individual ID, long/lat, and partitions
	for a series of K values, see example files for format
calcpartoutfile = path to output file, will contain a table of ranked partition distances between env and gen
example: ./calcpart ./sampledsites.txt ./genpartSTRUCTURE.txt ./cpout.txt
*/

 


/***************VARIABLES*****************/
double pi = 3.14159265358979323846;
double earthRadiusKm = 6371.0;


/***************SUBROUTINES*****************/
//reads a tab delimited table, returns a 2D string vector
vector<vector<string> > readmat(const char* PathToFile)
{
	//find number of samples in genetic partitions file
	int i, j, c;
	string foo;  //will contain current line
	ifstream infile(PathToFile); //set up input file stream
	
	if (infile.good()) //test whether file exists and is readable
	{
		c = 0;  //c will tell the number of samples
		while ( !infile.eof() ) //count total number of lines to size vector
		{	
			getline(infile, foo);
			c++;
		}
		infile.clear();
		infile.seekg(0); //go back to beginning of infile

		//read genetic partitions file into a 2D vector
		vector<vector<string> > PartitionMatrix(c-1);
		i=0;		
		string field;
		while( !infile.eof() )  //read through rest of file to eof
		{
			getline(infile, foo);
			stringstream ss(foo);

			j = 0;
			while (getline(ss, field, '\t')) //split stream ss into fields using tab as delimiter
			{
				//remove any hanging newline character from field
				field.erase(std::remove(field.begin(), field.end(), '\n'), field.end());
				field.erase(std::remove(field.begin(), field.end(), '\r'), field.end());
				PartitionMatrix[i].push_back(field);
				j++;
			}
			i++;
		}
	
		return PartitionMatrix;
	}
	
	else 
	{
		cout << "The file: " << PathToFile << " either does not exist or is not readable.  Quitting...\n\n";
		exit(-1);
	}
}

// This function converts decimal degrees to radians
double deg2rad(double deg) 
{
  	return (deg * pi / 180);
}

//  This function converts radians to decimal degrees
double rad2deg(double rad) 
{
  	return (rad * 180 / pi);
}

/*
* Returns the distance between two points on the Earth.
* Direct translation from http://en.wikipedia.org/wiki/Haversine_formula
* @param lat1d Latitude of the first point in degrees
* @param lon1d Longitude of the first point in degrees
* @param lat2d Latitude of the second point in degrees
* @param lon2d Longitude of the second point in degrees
* @return The distance between the two points in kilometers
*/
double distanceEarth(double lon1d, double lat1d, double lon2d, double lat2d) 
{
	double lat1r, lon1r, lat2r, lon2r, u, v;
	lat1r = deg2rad(lat1d);
	lon1r = deg2rad(lon1d);
	lat2r = deg2rad(lat2d);
	lon2r = deg2rad(lon2d);
	u = sin(lat2r - lat1r);
	v = sin(lon2r - lon1r);
	return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

//find nearest neighbors of genetic samples in location-weighted random sample from niche model
vector<int> FindNNR(vector<vector<string> > GenPart, vector<vector<string> > EnvPart)
{
	//perform naive nearest neighbor search (calculate all pairwise great circle distances)
	cout << "Finding nearest neighbors...\n";
	int i, j;
	int n; //nearest neighbor index in EnvPart
	double gc ; //great circle distance
	double gcmin; //great circle distance to current closest point
	vector<int> nnr; //will hold the index in EnvPart of the nearest neighbor for all samples in GenPart
	double lon1d, lat1d, lon2d, lat2d;
	stringstream fs;
	for (i=1;i<GenPart.size();i++) //GenPart includes header, skip that by using i=1
	{
		cout << "  Sample " << i << endl;
		gcmin = 999999999999.9; //set to arbitrary large value
		
		fs.str(GenPart[i][1]); //conversion to double using stringstream
		fs >> lon1d;
		fs.str("");
		fs.clear();
			
		fs.str(GenPart[i][2]);
		fs >> lat1d;
		fs.str("");
		fs.clear();

		for (j=1;j<EnvPart.size();j++)
		{
			fs.str(EnvPart[j][1]);
			fs >> lon2d;
			fs.str("");
			fs.clear();

			fs.str(EnvPart[j][2]);
			fs >> lat2d;
			fs.str("");
			fs.clear();

			gc = distanceEarth(lon1d, lat1d, lon2d, lat2d);	
			
			//cout << "i=" << i << " j=" << j << "  gc=" << gc << "   lon1d=" << lon1d << " lat1d=" << lat1d << " lon2d=" << lon2d << " lat2d=" << lat2d << endl;
			//double distanceEarth(double lon1d, double lat1d, double lon2d, double lat2d) 
			if (gc < gcmin) 
			{
				gcmin = gc;
				n = j; 
			}
		}
		nnr.push_back(n);
	}
	
	return nnr;
}

//calculate all pairwise partition distances between environmental and genetic partitions
vector<vector<int> > calcPDmat(vector<vector<string> > GenPart, vector<vector<string> > EnvCGen, int g, int e, vector<pair<int, int> >& pqmatrix, ofstream& outf)
{
	int i, j, c, m;
	vector<vector<int> > gen(g), env(e); //2D vector to contain all partitions for gen and env data
	vector<vector<int> > PDmatrix(e);
	
	//create vector of env partitions
	for (c=0;c<e;c++) //column index
	{
		for (i=1;i<EnvCGen.size();i++)  //start at index 1 because of header
		{
			m = atoi(EnvCGen[i][c+2].c_str()); //convert string to int
			env[c].push_back(m); //partitions start at column index 2 hence (c+2)
		}
	}		
			
	//create vector of gen partitions
	for (c=0;c<g;c++) //column index
	{
		for (i=1;i<GenPart.size();i++)  //start at index 1 because of header
		{
			m = atoi(GenPart[i][c+3].c_str()); //convert string to int
			gen[c].push_back(m); //partitions start at column index 3 hence (c+3)
		}
	}		
	
	//create ofstream to accept output for partdist, notably from Liu's program called therein
  	//only necessary if running as standalone, otherwise cout redirected to log file referenced by outf
  	/*
  	ofstream outf; //open file stream for log file
	outf.open ("./log.txt");
	outf.close(); //quick open close done to clear any existing file each time envclus is run
	outf.open ("./log.txt", ios::out | ios::app); //open file in append mode
	*/

	//set up status indicator
	c=0;
	m=0;
	double d=env.size()/6.0;
	vector<string> wordvec;
	string word;
	stringstream timer("Good things come to those who");
	while (timer.good()) //load vector with catch phrase
	{
    	timer >> word;
    	wordvec.push_back(word);
    }  
    
	//calculate distance between gen(j) and env(i)
	int pd;
	for (i=0;i<env.size();i++)
	{
		//display status of calculations
		c++;
		if ( c >= d )
		{
			cout << "  " << string(m, ' ') + wordvec[m] << endl;
			m++;
			d = d + env.size()/6.0;
		}
 
		//perform calculations
		for (j=0;j<gen.size();j++)
		{
			//redirect output from Liu's program and partdist to log.txt, served by the reference outf
    		std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    		std::cout.rdbuf(outf.rdbuf()); //redirect std::cout to log.txt!

			pd = partdist(env[i], gen[j], outf, pqmatrix);
			//redirect cout back to standard output
			std::cout.rdbuf(coutbuf);
			
			PDmatrix[i].push_back(pd);
		}

	}
	
	outf.close();

	return PDmatrix;
}

//normalize pairwise partition distance matrix using maximum partition distance possible
//uses formulas for maximum "transfer distance" from Charon et al. 2006, Journal of Classification 23:103
vector<vector<double> > normPDmat(vector<vector<int> > PDmatrix, vector<pair<int, int> >& pqmatrix, int n)
{
	int i, j, m, envmax, genmax;
	double p, q, c, pdmax;
	vector<vector<double> > nPDmatrix(PDmatrix.size());
		
	m=0;
	for (i=0;i<PDmatrix.size();i++)
	{
		for (j=0;j<PDmatrix[0].size();j++)
		{
			//determine p and q
			envmax = pqmatrix[m].first;
			genmax = pqmatrix[m].second;
			
			if (envmax < genmax)
			{
				p = double(envmax);
				q = double(genmax);
			}
			else
			{
				p = double(genmax);
				q = double(envmax);
			}
		
			//calculate max possible partition distance, 3 cases
			if (n <= (p + q - 2)) 
			{
				pdmax = ( (2 * n) - p - q); 
				//cout << "case 1, pdmax=" << pdmax << endl;
			}											 
			else if ( ((p + q - 1) <= n) && (n <= ((p - 1) * q)) )
			{
				pdmax = ( n - (ceil ((n + q - p) / q)) );
				//cout << "case 2, pdmax=" << pdmax << endl;
			}
			else if ( ((p - 1) * q) < n  )
			{
				pdmax = ( n - (ceil (n / q)) );
				//cout << "case 3, pdmax=" << pdmax << "n/q=" << n/q <<endl;
			}
			
			//add normalized value to nPDmatrix
			c = ( PDmatrix[i][j] / pdmax );
						
			nPDmatrix[i].push_back(c);
			m++;
		}
	}
	
	cout << "        wait.\n";

	return nPDmatrix;

}

//statistical analysis of pairwise partition distances
double dostats(vector<vector<double> > nPDmatrix, vector<vector<string> > EnvPart, vector<vector<string> > GenPart, const char* PathToCalcPartOutput)
{
	int i, j;
	
	//find best mapping of environmental and genetic partitions
	double c = 1.0; //minimum value
	pair<int, int> best; //in order env, gen
	
	for (i=0;i<nPDmatrix.size();i++)
	{
		for (j=0;j<nPDmatrix[0].size();j++)
		{
			if (nPDmatrix[i][j] < c) 
			{
				best = make_pair(i, j);
				c = nPDmatrix[i][j];
			}	
		}
	}
	
	//perform Jarque-Bera normality test on all pairwise partition distances, uses ALGLIB
	alglib::ae_int_t n;
	n = nPDmatrix.size() * nPDmatrix[0].size();
	int d;
	double e, f, p;
	alglib::real_1d_array x;
	x.setlength(n);
	vector<double> nPD1D; //convert 2D to 1D for later rank ordering
		
	d=0;
	for (i=0;i<nPDmatrix.size();i++)
	{
		for (j=0;j<nPDmatrix[0].size();j++)
		{
			//for Jarque-Bera
			e = nPDmatrix[i][j];
			x[d] = e;
			d++;
			//for ranking
			nPD1D.push_back(e);
		}
	}

	//run normality test
	alglib::jarqueberatest(x, n, p);
	
	//run additional statistics
	double mean, variance, skewness, kurtosis;
	samplemoments(x, mean, variance, skewness, kurtosis);

	//calculate probability of each partition distance, assuming normal
	double cdf;
	vector<double> nPDcdf; //make mean 0 and unit variance normal distribution like empirical distribution
	for (i=0;i<nPDmatrix.size();i++)
	{
		for (j=0;j<nPDmatrix[0].size();j++)
		{
			//for std normal probability estimation
			e = nPDmatrix[i][j];
			f = (e - mean)/(sqrt(variance)); //normalize to mean 0, unit variance normal distribution
			cdf = alglib::normaldistribution(f); //calculate area under the PDF, from -infinity to f
			nPDcdf.push_back(cdf);
		}
	}
	
	//rank order partition distances
	pair<int,int> yx;
	pair<double,pair<int,int> > PDyx;
	vector<pair<double,pair<int,int> > > ranker;
	//construct vector of nested pairs to retain index position before sort: PD (row, col)
	for (i=0;i<nPDmatrix.size();i++)
	{
		for (j=0;j<nPDmatrix[0].size();j++)
		{
			yx = make_pair(i,j);
			e = nPDmatrix[i][j];
			PDyx = make_pair(e,yx);
			ranker.push_back(PDyx);
		}
	}
	//sort vector of pairs by PDyx.first()
	std::sort(ranker.begin(), ranker.end());
			
	
	//print out results
	cout << "Statistics on matrix of pairwise partition distances:\n";
	cout << "  n = " << x.length() << endl;
	cout << "  mean = " << mean << endl;
	cout << "  sd = " << sqrt(variance) << endl;
	cout << "  skewness = " << skewness << endl;
	cout << "  kurtosis = " << kurtosis << endl;
	cout << "  Jarque-Bera normality test, p = " << p << endl << endl; 
	
	//optimal mapping
	cout << "The optimal mapping of genetic structure to environmental variables occurs for:\n";
	cout << "  env: " << EnvPart[0][(best.first + 6)] << "\n  gen: " << GenPart[0][(best.second + 3)] << endl; 
	
	//table of ranked mappings with p value, print to terminal and write to file
	cout << "\nRank order correspondence between genetic and environmental partitions\n(cdf assumes normality):\n";
	cout << "\nrank\tpd\tenv\tgen\tcdf\n";
	int pind; //pind is the index for a 1D vector calculated from 2D vector[i][j]
	ofstream outputf;
	outputf.open(PathToCalcPartOutput);
	outputf << "rank\tpd\tenv\tgen\tcdf\n";
	
	for (i=0;i<ranker.size();i++)
	{
		e = ranker[i].first;
		yx = ranker[i].second;
		j = yx.first; //env row index
		d = yx.second; //gen col index
		pind = j * nPDmatrix[0].size() + d;
		cdf = nPDcdf[pind]; 
		
		//print to terminal
		cout << i+1 << "\t" << e << "\t" << EnvPart[0][(j + 6)] ;
		cout << "\t" << GenPart[0][(d + 3)] << "\t" << cdf << endl;
		//write to outfile
		outputf << i+1 << "\t" << e << "\t" << EnvPart[0][(j + 6)] ;
		outputf << "\t" << GenPart[0][(d + 3)] << "\t" << cdf << endl;
	}
	outputf.close();

	
	return ranker[0].first;  //return smallest partition distance
}

/***************MAIN*****************/

double ibdibe (const char* PathToEnvPartition, const char* PathToGenPartition, const char* PathToCalcPartOutput, ofstream& outf)
{
	int i, j, n;
	/*
	const char* PathToEnvPartition = argv[1];
	const char* PathToGenPartition = argv[2];
	const char* PathToCalcPartOutput = argv[3];
	*/
		
	//read in genetic partitions file
	vector<vector<string> > GenPart;
	GenPart = readmat(PathToGenPartition);
	n = ( GenPart.size() - 1 );
	cout << "Number of samples in genetic partition = " << n << endl;

	//read in sampledsites file (with environmental partitions)
	vector<vector<string> > EnvPart;
	EnvPart = readmat(PathToEnvPartition);
	
	//find nearest neighbors to sites sampled for genetic data
	vector<int> nnr; //index to nearest neighbors in EnvPart
	nnr = FindNNR(GenPart, EnvPart);
	
	//calculate matrix of partition distances
	int e, g; //number of partitions to consider for environmental and genetic data
	g = ( GenPart[0].size() - 3 );  //number of items in header minus Indiv, Long, Lat columns.
									//gen partitions at column index 3 thru 3+g
	e = ( EnvPart[0].size() - 7 );  //number of items in header minus sample, long, lat, column, row, logisticSuitability, seed columns.
									//env partitions at column index 6 thru 6+e. "seed" at end
	
	//create env partition corresponding to genetic sample localities
	vector<vector<string> > EnvCGen( (nnr.size()+1) );
	EnvCGen[0].push_back("nnlong"); //load header info into EnvCGen
	EnvCGen[0].push_back("nnlat");
	i=0;
	while (i<e)
	{
		EnvCGen[0].push_back(EnvPart[0][6+i]);
		i++;
	}
	//extract env partition info for gen sites using index values in nnr
	int m;
	for (i=0;i<nnr.size();i++)
	{
		m = nnr[i]; //get the line index from nnr
		EnvCGen[i+1].push_back(EnvPart[m][1]); //add longitude for nearest site
		EnvCGen[i+1].push_back(EnvPart[m][2]); //add latitude
		for (j=6;j<6+e;j++) //add columns containing partition information
		{
			EnvCGen[i+1].push_back(EnvPart[m][j]);
		}
	}

	
	//calculate all pairwise partition distances between GenPart and EnvCGen
	cout << "Calculating partition distances...\n";
	vector<vector<int> > PDmatrix;
	vector<pair<int, int> > pqmatrix;
	PDmatrix = calcPDmat(GenPart, EnvCGen, g, e, pqmatrix, outf);	
	
	//normalize matrix of partition distances using maximum possible values
	vector<vector<double> > nPDmatrix;
	nPDmatrix = normPDmat(PDmatrix, pqmatrix, n);
	
	//compute rank order and other stats for distribution of partition distances nPDmatrix
	double PDmin; //minimum partition distance observed
	PDmin = dostats(nPDmatrix, EnvPart, GenPart, PathToCalcPartOutput);

	
	/*****IBD/IBE test******/
	//generate null distribution
	vector<double> IBDpds; //contains partition distances between randomly selected regions and genetic data
	nullib(EnvPart, GenPart, numresamples);
	
	//
	
	
	
			
	cout << "Done!\n\n";	
		
/*		
		//print out PDmatrix
	cout << "PDmatrix:\n";
	for (i=0;i<PDmatrix.size();i++)
	{
		for (j=0;j<PDmatrix[i].size();j++)
		{
			cout << PDmatrix[i][j] << " ";
		}
		cout << endl;
	}	

	cout << "nPDmatrix:\n";
	for (i=0;i<nPDmatrix.size();i++)
	{
		for (j=0;j<nPDmatrix[i].size();j++)
		{
			cout << nPDmatrix[i][j] << " ";
		}
		cout << endl;
	}	
*/
	return PDmin;	
}
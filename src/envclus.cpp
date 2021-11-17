#include <cstdlib>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip> //for defining floating point precision
#include <cmath>
#include <vector>
#include <numeric>
#include <boost/filesystem.hpp>
#include <iterator>
#include "envclus.hpp"
namespace fs = boost::filesystem;
using namespace std;



/*
to compile:  g++ envclus.cpp kmeans.cpp ktests.cpp -o ../bin/envclus -I/opt/local/include  -L/opt/local/lib -lboost_system-mt -lboost_filesystem-mt
usage:  envclus numsamples kmax gridfile envlayersdir
where, 
numsamples = number of random points to sample, location weighted
kmax = maximum number of clusters desired, program will consider K=2 to Kmax automatically
gridfile = path to asc format GIS layer containing location weights, e.g. a niche model 
envlayersdir = path to a directory containing environmental layers in ESRI ascii format

example: ./envclus 1000 4 ./pumilus_avg.asc ./asciifiles/
*/


/***************GLOBAL VARIABLES*****************/
bool debug = false;  //if true: seed is fixed value, ...



/***************SUBROUTINES*****************/
//check for unix style line breaks
/*void VerifyLineBreaks(const char* GridFile)
{
	string foo;
	ifstream infile(GridFile, std::ifstream::binary); 
	//test for windows or old mac style line breaks, they must be unix
	getline(infile, foo);
	string query = "ncols";
	size_t pos;
	pos = foo.find(query);
	
	cout << "pos=" << pos << "\n";
	getchar();
	
	if (pos != std::string::npos)
	{
		//do nothing	
	}
	else
	{
		cout << "File " << GridFile << " does not have Unix line breaks.  Please convert.  Exiting...\n";
		exit (EXIT_FAILURE);
	}
}
*/

//check that header of grid files is consistent with ascii standard
void VerifyHeader(const char* GridFile)
{
	string foo;
	ifstream infile(GridFile, std::ifstream::binary); //binary mode to ignore line endings 

	//test whether it complies with ascii grid standard, i.e. first line contains "ncols"
	getline(infile, foo);
	string query = "ncols";
	size_t pos;
	pos = foo.find(query);
	if (pos != std::string::npos)
	{
		//do nothing	
	}
	else
	{
		cout << "File " << GridFile << " does not appear to be an ascii grid file.  Exiting...\n";
		exit (EXIT_FAILURE);
	}
}

//get grid file specifications
vector<double> GetGridSpecs(const char* GridFile)
{
	//create a vector containing grid specs
	string foo;
	ifstream infile(GridFile, std::ifstream::binary); //set up input file stream, binary to ignore line endings
	//extract consecutive lines, get last word, convert to proper variable type
	
	//ncols
	getline(infile, foo);
	string sncols = foo.substr(foo.rfind(" ")+1);
	
	//nrows
	getline(infile, foo);
	string snrows = foo.substr(foo.rfind(" ")+1);
	
	//xllcorner
	getline(infile, foo);
	string sxllcorner = foo.substr(foo.rfind(" ")+1);
	
	//yllcorner
	getline(infile, foo);
	string syllcorner = foo.substr(foo.rfind(" ")+1);
	
	//cellsize
	getline(infile, foo);
	string scellsize = foo.substr(foo.rfind(" ")+1);
	
	//NODATA_value
	getline(infile, foo);
	string sNODATA_value = foo.substr(foo.rfind(" ")+1);
	
	//convert to number formats, use .c_str() --must be a C string to use atoi
	int ncols, nrows, NODATA_value;
	long double xllcorner, yllcorner, cellsize;

	ncols = atoi(sncols.c_str());
	nrows = atoi(snrows.c_str());
	xllcorner = strtold(sxllcorner.c_str(), NULL); //string to long double
	yllcorner = strtold(syllcorner.c_str(), NULL);
	cellsize = strtold(scellsize.c_str(), NULL);
	NODATA_value = atoi(sNODATA_value.c_str());
	
	//add to vector and return
	vector<double> GridSpecs(6);
	GridSpecs[0] = ncols;
	GridSpecs[1] = nrows;
	GridSpecs[2] = xllcorner;
	GridSpecs[3] = yllcorner;
	GridSpecs[4] = cellsize;
	GridSpecs[5] = NODATA_value;

	return GridSpecs;
}


//verify that grid files are properly formatted and all have the same specs
void CompareGridSpecs(const char* GridFile, char* EnvLayersDir)
{
	int i;
	string foo;
	
	//verify that line breaks are Unix
	//VerifyLineBreaks(GridFile);
	//verify header on niche model
	VerifyHeader(GridFile);
	
	//verify header, line breaks on environmental layers
	//get vector of paths from EnvLayersDir
	typedef vector<fs::path> pvec;  // create a vector data type that contains paths       
	pvec EnvLayers;                            
	
	//copy the contents of folder EnvLayersDir into pvec EnvLayers, uses Boost::Filesystem
	copy(fs::directory_iterator(EnvLayersDir), fs::directory_iterator(), back_inserter(EnvLayers)); // directory_iterator::value_type
	
	sort(EnvLayers.begin(), EnvLayers.end()); //sort paths
	
	//verify each environmental layer	
	for (i=0; i < EnvLayers.size(); i++)
	{
		foo = EnvLayers[i].string();  //get path out of vector as string
		const char* CurrLayer = foo.c_str(); //convert string to const char*
		//VerifyLineBreaks(CurrLayer);
		VerifyHeader(CurrLayer);
	}
	

	//confirm that all grid files have the same specs
	//get specs for niche model
	vector<double> NicheGridSpecs(6);
	vector<double> EnvGridSpecs(6);
	NicheGridSpecs = GetGridSpecs(GridFile);
	
	//compare to each environmental layer
	for (i=0; i < EnvLayers.size(); i++)
	{
		foo = EnvLayers[i].string();  //get path out of vector as string
		const char* CurrLayer = foo.c_str(); //convert string to const char*
		vector<double>().swap(EnvGridSpecs); //clear EnvGridSpecs
		EnvGridSpecs = GetGridSpecs(CurrLayer);
		if (EnvGridSpecs != NicheGridSpecs)
			{
      		  	cout << "\nSpecs of grid " << CurrLayer << " differ from niche model.  Exiting...\n\n";
      		  	exit (EXIT_FAILURE);
    	    }
	}
}



//read in an ESRI ascii format GIS grid file, GridFile
//return by reference grid file specs: ncols nrows xllcorner yllcorner cellsize NODATA_value
//returns a vector containing the entire grid loaded into memory
vector<vector<double> > ReadGrid(const char* GridFile, int& ncols, int& nrows, long double& xllcorner, long double& yllcorner, long double& cellsize, int& NODATA_value)
{
	//local variables
	string foo;  //contains current line of file
	
	ifstream infile(GridFile, std::ifstream::binary); //set up input file stream for Env model
	//extract consecutive lines, get last word, convert to proper variable type
	
	//ncols
	getline(infile, foo);
	string sncols = foo.substr(foo.rfind(" ")+1);
	
	//nrows
	getline(infile, foo);
	string snrows = foo.substr(foo.rfind(" ")+1);
	
	//xllcorner
	getline(infile, foo);
	string sxllcorner = foo.substr(foo.rfind(" ")+1);
	
	//yllcorner
	getline(infile, foo);
	string syllcorner = foo.substr(foo.rfind(" ")+1);
	
	//cellsize
	getline(infile, foo);
	string scellsize = foo.substr(foo.rfind(" ")+1);
	
	//NODATA_value
	getline(infile, foo);
	string sNODATA_value = foo.substr(foo.rfind(" ")+1);
	
	//convert to number formats, use .c_str() --must be a C string to use atoi
	ncols = atoi(sncols.c_str());
	nrows = atoi(snrows.c_str());
	xllcorner = strtold(sxllcorner.c_str(), NULL); //string to long double
	yllcorner = strtold(syllcorner.c_str(), NULL);
	cellsize = strtold(scellsize.c_str(), NULL);
	NODATA_value = atoi(sNODATA_value.c_str());

	//print out
	cout << "\nReading gridfile " << GridFile << "\n";
	cout << setprecision(15);
	cout << "Grid file specs:\n";
	cout << "  ncols = " << ncols << "\n";
	cout << "  nrows = " << nrows << "\n";
	cout << "  xllcorner = " << xllcorner << "\n";
	cout << "  yllcorner = " << yllcorner << "\n";
	cout << "  cellsize = " << cellsize << "\n";
	cout << "  NODATA_value = " << NODATA_value << "\n";
	
	//read rest of asc grid file (infile) into a 2D vector of size ncol X nrow
	vector<vector<double> > GISGrid(nrows, vector<double>(ncols));  //declare vector, define size
	
	int s = 0; //monitors status of GIS file load
	int l = 5; //describes % of GIS file loaded
	int i = 0; //iterates through rows
	double f;
	while( !infile.eof() )
	{
		getline(infile, foo); //remove header line
		stringstream ss(foo);
		string field;

		int j = 0; //iterates through cols

		while (getline(ss, field, ' ')) //split stream ss into fields using space as delimiter
		{
			stringstream ds(field); //conversion to double using stringstream
			ds >> f;
			GISGrid[i][j] = f;
			j++;
		}
		i++;
		s++;
		//display status of read
		if ( (s >= (nrows/20)) && (l <= 100) )
		{
			cout << "  " << l << "%\n";
			s = 0; //reset s
			l = l + 5;
		}
	}

	return GISGrid;
}





//modifies a reference to a 2D vector, LongLatVec, to contain a random sample from a 
//species distribution, location weighted using a habitat suitability model.
//Also, writes the sampled cells to OutFile
//Also return a vector of grid specifications
vector<double> ranlocwt(int NumSamples, const char* GridFile, vector<vector<int> >& LongLatVec, const char* OutFile)
{
	// & LongLatVec is a reference to memory allocation for 2D vector with NumSamplesX2, sub returns void
	
	//declare variables
	int ncols, nrows, NODATA_value;
	long double xllcorner, yllcorner, cellsize;
	vector<vector<double> > GISGrid;  //declare vector
	
	//Read in species distribution model
	GISGrid = ReadGrid(GridFile, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value);
	
	//create a vector containing grid specs
	vector<double> EnvGridSpecs(6);
	EnvGridSpecs[0] = ncols;
	EnvGridSpecs[1] = nrows;
	EnvGridSpecs[2] = xllcorner;
	EnvGridSpecs[3] = yllcorner;
	EnvGridSpecs[4] = cellsize;
	EnvGridSpecs[5] = NODATA_value;
	
	
	//generate list of random, location weighted samples, write to OutFile
	cout << "Sampling...\n";
	ofstream outputf;
	outputf << setprecision(10);
	outputf.open(OutFile);
	
	//set random number seed
	int tt;
	if (debug == false) tt = (time(NULL));
	else if (debug == true) tt = 12345;  //1383318854;
	srand ( tt ); //initialize
	
	//determine random samples, write to file
	int c = 0;
	outputf << "sample	long	lat	column	row	logisticSuitability	seed=" << tt << "\n";
	while (c<NumSamples)
	{
		//get random row and column
		int ranrow = rand() % nrows;
		int rancol = rand() % ncols;
		
		//extract suitability value at that random cell
		double suit = GISGrid[ranrow][rancol];
		
		//get random number 0..1, test whether it is smaller than suit, if so, retain ranrow,rancol
		//this has the effect of excluding any regions with NO_DATA values, assuming NO_DATA value < 0
		//it has a problem in that it is inefficient because the input niche model may not have
		//a maximum cell value of 1.  But the weighted sampling works anyway.
		double z = ((double) rand() / (RAND_MAX));
		if (z <= suit)
		{
		
		//calculate long/lat
		double longi = (xllcorner + (rancol*cellsize));
		double lat = ( yllcorner + ( (nrows-ranrow)*cellsize ) );
		
		//write to output file
		outputf << c << "\t" << longi << "\t" << lat << "\t";
		outputf << (rancol+1) << "\t" << (ranrow+1) << "\t" << suit << "\n"; //correct so column and row are indexed starting with 1 not 0
		
		//write to referenced 2D array memory location
		LongLatVec[c][1] = rancol;
		LongLatVec[c][2] = ranrow;
		
		c++;
		}
	}	
	outputf.close();
	cout << "Sampled sites written to: " << OutFile << endl;;
	
	return EnvGridSpecs;
}


//Takes as input a vector of location weighted, randomly sampled sites (LongLatVec),
//and a directory (EnvLayersDir) containing ESRI ascii formatted environmental layers.
//Returns environmental data for all sites in LongLatVec, for all layers, formatted
//for K means clustering.  Uses boost::filesystem libraries.
vector<vector<double> > exenv(int NumSamples, vector<vector<int> >& LongLatVec, char* EnvLayersDir)
{
	
	// create a vector data type to contain paths
	typedef vector<fs::path> pvec;             
	pvec EnvLayers;                            
	
	//copy the contents of folder EnvLayersDir into a vector, EnvLayers
	copy(fs::directory_iterator(EnvLayersDir), fs::directory_iterator(), back_inserter(EnvLayers)); // directory_iterator::value_type
	
	//sort paths
	sort(EnvLayers.begin(), EnvLayers.end());
	
	//extract data at sampled sites for all layers in EnvLayers vector of paths
	int ncols, nrows, NODATA_value;
	long double xllcorner, yllcorner, cellsize;
	vector<vector<double> > GISGrid;  //declare 2D vector
	cout << "Extracting data from environmental layers.\n";
	
	int NumLayers = EnvLayers.size();
	vector<vector<double> > EnvDataMatrix(NumSamples, vector<double>(NumLayers));  //declare vector, define size, 2D vector holds extracted data for environmental layers
	double ccol;
	double rrow;
	int i;
	int j;
	for (i=0; i < NumLayers; i++)
	{
		//read environmental layer into memory as GISGrid
		vector<vector<double> >().swap(GISGrid); //clear GISGrid
		string foo = EnvLayers[i].string();  //get path out of vector as string
		const char* CurrLayer = foo.c_str(); //convert string to const char*
		GISGrid = ReadGrid(CurrLayer, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value);
			
		//extract data at LongLatVec
		cout << "Extracting from " << CurrLayer << "\n";
		//cout << "ccol	rrow	EnvDataMatrix\n";
		for (j=0; j < NumSamples; j++)
		{
			ccol = LongLatVec[j][1];
			rrow = LongLatVec[j][2];
			EnvDataMatrix[j][i] = GISGrid[rrow][ccol]; //write vertically, column by column
			//cout << ccol << "\t" << rrow << "\t" << EnvDataMatrix[j][i] << "\n";
		}
	}

	//purge EnvDataMatrix of invariant variables
	double y, z;
	vector<int> ColToKeep; //list of columns to retain, i.e. that contain variable data
	//find invariant variables
	for (j=0;j<EnvDataMatrix[0].size();++j)
	{
		z = EnvDataMatrix[0][j]; //get first item of column j
		for (i=1;i<EnvDataMatrix.size();++i)
		{
			y = EnvDataMatrix[i][j];
			if (y != z) //keep the column if all items are not identical to the first item
			{
				ColToKeep.push_back(j);
				break;	
			}
		}
	
	}
	if (ColToKeep.size() != EnvDataMatrix.size())
	{
		cout << "\nPurging invariant variables...";
		
		//identify the columns to delete
		vector<int> ColToDelete;
		for (i=0;i<EnvDataMatrix[0].size();++i)
		{
			if ( std::find(ColToKeep.begin(), ColToKeep.end(), i) != ColToKeep.end() )
			{
				//do nothing, ColToKeep contains the index i
			}
			else
			{
				ColToDelete.push_back(i); //add missing column to ColToDelete
			}
		}
		//purge the invariant columns, do it backwards because the 2d vector shrinks
		for ( j=(ColToDelete.size() - 1);j>=0;j=j-1 )
		{
			z = ColToDelete[j];
			for (i=0;i<EnvDataMatrix.size();++i)
			{
				EnvDataMatrix[i].erase(EnvDataMatrix[i].begin() + z); 
			}
		}
		//redefine the number of environmental layers
		NumLayers = ColToKeep.size();
	}
	
	//standardize data to zero mean and unit variance
	cout << "\nStandardizing...\n";
	double sum;
	double ssum;
	double mean;
	double sq_sum;
	double stdev;
	vector<double> v(NumSamples); //temporary vector to hold values for one environmental variable
	vector<double> diff(v.size());
	vector<double> nv(NumSamples); //temporary vector to hold standardized values
	vector<vector<double> > StdEnvDataMatrix(NumSamples, vector<double>(NumLayers));  //declare vector, define size, 2D vector will hold standardized data for environmental layers
		
	for (i=0; i<NumLayers; i++)
	{
		//extract one column in EnvDataMatrix to a vector
		for (j=0; j<NumSamples; j++)
		{
			v[j] = EnvDataMatrix[j][i];
		}
		
		//calculate mean
		ssum = std::accumulate(v.begin(), v.end(), 0.0);
		sum = 0;
		for (j=0;j<v.size();++j) sum = sum + v[j];
		if (sum != ssum)
		{
			cout << "red alert! sum = "<<sum<<" while ssum = "<<ssum<<"\n";
			getchar();
		}
		mean = sum / v.size();
		//calculate sample stdev
		std::transform(v.begin(), v.end(), diff.begin(), std::bind2nd(std::minus<double>(), mean));
		sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
		stdev = std::sqrt(sq_sum / (v.size() - 1));
		
		//standardize v, add to 2D vector of standardized environmental values
		for (j=0; j<NumSamples; j++)
		{
			nv[j] = ((v[j] - mean) / stdev);
			StdEnvDataMatrix[j][i] = nv[j];
		}		
	}
	
	//write environmental data matrix to file
	ofstream outputs;
	outputs << setprecision(10);
	outputs.open("./stdenvdatamatrix.txt");
	
	ofstream outputf;
	outputf << setprecision(10);
	outputf.open("./envdatamatrix.txt");
	
	//write out header
	i=0;
	outputf << "#";
	outputs << "#";
	for (j=0;j<ColToKeep.size();++j)
	//while (i<NumLayers)
	{
		i = ColToKeep[j]; //include only those columns that are variable
		string foo = fs::basename(EnvLayers[i]); //get environmental layer file name
		outputf << foo << "\t";
		outputs << foo << "\t";
	}
	outputf << "\n";
	outputs << "\n";
	
	//write out data matrix
	i=0;
	while (i<NumSamples)
	{
		j=0;
		while (j<NumLayers)
		{
			outputf << EnvDataMatrix[i][j] << "\t";
			outputs << StdEnvDataMatrix[i][j] << "\t";
			j++;
		}
		outputf << "\n";	
		outputs << "\n";	
		i++;
	}
	cout << "Resampled environmental data written to: ./envdatamatrix.txt\n\n";

	return EnvDataMatrix;
}



//adds a column of cluster assignment values to output file OutFile
void addclusters(int k, string bestr, const char* OutFile)
{
	int j;
	string foo;
	string fig;
	string foonew;
	ifstream infile(OutFile);
	//deal with header
	getline(infile, foo);

	//separate header line into a vector of column titles
	stringstream ss(foo);
	string field;
	vector<string> ct;
	while (getline(ss, field, '\t'))  //split stream ss into fields using tab as delimiter
	{
		ct.push_back(field);
	}

	//write out new header to temporary output file
	ofstream outputf;
	outputf.open("./.sampledsites.tmp");
	for (j=0; j<(ct.size() - 1);j++)
	{
	outputf << ct[j] << "\t";
	}
	outputf << "K=" << k << "\t" << ct[(ct.size()-1)] << "\n"; //add new cluster column header plus old random # seed column header

	
	//write out rest of new file, line by line
	ifstream clusfile(bestr.c_str());
	string f;

	while( true )
	{
		getline(infile, foo);
		getline(clusfile, fig);
		if ( infile.eof() ) break;
		stringstream ss(fig); //remove spaces
		ss >> f;
		outputf << foo << "\t" << f << "\n";
	}
	outputf.close();

	
	//manage files
	rename("./.sampledsites.tmp", OutFile);
	if (debug == false || debug == true)
	{
		remove(bestr.c_str());
	}
}


//chooses the best clustering scheme returned by ktests() by minimizing energy (variance)
string bestclusters(int k, vector<vector<double> > ResultSummary, vector<int>& kvec, vector<double>& enervec)
{
	int i;
	double t;
	int j = 0;
	
	//refine the ResultSummary to include only those rows where the KTEST converged
	vector<vector<double> > RSrefined;
	for (i=0;i<ResultSummary.size();i++)
	{
		t = ResultSummary[i][6]; 
		if (t != -9999) // -9999 is value that represents no convergence
		{ 
			RSrefined.push_back (ResultSummary[i]);
			j++;
		}
	}
	
	
	
	double best = RSrefined[0][0]; //initialize best model value, contains numeric identifier for test name
	string bestr; //return variable, returns a filename as string
	double en;
	int bestindex = 0;  //index position for best KMEANS model
	double enmin = RSrefined[0][6]; //initialize min energy value
	
	for (i=0;i<RSrefined.size();i++)
	{
		en = RSrefined[i][6]; //get energy for each kmeans test run
		if (en < enmin)  //found new lowest energy 
		{
			best = RSrefined[i][0];
			enmin = en;
			bestindex = i;
		}
	}

	//add latest K and ce_total to summary vector
	kvec.push_back(k);
	enervec.push_back(RSrefined[bestindex][6]);
	
	
	//convert best to string containing file name of best kmeans test
	stringstream ss;
	ss << ".K" << k << "test0" << best << "_clusters.txt";
	bestr = ss.str();
	
	//clean up
	if (debug == true) cout << "bestr= " << bestr << "\n";
	//delete suboptimal files, and all centers files
	if (debug == false || debug == true)
	{
		double t;
		for (i=0;i<RSrefined.size();i++)
		{
			t = RSrefined[i][0]; //get number of KMEANS test
			//remove "centers" file
			string rr;
			ss.str(""); //clear the stringstream
			ss << ".K" << k << "test0" << t << "_centers.txt"; //make file name
			rr = ss.str(); //convert to string
			remove(rr.c_str()); //convert to c_string and remove
			
			//remove "clusters" file
			if (t != best) //if current KMEANS test number is not best
			{
				ss.str(""); //clear the stringstream
				ss << ".K" << k << "test0" << t << "_clusters.txt";
				rr = ss.str(); //convert to string
				remove(rr.c_str()); //convert to c_string and remove
			}
		}
	}
	
	return bestr;
}

/***************MAIN*****************/

int envclus (const int NumSamples, const int k, const char* GridFile, char* EnvLayersDir, const char* OutFile, ofstream& outf)
{
	/*
	//get command line arguments 
	const int NumSamples = atoi(argv[1]);
	const int k = atoi(argv[2]); //k is kmax, the maximum number of clusters
	const char* GridFile = argv[3];
	char* EnvLayersDir = argv[4];
	*/


	//print out input variables
	cout << "\nInput variables:\n  NumSamples = " << NumSamples << "\n";
	cout << "  SDM GridFile = " << GridFile << "\n";
	
	//confirm that all grid files have the same specs
	CompareGridSpecs(GridFile, EnvLayersDir);
	
	//generate list of location weighted, randomly sampled sites
	int i=0;
	vector<vector<int> > LongLatVec(NumSamples, vector<int>(2));  //declare vector, define size
	ranlocwt(NumSamples, GridFile, LongLatVec, OutFile); //receive sites, LongLatVec is sent as a reference to sub ranlocwt
	
	
	//extract environmental values at location weighted, randomly sampled sites.
	//standardize values to mean zero and unit variance
	vector<vector<double> > EnvDataMatrix;  //declare vector
	EnvDataMatrix = exenv(NumSamples, LongLatVec, EnvLayersDir);
	
		
	//perform Kmeans clustering for K=2 to Kmax
	const char* PathToData = "./stdenvdatamatrix.txt";
	vector<vector<double> > ResultSummary;
	string bestr; //will contain path to cluster file for best ktest
	vector<int> kvec; //contains current k value
	vector<double> enervec; //contains current total energy value
	cout << "Clustering...\n";
		
	for (i=2; i<=k; i++)
	{
		cout << "  K = " << i << "\n";
		ResultSummary = ktests(i, PathToData, outf);
		bestr = bestclusters(i, ResultSummary, kvec, enervec); //choose best clustering scheme by minimizing energy/variance
		                                                       //add total energy and k value to kvec, enervec references
		addclusters(i, bestr, OutFile); //add cluster assignment for current K to output
	}
	
	outf.close();
	
	
	
	//create output table describing relationship between K and Kmeans test energy
	double cetotal; //cetotal is at ResultSummary[i][6]
	ofstream outputf;
	outputf.open("./kenergy.txt");
	outputf << "K	total energy	mean energy per cluster\n";
	for (i=0; i<kvec.size(); i++)
	{
		cetotal = enervec[i];
		outputf << kvec[i] << "\t" << cetotal << "\t" <<  (cetotal/(i+2)) << "\n";	
	}
	outputf.close();
		
	cout << "Done!\n\n";



	return 0;
}




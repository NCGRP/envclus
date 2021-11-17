# include <cstdlib>
# include <iostream>
# include <sstream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <complex>
#include <vector>

using namespace std;

# include "kmeans.hpp"
# include "envclus.hpp"

vector<vector<double> > ktests (int k, const char* PathToData);

void test01 (int k, string kstr, const char* PathToData);
void test02 (int k, string kstr, const char* PathToData);
vector<double> test03 (int k, string kstr, int it_max, const char* PathToData); //Applied Statistics Algorithm #58
vector<double> test04 (int k, string kstr, int it_max, const char* PathToData); //Applied Statistics Algorithm #136
vector<double> test05 (int k, string kstr, int it_max, const char* PathToData); //algorithm of unknown origin
void test06 (int k, string kstr, const char* PathToData);
void test07 (int k, string kstr, const char* PathToData);
void test08 (int k, string kstr, const char* PathToData);

//****************************************************************************80

vector<vector<double> > ktests (int k, const char* PathToData, ofstream& outf)

//****************************************************************************80
//
//  Purpose:
//
//    ktests is the main program for KMEANS_PRB.
//
//  Discussion:
//
//    KMEANS_PRB calls the KMEANS tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    1 March 2013 by Pat Reeves
//
//  Author:
//
//    John Burkardt
//
{
  
  //redirect output from KMEANS program to a log file, served by the reference outf
  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(outf.rdbuf()); //redirect std::cout to log.txt!
  
  timestamp ( );
  cout << "\n";
  cout << "KMEANS\n";
  cout << "  C++ version\n";
  cout << "  Test the KMEANS library.\n";

  //convert int k to string kstr
  stringstream ss;
  ss << k;
  string kstr = ss.str();
  
  //set up for replicated analyses
  int i;
  int it_max_start = 20; //maximum number of iterations, starting value
  int it_max; //maximum number of iterations
  int it_num;  //actual number of iterations taken
  int it_cutoff = 1280; //point at which it has taken so many iterations it seems it will never converge
  vector<double> retvec; 	//vector of diagnostics, inserted into ResultSummary
  							
  vector<vector<double> > ResultSummary(3); //size of vector is number of tests
  											//retvec[0] = 3; //kmeans method 1 = test 3
  											//retvec[1] = dim_num;   = num environmental grids
  											//retvec[2] = point_num; = num samples
  											//retvec[3] = it_max;    = max iterations allowed
  											//retvec[4] = it_num;    = num iterations taken
  											//retvec[5] = k;         = current K value
  											//retvec[6] = ce_total;  = total energy
  
  
  
  
  //exclude tests that use unreliable HMEANS algorithm
  //test01 (k, kstr, PathToData);
  //test02 (k, kstr, PathToData);
  
  //************test03************
  it_max = it_max_start;
  while (true)
  {
	retvec = test03 (k, kstr, it_max, PathToData); //perform test03, retrieve vector of run diagnostics
  	it_max = retvec[3];
  	it_num = retvec[4];
	if (it_num < it_max ) break;  //get out of the loop if fewer than max iterations were taken, or more than cutoff
  	if (it_num > it_cutoff)   //too many iterations, not converging, quit trying
  	{
  		retvec[6] = -9999;
  		std::cout.rdbuf(coutbuf); //redirect cout to standard output
  		cout << "    Kmeans algorithm 3 is not converging.  I give up..." << endl;
  		std::cout.rdbuf(outf.rdbuf()); //redirect std::cout to log.txt!
  		break;
  	}
  	it_max = it_max*2; //double it_max, repeat
  	
  	std::cout.rdbuf(coutbuf); //redirect cout to standard output
  	cout << "    Kmeans algorithm 3 did not converge.  Resetting max iterations to..." << it_max << endl;
  	std::cout.rdbuf(outf.rdbuf()); //redirect std::cout to log.txt!

  }
    ResultSummary[0] = retvec;


  //************test04************
  it_max = it_max_start;
  vector<double>().swap(retvec); //clear retvec
  while (true)
  {
	retvec = test04 (k, kstr, it_max, PathToData); //perform test04, retrieve vector of run diagnostics
  	it_max = retvec[3];
  	it_num = retvec[4];
	if (it_num < it_max) break;  //get out of the loop if fewer than max iterations were taken
  	if (it_num > it_cutoff)   //too many iterations, not converging, quit trying
  	{
  		retvec[6] = -9999;
  		std::cout.rdbuf(coutbuf); //redirect cout to standard output
  		cout << "    Kmeans algorithm 4 is not converging.  I give up..." << endl;
  		std::cout.rdbuf(outf.rdbuf()); //redirect std::cout to log.txt!
  		break;
  	}
  	it_max = it_max*2; //double it_max, repeat
  	
  	std::cout.rdbuf(coutbuf); //redirect cout to standard output
  	cout << "    Kmeans algorithm 4 did not converge.  Resetting max iterations to..." << it_max << endl;
  	std::cout.rdbuf(outf.rdbuf()); //redirect std::cout to log.txt!

  }
    ResultSummary[1] = retvec;
  

  //************test05************
  it_max = it_max_start;
  vector<double>().swap(retvec); //clear retvec
  while (true)
  {
	retvec = test05 (k, kstr, it_max, PathToData); //perform test05, retrieve vector of run diagnostics
  	it_max = retvec[3];
  	it_num = retvec[4];
	if (it_num < it_max) break;  //get out of the loop if fewer than max iterations were taken
  	if (it_num > it_cutoff)   //too many iterations, not converging, quit trying
  	{
  		retvec[6] = -9999;
  		std::cout.rdbuf(coutbuf); //redirect cout to standard output
  		cout << "    Kmeans algorithm 5 is not converging.  I give up..." << endl;
  		std::cout.rdbuf(outf.rdbuf()); //redirect std::cout to log.txt!
  		break;
  	}
  	it_max = it_max*2; //double it_max, repeat
  	
  	std::cout.rdbuf(coutbuf); //redirect cout to standard output
  	cout << "    Kmeans algorithm 5 did not converge.  Resetting max iterations to..." << it_max << endl;
  	std::cout.rdbuf(outf.rdbuf()); //redirect std::cout to log.txt!

  }
    ResultSummary[2] = retvec;
  

	if (debug == true)
	{
		int j;
		cout << "\n\n\t**debug\nResultSummary:\n";
		for (i=0; i<ResultSummary.size(); i++)
		{
			for (j=0; j<7; j++)
			{
				cout << j << " " << ResultSummary[i][j] << "\n";
			}
			cout << "\n";
		}
	}
	cout << "\t **\n\n";
  
  //exclude tests that combine HMEANS and KMEANS algorithms
  //test06 (k, kstr, PathToData);
  //test07 (k, kstr, PathToData);
  //test08 (k, kstr, PathToData);
//
//  Terminate.
//
  cout << "\n";
  cout << "KMEANS\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );
  
  //redirect cout back to standard output
  std::cout.rdbuf(coutbuf); //redirect cout to standard output


  return ResultSummary;
}
//****************************************************************************80

void test01 (int k, string kstr, const char* PathToData)

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tries out the HMEANS_01 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "K" + kstr + "test01_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "K" + kstr + "test01_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int it_max;
  int it_num;
  double *point;
  string point_filename = PathToData;
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test the HMEANS_01 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = k;
  it_max = 20;
  seed = 123456789;

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];
  
  hmeans_01 ( dim_num, point_num, cluster_num, it_max, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, 
    cluster_num, point, cluster, cluster_center );

    //print cluster summary, meanwhile obtain total energy, returned from cluster_print_summary
  double ce_total;
  ce_total = cluster_print_summary ( point_num, cluster_num, 
    cluster_population, cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  return;
}
//****************************************************************************80

void test02 (int k, string kstr, const char* PathToData)

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tries out the HMEANS_02 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "K" + kstr + "test02_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "K" + kstr + "test02_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int it_max;
  int it_num;
  double *point;
  string point_filename = PathToData;
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Test the HMEANS_02 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = k;
  it_max = 20;
  seed = 123456789;

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];
  
  hmeans_02 ( dim_num, point_num, cluster_num, it_max, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy, &seed );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, 
    cluster_num, point, cluster, cluster_center );

    //print cluster summary, meanwhile obtain total energy, returned from cluster_print_summary
  double ce_total;
  ce_total = cluster_print_summary ( point_num, cluster_num, 
    cluster_population, cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  return;
}
//****************************************************************************80

vector<double> test03 (int k, string kstr, int it_max, const char* PathToData)

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tries out the KMEANS_01 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = ".K" + kstr + "test03_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = ".K" + kstr + "test03_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  //int it_max;
  int it_num;
  double *point;
  string point_filename = PathToData;
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Test the KMEANS_01 algorithm.\n";
  cout << "  (Applied Statistics Algorithm #58)\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = k;
  //it_max = 20;
  seed = 123456789;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, point, 
    &seed );

  kmeans_01 ( dim_num, point_num, cluster_num, it_max, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of KMEANS_01 iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  //print cluster summary, meanwhile obtain total energy, returned from cluster_print_summary
  double ce_total;
  ce_total = cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  //assemble the vector to return troubleshooting data
  vector<double> retvec(7);
  retvec[0] = 3; //kmeans method 1 = test 3
  retvec[1] = dim_num;
  retvec[2] = point_num;
  retvec[3] = it_max;
  retvec[4] = it_num;
  retvec[5] = k;
  retvec[6] = ce_total;
  
  return retvec;
}
//****************************************************************************80

vector<double> test04 (int k, string kstr, int it_max, const char* PathToData)

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tries out the KMEANS_02 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = ".K" + kstr + "test04_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = ".K" + kstr + "test04_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  //int it_max;
  int it_num;
  double *point;
  string point_filename = PathToData;
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Test the KMEANS_02 algorithm.\n";
  cout << "  (Applied Statistics Algorithm #136)\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = k;
  //it_max = 20;
  seed = 123456789;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_1 ( dim_num, point_num, cluster_num, point );

  kmeans_02 ( dim_num, point_num, cluster_num, it_max, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

    //print cluster summary, meanwhile obtain total energy, returned from cluster_print_summary
  double ce_total;
  ce_total = cluster_print_summary ( point_num, cluster_num, 
    cluster_population, cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  //assemble the vector to return troubleshooting data
  vector<double> retvec(7);
  retvec[0] = 4; //kmeans method 2 = test 4
  retvec[1] = dim_num;
  retvec[2] = point_num;
  retvec[3] = it_max;
  retvec[4] = it_num;
  retvec[5] = k;
  retvec[6] = ce_total;
  
  return retvec;
}
//****************************************************************************80

vector<double> test05 (int k, string kstr, int it_max, const char* PathToData)

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tries out the KMEANS_03 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = ".K" + kstr + "test05_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = ".K" + kstr + "test05_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  //int it_max;
  int it_num;
  double *point;
  string point_filename = PathToData;
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Test the KMEANS_03 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = k;
  //it_max = 20;
  seed = 123456789;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_1 ( dim_num, point_num, cluster_num,
    point );

  kmeans_03 ( dim_num, point_num, cluster_num, it_max, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

    //print cluster summary, meanwhile obtain total energy, returned from cluster_print_summary
  double ce_total;
  ce_total = cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  //assemble the vector to return troubleshooting data
  vector<double> retvec(7);
  retvec[0] = 5; //kmeans method 3 = test 5
  retvec[1] = dim_num;
  retvec[2] = point_num;
  retvec[3] = it_max;
  retvec[4] = it_num;
  retvec[5] = k;
  retvec[6] = ce_total;
  
  return retvec;
}
//****************************************************************************80

void test06 (int k, string kstr, const char* PathToData)

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tries out the HMEANS_01 + KMEANS_01 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "K" + kstr + "test06_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "K" + kstr + "test06_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int it_max;
  int it_max_h;
  int it_max_k;
  int it_num;
  double *point;
  string point_filename = PathToData;
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  Test the HMEANS_01 + KMEANS_01 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = k;
  it_max = 20;
  seed = 123456789;
  it_max_h = 3;
  it_max_k = it_max;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of HMEANS_01 iterations allowed is " << it_max_h << "\n";
  cout << "  Number of KMEANS_01 iterations allowed is " << it_max_k << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );

  hmeans_01 ( dim_num, point_num, cluster_num, it_max_h, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of HMEANS_01 iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  kmeans_01 ( dim_num, point_num, cluster_num, it_max_k, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of KMEANS_01 iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num,
    point, cluster, cluster_center );

    //print cluster summary, meanwhile obtain total energy, returned from cluster_print_summary
  double ce_total;
  ce_total = cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  return;
}
//****************************************************************************80

void test07 (int k, string kstr, const char* PathToData)

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tries out the HMEANS_01 + KMEANS_02 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "K" + kstr + "test07_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "K" + kstr + "test07_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int it_max;
  int it_max_h;
  int it_max_k;
  int it_num;
  double *point;
  string point_filename = PathToData;
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  Test the HMEANS_01 + KMEANS_02 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = k;
  it_max = 20;
  seed = 123456789;

  it_max_h = 3;
  it_max_k = it_max;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of HMEANS_01 iterations allowed is " << it_max_h << "\n";
  cout << "  Number of KMEANS_02 iterations allowed is " << it_max_k << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );

  hmeans_01 ( dim_num, point_num, cluster_num, it_max_h, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of HMEANS_01 iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  kmeans_02 ( dim_num, point_num, cluster_num, it_max_k, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of KMEANS_02 iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

    //print cluster summary, meanwhile obtain total energy, returned from cluster_print_summary
  double ce_total;
  ce_total = cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  return;
}
//****************************************************************************80

void test08 (int k, string kstr, const char* PathToData)

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tries out the HMEANS_01 + KMEANS_03 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "K" + kstr + "test08_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "K" + kstr + "test08_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int it_max;
  int it_max_h;
  int it_max_k;
  int it_num;
  double *point;
  string point_filename = PathToData;
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  Test the HMEANS_01 + KMEANS_03 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = k;
  it_max = 20;
  seed = 123456789;
  it_max_h = 3;
  it_max_k = it_max;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Initialize by using a few steps of HMEANS_02:\n";
  cout << "  Number of HMEANS_01 iterations allowed is " << it_max_h << "\n";
  cout << "  Number of KMEANS_03 iterations allowed is " << it_max_k << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );
//
//  Initialize the clusters.
//
  hmeans_01 ( dim_num, point_num, cluster_num, it_max_h, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of HMEANS_01 iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  kmeans_03 ( dim_num, point_num, cluster_num, it_max_k, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of KMEANS_03 iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

    //print cluster summary, meanwhile obtain total energy, returned from cluster_print_summary
  double ce_total;
  ce_total = cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  return;
}
//****************************************************************************80



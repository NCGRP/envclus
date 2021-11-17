#include <cstdlib>
#include <vector>
#include <iostream>
#include "calcpart.hpp"
using namespace std;

//provides manually specified partitions to partdist to test algorithms
//to compile: g++ testpartdist.cpp partdist.cpp Lmunkres.cpp BipartiteGraph.cpp Hungarian.cpp -o testpartdist
//usage: ./testpartdist

/***************MAIN*****************/

int main( int argc, char* argv[] )
{
	// make an array to get data into vector, for test purposes only, real data will come already in a vector
  	int aints[] = {3,3,3,1,1,1,2,2,2};
  	int bints[] = {0,1,0,2,0,0,1,3,1};
  	vector<int> env (aints, aints + sizeof(aints) / sizeof(int) );
  	vector<int> gen (bints, bints + sizeof(bints) / sizeof(int) );

	
	int cost;
	ofstream outf; //dummy
	vector<pair<int, int> > pqmatrix; 
	
	cost = partdist( env, gen, outf, pqmatrix );
	
	cout << "cost=" << cost << endl;

	//string as = "110101";
	//string bs = "002110";
	//vector<int> a = {1,1,0,1,0,1};
	//vector<int> b = {0,0,2,1,1,0};
	//vector<int> a(as.begin(), as.end()); //add string to vector
	//vector<int> b(bs.begin(), bs.end());




}
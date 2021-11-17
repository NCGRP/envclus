#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <numeric>
#include "calcpart.hpp"
using namespace std;

//to compile: g++ partdist.cpp Lmunkres.cpp BipartiteGraph.cpp Hungarian.cpp -o partdist


/***************GLOBAL VARIABLES*****************/


/***************SUBROUTINES*****************/
//change a partition that starts at >1 to a partition that starts at 0
//then make all class values consecutive, i.e. 0,1,2,3, not 0,1,3,4
vector<int> zeroInd(vector<int> gen)
{
	int i, j, correctionfactor;
	vector<int> genZ;

	int amin = ( *min_element(gen.begin(),gen.end()) );
	correctionfactor = (0 - amin);
	
	//zero index the class values 
	for (i=0;i<gen.size();i++)
	{
		genZ.push_back( gen[i]+correctionfactor );
	}
	
	//test whether values are consecutive
	vector<int> genZSort = genZ;
	sort( genZSort.begin(), genZSort.end() ); //sort copy of vector genZ

	vector<int> genZUniq = genZSort;
	genZUniq.erase( unique( genZUniq.begin(), genZUniq.end() ), genZUniq.end() ); //reduce to unique values
	
	int amax = ( *max_element(genZUniq.begin(), genZUniq.end()) ); //get the maximum value of the zero indexed vector

	vector<int> translateVec = genZ;  //holds zeroed partition with consecutive values
	if ( (amax + 1) != genZUniq.size() ) //if largest value +1 does not equal vector size, they can't be consecutive
	{
		//clear out genZ values from translateVec
		translateVec.clear();
		for (i=0;i<genZ.size();i++)
		{
			for (j=0;j<genZUniq.size();j++)
			{
				if (genZ[i] == genZUniq[j]) 
				{
					translateVec.push_back(j);
					break;
				}
			}
		
		}		
	
	}
	

	return translateVec;
}



/***************MAIN*****************/

int partdist( vector<int> env, vector<int> gen, ofstream& outf, vector<pair<int, int> >& pqmatrix )
{
	int i, j, k;

	//zero index partition coding in env, gen.  make values consecutive
	vector<int> genZ, envZ;
	genZ = zeroInd(gen);
	envZ = zeroInd(env);

	vector<int> a = envZ;
	vector<int> b = genZ;
	
	cout << "env  ";
	for (k=0;k<a.size();k++) cout << a[k];
	cout << endl << "gen  ";
	for (k=0;k<b.size();k++) cout << b[k];
	cout << endl;
		
	//identify maximum number of subsets
	int maxP;
	int amax = ( *max_element(a.begin(),a.end()) ); //get maximum value of vectors a, b.  must * (dereference) to get value instead of iterator
	int bmax = ( *max_element(b.begin(),b.end()) ); 
	if (amax >= bmax) maxP = amax;
	else maxP = bmax;
	
	//add paired max number of subsets to pqmatrix 
	//envmax = amax;
	//genmax = bmax;
	pqmatrix.push_back( make_pair(amax+1, bmax+1) );
	
 
		
	//convert single string partition coding to multiple string binary coding
	vector<vector<int> > A(maxP+1);
	vector<vector<int> > B(maxP+1);
	
	for (j=0; j<(maxP+1); j++)
	{
		for (i=0; i<a.size(); i++)
		{
			if (a[i] == j) 
			{
				A[j].push_back(1);
				//cout << "a[i]=" << a[i] << "   1	" << A[j][i] << "\n";
			}
			else 
			{
				A[j].push_back(0);
				//cout << "a[i]=" << a[i] << "   0	" << A[j][i] << "\n";
			}
				
			if (b[i] == j) 
			{
				B[j].push_back(1);
				//cout << "b[i]=" << b[i] << "   0	" << B[j][i] << "\n";
			}
			else 
			{
				B[j].push_back(0);
				//cout << "b[i]=" << b[i] << "   0	" << B[j][i] << "\n";
			}
		}
	}
	
	
	//get Bbar, the complement (inverse) of B
	vector<vector<int> > Bbar(maxP+1);
	for (i=0;i<Bbar.size();i++)
	{
		for (j=0;j<a.size();j++)
		{
			if (B[i][j] == 0) Bbar[i].push_back(1);
			else Bbar[i].push_back(0);
		}
	}
	
	
	//get the intersection of A and Bbar
/*	vector<vector<vector<int> > > intx(maxP+1, maxP+1); //3D vector
	for (i=0;i<A.size();i++)
	{
		for (j=0;j<Bbar.size();j++)
		{
			for (k=0;k<a.size();k++)
			{
				if (A[i][k] == Bbar[j][k] && A[i][k] == 1) 
				{
					intx[i][j].push_back(1);
				}
				else intx[i][j].push_back(0);
			}
		}
	}
*/
	//as above, but new vector initialization strategy for newer compiler and related logic for
	//adding values to 3d vector
	vector<vector<vector<int> > > intx(maxP+1,vector<vector<int> >(maxP+1,vector<int>(a.size()))); //3D vector
	for (i=0;i<A.size();i++)
	{
		for (j=0;j<Bbar.size();j++)
		{
			for (k=0;k<a.size();k++)
			{
				if (A[i][k] == Bbar[j][k] && A[i][k] == 1) 
				{
					intx[i][j].at(k) = 1;
				}
				else intx[i][j].at(k) = 0;
			}
		}
	}
	
	
	//sum the values in the intersection of A and Bbar
	int m;
	vector<vector<int> > sumx(maxP+1);
	for (i=0;i<A.size();i++)
	{
		for (j=0;j<Bbar.size();j++)
		{
			m=0;
			for (k=0;k<a.size();k++)
			{
				if (intx[i][j][k] == 1) m++;	
			}
			sumx[i].push_back(m);
		}
	}

	

	//Solve the minimum weighted assignment problem using Liu's program which
	//uses the bipartite graph approach to Kuhn-Munkres algorithm
	int cost;
	cost = Lmunkres(sumx, outf);
		


	return cost;
}

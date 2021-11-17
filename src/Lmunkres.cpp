
///////////////////////////////////////////////////////////////////////////////
// File name: Lmunkres.cpp
// This file calls the components and algorithms to complete Hungarian task.
// Lantao Liu, Nov 1, 2009
//
// Modified by Pat Reeves, Apr 1, 2013
// To compile:  g++ Lmunkres.cpp BipartiteGraph.cpp Hungarian.cpp -o ha
//
/////////////////////////////////////////////////////////////////////////////// 

#include "Hungarian.h"
#include "BipartiteGraph.h"


int Lmunkres(vector<vector<int> > sumx, ofstream& outf)
{
  //define a matrix;
  Matrix m;

	//load sumx onto 2D vector x
	unsigned int i, j;
	int v;
	std::vector< std::vector<int> > x ;
	x = sumx;

	//print out
	cout << "x:\n";
	for (i=0; i<x.size(); i++)
	{
 	   for (j=0; j<x[0].size(); j++)
	    {
			cout << x[i][j] << "  ";
   		}
   		cout << std::endl;
    }


	//get the dimensions of the 2D vector
	unsigned int nrows = x.size();
	unsigned int ncols = x[0].size();
	
	//transpose 2D vector so nrows <= ncols
	std::vector< std::vector<int> > xnew(ncols);
	if (ncols < nrows) 
	{
		cout << "WARNING: I AM TWISTING THE VECTOR\n";
		for (i=0; i<ncols; i++)
	    {
 	   	for (j=0; j<nrows; j++)
	    	{
				v = x[j][i];
				xnew[i].push_back(v);
   			}
    	}
		x.swap(xnew);

	//print out
	cout << "new x:\n";
	for (i=0; i<x.size(); i++)
		{
 		   for (j=0; j<x[0].size(); j++)
	 	   {
				cout << x[i][j] << "  ";
   			}
   			cout << std::endl;
   		}
	}


	//get the -new- dimensions of the 2D vector (after transpose)
	nrows = x.size();
	ncols = x[0].size();

	//convert the data table for calculation of minimum weighted assignment
	//for each cell, use conversion formula y = max(x) - x
	//get max value in x
	int xmaxtemp;
	int xmax = 0;
	for (i=0;i<x.size();i++)
	{	
		xmaxtemp = ( *max_element(x[i].begin(),x[i].end()) );
		if (xmaxtemp > xmax) xmax = xmaxtemp;
	}
	
	cout << "xmax= " << xmax << "\n";

	//add values to y
	std::vector< std::vector<int> > y(nrows);
	for (i=0;i<nrows;i++)
	{
		for (j=0;j<ncols;j++)
		{
			v = xmax - x[i][j];
			y[i].push_back(v);	
		}
	}
	
	//print out
	cout << "y:\n";
	for (i=0; i<y.size(); i++)
	{
 	   for (j=0; j<y[0].size(); j++)
	    {
			cout << y[i][j] << "  ";
   		}
   		cout << std::endl;
    }
	
	
	
	//resize Matrix to match vector y
	m.resize(nrows);
	for(i=0; i<nrows; i++)
      m[i].resize(ncols);
	
	//add values from vector to matrix
    for (i=0; i<nrows; i++)
    {
    	for (j=0; j<ncols; j++)
    	{
    		v = y[i][j];
    		m[i][j].SetWeight(v);
   		}
    }


  //define a bipartite graph
  BipartiteGraph bg(m);
    	


  //run Hungarian method, return Kuhn-Munkres Hungarian Perfect Matching table as referenced variable matchtable
  //structure of matchtable is pair<int, int>. like: (0,0) (1,2) (2,1)  shows positions in y to be added to get cost
  vector<EID> matchtable;
  Hungarian h(bg);
  h.HungarianAlgo(matchtable);

	
	//compute minimum weighted assignment problem cost (= partition distance) from original vector x
	int cost = 0;		
	for (i=0;i<y.size();i++)
	{
		v = x[matchtable[i].first][matchtable[i].second];
		cost = cost + v; 
	}
	
  return cost;
}


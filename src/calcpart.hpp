#include <vector>
#include <fstream>
#include <iostream>
#include <utility>






/***************STRUCTS*****************/
struct sort_pred {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
        return left.second < right.second;
    }
};



/***************FUNCTIONS*****************/
int partdist( std::vector<int> env, std::vector<int> gen, std::ofstream& outf, std::vector<std::pair<int, int> >& pqmatrix );
int Lmunkres( std::vector<std::vector<int> > sumx, std::ofstream& outf );


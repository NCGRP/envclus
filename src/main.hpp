#include <iostream>

/***************FUNCTIONS*****************/
int envclus (const int NumSamples, const int k, const char* GridFile, char* EnvLayersDir, const char* OutFile, std::ofstream& outf );
double calcpart (const char* PathToEnvPartition, const char* PathToGenPartition, const char* PathToCalcPartOutput, std::ofstream& outf );



// compile: g++ -I/usr/local/ generate_FvsD.cpp -o fvsd
#include <string>
// output
#include "stdio.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
// stl
#include <vector>
#include <set>
#include <map>

#include "database_utility.h"
//#define VERBOSE

using namespace SQL_interface;

int main(int narg, char* argv[]){
  //cout << narg << flush;
  if (narg < 3) {
    cout << "usage: generate_FvsD <db file> <cluster size> \n" << flush; 
    exit(1);
  }
  // parse command line
  string inname(argv[1]);
  int n_neighbs = atoi(argv[2]);
  int ncols = offset_+3*n_neighbs;

  // read in
  vector < MatrixXd > clusters;
  vector < vector < double > > forces;
  string line, buf;
  ifstream infile(inname.c_str());
  int count = 0;
  //cout << "# read " << flush;
  while(getline(infile, line)){ 
#ifdef VERBOSE
    cout << "line >>> " << line << "\n" << flush;
#endif
    stringstream ss(line); 
    vector<string> tokens; 
    while (ss >> buf){ tokens.push_back(buf);}
    vector<double> f; 
    f.push_back(atof(tokens[2].c_str()));
    f.push_back(atof(tokens[3].c_str()));
    f.push_back(atof(tokens[4].c_str()));
    forces.push_back(f);
    MatrixXd C(n_neighbs,3);
    int i = 0;
    for (int j = offset_; j< ncols; j+=3){
      C(i,0) = atof(tokens[j].c_str());
      C(i,1) = atof(tokens[j+1].c_str());
      C(i,2) = atof(tokens[j+2].c_str());
      i++;
    }
    clusters.push_back(C);
    //if (count++ % 1000 == 0) cout << " " << count << flush;
  }
  //cout << "\n";
  infile.close();

  // write pairwise distances
  int nClusters = clusters.size();
  cout << setprecision(16);
  VectorXd weights;
  ofstream outfile("fd.dat");
  int k = 0;
  for (int i = 0; i < nClusters; ++i) {
    for (int j= i+1; j < nClusters; ++j) {
      double d = pwrmsd(clusters[i],clusters[j]);
//    double d =  wrmsd(clusters[i],clusters[j]);
//    double d =   rmsd(clusters[i],clusters[j]);
      vector<double> & fi = forces[i];
      vector<double> & fj = forces[j];
      double df0 = fi[0]-fj[0];
      double df1 = fi[1]-fj[1];
      double df2 = fi[2]-fj[2];
      double df = sqrt(df0*df0+df1*df1+df2*df2);
      //cout << i << " " << j << " " ;
      //cout << d << " " << df << "\n"; 
      outfile << d << " " << df << "\n"; 
      if (++k % 1000 ==0) { cout << k << " " << flush; }
#if 0
#ifdef PRINT_JI
      cout << j << " " << i << " " << d << " "; 
      cout << fi[0] << " " << fj[0] << " ";
      cout << fi[1] << " " << fj[1] << " ";
      cout << fi[2] << " " << fj[2] << "\n";
#endif
#endif
    }
  }
  cout << "\n";
}

// compile: g++ -lsqlite3 -I/usr/local/ database_builder.cpp -o db_builder
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

#define MAX_DIST 100
//#define VERBOSE
#include "database_utility.h"
using namespace SQL_interface;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// To do:
// * cluster ids 
// * parallel
// * permutation check
// * optimally winnow + acct for permutations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
// *********************************************************************
// histogram  
// *********************************************************************
void histogram (double max_d, map < IJ, double > & ds, string name ) {
  const int nbins = 10;
  double bin_size = max_d/nbins;
  int counts[nbins];
  double partitions[nbins];
  for (int i = 0; i < nbins; i++) {
    counts[i] = 0;
    partitions[i] = (i+1)*bin_size;
  }
  int nd = ds.size();
  for( map< IJ, double >::const_iterator itr = ds.begin();
      itr != ds.end(); itr++) {
    double d = itr->second;
    for (int i = 0; i < nbins; i++) {
      if (d < partitions[i]) { counts[i]++; break; }
    }
  }
  ofstream hfile(name.c_str());
  hfile << "# histogram:\n";
  for (int i = 0; i < nbins; i++) {
    hfile << i << " " << partitions[i] << " " << counts[i] << "\n";
    printf("%2d %7.4f %10d ",i,partitions[i],counts[i]);
    for (int j = 0; j < int(4*nbins*counts[i]/nd); j++) cout << "*";
    cout << "\n";
  }
  hfile.close();
}

//=====================================================================+
// MAIN
//=====================================================================+
int main(int narg, char* argv[]){
  cout << setprecision(fp_precision_);
  if (narg < 4) {
    cout << "usage: db_builder <input DFT file> <output database name> <cluster size> [number of reference clusters]\n" << flush; 
    cout << "format: center_atom_type neighbor_types fx fy fz n1x n1y n1z ...\n" << flush; 
    exit(1);
  }
  // parse command line
  string inname(argv[1]);
  string outname(argv[2]); 
  int n_neighbs = atoi(argv[3]);
  cout << "read: " << inname << "\n";
  cout << "write:" << outname << " clusters of size " << n_neighbs; 
// maxClose_
  bool writeDistances = false;
  int n_ref_clusters = 0;
  if (narg > 4) {
    n_ref_clusters = atoi(argv[4]);
    writeDistances = true;
    cout << " and distances (" << n_ref_clusters << " refs)";
  }
  cout << "\n";
  bool writeClusters = true;
  vector < MatrixXd > clusters;
  vector < Vector3d > forces;
  int ncols = offset_+3*n_neighbs;
//.....................................................................
// workspace
//.....................................................................
  string cmd;
  sqlite3 *db;
  sqlite3_stmt *stmt;
//.....................................................................
// read & write clusterdata
//.....................................................................
  ifstream infile(inname.c_str()); 
  int retval = sqlite3_open(outname.c_str(), &db);
  // setup
  cmd = "DROP TABLE IF EXISTS clusters";// removes existing table
  query(db,cmd,stmt);
  cmd = "CREATE TABLE clusters (cid INTEGER, ctype INTEGER, cluster_type INTEGER, fx double, fy double, fz double ";
  for (int j = 1; j<n_neighbs+1; j++){
    cmd += ", c"+SQL_interface::to_string(j)+"x TEXT, "+"c"+SQL_interface::to_string(j)+"y TEXT, "+"c"+SQL_interface::to_string(j)+"z TEXT"+", c"+SQL_interface::to_string(j)+"r TEXT";}
//  cmd += ", c"+SQL_interface::to_string(j)+"x double, "+"c"+SQL_interface::to_string(j)+"y double, "+"c"+SQL_interface::to_string(j)+"z double"+", c"+SQL_interface::to_string(j)+"r double";}
  cmd += ");";
  query(db,cmd,stmt);
  int count =0;
  string line, buf;
  cout << "clusters: " << flush;
  while(getline(infile, line)){ 
#ifdef VERBOSE
    cout << "line >>> " << line << "\n" << flush;
#endif
    count++; // starts with 1
    if (count % 100 == 0) cout << "." << flush;
    stringstream ss(line); 
    vector<string> tokens; 
    while (ss >> buf){ tokens.push_back(buf);}
    string cid = SQL_interface::to_string(count); 
    string ctype        = tokens[0]; 
    string cluster_type = tokens[1]; 
    string fx           = tokens[2];
    string fy           = tokens[3];
    string fz           = tokens[4];
    cmd = "INSERT INTO clusters VALUES ( ";
    cmd += cid + ", "+ctype+", "+cluster_type+", "+fx+", "+fy+", "+fz;
    MatrixXd C(n_neighbs,3);
    int i = 0;
    for (int j = offset_; j< ncols; j+=3){ 
      cmd += ", "+tokens[j]+", "+tokens[j+1]+", "+tokens[j+2];
      double r = mag(atof(tokens[j].c_str()),
                     atof(tokens[j+1].c_str()),
                     atof(tokens[j+2].c_str()));
      cmd += ", "+SQL_interface::to_string(1./r);
      if (writeDistances) {
        C(i,0) = atof(tokens[j].c_str());
        C(i,1) = atof(tokens[j+1].c_str());
        C(i,2) = atof(tokens[j+2].c_str());
        i++;
      }
    }
    if (writeDistances) {
      clusters.push_back(C);
    }
    cmd += ");";
    query(db,cmd,stmt);
  }
  cout << " wrote:"<< count << "\n";
  infile.close();
  int n_clusters = count;
  if (writeClusters) {
    string outname = inname + ".dat";
    cout << "writing " << outname << "\n";
    ofstream outfile(outname.c_str()); 
    for (int i = 0; i < n_clusters; i++) {
      MatrixXd & C = clusters[i];
      for (int j = 0; j<n_neighbs; j++){
        for (int k = 0; k<3; k++){
          outfile << C(j,k) << " ";
        }
      }
      outfile << "\n"; 
    }
  }
  if ( ! writeDistances ) exit(0);

//.....................................................................
//   calculate distances
//.....................................................................
  int maxDistances = min(n_clusters-1,MAX_DIST);
  int n_pairs = n_clusters*(n_clusters-1)/2;
  int heartbeat = max(n_pairs/10, 1);
  cout << "calculating " << n_pairs << " distances: " << flush;
  map < IJ, double > distances;
  int p = 0;
  double max_d = 0;
  double min_d = 1.e20;
  for (int i = 0; i < n_clusters; i++) {
    for (int j = i+1; j < n_clusters; j++) {
      IJ ij(i,j);
      double d = SQL_interface::d(clusters[i],clusters[j]);
      min_d = min(d,min_d);
      max_d = max(d,max_d);
      distances[ij] = d;
      if (p % heartbeat == 0) { cout << "." << flush; }
      p++;
    }
  }
  //cout << "calculated " << n_pairs << " unique distances\n" << flush;
  printf(" %7.4f:%7.4f",min_d,max_d);
  cout << "\n" << flush;

  histogram(max_d,distances,"pair_distance_histogram.dat");

//.....................................................................
//   write distance database
//.....................................................................
  cmd = "DROP TABLE IF EXISTS cluster_distances";// removes existing table
#ifdef DEBUG
  cout << cmd << "\n" << flush;
#endif
  query(db,cmd,stmt);
  cmd = "CREATE TABLE cluster_distances (cid INTEGER ";
#ifdef DEBUG
  cout << cmd << "\n" << flush;
#endif
  for (int j = 1; j<maxDistances+1; j++){
    cmd += ", c"+SQL_interface::to_string(j)+" INTEGER, "+"d"+SQL_interface::to_string(j)+" TEXT";}
//  cmd += ", c"+SQL_interface::to_string(j)+" INTEGER, "+"d"+SQL_interface::to_string(j)+" double";}
  cmd += ");";
#ifdef DEBUG
  cout << cmd << "\n" << flush;
#endif
  query(db,cmd,stmt);
  cout << "distances: " ;
  heartbeat = max(n_clusters/10,1);
  p = 0;
  for (int i = 0; i < n_clusters; i++) {
    vector < ID >  row(n_clusters-1);
    int k = 0;
    for (int j = 0; j < n_clusters; j++) { 
      if (j==i) continue;
      IJ ij(min(i,j),max(i,j));
      double d = distances[ij]; // 0 bases
      ID dj(j+1,d); // 1 based
      row[k++] =  dj;
    }
    sort(row.begin(),row.end(),comp_id);
#ifdef VERBOSE
    cout << i << ": min " << row.front().second << " max " << row.back().second << "\n" << flush;
#endif
#ifdef DEBUG
    cout << "  ";
    for (int j = 0; j < maxDistances; j++) { 
      cout << row[j].first << ":" << row[j].second << " ";
    }
    cout << "\n" << flush;
#endif
    cmd = "INSERT INTO cluster_distances VALUES( " + SQL_interface::to_string(i+1); //1-based
    for (int j = 0; j<maxDistances; j++){
      ID dj = row[j];
      int id =  dj.first;
      double d =  dj.second;
      cmd += ", "+SQL_interface::to_string(id)+", "+SQL_interface::to_string(d)+" ";
    }
    cmd += ");";
    query(db,cmd,stmt);
    if (p % heartbeat == 0) cout << "." << flush;
    p++;
  }
  cout << " wrote "<< n_clusters << "x" << maxDistances << flush;
  // reference clusters
  cmd = "DROP TABLE IF EXISTS reference_clusters";// removes existing table
  query(db,cmd,stmt);
  cmd = "CREATE TABLE reference_clusters (cid INTEGER); ";
  query(db,cmd,stmt);
  for (int i = 0; i < n_ref_clusters; i++) {
    cmd = "INSERT INTO reference_clusters VALUES("+SQL_interface::to_string(i)+");";
    query(db,cmd,stmt);
  }
  cout << " and " << n_ref_clusters << " references\n" << flush;

  retval = sqlite3_close(db);
}

#ifndef LMP_DATABASE_UTILITY_H
#define LMP_DATABASE_UTILITY_H
#include <utility>
#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <set>
#include <map>
#include <sqlite3.h>
#include <Eigen/Core>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
using std::fstream;
using std::vector;
using std::pair;
using std::string;

// to do
// remove sql ret values externally and make all checks internal
// history / cache 1st read into cache then search cache distance and clusters?
// limit memory distances and resets

#define VERBOSE
//#define PRINT_DEBUG
#define FL __FILE__
#define LN __LINE__

namespace SQL_interface {
//======================================================================
static inline string SQL_error(int i) {
  switch (i) {
  case SQLITE_OK        : return "Successful result";
  case SQLITE_ERROR     : return "SQL error or missing database";
  case SQLITE_INTERNAL  : return "Internal logic error in SQLite";
  case SQLITE_PERM      : return "Access permission denied";
  case SQLITE_ABORT     : return "Callback routine requested an abort";
  case SQLITE_BUSY      : return "The database file is locked";
  case SQLITE_LOCKED    : return "A table in the database is locked";
  case SQLITE_NOMEM     : return "A malloc() failed";
  case SQLITE_READONLY  : return "Attempt to write a readonly database";
  case SQLITE_INTERRUPT : return "Operation terminated by sqlite3_interrupt";
  case SQLITE_IOERR     : return "Some kind of disk I/O error occurred";
  case SQLITE_CORRUPT   : return "The database disk image is malformed";
  case SQLITE_NOTFOUND  : return "Unknown opcode in sqlite3_file_control";
  case SQLITE_FULL      : return "Insertion failed because database is full";
  case SQLITE_CANTOPEN  : return "Unable to open the database file";
  case SQLITE_PROTOCOL  : return "Database lock protocol error";
  case SQLITE_EMPTY     : return "Database is empty";
  case SQLITE_SCHEMA    : return "The database schema changed";
  case SQLITE_TOOBIG    : return "String or BLOB exceeds size limit";
  case SQLITE_CONSTRAINT: return "Abort due to constraint violation";
  case SQLITE_MISMATCH  : return "Data type mismatch";
  case SQLITE_MISUSE    : return "Library used incorrectly";
  case SQLITE_NOLFS     : return "Uses OS features not supported on host";
  case SQLITE_AUTH      : return "Authorization denied";
  case SQLITE_FORMAT    : return "Auxiliary database format error";
  case SQLITE_RANGE     : return "2nd parameter to sqlite3_bind out of range";
  case SQLITE_NOTADB    : return "File opened that is not a database file";
  case SQLITE_ROW       : return "sqlite3_step() has another row ready";
  case SQLITE_DONE      : return "sqlite3_step() has finished executing";
  }
  return "unknown";
}
//======================================================================
#ifdef VERBOSE
#define DB_ERROR_CHECK(S,F,L) if( S != SQLITE_OK && S != SQLITE_ROW && S != SQLITE_DONE) { cout << "SQL error code: " << SQL_error(S) << " at " << F << ":" << L << "\n" << flush;  }
#else
// fatal
#define DB_ERROR_CHECK(S,F,L) if( S != SQLITE_OK && S != SQLITE_ROW && S != SQLITE_DONE) { cout << "SQL error code: " << SQL_error(S) << " at " << F << ":" << L << "\n" << flush; exit(S); }
// non-fatal
//#define DB_ERROR_CHECK(S,F,L) if( S != SQLITE_OK && S != SQLITE_ROW && S != SQLITE_DONE) { cout << "SQL error code: " << SQL_error(S) << " at " << F << ":" << L << "\n" << flush;}
#endif
// constants ============================================================
const double big_ = 1.e20; 
const int offset_ = 5; // for sql builder
const int dxSize_ = 4; 
const int fp_precision_ = 16;
const double permTol_ = 1.e-8;
const string clusterTableName_ = "clusters";
const string distanceTableName_ = "cluster_distances";
const string referenceTableName_ = "reference_clusters";
// typedefs =============================================================
struct index_type_distance {int idx; int type;  double r;} ;
typedef index_type_distance  INDEX;
typedef vector< INDEX > INDICES;
struct triple {double x1; double x2;  double x3;} ;
typedef vector< triple > DATA;
typedef pair<int,int> IJ;
typedef pair<int,double> ID;
//========================================================================
static inline double mag (double x, double y, double z) { return sqrt(x*x+y*y+z*z);}
static inline double mag (double * v) {return sqrt(v[0]*v[0]+v[1]*v[1]*v[2]*v[2]);}
//========================================================================
static inline bool comp_index (INDEX i,INDEX j) { return (i.r<j.r); }
static inline bool comp_id (ID i,ID j) { return (i.second<j.second); } 
//========================================================================
inline string to_string(int i) { stringstream ss; ss << i; return ss.str(); }
inline string to_string(double i) { stringstream ss; ss << setprecision(fp_precision_); ss << i; return ss.str(); }
//========================================================================
inline void print_config(const MatrixXd & C) { 
  for (int jj=0;jj<C.rows();jj++) {
    //printf("%20.17g %20.17g %20.17g\n", C(jj,0), C(jj,1), C(jj,2));
    printf("%8.4f %8.4f %8.4f\n", C(jj,0), C(jj,1), C(jj,2));
  }
}
//========================================================================
inline string type_string(vector<int> types) { 
  string typeString;
  for (unsigned int i=0; i<types.size(); i++) {
    int itype = types[i];
    if (itype > 9) {cout << "!!! LAMMPS type " << itype <<" too large to handle\n" << flush; throw;}
    typeString += to_string(itype);
  }
  return typeString;
}
//========================================================================
inline string type_string(INDICES indices, unsigned int len = 0) { 
  string typeString;
  if (len == 0) { len = indices.size(); }
  for (unsigned int i=0; i< len; i++) {
    int itype = indices[i].type;
    if (itype > 9) {cout << "!!! LAMMPS type " << itype <<" too large to handle\n" << flush; throw;}
    typeString += to_string(itype);
  }
  return typeString;
}
//========================================================================
// Gaussian overlap
//========================================================================
static inline double g(double r2) {
  return exp(-r2);
}
static inline double ogto(const MatrixXd & C1, const MatrixXd & C2, 
  Matrix3d & R, double sigma = 1.0)
{
  int n1 = C1.rows();
  int n2 = C2.rows();
  double d = 0.;
  // C1.C1
  for (int i1=0; i1 < n1; ++i1) {
    for (int i2=0; i2 < n1; ++i2) {
      Vector3d dx = C1.row(i1)-C1.row(i2);
      d += g(dx.dot(dx));
    }
  }
  // C2.C2
  for (int i1=0; i1 < n2; ++i1) {
    for (int i2=0; i2 < n2; ++i2) {
      Vector3d dx = C2.row(i1)-C2.row(i2);
      d += g(dx.dot(dx));
    }
  }
  // C1.C2
  for (int i1=0; i1 < n1; ++i1) {
    for (int i2=0; i2 < n2; ++i2) {
      Vector3d dx = C1.row(i1)-C2.row(i2);
      d += g(dx.dot(dx));
    }
  }
  // NOTE R???
  return d;
}
//========================================================================
// SVD = C1.C2^T 
//========================================================================
static inline double rmsd(const MatrixXd & C1, const MatrixXd & C2, 
  Matrix3d & R, double n = 0)
{
  Matrix3d Covar = C1.transpose()*C2;
  JacobiSVD<MatrixXd> svd(Covar, ComputeFullU|ComputeFullV);
  MatrixXd  U = svd.matrixU();
  MatrixXd  V = svd.matrixV();
  VectorXd  S = svd.singularValues();
  if (n == 0) n = S.size();
  // flip if a reflection
  R = V*U.transpose();
  if (R.determinant() < 0){ S[2] = -S[2]; }
  double rmsd_sq = (C1.squaredNorm()+C2.squaredNorm()-2*S.sum());
  rmsd_sq = fabs(rmsd_sq); // handle small negative numbers
#ifdef PRINT_DEBUG
  if (fabs(R(0,0)-1) < 1.e-7) {
  cout << "C1\n" << C1 << "\n";
  cout << "C2\n" << C2 << "\n";
  cout << "R\n" << R << "\n";
  cout << "Covar\n" << Covar << "\n";
  cout << "S\n" << S << "\n";
  cout << "rmsd_sq " << rmsd_sq 
       << " ||C1||^2 " << C1.squaredNorm()
       << " ||C2||^2 " << C2.squaredNorm()
       << " 2 trS " << 2*S.sum()
       <<  "\n";
  }
#endif
  return sqrt(rmsd_sq/n); // normalize 
}
//========================================================================
static inline double rmsd(const MatrixXd & C1, const MatrixXd & C2)
{
  Matrix3d R;
  return rmsd(C1,C2,R);
}
//========================================================================
// SVD = C1.C2^T 
//========================================================================
static inline double rmsd(const MatrixXd & C1, const MatrixXd & C2, Matrix3d & R)
{
  Matrix3d Covar = C1.transpose()*C2;
  JacobiSVD<MatrixXd> svd(Covar, ComputeFullU|ComputeFullV);
  MatrixXd  U = svd.matrixU();
  MatrixXd  V = svd.matrixV();
  VectorXd  S = svd.singularValues();
  int n = S.size();
  // flip if a reflection
  R = V*U.transpose();
  if (R.determinant() < 0){ S[2] = -S[2]; }
  double rmsd_sq = (C1.squaredNorm()+C2.squaredNorm()-2*S.sum());
  rmsd_sq = fabs(rmsd_sq); // handle small negative numbers
  return sqrt(rmsd_sq/n); // normalize 
}
//========================================================================
// W = diag (1/r_1+1/r_2)
//========================================================================
static inline void weight_vector(const MatrixXd & C1,  const MatrixXd & C2, 
  VectorXd & W) 
{
  VectorXd W1 = C1.rowwise().norm();
  VectorXd W2 = C2.rowwise().norm();
  W = W1.array().pow(-1) + W2.array().pow(-1);
}
//========================================================================
// SVD = C1.W.C2^T = (sqrtW.C1).(sqrtW.C2)^T
//========================================================================
static inline double wrmsd(const MatrixXd & C1, const MatrixXd & W, const MatrixXd & C2,Matrix3d & R)
{
  MatrixXd sqrtW = W.array().sqrt();
  MatrixXd wC1 = sqrtW*C1;
  MatrixXd wC2 = sqrtW*C2;
  double w = W.sum();
  return rmsd(wC1,wC2,R,w);
}
//========================================================================
static inline double wrmsd(const MatrixXd & C1, const MatrixXd & C2, Matrix3d & R)
{
  VectorXd W;
  weight_vector(C1, C2,  W);
  return wrmsd(C1,W.asDiagonal(),C2,R);
}
//========================================================================
static inline double wrmsd(const MatrixXd & C1, const MatrixXd & C2)
{
  VectorXd W;
  weight_vector(C1, C2,  W);
  Matrix3d R;
  return wrmsd(C1,W.asDiagonal(),C2,R);
}
//========================================================================
// for a secondary sort to speed up permutation
//========================================================================
static inline int equivalence_sets(const set<double> & rs, double tol, set<set<int> > & esets) {
  esets.clear();
  double r = 0;
  set<int> aset;
  int nsets = 0, i = 0;
  for (set<double>::iterator itr = rs.begin(); itr != rs.end(); itr++) {
    if (abs(r-*itr) < tol) { aset.insert(i); }
    else { esets.insert(aset); aset.clear(); aset.insert(i); nsets++; }
    i++;
  } 
  esets.insert(aset); nsets++;
  return nsets++;
}
//========================================================================
// SVD = C1.W.P.C2^T = (sqrtW.C1).(sqrtW.P.C2)^T
//========================================================================
static inline double pwrmsd(const MatrixXd & C1, const MatrixXd & C2, Matrix3d & R,
  int maxPermutations = 0)
{
  VectorXd W;
  weight_vector(C1, C2,  W);
  int n = W.size();
  MatrixXd sqrtW = W.array().sqrt();
#ifdef DEBUG_PWRMSD 
  cout << "C1\n" << C1 << "\n";
  cout << "C2\n" << C2 << "\n";
  cout << "W\n" << W << "\n";
  cout << "sqrtW\n" << sqrtW << "\n";
  cout << flush;
#endif
  MatrixXd wC1 = (sqrtW.asDiagonal())*C1; 
  MatrixXd wC2 = (sqrtW.asDiagonal())*C2;
  double w = W.sum();
  vector<int> ids(n);
  double min_d = big_; //  rmsd(wC1,wC2,R,w);
  for (unsigned int i = 0; i < ids.size(); i++) { ids[i] = i; }
  int p = 0;
  do {
    if (p > 0) {
      MatrixXd sqrtW_P = (sqrtW.asDiagonal());
      for (unsigned int i = 0; i < ids.size(); i++) { 
       sqrtW_P(i,ids[i]) = sqrtW_P(i,i);
       sqrtW_P(i,i) = 0; 
      }
      wC2 = sqrtW_P*C2;
    }
    double d = rmsd(wC1,wC2,R,w);
    if (d < min_d) {
      min_d = d;
    } 
#ifdef DEBUG_PWRMSD
    std::cout << p << ": ";
    for (int i = 0; i < n; i++) { std::cout << ids[i]  << " "; } 
    std::cout << "d " << d << "\n";
#endif
    p++;
    if (maxPermutations > 0 && p > maxPermutations) break;
  } while ( std::next_permutation(ids.begin(),ids.end()) );
#ifdef DEBUG_PWRMSD
  std::cout << "------------- " << min_d << " ---------------\n";
  throw;
#endif
  return min_d;
}
//========================================================================
static inline double pwrmsd(const MatrixXd & C1, const MatrixXd & C2, 
  int maxPermutations = 0)
{
  Matrix3d R;
  return pwrmsd(C1,C2,R,maxPermutations);
}
//========================================================================
// generic distance interface 
//========================================================================
static inline double d(const MatrixXd & C1, const MatrixXd & C2, Matrix3d & R)
{
  return wrmsd(C1,C2,R);
}
//========================================================================
static inline double d(const MatrixXd & C1, const MatrixXd & C2)
{
  Matrix3d R;
  return wrmsd(C1,C2,R);
}
//======================================================================
// basic SQL
//======================================================================
static inline int query(sqlite3 * db, string cmd, sqlite3_stmt * & stmt) {
#ifdef VERBOSE
  cout << cmd << "\n";
#endif
  int retval = sqlite3_prepare_v2(db,cmd.c_str(), -1, &stmt,0);
  DB_ERROR_CHECK(retval,FL,LN)
  return retval = sqlite3_step(stmt);
}
//======================================================================
static inline int finalize(sqlite3 * db, sqlite3_stmt * & stmt) {
  return sqlite3_finalize(stmt);
}
//======================================================================
// SQLITE wrappers
//======================================================================
class SQL_Interface 
{
  public:
    SQL_Interface(sqlite3 *db, string tablename)
      :name_("unknown"),tablename_(tablename),db_(db) {};//assumes open
    SQL_Interface(string name,string tablename)
      :name_(name),tablename_(tablename) { 
      int retval = sqlite3_open(name.c_str(), &db_);
      if ( retval )  {
        cout << "can't open database " << name << "\n" << flush;
        exit(1);
      }
    };
    virtual ~SQL_Interface(){ if (name_ != "unknown") sqlite3_close(db_); }
    // basic
    inline sqlite3 * database() { return db_;}
    inline string name() { return name_;}
    inline int size(string tablename){
      string cmd = "SELECT Count(*) as tablesize FROM "+tablename+";";  
      query(cmd);
      int size  = col_int(0);
      end_query();
      return size;
    }
    inline int size(){ return size(tablename_); } 
    inline bool has_table(string name) {
      string cmd = "SELECT count(*) FROM sqlite_master WHERE type='table' AND name='"+name+"';";
      int retval = query(cmd);
      DB_ERROR_CHECK(retval,FL,LN)
      bool has = col_int(0) > 0;
      end_query();
      return has;
    }
    inline int number_of_columns()
    {
      string cmd = "SELECT * FROM "+tablename_+" where cid = 1;";
      int retval = query(cmd);
      DB_ERROR_CHECK(retval,FL,LN)
      int ncols = sqlite3_column_count(stmt_);
      end_query();
      return ncols;
    }
    inline int query(const string cmd)
    {
#ifdef VERBOSE
      cout << cmd << "\n";
#endif
      int retval = prep(cmd);
      DB_ERROR_CHECK(retval,FL,LN)
      step();
      return retval; // NOTE retval or done bool?
    }
    int prep(string cmd) {
#ifdef VERBOSE
      cout << cmd << "\n";
#endif
      return sqlite3_prepare(db_, cmd.c_str(),-1, &stmt_, 0);
    }
    bool step() { 
      int retval = sqlite3_step(stmt_); 
      DB_ERROR_CHECK(retval,FL,LN)
      return (! ( retval == SQLITE_DONE ) );
    }
    int end_query() { return sqlite3_finalize(stmt_); }
  protected:
    inline int    col_int(const int idx) {
      return sqlite3_column_int(stmt_,idx); }
    inline double col_dbl(const int idx) {
      const unsigned char * text = sqlite3_column_text(stmt_,idx); 
      std::stringstream ss;
      ss << text;
      return std::atof(ss.str().c_str()); }
//    return sqlite3_column_double(stmt_,idx); }
    inline string col_string(const int idx) {
      //double num = std::atof(sqlite3_column_text(stmt_,idx));
      double num = sqlite3_column_double(stmt_,idx);
      const string text = to_string(num);
      //const string text = sqlite3_column_text(stmt_,idx);
      return text; }
    SQL_Interface();
    string name_;
    string tablename_;
    //class sqlite3 *db_;
    //class sqlite3_stmt *stmt_; 
    struct sqlite3 *db_;
    struct sqlite3_stmt *stmt_; 

};// class
class SQL_Cluster_Interface : public SQL_Interface
{
  public:
    static const int ndx_ = dxSize_;
    static const int configOffset_ = offset_+1;
    SQL_Cluster_Interface(string name):SQL_Interface(name,clusterTableName_) { 
      nn_ = (number_of_columns()-configOffset_)/ndx_;
    }
    virtual ~SQL_Cluster_Interface(){};
    // basic
    inline int number_of_clusters() { return size(); }
    inline int configuration_size() { return nn_; }
    int prep_same_type(int ctype, string cluster_type) {
      string cmd = "SELECT * FROM "+clusterTableName_+" WHERE ctype = "+to_string(ctype)+" AND cluster_type = " + cluster_type + " ;";
      int retval = sqlite3_prepare(db_, cmd.c_str(),-1, &stmt_, 0);
      DB_ERROR_CHECK(retval,FL,LN)
      return 0; // HACK MAKE RETURN COUNT
    }
    // specific
    int id()           { return col_int(0); }
    int type()         { return col_int(1); } 
    string cluster_type() { return col_string(2);}
    void cluster_force(double * f) {
      f[0] = col_dbl(3);
      f[1] = col_dbl(4);
      f[2] = col_dbl(5);
    };
    double dx(int i,double * dx) { // relative position
      int p = configOffset_+ndx_*i;  
      dx[0] = col_dbl(p++);
      dx[1] = col_dbl(p++);
      dx[2] = col_dbl(p++);
      return col_dbl(p++);
    };
    double r(int i) { // radial distance
      int p = configOffset_+ndx_*(i+1)-1;  
      return col_dbl(p);
    };
    int configuration(MatrixXd & C) { 
      int p = configOffset_;
      for (int i =0; i<nn_; i++){
        C(i,0) = col_dbl(p);
        C(i,1) = col_dbl(p+1);
        C(i,2) = col_dbl(p+2);
        p += ndx_;
      }
      return nn_;
    };
  private:
    SQL_Cluster_Interface();
    int nn_;
};// class
class SQL_Distance_Interface : public SQL_Interface
{
  public:
    SQL_Distance_Interface(sqlite3* db)
      :SQL_Interface(db,distanceTableName_),nDist_(0),nRef_(0){
      nDist_ = (number_of_columns()-1)/2; 
      nRef_ = nreference_configurations();
      if (nRef_ == 0) {
        cout << "database has no reference configurations\n" << flush;
        exit(1);
      }
    };
    SQL_Distance_Interface(string name)
      :SQL_Interface(name,distanceTableName_),nDist_(0),nRef_(0){
      nDist_ = (number_of_columns()-1)/2; 
      nRef_ = nreference_configurations();
      if (nRef_ == 0) {
        cout << "database " << name << " has no reference configurations\n" << flush;
        exit(1);
      }
    };
    virtual ~SQL_Distance_Interface(){};
    int near(int ref_id, double tol, set<int> & ids) { 
      string cmd = "SELECT * FROM "+distanceTableName_+" WHERE cid = "+to_string(ref_id)+";";
      int retval =prep(cmd);
      DB_ERROR_CHECK(retval,FL,LN)
      step(); // return value is a bool signifying at end of list
      int ncols = 2*nDist_+1;
      for (int p=1;p<ncols;p+=2){
        double  d = col_dbl(p+1);
        int id = col_int(p);
#ifdef PRINT_DEBUG
        cout << " distance " << ref_id << " to " << id  << " is " << d << " < " << tol << "\n";
#endif
        if (d<tol) {
          ids.insert(id);
        } else { break; }
      }
      end_query();
      return ids.size();
    }
    int distances(int ref_id, set<int> & ids, map<int,double> & dmap) { 
      string cmd = "SELECT * FROM "+distanceTableName_+" WHERE cid = "+to_string(ref_id)+";";
      int retval =prep(cmd);
      DB_ERROR_CHECK(retval,FL,LN)
      set<int>::iterator itr;
      step(); // return value is a bool signifying at end of list
      int ncols = 2*nDist_+1;
      for (int p=1;p<ncols;p+=2){
        double  d = col_dbl(p+1);
        int id = col_int(p);
        itr=ids.find(id);
        if (itr != ids.end() ) {
#ifdef PRINT_DEBUG
        cout << " distance " << ref_id << " to " << id  << " is " << d << " < " << tol << "\n";
#endif
          dmap[id] = d;
        } 
        else {
          cout << "!!! distance " << ref_id << " to " << id  << " NOT FOUND\n";
          dmap[id] = 1.e8; // HACK dummy value
        }
      }
      end_query();
      return dmap.size();
    }
    int nreference_configurations() { return size(referenceTableName_); }
    int reference_ids(set<int> & ids) {
      ids.clear();
      string cmd = "SELECT * FROM "+referenceTableName_+" ;";
      prep(cmd);
      while (step()) {
        ids.insert(col_int(0));
      }
      return ids.size();
    }
  private:
    SQL_Distance_Interface();
    int nDist_;
    int nRef_;
};// class
class ClusterDatabase : public SQL_Cluster_Interface
{
  public:
    ClusterDatabase(string name, int nn, int max_local, double dtol );
    virtual ~ClusterDatabase();
    // ==================================================================
    // macro interface
    // ==================================================================
    int force(const MatrixXd & C, const int ctype, const string C_type, double * f);
    // ==================================================================
    // micro interface
    // ==================================================================
    double dist(const MatrixXd & Q,const MatrixXd & C, Matrix3d & R) {
      return  d(Q, C, R);
    }
    int distances(int ref_id, set<int> & ids, map<int,double> & dmap) { 
      return distanceDB_->distances(ref_id,ids,dmap);
    }
    int exhaustive_search(const MatrixXd & C, const int ctype, const string C_type);
    int metric_search(const MatrixXd & C, const int ctype, const string C_type);
    bool check(const MatrixXd & Q, const MatrixXd & C);
    void interpolate_force(double * f);
    bool has_distances() { return distanceDB_; }
    void delete_distance_database() { delete distanceDB_; distanceDB_ = NULL;}
    void use_exhaustive_search() { useExhaustiveSearch_ = true; } 
    void report();
    // ==================================================================
    // query interface
    // ==================================================================
    int prep_near(int ref_id, int ctype, string cluster_type, double tol) {
      set<int> ids;
      int n = distanceDB_->near(ref_id,tol,ids);
      if (n==0) return n;
      return prep_ids(ids);
    }
    int prep_ref(int ctype, string cluster_type) {
      set<int> ids;
      int n = distanceDB_->reference_ids(ids);
      if (n==0) return n;
      return prep_ids(ids);
    }
    int prep_ids(set<int> & ids) {
      int n = ids.size();
      string cmd = "SELECT * FROM "+clusterTableName_+" WHERE cid = ";
      for(set<int>::iterator itr = ids.begin(); itr != ids.end(); itr++) {
         cmd += to_string(*itr)+" OR cid =";
      }
      cmd.erase(cmd.end()-8,cmd.end());
      int retval = prep(cmd);
      DB_ERROR_CHECK(retval,FL,LN)
      return n;
    }
    void error(string err);
  protected:
    SQL_Distance_Interface * distanceDB_;
    // parameters
    int n_neighbs_; // size of neighborhood
    int maxLocalClusters_; // max count in ball
    double dtol_; // radius of interpolation ball
    // workspace
    int nDistComps_; // current number of distance comparisons
    int iquery_; // current query count
    MatrixXd C_;
    int clusterCount_;
    int min_id_;
    double min_d_;
    double * clusterForces_;
    double * clusterDistances_;
    double * clusterWeights_;
    int * clusterIds_;
    bool useExhaustiveSearch_;
  private:
    ClusterDatabase();
};//class
};// namespace
#endif

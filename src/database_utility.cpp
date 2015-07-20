#include "database_utility.h"
using namespace SQL_interface;
using std::cout;
using namespace std;
using namespace Eigen;

#define VERBOSE
#define PRINT_DEBUG

/*===================================================================
  cluster-cluster distance
=====================================================================*/
ClusterDatabase::ClusterDatabase(string name, 
                                 int nn, 
                                 int max_local,
                                 double dtol)
 :SQL_Cluster_Interface(name),
  distanceDB_(NULL),
  n_neighbs_(nn),
  maxLocalClusters_(max_local),
  dtol_(dtol),
  nDistComps_(0),
  iquery_(1),
  useExhaustiveSearch_(false)
{
  if (has_table(distanceTableName_))  {
    distanceDB_ = new SQL_Distance_Interface(database());
    // NOTE load refs into memory
  }
  C_.resize(n_neighbs_,3);
  clusterForces_    = new double[maxLocalClusters_*3];
  clusterDistances_ = new double[maxLocalClusters_];
  clusterWeights_   = new double[maxLocalClusters_];
  clusterIds_   = new int[maxLocalClusters_];
}
ClusterDatabase::~ClusterDatabase() 
{ 
  if (distanceDB_)       delete distanceDB_;
  if (clusterForces_)    delete clusterForces_;
  if (clusterDistances_) delete clusterDistances_;
  if (clusterWeights_)   delete clusterWeights_;
  if (clusterIds_)       delete clusterIds_;
}
/*-------------------------------------------------------------------
  Search of the database and interpolate forces
---------------------------------------------------------------------*/
int ClusterDatabase::force(const MatrixXd & C, 
                           const int ctype, 
                           const string Ctype, 
                           double * force) 
{
#ifdef VERBOSE
  cout << "====START== query number " << iquery_ << " ======\n" << std::flush;
  print_config(C);
#endif
  int size = 0;
  if (has_distances() && ! useExhaustiveSearch_) 
    size = metric_search    (C, ctype, Ctype);
  else                 
    size = exhaustive_search(C, ctype, Ctype);
  if (size > 0) interpolate_force(force);
  else {
    printf("----------------- begin cluster -----------------------\n");
    printf("%d type %s \n", ctype,Ctype.c_str());
    print_config(C_);
    printf("------------------  end cluster -----------------------\n");
    printf("query: %d, closest: %i, distance: %g > %g\n",iquery_,min_id_,min_d_,dtol_);
    set<int> ids; ids.insert(min_id_); prep_ids(ids); step();
    configuration(C_);
    printf("--------------- begin best match cluster -------------------\n");
    print_config(C_);
cout << ">>>> retrieve last ROTATION\n"; // HACK
    printf("---------------  end best match cluster --------------------\n");
    cout << "====END== query number " << iquery_ << " ======\n" << std::flush;
    error("no close clusters, add to database");
  }
#ifdef VERBOSE
  printf(">> force: %22.16g %22.16g %22.16g\n", force[0],force[1],force[2]);
  cout << "====END== query number " << iquery_ << " min dist " << min_d_ << " min id " << min_id_ << " ======\n" << std::flush;
#endif
  iquery_ += 1;
  return size;
}
/*-------------------------------------------------------------------
  Perform an exhaustive search of the database
---------------------------------------------------------------------*/
int ClusterDatabase::exhaustive_search(const MatrixXd & Q, 
                                       const int qtype, 
                                       const string Q_type)
{
#ifdef VERBOSE
  printf("> exhaustive search\n");
#endif
  prep_same_type(qtype,Q_type);
  min_id_ = -1; 
  min_d_ = big_;
  clusterCount_ = 0; 
  while(step()){  
    configuration(C_);
    check(Q, C_);
    if (clusterCount_ == maxLocalClusters_) break;
  }
#ifdef PRINT_DEBUG
  printf(">>found %i close clusters\n",clusterCount_);
#endif
  end_query();
  return clusterCount_; 
}
/*-------------------------------------------------------------------
  Perform an metric search of the database
---------------------------------------------------------------------*/
int ClusterDatabase::metric_search(const MatrixXd & Q, 
                                   const int qtype, 
                                   const string Q_type)
{
#ifdef PRINT_DEBUG
  printf("> metric search\n");
#endif
  min_id_ = 1; // should be a initial set : loop over refs or lasts
  min_d_ = big_;
  clusterCount_ = 0; 
  // (A) select reference
  int n = prep_ref(qtype,Q_type);
  if (n==0) { error("no reference clusters"); } 
  while(step()){ 
    configuration(C_);
    check(Q, C_);
    if (clusterCount_ == maxLocalClusters_) break;
  }
#ifdef PRINT_DEBUG
  printf(">> best reference id: %i of %i, distance %g\n",min_id_,n,min_d_);
  printf("> found %4i/%4i\n",clusterCount_,maxLocalClusters_);
#endif
  int last_id = -2;
  double R = 2.*min_d_;
  while( min_id_ != last_id && clusterCount_ != maxLocalClusters_) {
    last_id = min_id_;
    R = 2.*min_d_;
    // retrieve all near min id
    int n = prep_near(min_id_,qtype,Q_type,R);
#ifdef PRINT_DEBUG
    printf(">> search ball radius %g %g size %i\n",R,R/dtol_,n+1);// include center
    printf("> found %4i/%4i\n",clusterCount_,maxLocalClusters_);
#endif
    if ( n==0 || R < dtol_ ) { 
      if (clusterCount_ > 0) { break;} // no additional close clusters
      else { error("no clusters near search id"); }
    }
    // find closer 
    while(step()){ 
      configuration(C_);
      if (check(Q, C_) ) { break; }
      if (clusterCount_ == maxLocalClusters_) break;
    }
    end_query();
#ifdef PRINT_DEBUG
    printf(">> min id %i prev %i best %g\n",min_id_,last_id,min_d_);
#endif
  }
#ifdef PRINT_DEBUG
  printf(">> found %i close clusters\n",clusterCount_);
#endif
  return clusterCount_; 
}

/*-------------------------------------------------------------------
  Perform an detailed search for a particular pair, returns true stop looking
---------------------------------------------------------------------*/
bool ClusterDatabase::check(const MatrixXd & Q, const MatrixXd & C)
{
  nDistComps_++; // for statistics
  Matrix3d R;
  double d=dist(Q, C, R);
  int ii = id();
#ifdef VERBOSE
  if (isnan(d)) error("distance is NaN");
  printf("- distance to id %4i: %10.6f | %10.6f\n",ii,d,min_d_);
#endif
  bool new_min = false;
  if (d < min_d_) { 
    min_d_ = d; 
    min_id_ = ii; 
    new_min = true;
#ifdef VERBOSE
    printf("* new min id %i distance %g\n",min_id_,min_d_);
#endif
  }
  // inside interpolation ball, add configuration to list
  if (d < dtol_){
    clusterIds_[clusterCount_] = ii;
    clusterDistances_[clusterCount_] = d;
    double * cf = &(clusterForces_[3*clusterCount_]); // NOTE use Eigen
    double f[3];
    cluster_force(f);
    // matrix multiply NOTE use Eigen
    for (int i = 0; i<3; i++){
      cf[i] = R(i,0)*f[0]+R(i,1)*f[1]+R(i,2)*f[2];
    }
    clusterCount_++;
#ifdef VERBOSE
    printf("close ids:\n");
    for (int k = 0; k < clusterCount_; ++k) { printf("%4i:%10.6f\n",clusterIds_[k],clusterDistances_[k]); };
    printf("cluster:\n");
    print_config(C);
#endif
  }
  return new_min;
}
//--------------------------------------------------------------------
// interpolate cluster forces
//--------------------------------------------------------------------
double rbf(double r, double sigma) {
  r *= 1./sigma; 
  return 1./(1.-r*r);
}
void ClusterDatabase::interpolate_force(double * f)
{
  set<int> ids; 
  for (int i = 0; i < clusterCount_ ; i++) {
    ids.insert(clusterIds_[i]); 
    cout << i << " id:" << clusterIds_[i] << "\n";
  }
  MatrixXd A(clusterCount_,clusterCount_);
  for (int i = 0; i < clusterCount_ ; i++) {
    map<int,double> dmap;
    set<int> iids = ids; iids.erase(clusterIds_[i]);
    distances(clusterIds_[i],iids,dmap);
    A(i,i) = rbf(0.,dtol_);
    for (int j = i; j < clusterCount_ ; j++) {
      double dij = dmap[clusterIds_[j]]; 
      A(i,j) = A(j,i) = rbf(dij,dtol_);
    }
  }
  MatrixXd F(clusterCount_,3);
  for (int i = 0; i < clusterCount_ ; i++) {
    F(i,0) = clusterForces_[3*i]; 
    F(i,1) = clusterForces_[3*i+1]; 
    F(i,2) = clusterForces_[3*i+2];
  }
  MatrixXd B(clusterCount_,3);
  B = A.ldlt().solve(F);
#ifdef PRINT_DEBUG
//cout << A << "\n";
#endif
  f[0] = 0;
  f[1] = 0;
  f[2] = 0;
  for (int i = 0; i < clusterCount_ ; i++) {
    double w = rbf(clusterDistances_[i],dtol_);
    f[0] += w*B(i,0);
    f[1] += w*B(i,1);
    f[2] += w*B(i,2);
#ifdef PRINT_DEBUG
    printf("   %2i %4i distance %8.6g weight %4.2g coef %9.6g  %9.6g  %9.6g\n",i+1,clusterIds_[i],clusterDistances_[i],w,B(i,0),B(i,1),B(i,2));
#endif
  }
// HACK
  f[0] = clusterForces_[0];
  f[1] = clusterForces_[1];
  f[2] = clusterForces_[2];
}

//--------------------------------------------------------------------
// reporting
//--------------------------------------------------------------------
void ClusterDatabase::report() {
  double ave = nDistComps_/iquery_;
  double frac = ave/size();
  printf(">> ave. dist. comps: %g frac. of db used: %g\n",ave,frac);
}

//--------------------------------------------------------------------
// error handling
//--------------------------------------------------------------------
void ClusterDatabase::error(string error) {
  cout << "!!! " << error << " !!!\n" << std::flush;
  exit (EXIT_FAILURE);
}

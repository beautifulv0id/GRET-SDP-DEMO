#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <Eigen/Dense>
#include <sdpa_call.h>
#include <gr/accelerators/kdtree.h>
#include "shared.h"
#include <atomic>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "fusion.h"
using namespace mosek::fusion;
using namespace monty;


static std::string path_prefix  = "";
namespace pt = boost::property_tree;

static int n = 1000;
static const int d = 3;
static const int m = 5;
static const double samp = 0.7; // sample probability
using PointType = gr::Point3D<double>;
using Scalar = PointType::Scalar;
using MatrixType = Eigen::Matrix<Scalar, d+1, d+1>;
using VectorType = Eigen::Matrix<Scalar, d, 1>;
using RigidTransformation = std::tuple<MatrixType, VectorType>;
using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
using PatchType = vector<pair<PointType, int>>;


template <typename InRange>
void CreateSpiral(InRange& in, const size_t n){
  using Scalar_ = typename InRange::value_type::Scalar;
  using VectorType_ = typename InRange::value_type::VectorType;
  using PointType_ = typename InRange::value_type;

  Scalar_ step = 4.0*M_PI/ (float) n;
  VectorType_ tmp(0.1,0.1,01.); // x, y, z
  for(int i = 0; i<n; ++i){
    in.push_back(PointType_(tmp));
    tmp(2) += step;
    tmp(1) = 2*cos(tmp(2)) + ((double) rand() / (RAND_MAX));
    tmp(0) = 2*sin(tmp(2)) + ((double) rand() / (RAND_MAX));
  }
}

MatrixType GenerateTransformationMatrix(){
  // generate rotation matrix
  VectorX eul(d);
  Scalar sign;
  Eigen::Matrix3d rm;
  for(int i = 0; i < d; i++){
    sign = rand() % 2 == 0? -1 : 1;
    eul(i) = sign * M_PI * ((double) rand() / (RAND_MAX));
  }

  rm = Eigen::AngleAxisd(eul[0], Eigen::Vector3d::UnitX())
      * Eigen::AngleAxisd(eul[1], Eigen::Vector3d::UnitY())
      * Eigen::AngleAxisd(eul[2], Eigen::Vector3d::UnitZ()); 


  // generate translation vector
  VectorX t(d);
  Scalar scale = 2;
  for(int i = 0; i < d; i++){
    sign = rand() % 2 == 0? -1 : 1;
    t(i) = sign * ((double) rand() / (RAND_MAX)) * scale;
  }

  MatrixType trafo = MatrixType::Zero();
  trafo.block(0,0,d,d) = rm;
  trafo.block(0,d,d,1) = t;
  trafo(d,d) = 1;

  return trafo;
}


template <typename TrRange>
void GenerateTrafos(TrRange& transformations){
  for(int i = 0; i < m; i++){
    transformations.push_back(GenerateTransformationMatrix());
  }
}

template <typename InRange, typename PatchRange, typename TrRange>
void GeneratePatches(const InRange& pointCloud, PatchRange& patches, TrRange& transformations, 
                      Eigen::Ref<MatrixX> L, Eigen::Ref<MatrixX> B, Eigen::Ref<MatrixX> D){

  L = Eigen::MatrixXd::Zero(n+m,n+m);
  B = Eigen::MatrixXd::Zero(m*d, n+m);
  D = Eigen::MatrixXd::Zero(m*d, m*d);


  GenerateTrafos(transformations);
  VectorType v_ki;


  // for every global point k
  for(int k = 0; k < n; k++){
    // for every patch i
    for(int i = 0; i < m; i++){
      if(((double) rand() / (RAND_MAX)) < samp){
        v_ki = transformations.at(i).block(0,0,d,d).inverse() * (pointCloud.at(k).pos() -  transformations.at(i).block(0,d,d,1));

        patches[i].emplace_back(PointType(v_ki), k);

        L(k,k)++;
        L(k,n+i)--;
        L(n+i,k)--;
        L(n+i,n+i)++;

        B.block(i*d,k,d,1) += v_ki;
        B.block(i*d,n+i,d,1) -= v_ki;

        D.block(i*d, i*d, d, d) += v_ki*v_ki.transpose();
      }
    }

    // ensure that every point is in one patch
    if(L(k,k) == 0){
      int i = rand() % m;
      //v_ki = transformations.at(i).block(0,0,d,d) * pointCloud.at(k).pos() + transformations.at(i).block(0,d,d,1);
      v_ki = transformations.at(i).block(0,0,d,d).inverse() * (pointCloud.at(k).pos() -  transformations.at(i).block(0,d,d,1));
      
      patches[i].emplace_back(PointType(v_ki), k);

      L(k,k)++;
      L(k,n+i)--;
      L(n+i,k)--;
      L(n+i,n+i)++;

      B.block(i*d,k,d,1) += v_ki;
      B.block(i*d,n+i,d,1) -= v_ki;

      D.block(i*d, i*d, d, d) += v_ki*v_ki.transpose();
    } 
  }

}

void printVector(double* ele, int dim, char* printFormat,
         FILE* fpout);
void printMatrix(double* ele, int dim, char* printFormat,
         FILE* fpout);
void printDimacsError(double dimacs_error[7],char* printFormat,
              FILE* fpout);


void SolveSDP_SDPA(Eigen::Ref<const MatrixX> C, Eigen::Ref<MatrixX> G){
  SDPA	Problem;

  // All parameteres are renewed
  Problem.setParameterType(SDPA::PARAMETER_STABLE_BUT_SLOW);

  Problem.printParameters(stdout);

  int mDIM   = d*(d+1)/2*m;
  int nBlock = 1;
  Problem.inputConstraintNumber(mDIM);
  Problem.inputBlockNumber(nBlock);
  Problem.inputBlockSize(1,d*m);
  Problem.inputBlockType(1,SDPA::SDP);

  Problem.initializeUpperTriangleSpace();

  // c vec
  int cnt = 1;
  for(int i = 0; i < m; i++){
    for(int j = 0; j < d; j++){
      Problem.inputCVec(cnt++,1);
      for(int k = j+1; k < d; k++)
        Problem.inputCVec(cnt++,0);
    }
  } 

  // F0 = -C
  for(int i = 0; i < m*d; i++)
    for(int j = i; j < m*d; j++)
      Problem.inputElement(0, 1, i+1, j+1, -C(i,j));

 // Fi
  cnt = 1;
  for(int k = 0; k < m; k++){
    for(int i = 0; i < d; i++){
      for(int j = i; j < d; j++){
        Problem.inputElement(cnt++, 1, k*d+i+1, k*d+j+1, 1);
      }
    }
  }

  Problem.initializeUpperTriangle();
  Problem.initializeSolve();

  // if necessary, dump input data and initial point
  // Problem1.writeInputSparse((char*)"tmp.dat-s",(char*)"%+8.3e");
  // Problem1.writeInitSparse((char*)"tmp.ini-s",(char*)"%+8.3e");

  Problem.solve();

  fprintf(stdout, "\nStop iteration = %d\n",
      Problem.getIteration());
  char phase_string[30];
  Problem.getPhaseString(phase_string);
  fprintf(stdout, "Phase          = %s\n", phase_string);
  fprintf(stdout, "objValPrimal   = %+10.6e\n",
      Problem.getPrimalObj());
  fprintf(stdout, "objValDual     = %+10.6e\n",
      Problem.getDualObj());
  fprintf(stdout, "p. feas. error = %+10.6e\n",
      Problem.getPrimalError());
  fprintf(stdout, "d. feas. error = %+10.6e\n\n",
      Problem.getDualError());


  fprintf(stdout, "xVec = \n");
  // Problem1.printResultXVec();
  printVector(Problem.getResultXVec(),
          Problem.getConstraintNumber(), (char*)"%+8.3e",
          stdout);

  fprintf(stdout, "xMat = \n");
  // Problem1.printResultXMat();
  for (int l=0; l<Problem.getBlockNumber(); ++l) {
    if (Problem.getBlockType(l+1) == SDPA::SDP) {
      printMatrix(Problem.getResultXMat(l+1),
          Problem.getBlockSize(l+1), (char*)"%+8.3e",
          stdout);
    }
    else if (Problem.getBlockType(l+1) == SDPA::SOCP) {
      printf("current version does not support SOCP\n");
    }
    if (Problem.getBlockType(l+1) == SDPA::LP) {
      printVector(Problem.getResultXMat(l+1),
          Problem.getBlockSize(l+1), (char*)"%+8.3e",
          stdout);
    }
  }

  fprintf(stdout, "yMat = \n");
  // Problem1.printResultYMat();
  for (int l=0; l<Problem.getBlockNumber(); ++l) {
    if (Problem.getBlockType(l+1) == SDPA::SDP) {
      printMatrix(Problem.getResultYMat(l+1),
          Problem.getBlockSize(l+1), (char*)"%+8.3e",
          stdout);
    }
    else if (Problem.getBlockType(l+1) == SDPA::SOCP) {
      printf("current version does not support SOCP\n");
    }
    if (Problem.getBlockType(l+1) == SDPA::LP) {
      printVector(Problem.getResultYMat(l+1),
          Problem.getBlockSize(l+1), (char*)"%+8.3e",
          stdout);
    }
  }

  double* yMat = Problem.getResultYMat(1);
  Eigen::Map<Eigen::MatrixXd> Gmap(yMat, d*m, d*m);
  G = Gmap;
}


void SolveSDP_MOSEK(Eigen::Ref<const MatrixX> C, Eigen::Ref<MatrixX> G){
  // Create a model with the name 'lo1'
  Model::t M = new Model("gret-sdp"); auto _M = finally([&]() { M->dispose(); });

  auto c_ptr = new_array_ptr<Scalar, 2>(shape_t<2>(m*d,m*d));
  std::copy(C.data(), C.data()+C.size(), c_ptr->begin());
  Matrix::t C_ = Matrix::dense(c_ptr);

  M->setLogHandler([ = ](const std::string & msg) { std::cout << msg << std::flush; } );

  auto G_ = M->variable(Domain::inPSDCone(m*d));

  for (int i = 0; i < m; i++){
    for (int j = 0; j < d; j++){
          M->constraint(G_->index(i*d+j, i*d+j),Domain::equalsTo(1.0));
      for (int k = j+1; k < d; k++)
          M->constraint(G_->index(i*d+j, i*d+k), Domain::equalsTo(0.0));
    }
  }  
  
  // Set the objective function to (Tr(CG)=sum(C⋅G))
  M->objective("obj", ObjectiveSense::Minimize, Expr::dot(C_, G_));

  // Solve the problem
  M->solve();
  
  // Get the solution values
  auto sol = G_->level();
  Eigen::Map<Eigen::MatrixXd> Gmap(sol->begin(), d*m, d*m);
  G = Gmap;
}



template <typename TrRange>
void ComputeRelativeTrafos(const TrRange& transformations, TrRange& relTransformations, bool inv = false){
  relTransformations.reserve(m-1);
  MatrixType Tr;
  Eigen::MatrixXd O1(transformations[0].block(0, 0, d, d));
  Eigen::VectorXd t1(transformations[0].block(0, d, d, 1));
  Eigen::MatrixXd O(d,d);
  Eigen::VectorXd t(d);


  if(!inv){
    O1 = O1.transpose();
    for(int i = 1; i < m; i++){
      Tr = MatrixType::Zero();
      O = transformations[i].block(0, 0, d, d);
      t = transformations[i].block(0, d, d, 1);
      Tr.block(0, 0, d, d) = O1*O;
      Tr.block(0, d, d, 1) = O1*(t-t1);
      Tr(d,d) = 1;
      relTransformations.push_back(Tr);
    }
  } else {
    for(int i = 1; i < m; i++){
      Tr = transformations.at(i);
      O = Tr.block(0, 0, d, d);
      O = O.transpose();
      t = Tr.block(0, d, d, 1);
      Tr.block(0, 0, d, d) = O1*O;
      Tr.block(0, d, d, 1) = -O1*t+t1;
      relTransformations.push_back(Tr);
    } 
  }

}

void writeMatrix(Eigen::Ref<const Eigen::MatrixXd> m, std::string filename){
  std::ofstream file(path_prefix + filename);
  if(file.is_open()){
    file << m << '\n';
  }
  file.close();
}


template <typename PointRange>
void writePoints(const PointRange& vec, std::string filename){
  std::ofstream file(path_prefix + filename);

  for(int i = 0; i < d; i++){
   for(int j = 0; j < vec.size()-1; j++){
      file << vec[j][i] << " ";
    }
    file << vec[vec.size()-1][i] << std::endl;
  }

  file.close();
}


template <typename PointRange>
Scalar compute_lcp( const gr::KdTree<Scalar>& P, const PointRange& Q){
  using RangeQuery = typename gr::KdTree<Scalar>::template RangeQuery<>;
  const Scalar epsilon = 0.01;
  std::atomic_uint good_points(0);
  const size_t number_of_points = Q.size();
  Scalar best_LCP_ = 0;
  const size_t terminate_value = best_LCP_ * number_of_points;
  const Scalar sq_eps = epsilon*epsilon;

  for (size_t i = 0; i < number_of_points; ++i) {
    RangeQuery query;
    query.queryPoint = Q[i].pos();
    query.sqdist     = sq_eps;

    auto result = P.doQueryRestrictedClosestIndex( query );

    if ( result.first != gr::KdTree<Scalar>::invalidIndex() ) {
        good_points++;
    }

    // We can terminate if there is no longer chance to get better than the
    // current best LCP.
    if (number_of_points - i + good_points < terminate_value) {
        break;
    }
  }
  return Scalar(good_points) / Scalar(number_of_points);
}


template <typename PointRange>
gr::KdTree<Scalar> constructKdTree(const PointRange& Q){
  size_t number_of_points = Q.size();
  // Build the kdtree.
  gr::KdTree<Scalar> kd_tree(number_of_points);

  for (size_t i = 0; i < number_of_points; ++i) {
      kd_tree.add(Q[i].pos());
  }
  kd_tree.finalize();
  return kd_tree;
}


// transforms all patches to frame 0 using the relTrafos and stores them in P (in)
template <typename CommonPointRange, typename PatchRange, typename TrRange>
void transformToCommonFrame(CommonPointRange& P, const PatchRange& patches, const TrRange& relTrafos){
  MatrixType mat;
  VectorType vec;
  P = patches[0];
  // transform and add remaining patches
  for( int i = 1; i < patches.size(); i++){
    int patch_size = patches[i].size();
    mat = relTrafos[i-1];
    for(int j = 0; j < patch_size; j++){
      vec = (mat * patches[i][j].pos().homogeneous()).template head<3>();
      P.push_back(PointType(vec));
    }
  }
}

template <typename PatchRange, typename TrRange>
void exportConfig(const PatchRange& patches, const TrRange& transformations){

  // export patches
  PointType point;
  int index;
  std::ofstream file;
  for(int i = 0; i < patches.size(); i++){
    file.open("data/patch" + std::to_string(i) + ".dat");
    if(file.is_open()){
      file << patches[i].size() << std::endl;
      for(const std::pair<PointType,int> & point_with_index : patches[i]){
        point = point_with_index.first;
        index = point_with_index.second;
        file << index << " ";
        file << point.x() << " ";
        file << point.y() << " ";
        file << point.z() << " ";
        file << std::endl;
      }
    }
    file.close();
  }

  // export config file
  file.open("data/config.json");
  file << "{" << std::endl;
  file << "\t" << "\"n\":" << "\t" << n << "," << std::endl;
  file << "\t" << "\"m\":" << "\t" << m << "," << std::endl;
  file << "\t" << "\"d\":" << "\t" << d << "," << std::endl;

  file << "\t" <<  "\"patches\":" << "[";
  for(int i = 0; i < patches.size(); i++){
    file << "\"patch" + std::to_string(i) + ".dat\"" << ((i == m-1)? "" : ", ");
  }
  file << "]," << std::endl;
  file << "\t" <<  "\"gt_trafos\":" << "[";
  for(int i = 0; i < patches.size(); i++){
    file << "\"gt_trafo" + std::to_string(i) + ".dat\"" << ((i == m-1)? "" : ", ");
  }
  file << "]" << std::endl;
  file << "}";
  file.close();

  // export ground truth transformations
  for (int i = 0; i < transformations.size(); i++){
    file.open("data/gt_trafo" + std::to_string(i) + ".dat");
    if(file.is_open()){
      file << d+1 << " " << d+1 << std::endl;
      file << transformations[i];
    }
    file.close();
  }
}


int main (int argc, char** argv)
{
  if(argc > 1)
    n = atoi(argv[1]);

  srand(time(NULL));

  vector<gr::Point3D<Scalar>> spiral;
  spiral.reserve(n);
  CreateSpiral(spiral, n);
  //writePoints(spiral, "spiral.dat");

  vector<MatrixType> original_transformations;
  vector<PatchType> patches(m);
  MatrixX L(n+m, n+m);
  MatrixX B(m*d, n+m);
  MatrixX D(m*d, m*d);
  MatrixX Linv(n+m, n+m);
  MatrixX C(m*d, m*d);
  GeneratePatches(spiral, patches, original_transformations, L, B, D);
  
  // compute Linv, the Moore-Penrose pseudoinverse of L
  Linv = L.completeOrthogonalDecomposition().pseudoInverse();

  // compute C
  C = D - B * Linv * B.transpose();

  // solve the SDP (P2) using C
  Eigen::Matrix<Scalar, m*d, m*d> G;

  //SolveSDP_SDPA(C, G);
  SolveSDP_MOSEK(C,G);

  // compute top d eigenvalues and eigenvectors
  Eigen::EigenSolver<Eigen::Matrix<Scalar, m*d, m*d>> s(G);
  Eigen::VectorXcd eigvals = s.eigenvalues();
  Eigen::MatrixXcd eigvecs = s.eigenvectors();
  std::vector<double> re_eigvals;
  std::vector<Eigen::VectorXd> re_eigvecs;

  for(int i = 0; i < eigvals.rows(); i++) {
      if(eigvals(i).imag() == 0){
          re_eigvals.push_back(eigvals(i).real());
          re_eigvecs.push_back(eigvecs.col(i).real());
      }
  }   

  std::vector<std::pair<double, Eigen::VectorXd>> eig_pairs;
  eig_pairs.reserve(d);
  std::transform(re_eigvals.begin(), re_eigvals.end(), re_eigvecs.begin(), std::back_inserter(eig_pairs),
               [](double a, const Eigen::VectorXd& b) { return std::make_pair(a, b); });

  sort(eig_pairs.begin(), eig_pairs.end(),
    [&](std::pair<double, Eigen::VectorXd>& a, std::pair<double, Eigen::VectorXd>& b) {
        return (a.first > b.first);
    }
  );

  // construct W
  Eigen::Matrix<Scalar, d, m*d> W;
  for(int i = 0; i < d; i++){
      W.row(i) = std::sqrt(eig_pairs[i].first) * eig_pairs[i].second.transpose();
  }

  // compute transformations O
  Eigen::Matrix<Scalar, d, m*d> O;
  for(int i = 0; i < m; i++){
      Eigen::Ref<Eigen::Matrix<Scalar, d, d>> w(W.block(0, i*d, d, d));
      Eigen::JacobiSVD<Eigen::Matrix<Scalar, d, d>> svd(w, Eigen::ComputeFullU | Eigen::ComputeFullV);
      O.block<d, d>(0, i*d) = svd.matrixU() * svd.matrixV().transpose();
  }

  // compute Z
  MatrixX Z(d, n+m);
  Z = O*B*Linv;

  // compute trafos
  std::vector<MatrixType> transformations;
  MatrixType trafo;
  transformations.reserve(m);
  for(int i = 0; i < m; i++){
    trafo = MatrixType::Zero();
    trafo.block(0, 0, d, d) = O.block(0, i*d, d, d);
    trafo.block(0, d, d, 1) = Z.block(0, n+i, d, 1);
    trafo(d,d) = 1;
    transformations.push_back(trafo);
  }

  vector<PointType> reg_transformed_patches;
  vector<PointType> ori_transformed_patches;
  int accum_patch_size = 0;
  for(int i = 0; i < m; i++) accum_patch_size += patches[i].size();

  reg_transformed_patches.reserve(accum_patch_size);
  ori_transformed_patches.reserve(accum_patch_size);

  VectorType tmp;
  Eigen::Matrix3d reg_O_0(transformations[0].block(0,0,d,d).transpose());
  Eigen::Vector3d reg_t_0(transformations[0].block(0,d,d,1));
  Eigen::Matrix3d ori_O_0(original_transformations[0].block(0,0,d,d).transpose());
  Eigen::Vector3d ori_t_0(original_transformations[0].block(0,d,d,1));


  // computing the transformed patches for reference frame of patch one
  PointType point;
  int index;
  for(int i = 0; i < m; i++){
    for (const pair<PointType, int>& point_with_index : patches[i]){
        point = point_with_index.first;
        index = point_with_index.second;
        // using computed transformations
        tmp = (transformations[i]*point.pos().homogeneous()).template head<3>();
        tmp = reg_O_0 * (tmp - reg_t_0);
        reg_transformed_patches.emplace_back(PointType(tmp));

        // using ground truth transformations
        tmp = (original_transformations[i]*point.pos().homogeneous()).template head<3>();
        tmp = ori_O_0 * (tmp - ori_t_0);
        ori_transformed_patches.emplace_back(PointType(tmp));
    }
  }

  // construct kd_tree
  gr::KdTree<Scalar> kd_tree(constructKdTree(ori_transformed_patches));
  // compute lcp
  Scalar lcp = compute_lcp(kd_tree, reg_transformed_patches);
  std::cout << "lcp = " << lcp << std::endl;


  exportConfig(patches, original_transformations);
  // export matrices for usage in matlab
  // writePoints(spiral, "SpiralMat.dat");
  // for(int i = 0; i < m; i++){
  //   writePoints(patches[i], "patch" + std::to_string(i+1) + ".dat");
  // }

  // writeMatrix(L, "LMat.dat");
  // writeMatrix(B, "BMat.dat");
  // writeMatrix(Linv, "LinvMat.dat");
  // writeMatrix(D, "DMat.dat");
  // writeMatrix(C, "CMat.dat");
  // writeMatrix(G, "GMat.dat");
  // writeMatrix(Z, "ZMat.dat");
  // writeMatrix(W, "WMat.dat");
  // writeMatrix(O, "OMat.dat");

  // writePoints(reg_transformed_patches, "reg_transformed_patches.dat");
  // writePoints(ori_transformed_patches, "ori_transformed_patches.dat");
  // writePoints(reg_transformed_patches_frame0, "reg_transformed_patches_frame0.dat");
  // writePoints(ori_transformed_patches_frame0, "ori_transformed_patches_frame0.dat");
}



void printVector(double* ele, int dim, char* printFormat, FILE* fpout)
{
  fprintf(fpout,"[ ");
  for (int k=0; k<dim-1; ++k) {
    fprintf(fpout,printFormat,ele[k]);
    fprintf(fpout," ");
  }
  fprintf(fpout,printFormat,ele[dim-1]);
  fprintf(fpout,"]; \n");
}

void printMatrix(double* ele, int dim, char* printFormat, FILE* fpout)
{
  fprintf(fpout,"[\n");
  for (int i=0; i<dim; ++i) {
    fprintf(fpout,"[ ");
    for (int j=0; j<dim-1; ++j) {
      fprintf(fpout,printFormat,ele[i+dim*j]);
      fprintf(fpout," ");
    }
    fprintf(fpout,printFormat,ele[i+dim*(dim-1)]);
    fprintf(fpout,"]; \n");
  }
  fprintf(fpout,"]; \n");
}

void printDimacsError(double dimacs_error[7],char* printFormat,
              FILE* fpout)
{
  fprintf(fpout,  "\n");
  fprintf(fpout,  "* DIMACS_ERRORS * \n");
  fprintf(fpout,  "err1 = ");
  fprintf(fpout,  printFormat, dimacs_error[1]);
  fprintf(fpout, "  [||Ax-b|| / (1+||b||_1)]\n");
  fprintf(fpout,  "err2 = ");
  fprintf(fpout,  printFormat, dimacs_error[2]);
  fprintf(fpout, "  [max(0, -lambda(x)/(1+||b||_1))]\n");
  fprintf(fpout,  "err3 = ");
  fprintf(fpout,  printFormat, dimacs_error[3]);
  fprintf(fpout, "  [||A^Ty + z - c || / (1+||c||_1)]\n");
  fprintf(fpout,  "err4 = ");
  fprintf(fpout,  printFormat, dimacs_error[4]);
  fprintf(fpout, "  [max(0, -lambda(z)/(1+||c||_1))]\n");
  fprintf(fpout,  "err5 = ");
  fprintf(fpout,  printFormat, dimacs_error[5]);
  fprintf(fpout, "  [(<c,x> - <b,y>) / (1 + |<c,x>| + |<b,y>|)]\n");
  fprintf(fpout,  "err6 = ");
  fprintf(fpout,  printFormat, dimacs_error[6]);
  fprintf(fpout, "  [<x,z> / (1 + |<c,x>| + |<b,y>|)]\n");
  fprintf(fpout,  "\n");
}
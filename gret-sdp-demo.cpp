#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <Eigen/Dense>
#include <sdpa_call.h>
#include "shared.h"

static const int n = 1000;
static const int d = 3;
static const int m = 5;
static const double samp = 0.7; // sample probability
using Scalar = gr::Point3D<double>::Scalar;
using MatrixType = Eigen::Matrix<Scalar, d+1, d+1>;
using VectorType = Eigen::Matrix<Scalar, d, 1>;
using RigidTransformation = std::tuple<MatrixType, VectorType>;
using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

template <typename InRange>
void CreateSpiral(InRange& in, const size_t n){
  using Scalar_ = typename InRange::value_type::Scalar;
  using VectorType_ = typename InRange::value_type::VectorType;
  using PointType_ = typename InRange::value_type;

  Scalar_ step = 4.0*M_PI/ (double) n;
  VectorType_ tmp(0,0,0); // x, y, z
  for(int i = 0; i<n; ++i){
    in.push_back(PointType_(tmp));
    tmp(2) += step;
    tmp(1) = 2*cos(tmp(2)) + ((double) rand() / (RAND_MAX));
    tmp(0) = 2*sin(tmp(2)) + ((double) rand() / (RAND_MAX));
  }
}

MatrixType GenerateTransformationMatrix(){
  MatrixType rm = MatrixType::Zero();
  rm(d,d) = 1;

  // generate rotation matrix
  VectorX eul(d);
  Scalar sign;
  for(int i = 0; i < d; i++){
    sign = rand() % 2 == 0? -1 : 1;
    eul(i) = sign * M_PI * ((double) rand() / (RAND_MAX));
  }

  Scalar a = eul(0);
  Scalar b = eul(1);
  Scalar y = eul(2);

  rm(0,0) = cos(a)*cos(b);
  rm(0,1) = cos(a)*sin(b)-sin(a)*cos(y);
  rm(0,2) = cos(a)*sin(b)*cos(y)+sin(a)*sin(y);
  rm(1,0) = sin(a)*cos(b);
  rm(1,1) = sin(a)*sin(b)*sin(y)+cos(a)*cos(y);
  rm(1,2) = sin(a)*sin(b)*cos(y)-cos(a)*sin(y);
  rm(2,0) = -sin(b);
  rm(2,1) = cos(b)*sin(y);
  rm(2,2) = cos(b)*cos(y);

  // generate translation vector
  Scalar scale = 2;
  for(int i = 0; i < d; i++){
    sign = rand() % 2 == 0? -1 : 1;
    rm(i,3) = sign * ((double) rand() / (RAND_MAX)) * scale;
  }
  return rm;
}


template <typename TrRange>
void GenerateTrafos(TrRange& transformations){
  for(int i = 0; i < m; i++){
    transformations.push_back(GenerateTransformationMatrix());
  }
}

template <typename InRange, typename TrRange>
void GeneratePatches(const InRange& pointCloud, TrRange& transformations,
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
        v_ki = transformations.at(i).block(0,0,d,d) * pointCloud.at(k).pos() + transformations.at(i).block(0,d,d,1);

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
      v_ki = transformations.at(i).block(0,0,d,d) * pointCloud.at(k).pos() + transformations.at(i).block(0,d,d,1);

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


void SolveSDP(Eigen::Ref<const MatrixX> C, Eigen::Ref<MatrixX> G){
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


template <typename TrRange>
void ComputeRelativeTrafos(const TrRange& transformations, TrRange& relTransformations, bool inv = false){
  relTransformations.reserve(m-1);
  MatrixType Tr;
  Eigen::MatrixXd O1(transformations.at(0).block(0, 0, d, d));
  if(!inv)
    O1 = O1.inverse();
  Eigen::VectorXd t1(transformations.at(0).block(0, d, d, 1));
  Eigen::MatrixXd O(d,d);
  Eigen::VectorXd t(d);

  for(int i = 1; i < m; i++){
    Tr = transformations.at(i);
    O = Tr.block(0, 0, d, d);
    if(inv) O = O.inverse();
    t = Tr.block(0, d, d, 1);
    Tr.block(0, 0, d, d) = O * O1;
    Tr.block(0, d, d, 1) = (t-t1)*O1;
    relTransformations.push_back(Tr);
  }
}

void writeMatrix(Eigen::Ref<Eigen::MatrixXd> m, std::string filename){
  std::ofstream file(filename);
  if(file.is_open()){
    file << m << '\n';
  }
  file.close();
}

template <typename PointRange>
void writePoints(const PointRange& vec, std::string filename){
  std::ofstream file(filename);
  if(file.is_open()){
    for(auto& p : vec){
      file << p.pos().transpose() << std::endl;
    }
  }
  file.close();
}

int main ()
{
  
  vector<gr::Point3D<Scalar>> spiral;
  spiral.reserve(n);
  CreateSpiral(spiral, n);
  writePoints(spiral, "spiral.dat");

  vector<MatrixType> originalTransformations;
  originalTransformations.reserve(m);
  MatrixX L(n+m, n+m);
  MatrixX B(m*d, n+m);
  MatrixX D(m*d, m*d);
  MatrixX Linv(n+m, n+m);
  MatrixX C(m*d, m*d);
  GeneratePatches(spiral, originalTransformations, L, B, D);
  
  // compute Linv, the Moore-Penrose pseudoinverse of L
  Linv = L.completeOrthogonalDecomposition().pseudoInverse();

  // compute C
  C = D - B * Linv * B.transpose();

  // solve the SDP (P2) using C
  Eigen::Matrix<Scalar, m*d, m*d> G;

  SolveSDP(C, G);

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
               [](double a, Eigen::VectorXd b) { return std::make_pair(a, b); });

  sort(eig_pairs.begin(), eig_pairs.end(),
    [&](std::pair<double, Eigen::VectorXd>& a, std::pair<double, Eigen::VectorXd>& b) {
        return (a.first > b.first);
    }
  );


  // construct W
  Eigen::Matrix<Scalar, d, m*d> W;
  for(int i = 0; i < d; i++){
      W.row(i) = sqrt(eig_pairs[i].first) * eig_pairs[i].second.transpose();
  }

  // compute transformations O
  Eigen::Matrix<Scalar, d, m*d> O;
  for(int i = 0; i < m; i++){
      Eigen::Ref<Eigen::Matrix<Scalar, d, d>> w(W.block(0, i*d, d, d));
      Eigen::JacobiSVD<Eigen::Matrix<Scalar, d, d>> svd(w, Eigen::ComputeFullU | Eigen::ComputeFullV);
      O.block<d, d>(0, i*d) = svd.matrixU() * svd.matrixV().transpose();
  }

  // compute Z
  Eigen::Matrix<Scalar, d, n+m> Z;
  Z = O*B*Linv;

  // writeMatrix(L, "LMat.dat");
  // writeMatrix(B, "BMat.dat");
  // writeMatrix(Linv, "LinvMat.dat");
  // writeMatrix(D, "DMat.dat");
  // writeMatrix(C, "CMat.dat");
  // writeMatrix(G, "GMat.dat");
  // writeMatrix(Z, "ZMat.dat");
  // writeMatrix(W, "WMat.dat");
  // writeMatrix(O, "OMat.dat");

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

  std::vector<MatrixType> relTransformations;
  ComputeRelativeTrafos(transformations, relTransformations);

  std::vector<MatrixType> origRelTransformations;
  ComputeRelativeTrafos(originalTransformations, origRelTransformations, true);

  for(int i = 0; i < m-1; i++){
    std::cout << "Tr*[" << i << "] " << endl << relTransformations.at(i) << endl;
    std::cout << "Tr**[" << i << "] " << endl << origRelTransformations.at(i) << endl;
  }
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
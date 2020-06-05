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
using Scalar = gr::Point3D<float>::Scalar;
using MatrixType = Eigen::Matrix<Scalar, d+1, d+1>;
using VectorType = Eigen::Matrix<Scalar, d, 1>;
using RigidTransformation = std::tuple<MatrixType, VectorType>;
using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

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

template <typename InRange, typename TrRange>
void GeneratePatches(const InRange& pointCloud, TrRange& transformations, 
                      Eigen::Ref<MatrixX> L, Eigen::Ref<MatrixX> B, Eigen::Ref<MatrixX> D){

  for(int i = 0; i < m; i++){
    
  }
}


int main ()
{
  
  vector<gr::Point3D<Scalar>> spiral;
  spiral.reserve(n);
  CreateSpiral(spiral, n);

  vector<RigidTransformation> originalTransformations;
  originalTransformations.reserve(m);
  MatrixX L(n+m, n+m);
  MatrixX B(m*d, n+m);
  MatrixX D(m*d, m*d);
  MatrixX Linv(n+m, n+m);
  MatrixX C(m*d, m*d);
  GeneratePatches(spiral, originalTransformations, L, B, D);

}




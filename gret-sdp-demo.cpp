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
  GenerateTrafos(transformations);
  for(int i = 0; i < m; i++){
  }
}


int main ()
{
  
  vector<gr::Point3D<Scalar>> spiral;
  spiral.reserve(n);
  CreateSpiral(spiral, n);

  vector<MatrixType> originalTransformations;
  originalTransformations.reserve(m);
  MatrixX L(n+m, n+m);
  MatrixX B(m*d, n+m);
  MatrixX D(m*d, m*d);
  MatrixX Linv(n+m, n+m);
  MatrixX C(m*d, m*d);
  GeneratePatches(spiral, originalTransformations, L, B, D);

}




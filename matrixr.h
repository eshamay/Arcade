#ifndef MATRIXR_H_
#define MATRIXR_H_

#include "vecr.h"
#include <list>

//#define EIGEN_MATRIXBASE_PLUGIN "/home/src/Arcade/EigenMatrixAddon.h"
//#include <Eigen/Core>
//#include <Eigen/Dense>
#include <Eigen/LU>
//USING_PART_OF_NAMESPACE_EIGEN

typedef enum {xx=0, yx=1, zx=2, xy=3, yy=4, zy=5, xz=6, yz=7, zz=8} element;

typedef Matrix3d	MatR;
typedef std::list<MatR, Eigen::aligned_allocator<MatR> > 	MatR_list;
typedef MatR_list::const_iterator MatR_it;

#endif

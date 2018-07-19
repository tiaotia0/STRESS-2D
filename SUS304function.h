#ifndef _SUS304FUNCTION_
#define _SUS304FUNCTION_
#include<stdafx.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<interpolation.h>
#include<Eigen/Dense>
#include"SUS304properties.h"
using namespace Eigen;
using namespace alglib;
/*-----------------------------------*/
/*       SUS304性能及插值曲面函数    */
/*-----------------------------------*/
const spline1dinterpolant& Interpolant_Yongmodulus();
double Youngmodulus_vs_T(const double &Temperature);
double Youngmodulus_In_T(const double &Temperature);
const spline1dinterpolant& Interpolant_PoissonRatio();
double PoissonRatio_vs_T(const double &Temperature);
double PoissonRatio_In_T(const double &Temperature);
const spline1dinterpolant& Interpolant_Expantion();
double Expansion_In_T(const double &Temperature);
double T300K_TensileStress_vs_Strain(const double &Strain);
const spline2dinterpolant& Interpolant_TensileStress();
double TensileStress(const double &T, const double &Strain);
double TensileStress_vs_Strain(const double &T, const double &Strain);
double TensileStress_vs_Temper(const double &T, const double &Strain);
Matrix<double, 3, 3> De_In_T(const double &T);
#endif //!_SUS304FUNCTION_


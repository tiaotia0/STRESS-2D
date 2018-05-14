#ifndef _STRESS2D_
#define _STRESS2D_
#include<stdafx.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<interpolation.h>
#include<Eigen/Dense>

using namespace Eigen;
using namespace alglib;
struct Dep_and_ST
{
	//由于计算Dep与温度变化引起的应力增量StressInrease_T(ST)式子几乎一样，故可同时计算，封装于此结构体内
	Matrix<double, 3, 3> Dep;
	Matrix<double, 3, 1> StressIncrease_T;
};
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
Matrix<double, 3, 1> Stress_Derivative(const Matrix<double, 3, 1> &stress);
double Equal_PlasticStrain(const Matrix<double, 3, 1> &Stress, const Matrix<double, 3, 1> &Strain, const double &Temperature);
double Equal_Stress(const Matrix<double, 3, 1> &Stress);
Dep_and_ST Dep_and_StressInrease_T(const Matrix<double, 3, 1> &Stress, const Matrix<double, 3, 1> &Strain,const double &Temperature,const double &dt);
Matrix<double, 3, 1> StrainIncrease_Temperature(const Matrix<double, 3, 1> &stress, const double &Temperature,const double &dt);
Matrix<double, 3, 1> Stress_Process(const Matrix<double, 3, 1> &sumStress, const Matrix<double, 3, 1> &sumStrain, const Matrix<double, 3, 1> &dStrain, const double &Temperature,const double &dt);
#endif // !_STRESS_


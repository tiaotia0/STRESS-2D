#include"stressfunction.h"

/*-----------------------------------*/
/*    Euqal_Strain_Ele_Property类    */
/*-----------------------------------*/
double Euqal_Strain_Ele_Property::eq_plasticstrain = 0;
double Euqal_Strain_Ele_Property::preTemperature = 0;
Matrix<double, 3, 1> Euqal_Strain_Ele_Property::preStress = { 0,0,0 };
Matrix<double, 3, 1> Euqal_Strain_Ele_Property::preStrain = { 0,0,0 };
void Euqal_Strain_Ele_Property::Change_eq_plasticstrain(const double& c)
{
	eq_plasticstrain = c;
}
void Euqal_Strain_Ele_Property::Change_preTemperatrue(const double& c)
{
	preTemperature = c;
}
void Euqal_Strain_Ele_Property::Change_preStress(const Matrix<double, 3, 1>& c)
{
	preStress = c;
}
void Euqal_Strain_Ele_Property::Change_preStrain(const Matrix<double, 3, 1>& c)
{
	preStrain = c;
}
double Euqal_Strain_Ele_Property::Query_eq_plasticstrain()
{
	return eq_plasticstrain;
}
double Euqal_Strain_Ele_Property::Query_preTemperatrue()
{
	return preTemperature;
}
Matrix<double, 3, 1> Euqal_Strain_Ele_Property::Query_preStress()
{
	return preStress;
}
Matrix<double, 3, 1> Euqal_Strain_Ele_Property::Query_preStrain()
{
	return preStrain;
}
void Euqal_Strain_Ele_Property::Reset_all()
{
	eq_plasticstrain = 0;
	preTemperature = 0;
	preStress = { 0,0,0 };
	preStrain = { 0,0,0 };
}

/*-----------------------------------*/
/*          应力相关计算函数         */
/*-----------------------------------*/
Matrix<double, 3, 1> Stress_Derivative(const Matrix<double, 3, 1> &stress)
{
	//计算等效应力对应力矢量的偏导
	Matrix<double, 3, 1> result;
	result.setZero();
	double a1 = stress(0, 0), a2 = stress(1, 0), a3 = stress(2, 0);
	double amid = (a1 + a2) / 3;
	double s1 = a1 - amid, s2 = a2 - amid, s3 = a3;
	double equal_stress = Equal_Stress(stress);
	if (equal_stress == 0) return result; //直接返回0矩阵
	else
	{
		result(0, 0) = s1; result(1, 0) = s2; result(2, 0) = 2 * s3;
		result *= (1.5 / equal_stress);
		return result;
	}
}
double Equal_PlasticStrain(const Matrix<double, 3, 1> &Stress, const Matrix<double, 3, 1> &Strain, const double &Temperature)
{
	//计算塑性等效应变，Strain为全应变,若新的应力会产生新的塑性变形，则更新并返回，否则返回原值
	//staus为试算与否，若为试算，则不更新，若不为试算，则更新。
	double eq_plasticstrain = Euqal_Strain_Ele_Property::Query_eq_plasticstrain(),
			preTemperature = Euqal_Strain_Ele_Property::Query_preTemperatrue();
	Matrix<double, 3, 1> preStress = Euqal_Strain_Ele_Property::Query_preStress(),
						preStrain = Euqal_Strain_Ele_Property::Query_preStrain();
	double equal_stress = Equal_Stress(Stress);
	double y_stress = TensileStress(Temperature, eq_plasticstrain); //T温度下,等效塑性应变为eq_plasticstrain时的屈服应力
	if (equal_stress> y_stress)  //大于屈服应力时才会计算新等效塑性应变
	{
		Matrix<double, 3, 3> De = De_In_T(Temperature), pre_De = De_In_T(preTemperature);
		Matrix<double, 3, 1> eStrain, pre_eStrain, d_pStrain;
		eStrain = De.inverse()*Stress;
		pre_eStrain = pre_De.inverse()*preStress;
		double expansion_value = Expansion_In_T(Temperature),
			pre_expansion_value = Expansion_In_T(preTemperature);
		Matrix<double, 3, 1> expansion(expansion_value, expansion_value, 0),
			pre_expansion(pre_expansion_value, pre_expansion_value, 0);
		d_pStrain = (Strain - preStrain) - (eStrain - pre_eStrain) - (expansion - pre_expansion);
		double dp1, dp2, dp3;
		dp1 = d_pStrain(0, 0); dp2 = d_pStrain(1, 0); dp3 = d_pStrain(2, 0);
		//是否需计算Z向自由变形的塑性应变？
		//double E = Youngmodulus_In_T(Temperature);
		//double Zpstrain = -0.5 *(Stress(0, 0) + Stress(1, 0)) / E;
		double d_eq_plasticstrain = sqrt(2) / sqrt(3)*sqrt(dp1*dp1 + dp2*dp2 + 0.5*dp3*dp3);
		eq_plasticstrain += d_eq_plasticstrain;
		//更新eq_plasticstrain、preStress及preStrain，preTemprature
		Euqal_Strain_Ele_Property::Change_eq_plasticstrain(eq_plasticstrain);
		Euqal_Strain_Ele_Property::Change_preStress(Stress);
		Euqal_Strain_Ele_Property::Change_preStrain(Strain);
		Euqal_Strain_Ele_Property::Change_preTemperatrue(Temperature);
		return eq_plasticstrain;
	}
	else  //若等效应力未增加，则不重新计算，返回上一个状态的等效塑性应变
	{
		//更新eq_plasticstrain、preStress及preStrain，preTemprature
		Euqal_Strain_Ele_Property::Change_preStress(Stress);
		Euqal_Strain_Ele_Property::Change_preStrain(Strain);
		Euqal_Strain_Ele_Property::Change_preTemperatrue(Temperature);
		return eq_plasticstrain;
	}
}
double Equal_Stress(const Matrix<double, 3, 1> &Stress)
{
	//计算等效应力
	double s1 = Stress(0, 0), s2 = Stress(1, 0), s3 = Stress(2, 0);
	return sqrt(s1*s1 + s2*s2 - s1*s2 + 3 * s3*s3);
}
Dep_and_ST Dep_and_StressInrease_T(const Matrix<double, 3, 1> &Stress, const Matrix<double, 3, 1> &Strain, const double &Temperature, const double &dt)
{
	//由于计算Dp与温度变化引起的应力增量StressInrease_T(ST)式子几乎一样，故可同时计算，封装于结构体内Dep_and_ST内
	//Stress为当前应力状态，Strain为当前___全应变___ ,dt为温度增量
	double stress_vs_strain, stress_vs_temper, equal_plasticstrain;
	Matrix<double, 3, 1> ST, stress_derivative;//stress_derivative为应力偏导
	Matrix<double, 3, 3> De, Dp, Dep;
	Dep_and_ST result;
	//先计算公共量De,stress_vs_strain(屈服应力对等效塑性应变导数)、stress_vs_Temper(屈服应力对温度导数)、stress_derivative（等效应力对应力矢量偏导）
	//De
	De = De_In_T(Temperature);
	//stress_derivative
	stress_derivative = Stress_Derivative(Stress);
	//stress_vs_strain
	equal_plasticstrain = Equal_PlasticStrain(Stress, Strain, Temperature); //塑性等效应变
	stress_vs_strain = TensileStress_vs_Strain(Temperature, equal_plasticstrain);
	//stress_vs_Temper
	stress_vs_temper = TensileStress_vs_Temper(Temperature, equal_plasticstrain);
	//判断是否进入屈服阶段
	if (Equal_Stress(Stress) >= TensileStress(Temperature, equal_plasticstrain))
	{
		//屈服阶段(包括强化)
		auto temp = stress_derivative.transpose()*De*stress_derivative;
		Dp = (De*stress_derivative*stress_derivative.transpose()*De) / (stress_vs_strain + temp(0, 0));
		ST = (De*stress_derivative*stress_vs_temper*dt) / (stress_vs_strain + temp(0, 0));
		Dep = De - Dp;
	}
	else
	{
		//弹性阶段,(包括强化)
		Dep = De;
		ST.setZero();
	}
	result.Dep = Dep;
	result.StressIncrease_T = ST;
	return result;
}
Matrix<double, 3, 1> StrainIncrease_Temperature(const Matrix<double, 3, 1> &stress, const double &Temperature, const double &dt)
{
	//计算T温度下由温度变化引起的应变增量,注：此处由于dt可能非常大，会造成误差
	//注：此处不包含热膨胀增量
	Matrix<double, 3, 1> result;
	double E = Youngmodulus_In_T(Temperature);
	double Et = Youngmodulus_vs_T(Temperature);
	double P = PoissonRatio_In_T(Temperature);
	double Pt = PoissonRatio_vs_T(Temperature);
	//De的逆矩阵对T求导
	double temp1 = -Et;
	double temp2 = P*Et - Pt*E;
	double temp3 = 2 * (Pt*E - (1 + P)*Et);
	Matrix < double, 3, 3> TEMP;
	TEMP.setZero();
	TEMP(0, 0) = temp1;	TEMP(1, 1) = temp1; TEMP(2, 2) = temp3;
	TEMP(0, 1) = temp2;	TEMP(1, 0) = temp2;
	result = TEMP*stress*dt / (E*E);
	return result;
}
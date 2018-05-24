#include"stress2d.h"
#include"stress2dproperty.h"
#include<limits>
const spline1dinterpolant& Interpolant_Yongmodulus()
{
	//插值杨氏模量-T的曲线
	static bool Is_caculated = false;
	static spline1dinterpolant s;
	if (Is_caculated == false)  //只需要第一次计算保存，后边查询即可
	{
		extern const std::vector<double> YoungModulus;
		extern const std::vector<double> TemperInPandY;
		auto i = YoungModulus.size(), j = TemperInPandY.size();
		double *p1 = new double[i];
		double *p2 = new double[j];
		memcpy(p1, &YoungModulus[0], i * sizeof(double));
		memcpy(p2, &TemperInPandY[0], j * sizeof(double));
		real_1d_array y, T;
		y.setcontent(i, p1);
		T.setcontent(j, p2);
		spline1dbuildcubic(T, y, s);
		Is_caculated = true;
		delete p1;
		delete p2;
	}
	return s;
}
double Youngmodulus_vs_T(const double &Temperature)
{
	//计算杨氏模量对温度的导数
	double r,ds,d2s; //r是T温度下的杨氏模量，ds是一阶导数，d2s是二阶导
	if (Temperature <= 300) return 0; //假定小于300K时杨氏模量无变化
	else 
	{
		const spline1dinterpolant s(Interpolant_Yongmodulus());
		spline1ddiff(s, Temperature, r, ds, d2s);
		return ds;
	}
}
double Youngmodulus_In_T(const double &Temperature)
{
	//插值得到T温度下的杨氏模量
	const spline1dinterpolant s(Interpolant_Yongmodulus());
	if (Temperature <= 300) return spline1dcalc(s, 300); //假定小于300K时杨氏模量为300K时的值
	else return spline1dcalc(s, Temperature);
}
const spline1dinterpolant& Interpolant_PoissonRatio()
{
	//插值得到T温度下的泊松比-T曲线
	static bool Is_caculated = false;
	static spline1dinterpolant s;
	if (Is_caculated == false)  //只需要第一次计算保存，后边查询即可
	{
		extern const std::vector<double> PoissonRatio;
		extern const std::vector<double> TemperInPandY;
		auto i = PoissonRatio.size(), j = TemperInPandY.size();
		double *p1 = new double[i];
		double *p2 = new double[j];
		memcpy(p1, &PoissonRatio[0], i * sizeof(double));
		memcpy(p2, &TemperInPandY[0], j * sizeof(double));
		real_1d_array y, T;
		y.setcontent(i, p1);
		T.setcontent(j, p2);
		spline1dbuildcubic(T, y, s);
		Is_caculated = true;
		delete p1;
		delete p2;
	}
	return s;
}
double PoissonRatio_vs_T(const double &Temperature)
{
	//计算泊松比对温度的导数
	double r, ds, d2s; //r是T温度下的泊松比，ds是一阶导数，d2s是二阶导
	if (Temperature <= 300) return 0;  //假定小于300K时泊松比无变化
	else
	{
		const spline1dinterpolant s(Interpolant_PoissonRatio());
		spline1ddiff(s, Temperature, r, ds, d2s);
		return ds;
	}
}
double PoissonRatio_In_T(const double &Temperature)
{
	//插值得到T温度下的泊松比
	const spline1dinterpolant s(Interpolant_PoissonRatio());
	if (Temperature<=300) return spline1dcalc(s, 300); //假定小于300K时泊松比为300K时的值
	else return spline1dcalc(s, Temperature);
}
const spline1dinterpolant& Interpolant_Expantion()
{
	//插值热膨胀曲线
	static bool Is_caculated = false;
	static spline1dinterpolant s;
	if (Is_caculated == false)  //只需要第一次计算保存，后边查询即可
	{
		extern const std::vector<double> Expansion;
		extern const std::vector<double> TemperInExpan;
		auto i = Expansion.size(), j = TemperInExpan.size();
		double *p1 = new double[i];
		double *p2 = new double[j];
		memcpy(p1, &Expansion[0], i * sizeof(double));
		memcpy(p2, &TemperInExpan[0], j * sizeof(double));
		real_1d_array y, T;
		y.setcontent(i, p1);
		T.setcontent(j, p2);
		spline1dbuildcubic(T, y, s);
		Is_caculated = true;
	}
	return s;
}
double Expansion_In_T(const double &Temperature)
{
	//插值T温度下的膨胀量
	const spline1dinterpolant s(Interpolant_Expantion());
	if (Temperature <= 300) return 0;  //基准zero=300k
	else return spline1dcalc(s, Temperature)*(Temperature-300);
}
double T300K_TensileStress_vs_Strain(const double &Strain)
{
	//因为插值的屈服应力（关于T和等效塑性应变的）曲面在T<=300K时无导数，故设定当T<=300K时，屈服应力对等效塑性应变偏导等于300K时的屈服应力对塑性等效应变的偏导，屈服应力对T导数为0
	//当等效塑性应变为0时，设定屈服应力对等效塑性应变导数为无穷
	//注：此处是对300K时屈服应力对塑性等效应变求导，Strain是等效塑性应变
	static bool Is_caculated = false;
	static spline1dinterpolant s;
	//插值得到300k时屈服应力-塑性应变曲线
	if (Is_caculated == false)  //只需要第一次计算保存，后边查询即可
	{
		extern const std::vector<double> STRESS300K;
		extern const std::vector<double> STRAIN;
		auto i = STRESS300K.size(), j = STRAIN.size();
		double *p1 = new double[i];
		double *p2 = new double[j];
		memcpy(p1, &STRESS300K[0], i * sizeof(double));
		memcpy(p2, &STRAIN[0], j * sizeof(double));
		real_1d_array y, S; //S是塑性应变
		y.setcontent(i, p1);
		S.setcontent(j, p2);
		spline1dbuildcubic(S, y, s);  //splinecubic插值
		Is_caculated = true;
		delete p1;
		delete p2;
	}
	double r, ds, d2s; //r是S时的屈服应力值，ds是一阶导数，d2s是二阶导
	if (Strain < 0) std::cout << "erro: equal plastic strain less than zero! (in function T300K_TensileStress_vs_Strain)" << std::endl;
	if (Strain == 0) return std::numeric_limits<double>::max(); //Strain是等效塑性应变
	else
	{
		spline1ddiff(s, Strain, r, ds, d2s);
		return ds;
	}
}
const spline2dinterpolant& Interpolant_TensileStress()
{
	//插值 温度与塑性应变为变量的实验应力曲面，该应力曲面即为屈服应力曲面
	static bool Is_caculated = false;
	static spline2dinterpolant s;
	if (Is_caculated == false)  //只需要第一次计算保存，后边查询即可
	{
		extern const std::vector<double> TEMPERATURE;
		extern const std::vector<double> STRAIN;
		extern const std::vector<double> STRESS;
		auto i = TEMPERATURE.size(), j = STRAIN.size(), k = STRESS.size();
		double *pi = new double[i];
		double *pj = new double[j];
		double *pk = new double[k];
		memcpy(pi, &TEMPERATURE[0], i * sizeof(double));
		memcpy(pj, &STRAIN[0], j * sizeof(double));
		memcpy(pk, &STRESS[0], k * sizeof(double));
		real_1d_array T,strain,stress;
		T.setcontent(i, pi);
		strain.setcontent(j, pj);
		stress.setcontent(k, pk);
		spline2dbuildbicubicv(strain, j, T, i, stress, 1, s);
		Is_caculated = true;
		delete pi;
		delete pj;
		delete pk;
	}
	return s;
}
double TensileStress(const double &T, const double &Strain)
{
	//计算T温度下等效塑性应变为Strain时的屈服应力
	if (Strain < 0) std::cout << "erro:equal plastic strain less than zero!(in function TensileStress)" << std::endl;
	const spline2dinterpolant s(Interpolant_TensileStress());
	//传入的塑性等效应变总大于等于0,且把温度低于300k时的状态都认为是300k时状态
	if (T <= 300) return spline2dcalc(s, Strain, 300);
	else return spline2dcalc(s, Strain, T);
}
double TensileStress_vs_Strain(const double &T, const double &Strain)
{
	//计算T温度下等效塑性应变为Strain时的实验应力（屈服应力）对于塑性应变的偏导。
	//因为插值的屈服应力曲面在T<=300K时无导数，故设定当T<=300K时，屈服应力对等效塑性应变等于300K时的屈服应力对塑性等效应变的偏导，屈服应力对T导数为0。
	//当等效塑性应变为0时，设定屈服应力对等效塑性应变导数为无穷。
	if (Strain < 0) std::cout << "erro:equal plastic strain less than zero!(in function TensileStress_vs_Strain)" << std::endl;
	if (Strain == 0) { return std::numeric_limits<double>::max(); } //未屈服时，导数为无穷
	else
	{
		if (T <= 300) return T300K_TensileStress_vs_Strain(Strain);
		else {
			double Dtemper, Dstrain, TensileStress, Dtemper_and_strain;
			const spline2dinterpolant s(Interpolant_TensileStress());
			spline2ddiff(s, Strain, T, TensileStress, Dstrain, Dtemper, Dtemper_and_strain);
			return Dstrain;
		}
	}
}
double TensileStress_vs_Temper(const double &T, const double &Strain)
{
	//计算T温度下 /等效塑性应变/ 为Strain时的  /实验应力（屈服应力）/  对于温度的偏导
	//因为插值的屈服应力于T和等效塑性应变的曲面在T<=300K时无导数，故设定当T<=300K时，屈服应力对等效塑性应变等于300K时的屈服应力对塑性等效应变的偏导，屈服应力对T导数为0
	//当等效塑性应变为0时，设定屈服应力对等效塑性应变导数为无穷
	if (Strain < 0) std::cout << "erro:equal plastic strain less than zero!(in function TensileStress_vs_Temper)" << std::endl;
	if (T < 300) return 0;
	else
	{
		double Dtemper, Dstrain, TensileStress, Dtemper_and_strain;
		const spline2dinterpolant s(Interpolant_TensileStress());
		spline2ddiff(s, Strain, T, TensileStress, Dstrain, Dtemper, Dtemper_and_strain);
		return Dtemper;
	}
}
Matrix<double, 3, 3> De_In_T(const double &T)
{
	//计算T温度时的弹性矩阵De，u为泊松比，E为杨氏模量
	double E, u;
	E = Youngmodulus_In_T(T);
	u = PoissonRatio_In_T(T);
	Matrix<double, 3, 3> De;
	De.setZero();
	double temp = E/(1-u*u);
	De(0, 0) = 1; De(1, 1) = 1; De(2, 2) = (1-u)/2;
	De(0, 1) = u; De(1, 0) = u;
	return De*temp;
}
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
		result(0, 0) = s1; result(1, 0) = s2; result(2, 0) = 2*s3;
		result *= (1.5 / equal_stress);
		return result;
	}
}
double Equal_PlasticStrain(const Matrix<double, 3, 1> &Stress, const Matrix<double, 3, 1> &Strain, const double &Temperature)
{
	//计算塑性等效应变，Strain为全应变,若新的应力会产生新的塑性变形，则更新并返回，否则返回原值
	//staus为试算与否，若为试算，则不更新，若不为试算，则更新。
	static double eq_plasticstrain = 0, preTemperature = 0;
	static Matrix<double, 3, 1> preStress = { 0, 0, 0 }, preStrain = { 0, 0, 0 };
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
		//更新preStress及preStrain，preTemprature
		preStress = Stress;
		preStrain = Strain; std::cout << 282 << std::endl;
		preTemperature = Temperature; std::cout << 283 << std::endl;
		return eq_plasticstrain;
	}
	else  //若等效应力未增加，则不重新计算，返回上一个状态的等效塑性应变
	{
		//更新preStress及preStrain，preTemprature
		preStress = Stress;
		preStrain = Strain;
		preTemperature = Temperature;
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
	//Stress为当前应力状态，Strain为当前//全应变// ,dt为温度增量
	Matrix<double, 3, 3> De,Dp,Dep;
	Matrix<double, 3, 1> ST,stress_derivative;//stress_derivative为应力偏导
	double stress_vs_strain, stress_vs_temper, equal_plasticstrain;
	Dep_and_ST result;
	//先计算公共量De,stress_vs_strain(屈服应力对等效塑性应变导数),stress_derivative（等效应力对应力矢量偏导）
	//De
	De = De_In_T(Temperature);
	//stress_vs_strain
	equal_plasticstrain = Equal_PlasticStrain(Stress, Strain, Temperature); //计算塑性等效应变
	if(Equal_Stress(Stress)>=TensileStress(Temperature,equal_plasticstrain))
		stress_vs_strain = TensileStress_vs_Strain(Temperature, equal_plasticstrain);  //屈服应力对塑性等效应变偏导
	else stress_vs_strain = std::numeric_limits<double>::max(); //屈服后卸载，再加载时按弹性加载，则此时屈服应力对塑性应变偏导为无穷
	//stress_derivative
	stress_derivative = Stress_Derivative(Stress);
	//计算Dep与计算ST
	if (stress_vs_strain==std::numeric_limits<double>::max())
	{
		//弹性阶段，此时Dp为0
		Dep=De;
		//由于stress_vs_strain=无穷，stress_vs_temper为小于0的一个有极限的数，此时ST为0
		ST.setZero();
	}
	else 
	{
		auto temp = stress_derivative.transpose()*De*stress_derivative;
		stress_vs_temper = TensileStress_vs_Temper(Temperature, equal_plasticstrain);
		Dp = (De*stress_derivative*stress_derivative.transpose()*De) / (stress_vs_strain + temp(0, 0));
		ST = (De*stress_derivative*stress_vs_temper*dt) / (stress_vs_strain + temp(0, 0));
		Dep = De - Dp;
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
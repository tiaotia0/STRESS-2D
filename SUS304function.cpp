#include"SUS304function.h"
/*-----------------------------------*/
/*       SUS304性能及插值曲面函数    */
/*-----------------------------------*/
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
	double r, ds, d2s; //r是T温度下的杨氏模量，ds是一阶导数，d2s是二阶导
	if (Temperature < 300) return 0; //假定小于300K时杨氏模量无变化
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
	if (Temperature < 300) return 0;  //假定小于300K时泊松比无变化
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
	if (Temperature <= 300) return spline1dcalc(s, 300); //假定小于300K时泊松比为300K时的值
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
	else return spline1dcalc(s, Temperature)*(Temperature - 300);
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
		real_1d_array T, strain, stress;
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
	if (Strain == 0) return 0;  //未屈服时，屈服应力对温度导数为0
	else
	{
		if (T < 300) return 0;
		else
		{
			double Dtemper, Dstrain, TensileStress, Dtemper_and_strain;
			const spline2dinterpolant s(Interpolant_TensileStress());
			spline2ddiff(s, Strain, T, TensileStress, Dstrain, Dtemper, Dtemper_and_strain);
			return Dtemper;
		}
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
	double temp = E / (1 - u*u);
	De(0, 0) = 1; De(1, 1) = 1; De(2, 2) = (1 - u) / 2;
	De(0, 1) = u; De(1, 0) = u;
	return De*temp;
}
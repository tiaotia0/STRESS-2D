#include"stress2d.h"
#include"stress2dproperty.h"
#include<limits>
const spline1dinterpolant& Interpolant_Yongmodulus()
{
	//��ֵ����ģ��-T������
	static bool Is_caculated = false;
	static spline1dinterpolant s;
	if (Is_caculated == false)  //ֻ��Ҫ��һ�μ��㱣�棬��߲�ѯ����
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
	//��������ģ�����¶ȵĵ���
	double r,ds,d2s; //r��T�¶��µ�����ģ����ds��һ�׵�����d2s�Ƕ��׵�
	if (Temperature <= 300) return 0; //�ٶ�С��300Kʱ����ģ���ޱ仯
	else 
	{
		const spline1dinterpolant s(Interpolant_Yongmodulus());
		spline1ddiff(s, Temperature, r, ds, d2s);
		return ds;
	}
}
double Youngmodulus_In_T(const double &Temperature)
{
	//��ֵ�õ�T�¶��µ�����ģ��
	const spline1dinterpolant s(Interpolant_Yongmodulus());
	if (Temperature <= 300) return spline1dcalc(s, 300); //�ٶ�С��300Kʱ����ģ��Ϊ300Kʱ��ֵ
	else return spline1dcalc(s, Temperature);
}
const spline1dinterpolant& Interpolant_PoissonRatio()
{
	//��ֵ�õ�T�¶��µĲ��ɱ�-T����
	static bool Is_caculated = false;
	static spline1dinterpolant s;
	if (Is_caculated == false)  //ֻ��Ҫ��һ�μ��㱣�棬��߲�ѯ����
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
	//���㲴�ɱȶ��¶ȵĵ���
	double r, ds, d2s; //r��T�¶��µĲ��ɱȣ�ds��һ�׵�����d2s�Ƕ��׵�
	if (Temperature <= 300) return 0;  //�ٶ�С��300Kʱ���ɱ��ޱ仯
	else
	{
		const spline1dinterpolant s(Interpolant_PoissonRatio());
		spline1ddiff(s, Temperature, r, ds, d2s);
		return ds;
	}
}
double PoissonRatio_In_T(const double &Temperature)
{
	//��ֵ�õ�T�¶��µĲ��ɱ�
	const spline1dinterpolant s(Interpolant_PoissonRatio());
	if (Temperature<=300) return spline1dcalc(s, 300); //�ٶ�С��300Kʱ���ɱ�Ϊ300Kʱ��ֵ
	else return spline1dcalc(s, Temperature);
}
const spline1dinterpolant& Interpolant_Expantion()
{
	//��ֵ����������
	static bool Is_caculated = false;
	static spline1dinterpolant s;
	if (Is_caculated == false)  //ֻ��Ҫ��һ�μ��㱣�棬��߲�ѯ����
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
	//��ֵT�¶��µ�������
	const spline1dinterpolant s(Interpolant_Expantion());
	if (Temperature <= 300) return 0;  //��׼zero=300k
	else return spline1dcalc(s, Temperature)*(Temperature-300);
}
double T300K_TensileStress_vs_Strain(const double &Strain)
{
	//��Ϊ��ֵ������Ӧ��������T�͵�Ч����Ӧ��ģ�������T<=300Kʱ�޵��������趨��T<=300Kʱ������Ӧ���Ե�Ч����Ӧ��ƫ������300Kʱ������Ӧ�������Ե�ЧӦ���ƫ��������Ӧ����T����Ϊ0
	//����Ч����Ӧ��Ϊ0ʱ���趨����Ӧ���Ե�Ч����Ӧ�䵼��Ϊ����
	//ע���˴��Ƕ�300Kʱ����Ӧ�������Ե�ЧӦ���󵼣�Strain�ǵ�Ч����Ӧ��
	static bool Is_caculated = false;
	static spline1dinterpolant s;
	//��ֵ�õ�300kʱ����Ӧ��-����Ӧ������
	if (Is_caculated == false)  //ֻ��Ҫ��һ�μ��㱣�棬��߲�ѯ����
	{
		extern const std::vector<double> STRESS300K;
		extern const std::vector<double> STRAIN;
		auto i = STRESS300K.size(), j = STRAIN.size();
		double *p1 = new double[i];
		double *p2 = new double[j];
		memcpy(p1, &STRESS300K[0], i * sizeof(double));
		memcpy(p2, &STRAIN[0], j * sizeof(double));
		real_1d_array y, S; //S������Ӧ��
		y.setcontent(i, p1);
		S.setcontent(j, p2);
		spline1dbuildcubic(S, y, s);  //splinecubic��ֵ
		Is_caculated = true;
		delete p1;
		delete p2;
	}
	double r, ds, d2s; //r��Sʱ������Ӧ��ֵ��ds��һ�׵�����d2s�Ƕ��׵�
	if (Strain < 0) std::cout << "erro: equal plastic strain less than zero! (in function T300K_TensileStress_vs_Strain)" << std::endl;
	if (Strain == 0) return std::numeric_limits<double>::max(); //Strain�ǵ�Ч����Ӧ��
	else
	{
		spline1ddiff(s, Strain, r, ds, d2s);
		return ds;
	}
}
const spline2dinterpolant& Interpolant_TensileStress()
{
	//��ֵ �¶�������Ӧ��Ϊ������ʵ��Ӧ�����棬��Ӧ�����漴Ϊ����Ӧ������
	static bool Is_caculated = false;
	static spline2dinterpolant s;
	if (Is_caculated == false)  //ֻ��Ҫ��һ�μ��㱣�棬��߲�ѯ����
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
	//����T�¶��µ�Ч����Ӧ��ΪStrainʱ������Ӧ��
	if (Strain < 0) std::cout << "erro:equal plastic strain less than zero!(in function TensileStress)" << std::endl;
	const spline2dinterpolant s(Interpolant_TensileStress());
	//��������Ե�ЧӦ���ܴ��ڵ���0,�Ұ��¶ȵ���300kʱ��״̬����Ϊ��300kʱ״̬
	if (T <= 300) return spline2dcalc(s, Strain, 300);
	else return spline2dcalc(s, Strain, T);
}
double TensileStress_vs_Strain(const double &T, const double &Strain)
{
	//����T�¶��µ�Ч����Ӧ��ΪStrainʱ��ʵ��Ӧ��������Ӧ������������Ӧ���ƫ����
	//��Ϊ��ֵ������Ӧ��������T<=300Kʱ�޵��������趨��T<=300Kʱ������Ӧ���Ե�Ч����Ӧ�����300Kʱ������Ӧ�������Ե�ЧӦ���ƫ��������Ӧ����T����Ϊ0��
	//����Ч����Ӧ��Ϊ0ʱ���趨����Ӧ���Ե�Ч����Ӧ�䵼��Ϊ���
	if (Strain < 0) std::cout << "erro:equal plastic strain less than zero!(in function TensileStress_vs_Strain)" << std::endl;
	if (Strain == 0) { return std::numeric_limits<double>::max(); } //δ����ʱ������Ϊ����
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
	//����T�¶��� /��Ч����Ӧ��/ ΪStrainʱ��  /ʵ��Ӧ��������Ӧ����/  �����¶ȵ�ƫ��
	//��Ϊ��ֵ������Ӧ����T�͵�Ч����Ӧ���������T<=300Kʱ�޵��������趨��T<=300Kʱ������Ӧ���Ե�Ч����Ӧ�����300Kʱ������Ӧ�������Ե�ЧӦ���ƫ��������Ӧ����T����Ϊ0
	//����Ч����Ӧ��Ϊ0ʱ���趨����Ӧ���Ե�Ч����Ӧ�䵼��Ϊ����
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
	//����T�¶�ʱ�ĵ��Ծ���De��uΪ���ɱȣ�EΪ����ģ��
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
	//�����ЧӦ����Ӧ��ʸ����ƫ��
	Matrix<double, 3, 1> result;
	result.setZero();
	double a1 = stress(0, 0), a2 = stress(1, 0), a3 = stress(2, 0);
	double amid = (a1 + a2) / 3;
	double s1 = a1 - amid, s2 = a2 - amid, s3 = a3;
	double equal_stress = Equal_Stress(stress);
	if (equal_stress == 0) return result; //ֱ�ӷ���0����
	else 
	{
		result(0, 0) = s1; result(1, 0) = s2; result(2, 0) = 2*s3;
		result *= (1.5 / equal_stress);
		return result;
	}
}
double Equal_PlasticStrain(const Matrix<double, 3, 1> &Stress, const Matrix<double, 3, 1> &Strain, const double &Temperature)
{
	//�������Ե�ЧӦ�䣬StrainΪȫӦ��,���µ�Ӧ��������µ����Ա��Σ�����²����أ����򷵻�ԭֵ
	//stausΪ���������Ϊ���㣬�򲻸��£�����Ϊ���㣬����¡�
	static double eq_plasticstrain = 0, preTemperature = 0;
	static Matrix<double, 3, 1> preStress = { 0, 0, 0 }, preStrain = { 0, 0, 0 };
	double equal_stress = Equal_Stress(Stress);
	double y_stress = TensileStress(Temperature, eq_plasticstrain); //T�¶���,��Ч����Ӧ��Ϊeq_plasticstrainʱ������Ӧ��
	if (equal_stress> y_stress)  //��������Ӧ��ʱ�Ż�����µ�Ч����Ӧ��
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
		//�Ƿ������Z�����ɱ��ε�����Ӧ�䣿
		//double E = Youngmodulus_In_T(Temperature);
		//double Zpstrain = -0.5 *(Stress(0, 0) + Stress(1, 0)) / E;
		double d_eq_plasticstrain = sqrt(2) / sqrt(3)*sqrt(dp1*dp1 + dp2*dp2 + 0.5*dp3*dp3);
		eq_plasticstrain += d_eq_plasticstrain;
		//����preStress��preStrain��preTemprature
		preStress = Stress;
		preStrain = Strain; std::cout << 282 << std::endl;
		preTemperature = Temperature; std::cout << 283 << std::endl;
		return eq_plasticstrain;
	}
	else  //����ЧӦ��δ���ӣ������¼��㣬������һ��״̬�ĵ�Ч����Ӧ��
	{
		//����preStress��preStrain��preTemprature
		preStress = Stress;
		preStrain = Strain;
		preTemperature = Temperature;
		return eq_plasticstrain;
	} 
}
double Equal_Stress(const Matrix<double, 3, 1> &Stress)
{
	//�����ЧӦ��
	double s1 = Stress(0, 0), s2 = Stress(1, 0), s3 = Stress(2, 0);
	return sqrt(s1*s1 + s2*s2 - s1*s2 + 3 * s3*s3);
}
Dep_and_ST Dep_and_StressInrease_T(const Matrix<double, 3, 1> &Stress, const Matrix<double, 3, 1> &Strain, const double &Temperature, const double &dt)
{
	//���ڼ���Dp���¶ȱ仯�����Ӧ������StressInrease_T(ST)ʽ�Ӽ���һ�����ʿ�ͬʱ���㣬��װ�ڽṹ����Dep_and_ST��
	//StressΪ��ǰӦ��״̬��StrainΪ��ǰ//ȫӦ��// ,dtΪ�¶�����
	Matrix<double, 3, 3> De,Dp,Dep;
	Matrix<double, 3, 1> ST,stress_derivative;//stress_derivativeΪӦ��ƫ��
	double stress_vs_strain, stress_vs_temper, equal_plasticstrain;
	Dep_and_ST result;
	//�ȼ��㹫����De,stress_vs_strain(����Ӧ���Ե�Ч����Ӧ�䵼��),stress_derivative����ЧӦ����Ӧ��ʸ��ƫ����
	//De
	De = De_In_T(Temperature);
	//stress_vs_strain
	equal_plasticstrain = Equal_PlasticStrain(Stress, Strain, Temperature); //�������Ե�ЧӦ��
	if(Equal_Stress(Stress)>=TensileStress(Temperature,equal_plasticstrain))
		stress_vs_strain = TensileStress_vs_Strain(Temperature, equal_plasticstrain);  //����Ӧ�������Ե�ЧӦ��ƫ��
	else stress_vs_strain = std::numeric_limits<double>::max(); //������ж�أ��ټ���ʱ�����Լ��أ����ʱ����Ӧ��������Ӧ��ƫ��Ϊ����
	//stress_derivative
	stress_derivative = Stress_Derivative(Stress);
	//����Dep�����ST
	if (stress_vs_strain==std::numeric_limits<double>::max())
	{
		//���Խ׶Σ���ʱDpΪ0
		Dep=De;
		//����stress_vs_strain=���stress_vs_temperΪС��0��һ���м��޵�������ʱSTΪ0
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
	//����T�¶������¶ȱ仯�����Ӧ������,ע���˴�����dt���ܷǳ��󣬻�������
	//ע���˴�����������������
	Matrix<double, 3, 1> result;
	double E = Youngmodulus_In_T(Temperature);  
	double Et = Youngmodulus_vs_T(Temperature);
	double P = PoissonRatio_In_T(Temperature);
	double Pt = PoissonRatio_vs_T(Temperature);
	//De��������T��
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
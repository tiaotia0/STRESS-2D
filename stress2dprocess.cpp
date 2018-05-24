#include"stress2d.h"
#include<limits>

Matrix<double, 3, 1> Stress_Process(const Matrix<double, 3, 1> &sumStress, const Matrix<double, 3, 1> &sumStrain, const Matrix<double, 3, 1> &dStrain, const double &Temperature, const double &dt)
{
	//���Ʊ������ȣ���С���ݹ���ü���dStress��
	Matrix<double, 3, 1> dStress;
	dStress.setZero();
	//����¶���Ӧ���Ƿ����仯��������仯�����ֺ���㲢����
	int n = abs(dt) / 20, n0 = abs(dStrain(0, 0) / 1e-5), n1 = abs(dStrain(1, 0) / 1e-5), n2 = abs(dStrain(2, 0) / 1e-5);
	if (n0 > n) n = n0; if (n1 > n) n = n1; if (n2 > n) n = n2;
	if (n != 0)
	{ //����仯
		Matrix<double, 3, 1> strainstep = dStrain / (n + 1);
		double dt_new = dt / (n + 1);
		Matrix<double, 3, 1> sumStress_temp = sumStress, sumStrain_temp = sumStrain;
		double Temperature_temp = Temperature;
		for (int i = 0; i < n + 1; i++)
		{
			Matrix<double, 3, 1> dStress_temp = Stress_Process(sumStress_temp, sumStrain_temp, strainstep, Temperature_temp, dt_new);
			sumStress_temp += dStress_temp;
			sumStrain_temp += strainstep;
			Temperature_temp += dt_new;
			dStress += dStress_temp;
		}
		return dStress;
	}
	else
	{ //С���仯
		double e = Expansion_In_T(Temperature), enext = Expansion_In_T(Temperature + dt);
		Matrix<double, 3, 1> dexpansion(enext - e, enext - e, 0);
		//�Լ�������
		Matrix<double, 3, 1> de0 = StrainIncrease_Temperature(sumStress, Temperature, dt);
		Dep_and_ST Dep_ST = Dep_and_StressInrease_T(sumStress, sumStrain, Temperature, dt);
		dStress = Dep_ST.Dep*(dStrain - de0 - dexpansion) + Dep_ST.StressIncrease_T;
		//������ء�ж��״̬���ж�����
		double s1 = sumStress(0, 0), s2 = sumStress(1, 0), s3 = sumStress(2, 0),
			ds1 = dStress(0, 0), ds2 = dStress(1, 0), ds3 = dStress(2, 0);
		double CRI = (2 * s1*ds1 + 2 * s2*ds2 - s1*ds2 - s2*ds1 + 6 * s3*ds3) / 3;  //dJ2
		if (CRI >= 0)
		{//���ع���
			return dStress;
		}
		else
		{//ж�ع��̣�������������뵥��ж�أ�
			Matrix<double, 3, 3> De = De_In_T(Temperature);
			dStress = De*(dStrain - de0 - dexpansion); 
			double ds1_unload = dStress(0, 0), ds2_unload = dStress(1, 0), ds3_unload = dStress(2, 0);
			double CRI_unload = (2 * s1*ds1_unload + 2 * s2*ds2_unload - s1*ds2_unload - s2*ds1_unload + 6 * s3*ds3_unload) / 3;
			if (CRI_unload<0)   return dStress;//����ж��
			else
			{ //������ع��̣��ٶ���ж����Ӧ����ȫΪ0״̬���Ҵ�ʱ�¶ȱ仯�������̵�һ��
				Matrix<double, 3, 1> zero, d_elastic_strain, mid_Strain, dStress1, dStress2;
				zero.setZero(); d_elastic_strain.setZero();
				dStress1 = zero - sumStress;
				d_elastic_strain = De.inverse()*dStress1;
				mid_Strain = sumStrain + d_elastic_strain + de0 / 2 + dexpansion / 2;
				Matrix<double, 3, 1> dStrain_now = sumStrain + dStrain - mid_Strain;
				dStress2 = Stress_Process(zero, mid_Strain, dStrain_now, Temperature + dt / 2, dt / 2);
				return dStress1 + dStress2;
			}
		}    
	}
}
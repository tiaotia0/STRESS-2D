#include<iostream>
#include<fstream>
#include"stress2d.h"

int main()
{
	std::ifstream filein("C:\\Users\\yl_so\\Desktop\\test.txt"); //�ļ���ʽ��xxӦ�䣬yyӦ�䣬xyӦ�䣬�¶�T����λΪ���϶ȣ�
	std::ofstream fileout("C:\\Users\\yl_so\\Desktop\\out.txt");//�ļ���ʽ��xxӦ����yyӦ����xyӦ��
	fileout << "Stress xx" << '\t' << "Stress yy" << '\t' << "Stress xy" << std::endl<<0<<'\t'<<0<<'\t'<<0<<std::endl;
	//��ʼ״̬����
	Matrix<double, 3, 1> sumStress, sumStrain;
	sumStress.setZero(); sumStrain.setZero();
	double Temperature=300;
	while (1)
	{
	//����Ӧ�����¶�
		double dt,sxx,syy,sxy,T; //sxxΪxx����Ӧ��
		filein >> sxx >> syy >> sxy >> T;
		if (filein.good())
		{
			//����Ӧ��/�¶�����
			Matrix<double, 3, 1> dStrain, dStress;
			dStrain.setZero(); dStress.setZero();
			dStrain(0, 0) = sxx - sumStrain(0, 0);
			dStrain(1, 0) = syy - sumStrain(1, 0);
			dStrain(2, 0) = sxy - sumStrain(2, 0);
			if (T < 27) T = 27;  //TС��27��ʱ�¶�����Ϊ27�ȣ�300K
			dt = T+273 - Temperature;
			//����Ӧ������
			dStress = Stress_Process(sumStress, sumStrain, dStrain, Temperature, dt);
			//������Ӧ��/Ӧ����
			sumStress += dStress;
			sumStrain += dStrain;
			Temperature += dt;
	//���
			fileout << sumStress(0, 0) << '\t' << sumStress(1, 0) << '\t' << sumStress(2, 0) <<std::endl;
		}
		else break;
	}
	filein.close();
	fileout.close();
	return 0;
}

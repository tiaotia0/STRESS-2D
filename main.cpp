#include<iostream>
#include<fstream>
#include"stress2d.h"

int main()
{
	std::ifstream filein("C:\\Users\\yl_so\\Desktop\\test.txt"); //文件格式：xx应变，yy应变，xy应变，温度T（单位为摄氏度）
	std::ofstream fileout("C:\\Users\\yl_so\\Desktop\\out.txt");//文件格式：xx应力，yy应力，xy应力
	fileout << "Stress xx" << '\t' << "Stress yy" << '\t' << "Stress xy" << std::endl<<0<<'\t'<<0<<'\t'<<0<<std::endl;
	//初始状态设置
	Matrix<double, 3, 1> sumStress, sumStrain;
	sumStress.setZero(); sumStrain.setZero();
	double Temperature=300;
	while (1)
	{
	//读入应变与温度
		double dt,sxx,syy,sxy,T; //sxx为xx方向应变
		filein >> sxx >> syy >> sxy >> T;
		if (filein.good())
		{
			//构造应变/温度增量
			Matrix<double, 3, 1> dStrain, dStress;
			dStrain.setZero(); dStress.setZero();
			dStrain(0, 0) = sxx - sumStrain(0, 0);
			dStrain(1, 0) = syy - sumStrain(1, 0);
			dStrain(2, 0) = sxy - sumStrain(2, 0);
			if (T < 27) T = 27;  //T小于27度时温度设置为27度，300K
			dt = T+273 - Temperature;
			//计算应力增量
			dStress = Stress_Process(sumStress, sumStrain, dStrain, Temperature, dt);
			//更新总应变/应力量
			sumStress += dStress;
			sumStrain += dStrain;
			Temperature += dt;
	//输出
			fileout << sumStress(0, 0) << '\t' << sumStress(1, 0) << '\t' << sumStress(2, 0) <<std::endl;
		}
		else break;
	}
	filein.close();
	fileout.close();
	return 0;
}

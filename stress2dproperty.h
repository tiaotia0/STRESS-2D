#ifndef _STRESS2DPROPERTY_
#define _STRESS2DPROPERTY_
#include<vector>
static const std::vector<double> PoissonRatio = { 0.28,0.29,0.3,0.285,0.25,0.25,0.25,0.25,0.25 };
static const std::vector<double> YoungModulus = { 207e3,186e3,167e3,145e3,114e3,62100,27600,10e3,2e3 };
static const std::vector<double> TemperInPandY = { 300,477,700,922,1144,1366,1477,1550,1728 };// 泊松比与杨氏模量对应的温度
static const std::vector<double> Expansion = { 2.080e-5,2.186e-5,2.226e-5,2.151e-5,2.151e-5,2.151e-5,2.048e-5,1.920e-5,1.807e-5 };
static const std::vector<double> TemperInExpan = { 300,600,1000,1300,1672,1728,1800,1900,2000 };
static const std::vector<double> TEMPERATURE = {300,693,923,1089,1366,1550,1728}; //拉伸实验温度
static const std::vector<double> STRAIN = { 0,0.01,0.05,0.1,0.2,0.4,0.8};// 实验塑性应变
static const std::vector<double> STRESS = {
	//实验应力值
	//2d插值格式参见alglib的说明
	//T300
	241.3,289.6,393,482.6,620.5,827.4,1276,
	//T693
	117.2,193.1,275.8,358.5,482.6,620.5,724,
	//T923
	96.5,144.8,213.7,282.7,379.2,448.2,482.6,
	//T1089
	68.9,110.3,144.8,172.4,206.8,241.3,241.3,
	//T1366
	27.6,38.6,41.4,41.4,41.4,41.4,41.4,
	//T1550
	5,6,6.3,6.3,6.3,6.3,6.3,
	//T1728
	1,1,1,1,1,1,1
};
static const std::vector<double> STRESS300K = { 241.3,289.6,393,482.6,620.5,827.4,1276 };//300k时的屈服应力
static const std::vector<double> YEILEDSTRESS = { 241.3,117.2,96.5,68.9,27.6,5,1 };//不同温度下的首次加载屈服应力
#endif // !_STRESS2DPROPERTY_

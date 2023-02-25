#include"Model.h"
#include"Field.h"
#include <math.h> 

#include <time.h>
#include <cstdlib>
#include <random>
#include <chrono>
//double floor(double x); 
//double ceil(double x);

/*---------------------------------------------*/
/*--------------円柱Mie散乱--------------------*/
/*---------------------------------------------*/
FazzyMieModel::FazzyMieModel(Field *f) :
	FazzyModel(f)
{
	ep = 1.6*1.6*EPSILON_0_S;			//誘電率 = (屈折率)^2
	
}

string FazzyMieModel::mkdir(string root) {
	_mkdir((root + "Mie").c_str());

	string name = "Mie/" + to_s((int)(mField->cellToNano(r))) + "nm," + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}

double FazzyMieModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	//r = mField->nanoToCell(60);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;

	double _x = mx - 0.5*mField->getNx();//(N_X/2, N_Y/2)を原点にシフト
	double _y = my - 0.5*mField->getNy();

	//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
	if ((_x*_x) + (_y * _y) >= pow(r + 1, 2.0))

		return EPSILON_0_S;

	//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外
	if ((_x*_x)+ (_y * _y)<= pow(r - 1, 2.0))
		return ep;

	double s = 0;

	double a = 1.0, b = 1.0;
	if (f == D_X) b = 0;
	if (f == D_Y) a = 0;
	for (double i = -16 + 0.5; i < 16; i += 1)
		for (double j = -16 + 0.5; j < 16; j += 1)
			if (pow(_x + a * i / 32.0, 2.0) + pow(_y + b * j / 32.0, 2.0) <= r * r)
				s += 1;
	s /= 32.0*32.0;
	return ep * s + EPSILON_0_S * (1 - s);

}


SapphirinidModelSmooth::SapphirinidModelSmooth(Field* f, double lambda) :FazzyModel(f), epOfCrystal(1.73), epOfCytoplasm(1.40)

{
	setLambda(lambda);
	cout << "ep= : " << endl;
}



string SapphirinidModelSmooth::mkdir(string root) {
	_mkdir((root + "SapphirinidModel").c_str());

	string name = "SapphirinidModel\\" + to_s(int(lambda));
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "\\";
}

double SapphirinidModelSmooth::calcEPS(const double& x, const double& y, enum INTEG f) {
	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;

	double startLine = 525.0;
	double weight_w = calcMultilayer((double)(mField->cellToNano(mx)), (double)(mField->cellToNano(my)), startLine);
	if (weight_w < 0)
		return EPSILON_0_S;
	return weight_w * epOfCrystal*epOfCrystal + (1 - weight_w)*epOfCytoplasm*epOfCytoplasm;//余数为1 就是cystal;余数为0就是plasm,Llength
}

double SapphirinidModelSmooth::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;

	
	double startLine = 0;
	double weight_w = calcMultilayer(mx, my, startLine);
	if (weight_w < 0)
		return 0;
	return weight_w * 0.0 + (1 - weight_w)*0.0;


}
double SapphirinidModelSmooth::calcMultilayer(double x, double y, double start) {
	double layer = selectLayer(y, start);
	if (layer < 0)
		return -1;
	return (double)((int)layer % 2); /*1white,0black*/
}

double SapphirinidModelSmooth::selectLayer(double y, double start) {
	int length = y - start;
	
	double all = 145.0;
	double Llength =30.0;
	double Hlength = 115.0;
	if (length < 0)
		return -1;
	if (length >725.0) //90+36;270+14*9/100is 8layer for 50/70    yellow:285 is 150cy 80cr 13layers       miras 97 8layer 50_70   red 160 70 11*7+16*10    8*11+4*10 180*2+60*3   20*2+7*3 21*2+6*3
		//整个形状高度
		return -1;
	

	if (fmod(abs(length), all) > Llength) // 160/80 //180/60  //200 70 //210/60 单数是Hlength，厚的
		return 3;
	return 4;
}



//(int)(mField->cellToNano(mx)
double * getRandom(int n)
{
	static double  r[100];
	mt19937 gen((unsigned int)time(nullptr)); // 定义随机数生成器对象gen，使用time(nullptr)作为随机数生成器的种子
	uniform_real_distribution<double> dis(100,400);
	
	// 设置种子
	//srand((unsigned)time(NULL));
	for (int i = 0; i < 100; ++i)
	{
		//r[i] = rand() % n + 1;//生成区间r~l的随机数 
		r[i]=dis(gen);
		//cout << dis(gen) << endl;
	}

	return r;
}

double * getRandom2(int n)
{
	static double  r[100];

	srand((unsigned)time(NULL));
	for (int i = 0; i < 100; ++i)
	{
		r[i] = rand() % n + 1;//生成区间r~l的随机数 
	
	}

	return r;
}

Fazzy_amber_noise::Fazzy_amber_noise(Field *f) :
	FazzyModel(f)
{

	ep1 = 1.73*1.73*EPSILON_0_S;			//誘電率 = (屈折率)^2
	ep2 = 1.4*1.4*EPSILON_0_S;
	r1 = mField->nanoToCell(1000);
}

string Fazzy_amber_noise::mkdir(string root) {
	_mkdir((root + "Fazzy_amber_noise").c_str());

	string name = "Fazzy_amber_noise/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}


double Fazzy_amber_noise::calcEPS(const double& x, const double& y, enum INTEG f) {
	
	
	double *p;
	p = getRandom2(100);//生成随机数的通常范围为0~32767，这里通过取模控制取值为0~100 
	//r = mField->nanoToCell(1);
	double a = mField->nanoToCell(500);//半径
	double b = mField->nanoToCell(57.5);
	double La = mField->nanoToCell(50);
	double Lb = mField->nanoToCell(14.5);
	
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;
	double Llength = mField->nanoToCell(29);
	double Hlength = mField->nanoToCell(115);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml(); //- (Hlength + Llength) * 5;
	if (mx < mField->nanoToCell(0.0) || my < 0 || mx >= mField->getNx() - mField->nanoToCell(0.0) || my >= mField->getNy()) return EPSILON_0_S;

	double _x = mx - 0.5*mField->getNx();//(N_X/2, N_Y/2)を原点にシフト
	double _y = my - 0.5*mField->getNy();

	int numout = rand() % 100 + 1;
	if (my > 0.5 *mField->getNy() - (Hlength + Llength) * 4 && my < 0.5*mField->getNy()) {


		for (int j = 0; j < 5; j++) {
			for (int i = 0; i < 11; i++) {
				double index = 0.1 * i;
				//cout << numx << endl;

				//cout << p[(i + 1)*(j + 1)] << endl;
				int numx = rand() % 100 + 1;
				
				double _x = mx - index * mField->getNx();//(N_X/2, N_Y/2)を原点にシフト mx:left side; mx-1:right shift
				
				
				double _y = my - 0.5 *mField->getNy() + (Hlength + Llength)*j;

				//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
				//if ((_x*_x) / a + (_y * _y) / b >= pow(r, 2.0))

				//	return EPSILON_0_S;

				//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外

				//cout << a + mField->nanoToCell(p[numx]) << endl;
				
				if ((_x*_x) / pow(a - mField->nanoToCell(p[numx]), 2.0) + (_y * _y) / pow(b, 2.0) <= pow(1, 2.0)) { // 

					return ep1;
				}
			}
		}return ep2;

	}

	for (int i = 0; i < 11; i++) {
		double index = 0.1 * i;
		//cout << numx << endl;

		//cout << p[(i + 1)*(j + 1)] << endl;
		int numx = rand() % 100 + 1;
		double _x = mx - index * mField->getNx();//(N_X/2, N_Y/2)を原点にシフト mx:left side; mx-1:right shift
		double _y = my - 0.5 *mField->getNy();

		//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
		//if ((_x*_x) / a + (_y * _y) / b >= pow(r, 2.0))

		//	return EPSILON_0_S;

		//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外

		
		if ((_x*_x) / pow(a - mField->nanoToCell(p[numx]), 2.0) + (_y * _y) / pow(b, 2.0) <= pow(1, 2.0)) { // 

			return ep1;
		}
	}
	for (int i = 0; i < 11; i++) {
		double index = 0.1 * i;
		//cout << numx << endl;

		//cout << p[(i + 1)*(j + 1)] << endl;
		int numx = rand() % 100 + 1;
		double _x = mx - index * mField->getNx();//(N_X/2, N_Y/2)を原点にシフト mx:left side; mx-1:right shift
		double _y = my - 0.5 *mField->getNy() + (Hlength + Llength) * 4;

		//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
		//if ((_x*_x) / a + (_y * _y) / b >= pow(r, 2.0))

		//	return EPSILON_0_S;

		//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外


		if ((_x*_x) / pow(a - mField->nanoToCell(p[numx]), 2.0) + (_y * _y) / pow(b, 2.0) <= pow(1, 2.0)) { // 

			return ep1;
		}
	}
	for (int i = 0; i < 11; i++) {
		double index = 0.1 * i;
		//cout << numx << endl;

		//cout << p[(i + 1)*(j + 1)] << endl;
		int numx = rand() % 50 + 1;
		double _x = mx - index * mField->getNx();//(N_X/2, N_Y/2)を原点にシフト mx:left side; mx-1:right shift
		double _y = my - 0.5 *mField->getNy() + (Hlength + Llength) * 4 + Hlength * 0.5;

		//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
		//if ((_x*_x) / a + (_y * _y) / b >= pow(r, 2.0))

		//	return EPSILON_0_S;

		//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外


		//if ((_x*_x) / pow(La + p[numx], 2.0) + (_y * _y) / pow(Lb, 2.0) <= pow(1, 2.0)) { // 
		if ((_x*_x) / pow(La + p[numx], 2.0) + (_y * _y) / pow(Lb, 2.0) <= pow(1, 2.0)) { // 
			return ep2;
		}
	}

	return EPSILON_0_S; 



	/*
	if (my > 0.3*mField->getNy()&& my < 0.5*mField->getNy()) {
	//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
	if ((_x*_x) + (_y * _y) >= pow(r + 1, 2.0))

		return ep2;

	//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外
	if ((_x*_x) + (_y * _y) <= pow(r - 1, 2.0))
		return ep1;
	
	double s = 0;

	double a = 1.0, b = 1.0;
	if (f == D_X) b = 0;
	if (f == D_Y) a = 0;
	for (double i = -16 + 0.5; i < 16; i += 1)
		for (double j = -16 + 0.5; j < 16; j += 1)
			if (pow(_x + a * i / 32.0, 2.0) + pow(_y + b * j / 32.0, 2.0) <= r * r)
				s += 1;
	s /= 32.0*32.0;
	return ep1 * s + EPSILON_0_S * (1 - s);
	}

	return EPSILON_0_S;
*/
/*
	int numout = rand() % 100 + 1;
	if (my > 0.5 *mField->getNy() - (Hlength + Llength)*4 && my < 0.5*mField->getNy()) {
		
		
		for (int j = 0; j < 5; j++) {
			for (int i = 0; i < 11; i++) {
				double index = 0.1 * i;
				//cout << numx << endl;
				
				//cout << p[(i + 1)*(j + 1)] << endl;
				int numx = rand() % 100 + 1;
				double _x = mx - index * mField->getNx();//(N_X/2, N_Y/2)を原点にシフト mx:left side; mx-1:right shift
				double _y = my - 0.5 *mField->getNy()+ (Hlength+Llength)*j;

				//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
				//if ((_x*_x) / a + (_y * _y) / b >= pow(r, 2.0))

				//	return EPSILON_0_S;

				//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外
				

				if ((_x*_x) / pow(a + p[numx], 2.0) + (_y * _y) / pow(b, 2.0) <= pow(1, 2.0)) { // 

					return ep1;
				}
			}
		}return ep2;
	
	}
	
		for (int i = 0; i < 11; i++) {
			double index = 0.1 * i;
			//cout << numx << endl;

			//cout << p[(i + 1)*(j + 1)] << endl;
			int numx = rand() % 100 + 1;
			double _x = mx - index * mField->getNx();//(N_X/2, N_Y/2)を原点にシフト mx:left side; mx-1:right shift
			double _y = my - 0.5 *mField->getNy();

			//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
			//if ((_x*_x) / a + (_y * _y) / b >= pow(r, 2.0))

			//	return EPSILON_0_S;

			//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外


			if ((_x*_x) / pow(a + p[numx], 2.0) + (_y * _y) / pow(b, 2.0) <= pow(1, 2.0)) { // 

				return ep1;
			}
		}
		for (int i = 0; i < 11; i++) {
			double index = 0.1 * i;
			//cout << numx << endl;

			//cout << p[(i + 1)*(j + 1)] << endl;
			int numx = rand() % 100 + 1;
			double _x = mx - index * mField->getNx();//(N_X/2, N_Y/2)を原点にシフト mx:left side; mx-1:right shift
			double _y = my - 0.5 *mField->getNy() + (Hlength + Llength)*4;

			//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
			//if ((_x*_x) / a + (_y * _y) / b >= pow(r, 2.0))

			//	return EPSILON_0_S;

			//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外


			if ((_x*_x) / pow(a + p[numx], 2.0) + (_y * _y) / pow(b, 2.0) <= pow(1, 2.0)) { // 

				return ep1;
			}
		}
		for (int i = 0; i < 11; i++) {
			double index = 0.1 * i;
			//cout << numx << endl;

			//cout << p[(i + 1)*(j + 1)] << endl;
			int numx = rand() % 50 + 1;
			double _x = mx - index * mField->getNx();//(N_X/2, N_Y/2)を原点にシフト mx:left side; mx-1:right shift
			double _y = my - 0.5 *mField->getNy() + (Hlength + Llength) * 4 + Hlength*0.5;

			//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
			//if ((_x*_x) / a + (_y * _y) / b >= pow(r, 2.0))

			//	return EPSILON_0_S;

			//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外


			if ((_x*_x) / pow(La + p[numx], 2.0) + (_y * _y) / pow(Lb, 2.0) <= pow(1, 2.0)) { // 

				return ep2;
			}
		}
	
	return EPSILON_0_S;*/
}


Fazzy_amber_meso::Fazzy_amber_meso(Field *f) :
	FazzyModel(f)
{
	
	ep1 = 1.73*1.73*EPSILON_0_S;			//誘電率 = (屈折率)^2
	ep2 = 1.4*1.4*EPSILON_0_S;
	r1 = mField->nanoToCell(1000);
}

string Fazzy_amber_meso::mkdir(string root) {
	_mkdir((root + "Fazzy_amber_meso").c_str());

	string name = "Fazzy_amber_meso/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}

double Fazzy_amber_meso::calcEPS(const double& x, const double& y, enum INTEG f) {
	//r = mField->nanoToCell(320);
	double Llength = mField->nanoToCell(29);
	double Hlength = mField->nanoToCell(115);
	double mx = x - mField->getNpml();
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double _x = mx - 0.5*mField->getNx();
	double _y = my - 0.0001*mField->getNy()- mField->nanoToCell(250);
/*
	if(my> 0.5*mField->getNy()){
		if (_x*_x + _y * _y >= pow(r1 , 2.0))

			return EPSILON_0_S;

		if (_x*_x + _y * _y <= pow(r1 , 2.0) && _x*_x + _y * _y >= pow(r1-Hlength, 2.0))
			return ep1;
		if (_x*_x + _y * _y <= pow(r1 - Hlength, 2.0) && _x*_x + _y * _y >= pow(r1 - Hlength - Llength, 2.0))
			return ep2;
		if (_x*_x + _y * _y <= pow(r1 - Hlength - Llength, 2.0) && _x*_x + _y * _y >= pow(r1 - Hlength - Llength - Hlength, 2.0))
			return ep1;
		if (_x*_x + _y * _y <= pow(r1 - Hlength - Llength - Hlength, 2.0) && _x*_x + _y * _y >= pow(r1 - Hlength - Llength - Hlength - Llength, 2.0))
			return ep2;
		if (_x*_x + _y * _y <= pow(r1 - Hlength - Llength - Hlength - Llength, 2.0) && _x*_x + _y * _y >= pow(r1 - Hlength - Llength - Hlength - Llength - Hlength, 2.0))
			return ep1;
		if (_x*_x + _y * _y <= pow(r1 - Hlength * 3 - Llength *2 , 2.0) && _x*_x + _y * _y >= pow(r1 - Hlength * 3 - Llength * 3, 2.0))
			return ep2;
		if (_x*_x + _y * _y <= pow(r1 - Hlength * 3 - Llength * 3, 2.0) && _x*_x + _y * _y >= pow(r1 - Hlength * 4 - Llength * 3, 2.0))
			return ep1;
		if (_x*_x + _y * _y <= pow(r1 - Hlength * 4 - Llength * 3, 2.0) && _x*_x + _y * _y >= pow(r1 - Hlength * 4 - Llength * 4, 2.0))
			return ep2;
		if (_x*_x + _y * _y <= pow(r1 - Hlength * 4 - Llength * 4, 2.0) && _x*_x + _y * _y >= pow(r1 - Hlength * 5 - Llength * 4, 2.0))
			return ep1;
		if (_x*_x + _y * _y <= pow(r1 - Hlength * 5 - Llength * 4, 2.0) && _x*_x + _y * _y >= pow(r1 - Hlength * 5 - Llength * 5, 2.0))
			return ep2;
		if (_x*_x + _y * _y <= pow(r1 - Hlength * 5 - Llength * 5, 2.0))
			return EPSILON_0_S;
	}
}
*/

	mx=mx- 0.5*mField->getNx();
	my=my - mField->nanoToCell(250);
	if(my> mField->nanoToCell(250)){
	if (mx*mx + my * my >= pow(r1, 2.0))

		return EPSILON_0_S;

	if (mx*mx + my * my <= pow(r1, 2.0) && mx*mx + my * my >= pow(r1 - Hlength, 2.0))
		return ep1;
	if (mx*mx + my * my <= pow(r1 - Hlength, 2.0) && mx*mx + my * my >= pow(r1 - Hlength - Llength, 2.0))
		return ep2;
	if (mx*mx + my * my <= pow(r1 - Hlength - Llength, 2.0) && mx*mx + my * my >= pow(r1 - Hlength - Llength - Hlength, 2.0))
		return ep1;
	if (mx*mx + my * my <= pow(r1 - Hlength - Llength - Hlength, 2.0) && mx*mx + my * my >= pow(r1 - Hlength - Llength - Hlength - Llength, 2.0))
		return ep2;
	if (mx*mx + my * my <= pow(r1 - Hlength - Llength - Hlength - Llength, 2.0) && mx*mx + my * my >= pow(r1 - Hlength - Llength - Hlength - Llength - Hlength, 2.0))
		return ep1;
	if (mx*mx + my * my <= pow(r1 - Hlength * 3 - Llength * 2, 2.0) && mx*mx + my * my >= pow(r1 - Hlength * 3 - Llength * 3, 2.0))
		return ep2;
	if (mx*mx + my * my <= pow(r1 - Hlength * 3 - Llength * 3, 2.0) && mx*mx + my * my >= pow(r1 - Hlength * 4 - Llength * 3, 2.0))
		return ep1;
	if (mx*mx + my * my <= pow(r1 - Hlength * 4 - Llength * 3, 2.0) && mx*mx + my * my >= pow(r1 - Hlength * 4 - Llength * 4, 2.0))
		return ep2;
	if (mx*mx + my * my <= pow(r1 - Hlength * 4 - Llength * 4, 2.0) && mx*mx + my * my >= pow(r1 - Hlength * 5 - Llength * 4, 2.0))
		return ep1;
	if (mx*mx + my * my <= pow(r1 - Hlength * 5 - Llength * 4, 2.0) && mx*mx + my * my >= pow(r1 - Hlength * 5 - Llength * 5, 2.0))
		return ep2;
	if (mx*mx + my * my <= pow(r1 - Hlength * 5 - Llength * 5, 2.0))
		return EPSILON_0_S;
}
	else {
		return EPSILON_0_S;
	}
}

/*---------------------------------------------*/
/*--------------多層膜-------------------------*/
/*---------------------------------------------*/
FazzySlabModel::FazzySlabModel(Field* f) :
	FazzyModel(f), ep1(2.0*2.0*EPSILON_0_S), ep2(EPSILON_0_S), width1(250), width2(50)
{
}

double FazzySlabModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	//左100nmから,250nm間隔で50nmのスラブを入れていく  **左250nmから(L70.71)10nmスラブに変更(L73)
	//多層膜

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;

	int k = (int)(mField->cellToNano(mx) - 250) % 250;
	double l = (mField->cellToNano(mx) - 250) / 250;

	if (k > 0 && k <= 10 && l < 5)
		return ep1;
	else
		return ep2;

}

string FazzySlabModel::mkdir(string root) {
	_mkdir((root + "SlabModel").c_str());

	string name = "SlabModel/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}


Fazzy_amber::Fazzy_amber(Field* f) :
	FazzyModel(f), ep1(1.73*1.73*EPSILON_0_S), ep2(1.4*1.4*EPSILON_0_S), width1(250), width2(50) {
}

double Fazzy_amber::calcEPS(const double& x, const double& y, enum INTEG f) {
	//左100nmから,250nm間隔で50nmのスラブを入れていく  **左250nmから(L70.71)10nmスラブに変更(L73)
	//多層膜

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();

	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;


	double length = mField->nanoToCell(145.0);
	double Llength = mField->nanoToCell(30);
	double Hlength = mField->nanoToCell(115);
	double firstlayer = mField->nanoToCell(50);


	if (mx < mField->nanoToCell(0) ||  my < 0 || mx >= mField->getNx() - mField->nanoToCell(0) || my >= mField->getNy()) return EPSILON_0_S;

	//double  k = int(my)  % length;//13%5=3
	//double k = fmod(abs(my - cy), length);
	//int l = abs(my - cy) / length; //13/5=2	
	double k = fmod(abs(my - cy), length);
	int l = abs(my - cy) / length; //13/5=2	
//cout << "length : " + to_s(length) << endl;
	//cout << "my : " + to_s(my) << endl;
	//cout << "K : " + to_s(k)<< endl;


/*
if (my>cy &&my-cy<firstlayer )
		return ep1;
	if (k > 0 && k <= Llength && l < 5)
		return ep2;
	if (k >Llength && l < 5)
		return ep1;
	else
		return EPSILON_0_S;

*/


//if (my>cy &&my-cy<firstlayer )
//	return ep2;

	if (my < cy && k > 0 && k <= Hlength && l < 5)
		return ep1;
	if (my<cy &&k > Hlength && l < 5)
		return ep2;
	else
		return EPSILON_0_S;

	/*

	if ( k > 0 && k <= Llength && l < 5)
		return ep2;
	if ( k > Llength && l < 5)
		return ep1;
	else
		return EPSILON_0_S;
*/

}

string Fazzy_amber::mkdir(string root) {
	_mkdir((root + "SlabModel").c_str());

	string name = "SlabModel/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}


/*---------------------------------------------*/
Fazzy_hairslab::Fazzy_hairslab(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(1.4*1.4*EPSILON_0_S), width1(250), width2(50) {
}

double Fazzy_hairslab::calcEPS(const double& x, const double& y, enum INTEG f) {
	//左100nmから,250nm間隔で50nmのスラブを入れていく  **左250nmから(L70.71)10nmスラブに変更(L73)
	//多層膜

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();

	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;


	double length = mField->nanoToCell(580.0);
	double Llength = mField->nanoToCell(80);
	double Hlength = mField->nanoToCell(500);
	double firstlayer = mField->nanoToCell(50);


	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;

	//double  k = int(my)  % length;//13%5=3
	//double k = fmod(abs(my - cy), length);
	//int l = abs(my - cy) / length; //13/5=2	
	double k = fmod(abs(my - cy), length);
	int l = abs(my - cy) / length; //13/5=2	
//cout << "length : " + to_s(length) << endl;
	//cout << "my : " + to_s(my) << endl;
	//cout << "K : " + to_s(k)<< endl;


/*
if (my>cy &&my-cy<firstlayer )
		return ep1;
	if (k > 0 && k <= Llength && l < 5)
		return ep2;
	if (k >Llength && l < 5)
		return ep1;
	else
		return EPSILON_0_S;

*/


//if (my>cy &&my-cy<firstlayer )
//	return ep2;

	if (my < cy && k > 0 && k <= Hlength && l < 5)
		return ep1;
	if (my<cy &&k > Hlength && l < 5)
		return ep2;
	else
		return EPSILON_0_S;



}

string Fazzy_hairslab::mkdir(string root) {
	_mkdir((root + "Fazzy_hairslab").c_str());

	string name = "Fazzy_hairslab/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}



FazzyHair_cuti_only::FazzyHair_cuti_only(Field* f) :
	FazzyModel(f), ep1(1.5*1.5*EPSILON_0_S), ep2(EPSILON_0_S), alpha(3), cwidth(0.5), length(60), r(64), cmc(0.05), e(0.6), beta(3)

{
	alphaR = alpha * PI / 180;

	int n = length * sin(alphaR) / (cwidth*cos(alpha));


}


struct Point
{
	float x;
	float y;
	Point(float x, float y)
	{
		this->x = x;
		this->y = y;
	}
};


float GetCross(Point& p1, Point& p2, Point& p)
{
	return (p2.x - p1.x) * (p.y - p1.y) - (p.x - p1.x) * (p2.y - p1.y);
}

bool IsPointInMatrix(Point& p, Point&p1, Point&p2, Point&p3, Point& p4)
{


	return GetCross(p1, p2, p) * GetCross(p3, p4, p) >= 0 && GetCross(p2, p3, p) * GetCross(p4, p1, p) >= 0;
	//return false;
}

double FazzyHair_cuti_only::calcEPS(const double& x, const double& y, enum INTEG f) {

	ax = 6;
	by = 2.5;

	alphaR = alpha * PI / 180;
	ln = mField->nanoToCell(length * 1000);
	lx = ln * cos(alphaR);
	ly = ln * sin(alphaR);
	rn = mField->nanoToCell(r * 1000); //cortex
	cn = mField->nanoToCell(cwidth * 1000);
	mn = mField->nanoToCell(cmc * 1000);
	h_cn = abs(cn * cos(alphaR));;
	h_cnmn = abs((cn + mn)*cos(alphaR));
	//h_cnmn_up = abs((cn + mn)*sin(beta));


	h_out = 0.1*ln*sin(alphaR);
	l_out = abs(h_out / tan(alphaR));



	double mx = x - mField->getNpml();
	int bound_x = mField->getNx();
	int bound_y = mField->getNy();
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2; //zhongdian
	//double cx = 0;
	double cy = mField->getNy() / 2;

	//cout << "i : " + to_s(floor(cx*2/l_out)) + "micro" << endl;

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;


	//cout << "i : " + to_s(mField->getNpml()) + "micro" << endl;

	//cout << "ny : " + to_s(bound_y) + "micro" << endl;
	if ((mx > cx - abs(ln*sin(alpha)) - h_cnmn - abs((cn + mn)*cos(alphaR)) * 5)) {

		for (int i = 0; i < 15; i++) {

			double p4_x, p4_y, p44_x, p44_y, p3_x, p3_y, p33_x, p33_y, p2_x, p22_x, p22_y, p1_x, p1_y, p2_y, p6_x, p6_y, p5_x, p5_y;
			double upcuti = abs((cn + mn) / tan(alpha));


			
			
			//Field(6000,6000,20, 20);
			p2_y = (i - 10) * abs(((cn + mn) / sin(alphaR))) - 900;
			//p2_y = (i - 2) * 2 * abs(((cn + mn) / sin(alphaR))) - 900;//Field(6000,6000,20, 20);
			//p2_y = (i - 2) * 2 * abs(((cn + mn) / sin(alphaR))) - 950;
			p2_x = cx - h_cnmn * 0;
			p22_y = p2_y - abs(mn * sin(alphaR));
			p22_x = p2_x + abs(mn * cos(alphaR));
			p1_y = p2_y - abs((cn + mn)*sin(alphaR));
			p1_x = cx + h_cnmn - h_cnmn * 0;

			p3_y = p2_y + abs(ln*cos(alphaR));
			p3_x = cx + abs(ln*sin(alphaR)) - h_cnmn * 0;
			p33_y = p3_y - abs(mn*sin(alphaR));
			p33_x = p3_x + abs(mn * cos(alphaR));
			p4_y = p3_y - abs((cn + mn)*sin(alphaR));
			p4_x = p3_x + abs((cn + mn)*cos(alphaR));


			
				/*
			//p2_y = (i - 2) * 2 * abs(((cn + mn) / sin(alphaR))) - 1800;//Field(6000,6000,10, 20);
			//p2_y = (i - 2) * 2 * abs(((cn + mn) / sin(alphaR))) - 1750;//Field(6000,6000,10, 20);
			p2_y = (i - 10) * abs(((cn + mn) / sin(alphaR))) - 950;//(4000,4000,20, 20);
			p2_x = cx - h_cnmn * 0;
			p22_y = p2_y - abs(mn * sin(alphaR));
			p22_x = p2_x + abs(mn * cos(alphaR));
			p1_y = p2_y - abs((cn + mn)*sin(alphaR));
			p1_x = cx + h_cnmn - h_cnmn * 0;

			p3_y = p2_y + abs(ln*cos(alphaR));
			p3_x = cx + abs(ln*sin(alphaR)) - h_cnmn * 0;
			p33_y = p3_y - abs(mn*sin(alphaR));
			p33_x = p3_x + abs(mn * cos(alphaR));
			p4_y = p3_y - abs((cn + mn)*sin(alphaR));
			p4_x = p3_x + abs((cn + mn)*cos(alphaR));
			*/
			
			/*	//(6000,6000,10, 20);
			//p2_y = (i - 10) * abs(((cn + mn) / sin(alphaR))) - 650;
			p2_y = (i - 10) * abs(((cn + mn) / sin(alphaR))) - 1800;
			p2_x = cx - h_cnmn * 0;
			p22_y = p2_y - abs(mn * sin(alphaR));
			p22_x = p2_x + abs(mn * cos(alphaR));
			p1_y = p2_y - abs((cn + mn)*sin(alphaR));
			p1_x = cx + h_cnmn - h_cnmn * 0;

			p3_y = p2_y + abs(ln*cos(alphaR));
			p3_x = cx + abs(ln*sin(alphaR)) - h_cnmn * 0;
			p33_y = p3_y - abs(mn*sin(alphaR));
			p33_x = p3_x + abs(mn * cos(alphaR));
			p4_y = p3_y - abs((cn + mn)*sin(alphaR));
			p4_x = p3_x + abs((cn + mn)*cos(alphaR));
			*/
			/*

			//8000 10cell 
			p2_y = (i -2) * 2 * abs(((cn + mn) / sin(alphaR))) - 1700;
			p2_x = cx - h_cnmn * 0;
			p22_y = p2_y - abs(mn * sin(alphaR));
			p22_x = p2_x + abs(mn * cos(alphaR));
			p1_y = p2_y - abs((cn + mn)*sin(alphaR));
			p1_x = cx + h_cnmn - h_cnmn * 0;

			p3_y = p2_y + abs(ln*cos(alphaR));
			p3_x = cx + abs(ln*sin(alphaR)) - h_cnmn * 0;
			p33_y = p3_y - abs(mn*sin(alphaR));
			p33_x = p3_x + abs(mn * cos(alphaR));
			p4_y = p3_y - abs((cn + mn)*sin(alphaR));
			p4_x = p3_x + abs((cn + mn)*cos(alphaR));
			*/
			/*
			p2_y = (i-10) *2* abs(((cn + mn) / sin(alphaR)))  -1700;
			p2_x = cx - h_cnmn * 0;
			p22_y = p2_y - abs(mn * sin(alphaR));
			p22_x = p2_x + abs(mn * cos(alphaR));
			p1_y = p2_y - abs((cn + mn)*sin(alphaR));
			p1_x = cx + h_cnmn - h_cnmn *0;

			p3_y = p2_y + abs(ln*cos(alphaR));
			p3_x = cx + abs(ln*sin(alphaR)) - h_cnmn *0;
			p33_y = p3_y - abs(mn*sin(alphaR));
			p33_x = p3_x + abs(mn * cos(alphaR));
			p4_y = p3_y - abs((cn + mn)*sin(alphaR));
			p4_x = p3_x + abs((cn + mn)*cos(alphaR));
			
			*/

			/*精度20
			p2_y = (i - 10) * abs(((cn + mn) / sin(alphaR))) - 430;
			p2_x = cx - h_cnmn * 2;
			p22_y = p2_y - abs(mn * sin(alphaR));
			p22_x = p2_x + abs(mn * cos(alphaR));
			p1_y = p2_y - abs((cn + mn)*sin(alphaR));
			p1_x = cx + h_cnmn - h_cnmn * 2;

			p3_y = p2_y + abs(ln*cos(alphaR));
			p3_x = cx + abs(ln*sin(alphaR)) - h_cnmn * 2;
			p33_y = p3_y - abs(mn*sin(alphaR));
			p33_x = p3_x + abs(mn * cos(alphaR));
			p4_y = p3_y - abs((cn + mn)*sin(alphaR));
			p4_x = p3_x + abs((cn + mn)*cos(alphaR));
			
			/*精度10
			p2_y = (i - 10) * abs(((cn + mn) / sin(alphaR))) -850;
			p2_x = cx - h_cnmn *2.5;
			p22_y = p2_y - abs(mn * sin(alphaR));
			p22_x = p2_x + abs(mn * cos(alphaR));
			p1_y = p2_y - abs((cn + mn)*sin(alphaR));
			p1_x = cx + h_cnmn - h_cnmn * 2.5;

			p3_y = p2_y + abs(ln*cos(alphaR));
			p3_x = cx + abs(ln*sin(alphaR)) - h_cnmn *2.5;
			p33_y = p3_y - abs(mn*sin(alphaR));
			p33_x = p3_x + abs(mn * cos(alphaR));
			p4_y = p3_y - abs((cn + mn)*sin(alphaR));
			p4_x= p3_x + abs((cn + mn)*cos(alphaR));
			*/

			Point testPoint(mx, my);
			Point p1(p1_x, p1_y);
			Point p2(p2_x, p2_y);
			Point p3(p3_x, p3_y);
			Point p22(p22_x, p22_y);
			Point p33(p33_x, p33_y);
			Point p4(p4_x, p4_y);

			if (IsPointInMatrix(testPoint, p22, p2, p3, p33)) return 1.42*1.42*EPSILON_0_S;
			if (IsPointInMatrix(testPoint, p1, p22, p33, p4))	return ep1;
			else
				continue;
			//break;

		}
		return ep2;

	}

	return ep2;


}

double FazzyHair_cuti_only::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {

	alphaR = alpha * PI / 180;
	ln = mField->nanoToCell(length * 1000);
	lx = ln * cos(alphaR);
	ly = ln * sin(alphaR);
	rn = mField->nanoToCell(r * 1000); //cortex
	cn = mField->nanoToCell(cwidth * 1000);
	mn = mField->nanoToCell(cmc * 1000);
	h_cn = abs(cn * cos(alpha));;
	h_cnmn = abs((cn + mn)*cos(alpha));

	int bound = mField->getNx();
	h_out = 0.1*ln*sin(alpha);
	l_out = abs(h_out / tan(alpha));
	double mx = x - mField->getNpml();
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;


	return 0;

}

string FazzyHair_cuti_only::mkdir(string root) {
	_mkdir((root + "FazzyHair_cuti_only").c_str());

	string name = "FazzyHair_cuti_only/" + to_s((int)(mField->cellToNano(r))) + "nm," + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//僨傿儗僋僩儕偺嶌惉
	return name + "/";
}


/*---------------------------------------------*/
Fazzy_mesoscutum::Fazzy_mesoscutum(Field* f) :
	FazzyModel(f), ep1(1.73*1.73*EPSILON_0_S), ep2(1.4*1.4*EPSILON_0_S), width1(250), width2(50) {
}

double Fazzy_mesoscutum::calcEPS(const double& x, const double& y, enum INTEG f) {

	double length = mField->nanoToCell(144.0);
	double Llength = mField->nanoToCell(29);
	double Hlength = mField->nanoToCell(115);
	double Tlength = mField->nanoToCell(720);
	r = mField->nanoToCell(1000);

	double mx = x - mField->getNpml(); //計算領域内へ写像 
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;



	for (int m = 11; m < 30; m++) {
		for (int n = 6; n <= 15; n++) {

		}
	}
	double _lx = mx - 0.25*mField->getNx();//(N_X/2, N_Y/2)を原点にシフト
	double _rx = mx - 0.75*mField->getNx();//(N_X/2, N_Y/2)を原点にシフト
	double _x = mx - 0.5*mField->getNx();//(N_X/2, N_Y/2)を原点にシフト
	double _y = my - 0.5*mField->getNy();

	//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
	if (my< 0.5*mField->getNy() && _x*_x + _y * _y <= pow(r, 2.0) && _x*_x + _y * _y > pow(r - Llength, 2.0))
		return ep2;

	if (my < 0.5*mField->getNy() && _x*_x + _y * _y <= pow(r - Llength, 2.0) && _x*_x + _y * _y > pow(r - Llength - Hlength, 2.0))
		return ep1;
	if (my < 0.5*mField->getNy() && _x*_x + _y * _y <= pow(r - Llength - Hlength, 2.0) && _x*_x + _y * _y > pow(r - Llength - Hlength - Llength, 2.0))
		return ep2;
	if (my < 0.5*mField->getNy() && _x*_x + _y * _y <= pow(r - Llength - Hlength - Llength, 2.0) && _x*_x + _y * _y > pow(r - Llength - Hlength - Hlength - Llength, 2.0))
		return ep1;

	if (my < 0.5*mField->getNy() && _x*_x + _y * _y <= pow(r - 2 * Llength - 2 * Hlength, 2.0) && _x*_x + _y * _y > pow(r - 3 * Llength - 2 * Hlength, 2.0))
		return ep2;
	if (my < 0.5*mField->getNy() && _x*_x + _y * _y <= pow(r - 3 * Llength - 2 * Hlength, 2.0) && _x*_x + _y * _y > pow(r - 3 * Llength - 3 * Hlength, 2.0))
		return ep1;
	if (my < 0.5*mField->getNy() && _x*_x + _y * _y <= pow(r - 3 * Llength - 3 * Hlength, 2.0) && _x*_x + _y * _y > pow(r - 4 * Llength - 3 * Hlength, 2.0))
		return ep2;
	if (my < 0.5*mField->getNy() && _x*_x + _y * _y <= pow(r - 4 * Llength - 3 * Hlength, 2.0) && _x*_x + _y * _y > pow(r - 4 * Llength - 4 * Hlength, 2.0))
		return ep1;
	if (my < 0.5*mField->getNy() && _x*_x + _y * _y <= pow(r - 4 * Llength - 4 * Hlength, 2.0) && _x*_x + _y * _y > pow(r - 5 * Llength - 4 * Hlength, 2.0))
		return ep2;
	if (my < 0.5*mField->getNy() && _x*_x + _y * _y <= pow(r - 5 * Llength - 4 * Hlength, 2.0) && _x*_x + _y * _y > pow(r - 5 * Llength - 5 * Hlength, 2.0))
		return ep1;
	//double k = fmod((_x*_x + _y * _y- pow(r - Tlength, 2.0)), length);
	//int l = (x*_x + _y * _y - pow(r - Tlength, 2.0)) / length;

	//cout << "length : " + to_s(k) << endl;

	//if (_x*_x + _y * _y - pow(r - Tlength, 2.0) > 0 && _x*_x + _y * _y - pow(r - Tlength, 2.0) < Tlength)
	//	return ep1;
	else return EPSILON_0_S;

}

string Fazzy_mesoscutum::mkdir(string root) {
	_mkdir((root + "Fazzy_mesoscutum").c_str());

	string name = "Fazzy_mesoscutum/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}


FazzyMorphoModel::FazzyMorphoModel(Field* f, double _h0, double _h1, enum STRUCTURE kind) :
	FazzyModel(f), shelf(kind)
{
	num = 8;				//愊傒廳偹傞悢
	// ep[1] = 3.5*3.5*EPSILON_0_S;	//桿揹棪
	// ep[0] = 1.45*1.45*EPSILON_0_S;
	ep[1] = 1.56*1.56*EPSILON_0_S;	//桿揹棪
	ep[0] = 1.0*1.0*EPSILON_0_S;
	width = mField->nanoToCell(150);	//横幅は150で固定

	min = mField->nanoToCell(120);
	max = mField->nanoToCell(120);
	height[1] = min;
	height[0] = min;
	cout << min << endl;
	cout << max << endl;
}

double FazzyMorphoModel::calcEPS(const double &x, const double &y, enum INTEG f) {
	double mx = x - mField->getNpml(); //寁嶼椞堟撪傊幨憸
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;

	//N_X/2を中心軸に,長方形を横に互い違いに配置
	int dis = height[0] + height[1];
	int oy = (mField->getNy() - num * dis) / 2.0;
	int ox = mField->getNx() / 2.0;
	double _x = mx - ox;	//ox,oy傪嵗昗偺尨揰偵
	double _y = my - oy;

	//モデルの左右か上下に離れていれば媒質の外
	if (abs(_x) > width + 1 || abs(_y - 1.0*num*dis / 2.0) > num*dis / 2.0 + 1)
		return EPSILON_0_S;

	double s[2] = { 0,0 };
	double a = 1.0, b = 1.0;
	if (f == D_X) b = 0;
	if (f == D_Y) a = 0;
	for (double i = -16 + 0.5; i < 16; i += 1) {
		for (double j = -16 + 0.5; j < 16; j += 1) {
			double sx = _x + a * i / 32.0;
			double sy = _y + b * j / 32.0;
			if (abs(sx) > width || abs(sy - 1.0*num*dis / 2.0) > 1.0*num*dis / 2.0) continue;

			bool k = ((int)sy%dis + 0.5) > height[0];	//境界上で比べないように0.5足している(0より大きく1未満なら何でもいい)
			//bool k =  (floor(sy/ dis)*dis < sy) && ( sy < floor(sy/ dis)*dis + height[0]);

			if (sx < 0 && shelf)
				k = !k;		//左右で反転, 互い違いでなかったら反転しない

			s[k] += 1;
		}
	}
	s[0] /= 32.0*32.0;
	s[1] /= 32.0*32.0;
	return EPSILON_0_S * (1 - s[0] - s[1]) + ep[0] * s[0] + ep[1] * s[1];
}

string FazzyMorphoModel::mkdir(string root) {
	string label = "Morpho(" + to_s(sqrt(ep[0])) + "," + to_s(sqrt(ep[1])) + ")M=" + to_s(num);
	//  string label = "Morpho";
	_mkdir((root + label).c_str());
	string name;

	if (shelf)
		name = label + "/" + to_s((int)(mField->cellToNano(height[0]))) + "nm" + mField->getStringCellInfo();
	else
		name = label + "/" + to_s((int)(mField->cellToNano(height[0]))) + "nm(nonShelf)" + mField->getStringCellInfo();

	_mkdir((root + name).c_str());	//僨傿儗僋僩儕偺嶌惉
	return name + "/";
}


/*--------------------------------*/
/*-----------儌僨儖側偟-----------*/
/*--------------------------------*/
bool FazzyMorphoModel::update(int dh) {
	height[0] += (int)mField->nanoToCell(dh);
	height[1] += (int)mField->nanoToCell(dh);

	if (height[0] > max)
		return false;

	return true;
}

string FazzyNoModel::mkdir(string root) {
	_mkdir((root + "NoModel").c_str());

	string name = "NoModel/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//僨傿儗僋僩儕偺嶌惉
	return name + "/";
}

FazzyHair_CrossSection::FazzyHair_CrossSection(Field* f) :
	FazzyModel(f), ep_cor(1.51*1.51*EPSILON_0_S), ep_air(EPSILON_0_S), e(0.6), r(32), ep_cuti(1.6*1.6*EPSILON_0_S), ep_med(1.2*1.2*EPSILON_0_S)
	//a:離心率  r:毛の半径(μm)
{
	cout << "楕円の離心率 = " + to_s((double)e) << endl;
	//ep = 1.6*1.6*EPSILON_0_S;			//誘電率 = (屈折率)^2
}

double FazzyHair_CrossSection::calcEPS(const double& x, const double& y, enum INTEG f) {
	//rn = mField->nanoToCell(r * 1000);
	ax = r;
	by = ax * sqrt(1 - e * e);


	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep_air;	//PML層

	double _x = mx - 0.5*mField->getNx();//(N_X/2, N_Y/2)を原点にシフト
	double _y = my - 0.5*mField->getNy();

	//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
	//if (_x*_x + _y * _y >= pow(rn + 1, 2.0))

		//return EPSILON_0_S;


	double _ax = ax + 1, _by = by + 1;
	if (((_x*_x) / (_ax*_ax)) + ((_y*_y) / (_by*_by)) >= 1)
		return EPSILON_0_S;

	//	_ax = ax - 1;
	//	_by = by - 1;



		//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外
	//if (_x*_x + _y * _y <= pow(r - 1, 2.0))
		//return ep1;


	if (((_x*_x) / (_ax*_ax)) + ((_y*_y) / (_by*_by)) <= 1 && ((_x*_x) / (_ax*_ax)) + ((_y*_y) / (_by*_by)) > 0.80)
		return ep_med;
	if (((_x*_x) / (_ax*_ax)) + ((_y*_y) / (_by*_by)) >= 0.05 && ((_x*_x) / (_ax*_ax)) + ((_y*_y) / (_by*_by)) < 0.80)
		return ep_cor;

	//double s = 0;
	if (((_x*_x) / (_ax*_ax)) + ((_y*_y) / (_by*_by)) <= 0.05)
		return ep_med;


	double s = 0;
	double m = 0;
	double n = 0;
	double a = 1.0, b = 1.0;
	if (f == D_X) b = 0;
	if (f == D_Y) a = 0;
	for (double i = -16 + 0.5; i < 16; i += 1)
		for (double j = -16 + 0.5; j < 16; j += 1)
			if (pow(_x + a * i / 32.0, 2.0) / (ax*ax) + pow(_y + b * j / 32.0, 2.0) / (by*by) <= 1 && pow(_x + a * i / 32.0, 2.0) / (ax*ax) + pow(_y + b * j / 32.0, 2.0) / (by*by) >= 0.8)
				s += 1;
	//if (pow(_x + a * i / 32.0, 2.0) + pow(_y + b * j / 32.0, 2.0) <= r * r)
			else if (pow(_x + a * i / 32.0, 2.0) / (ax*ax) + pow(_y + b * j / 32.0, 2.0) / (by*by) < 0.8 && pow(_x + a * i / 32.0, 2.0) / (ax*ax) + pow(_y + b * j / 32.0, 2.0) / (by*by) >= 0.05)
				m += 1;
			else if (pow(_x + a * i / 32.0, 2.0) / (ax*ax) + pow(_y + b * j / 32.0, 2.0) / (by*by) < 0.05)
				n += 1;
	//cout << setprecision(20) << ld << endl;
	printf("%d", s);
	s /= 32.0*32.0;
	m /= 32.0*32.0;
	n /= 32.0*32.0;
	//return ep1 * s + ep2 * (1 - s);
	return ep_cuti * s + ep_cor * m + ep_med * n + EPSILON_0_S * (1 - s - m - n);

}
string FazzyHair_CrossSection::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	_mkdir((root + "HairModel/CrossSection").c_str());

	string name = "HairModel/CrossSection/e=" + to_s((double)e);
	_mkdir((root + name).c_str());	//ディレクトリの作成

	name = "HairModel/CrossSection/e=" + to_s((double)e) + "/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}



/*---------------------------------------------*/
/*---------------毛髪--------------------------*/
/*---------------------------------------------*/

/*----------------------------*/
/*-----------縦断面-----------*/
/*----------------------------*/
FazzyHair_incidenceModel::FazzyHair_incidenceModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), alpha(5), cwidth(0.5), r(32 + (4.36 - cwidth))
	//alpha:キューティクルの角度(deg)  cwidth:キューティクルの幅(μm)  r:毛の半径(μm)(半径+キューティクルが重なる領域)
{
	alphaR = alpha * PI / 180;
	length = cwidth / sin(alphaR);
	cout << "キューティクルの角度 : " + to_s(alpha) + "deg" << endl;
	cout << "キューティクル幅 : " + to_s(cwidth) + "micro" << endl;
	cout << "キューティクル1枚の露出幅 : " + to_s(length) + "micro" << endl;
}

double FazzyHair_incidenceModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	alphaR = alpha * PI / 180;
	ln = mField->nanoToCell(length * 1000);
	lx = ln * cos(alphaR);
	rn = mField->nanoToCell(r * 1000);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML層
	if (my < cy)	my = 2 * cy - my;		//x軸に対して線対称


	/***************************************************/
	//上面だけキューティクルなし
//	if (y - mField->getNpml() >= mField->nanoToCell(32 * 1000) + cy)	return ep2;
	/***************************************************/


	int c = mField->getNx() / lx + 1;		//計算範囲内のキューティクルの数
	for (int i = 0; i < c; i++) {
		if (mx > i * lx + h && mx < (i + 1) * lx + h && mx < mField->getNx() - h) {
			//			if (my > tan(alphaR) * (mx - lx*i) + cy + rn)	return ep2;
			//			else return ep1;		//Fuzzyなし(Staircaseモデル)

			double dy1 = my - (tan(alphaR) * (mx - lx * i - h) + cy + rn);
			double dy2 = my - (tan(alphaR) * ((mx - lx * i - h) + 1) + cy + rn);
			double s;
			if (dy1 > 0 && dy2 > 0) return ep2;		//キューティクル直線の外側 (1)
			if (fabs(dy1) > 1 && fabs(dy2) > 1) return ep1;		//キューティクル直線の内側 (2)
			/*
			if (dy1 <= 0 && dy2 <= 0) {
				if (fabs(dy1) <= 1 && fabs(dy2) <= 1) {
					s = (fabs(dy1) + fabs(dy2)) * 1.0 / 2.0;
					return ep1 * s + ep2 * (1 - s);		// (3)
				}
				if (fabs(dy1) < 1 && fabs(dy2) > 1) {
					s = (1 - fabs(dy1)) * ((my - cy - rn) / tan(alphaR) - (mx - lx*i - h)) / 2;
					return ep2 * s + ep1 * (1 - s);		// (4)
				}
			}
			if (dy1 > 0 && dy2 < 0) {
				s = fabs(dy2) * (((mx - lx*i - h) + 1) - (my - cy - rn) / tan(alphaR)) / 2;
				return ep1 * s + ep2 * (1 - s);		// (5)
			}*/
		}

		else
			continue;
		//break;
	}

	return ep2;
}

double FazzyHair_incidenceModel::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	rn = mField->nanoToCell(32 * 1000);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;	//PML層
	if (my < cy)	my = 2 * cy - my;		//x軸に対して線対称

	if (my <= rn + cy) {
		int k = (int)(mField->cellToNano(my) - mField->cellToNano(cy) - 1500) % 4000;
		double l = (mField->cellToNano(my) - mField->cellToNano(cy) - 1500) / 4000;

		if (k > 0 && k <= 1000 && l < 8)	return 1.0;		//黒色の色素
		else	return 0;
	}
	else  return 0;
}

string FazzyHair_incidenceModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	string name;

	if (mField->sig == false) {
		_mkdir((root + "HairModel/incidenceplane").c_str());				//吸収係数なしの場合
		name = "HairModel/incidenceplane/" + mField->getStringCellInfo();
	}
	else if (mField->sig == true) {
		_mkdir((root + "HairModel/incidenceplane_withSig").c_str());		//吸収係数ありの場合
		name = "HairModel/incidenceplane_withSig/" + mField->getStringCellInfo();
	}

	_mkdir((root + name).c_str());	//ディレクトリの作成

	return name + "/";
}

//simple cuticle model
FazzyHair_simple_cuti_Model::FazzyHair_simple_cuti_Model(Field* f) :
	FazzyModel(f), ep1(1.6*1.6*EPSILON_0_S), ep2(EPSILON_0_S), alpha(13), cwidth(500), length(100), r(32), cmc(0.06)
	//alpha:キューティクルの角度(deg)  cwidth:キューティクルの厚さ(μm)  length:キューティクルの長さ(μm)	r:毛皮質範囲の半径(μm)  cmc:CMC範囲(μm)
{
	alphaR = alpha * PI / 180;

	//int n = length * sin(alphaR) / (cmc + cwidth);


	cout << "キューティクルの角度 : " + to_s(alpha) + "deg" << endl;
	cout << "キューティクル厚さ : " + to_s(cwidth) + "micro" << endl;
	cout << "キューティクル長さ : " + to_s(length) + "micro" << endl;
	cout << "キューティクル1枚の露出幅 : " + to_s((cmc + cwidth) / sin(alphaR)) + "micro" << endl;
	cout << "キューティクル範囲幅 : " + to_s(length*sin(alphaR)) + "micro" << endl;
	//cout << "キューティクル重なり枚数 : " + to_s(n) + "枚" << endl;
	cout << "the height of cuticle inside : " + to_s(length * sin(alphaR) - cmc - cwidth) + "mu_m" << endl;
}
double FazzyHair_simple_cuti_Model::calcEPS(const double& x, const double& y, enum INTEG f) {
	alphaR = alpha * PI / 180;
	ln = mField->nanoToCell(length * 100);
	lx = ln * cos(alphaR);
	ly = ln * sin(alphaR);
	rn = mField->nanoToCell(r * 100); //cortex
	cn = mField->nanoToCell(cwidth * 100);
	mn = mField->nanoToCell(cmc * 100);

	//h_in = ln * sin(alphaR) - mn - cn ;
	//h_out = mn + cn;
	//l_out = h_out / tan(alpha);

	h_out = 0.1*ln*sin(alpha);
	l_out = abs(h_out / tan(alpha));
	//cout << "h_in範囲幅 : " + to_s(h_in) + "micro" << endl;
	//cout << "h_out範囲幅 : " + to_s(h_out) + "micro" << endl;

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2; //zhongdian
	//double cx = 0;
	double cy = mField->getNy() / 2;
	cout << "i : " + to_s(mField->getNy()) + "micro" << endl;
	cout << "py : " + to_s(mField->getNpy()) + "micro" << endl;

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML層

	/**** 毛髪全体 ****/
	//if (my < cy)	my = 2 * cy - my;		//x軸に対して線対称  //对称到上边
	//if (my <= rn + cy)		  ep1;		//毛皮質部分
	//my = my - rn - cy;		//軸移動 减去cortex？，只留cuti


	//if (my <= cy-50) return ep1;
	//if (my >= cy + 30) return ep1;
	//if (my <= cy - 30) return ep1;

	if (my < cy + rn && my > cy) {
		//if (mx > cx + 50) 
		return 0.5*ep1;
	}
	else if (my > cy + rn && my < cy + rn + h_out) {
		for (int i = 0; i <= floor(cx * 2 / l_out); i++) {

			if (mx < (i + 1)*l_out && mx >= i * l_out) {
				//cout << "i的值：" << i << endl;
				//cout << "lout：" << (l_out) << endl;
				//cout << "mx-：" << ((mx - i * l_out) ) << endl;
				//cout << "mx：" << ((mx )) << endl;
				//cout << "mx：" << ((i*l_out)) << endl;
				//cout << "substract：" << (  (my - cy - rn)) << endl;
				//cout << "tan：" << ((l_out - (mx - i * l_out))*tan(alpha)) << endl;
				//cout << "tan：" << (my - cy - rn) << endl;
				//if (((mx - i * l_out) / l_out) < 0.5) 
				if (abs((l_out - (mx - i * l_out))*tan(alpha)) > (my - cy - rn)) return ep1;

				else return ep2;
				//if ((l_out - (mx - i * l_out))*tan(alpha) > (my - cy  - rn)) return ep1;
				//else return ep2;

			}
			//else return ep2 * 0.2;
		}

	}
	else if (my > cy - rn && my < cy) {
		//if (mx > cx + 50) 
		return 0.5*ep1;
	}
	else if (my < cy - rn && my > cy - rn - h_out) {
		for (int i = 0; i <= floor(cx * 2 / l_out); i++) {

			if (mx < (i + 1)*l_out && mx >= i * l_out) {

				if (abs((l_out - (mx - i * l_out))*tan(alpha)) > (cy - rn - my)) return ep1;

				else return ep2;
				//if ((l_out - (mx - i * l_out))*tan(alpha) > (my - cy  - rn)) return ep1;
				//else return ep2;

			}
			//else return ep2 * 0.2;
		}

	}
	else return ep2;

	//else if (my < (cy - rn)) return ep1;

/*
	if (my > ly)	return ep2;
	for (int i = -(ly / (cn + mn)); (i*mn + i * cn) / tan(alphaR) < mField->getNx(); i++) {//这时候半径都斜着。还能除tan？
		if (my <= ly - cn) {
			if (my >= tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) + 1 && my <= tan(alphaR) * (mx - (i*mn + i * cn) / tan(alphaR)) - 1)
				return ep1;

			else if (my >= tan(alphaR) * (mx - ((i + 1)*mn + (i + 1)*cn) / tan(alphaR)) + 1 && my <= tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) - 1)
				return ep2;

			else if ((my > tan(alphaR) * (mx - (i*mn + i * cn) / tan(alphaR)) - 1 && my < tan(alphaR) * (mx - (i*mn + i * cn) / tan(alphaR)) + 1)
				|| (my > tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) - 1 && my < tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) + 1)) {

				double s = 0;
				for (double a = -16; a < 16; a += 1)
					for (double b = -16; b < 16; b += 1)
						if (my + b / 32.0 >= tan(alphaR) * (mx + a / 32.0 - (i*mn + (i + 1)*cn) / tan(alphaR)) && my + b / 32.0 <= tan(alphaR) * (mx + a / 32.0 - (i*mn + i * cn) / tan(alphaR)) && mx <= (ly + i * mn + i * cn) / tan(alphaR))
							s += 1;
				s /= 32.0*32.0;
				return ep1 * s + ep2 * (1 - s);

			}

		}
		else if (my > ly - cn) {
			if (mx >= (my + i * mn + i * cn + 1) / tan(alphaR) && mx <= (ly + i * mn + i * cn) / tan(alphaR))
				return ep1;

			else if (mx < (my + i * mn + i * cn + 1) / tan(alphaR) && mx >(my + i * mn + i * cn - 1) / tan(alphaR)) {

				double s = 0;
				for (double a = -16; a < 16; a += 1)
					for (double b = -16; b < 16; b += 1)
						if (mx + a / 32.0 >= (my + b / 32.0 + i * mn + i * cn) / tan(alphaR) && mx <= (ly + i * mn + i * cn) / tan(alphaR))
							s += 1;
				s /= 32.0*32.0;
				return ep1 * s + ep2 * (1 - s);
			}
		}
	}
	*/

	//else return ep2;
}

double FazzyHair_simple_cuti_Model::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	rn = mField->nanoToCell(32 * 1000);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;	//PML層
	//if (my < cy)	my = 2 * cy - my;		//x軸に対して線対称


	if (my <= rn + cy) {
		int k = (int)(mField->cellToNano(my) - mField->cellToNano(cy) - 1500) % 4000;
		double l = (mField->cellToNano(my) - mField->cellToNano(cy) - 1500) / 4000;

		if (k > 0 && k <= 1000 && l < 8)	return 1.0;		//黒色の色素
		else	return 1;
	}
	else  return 1;
	;

}

string FazzyHair_simple_cuti_Model::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	string name;

	if (mField->sig == false) {
		_mkdir((root + "HairModel/incidenceLayer").c_str());				//吸収係数なしの場合
		name = "HairModel/incidenceLayer/" + mField->getStringCellInfo();
	}
	else if (mField->sig == true) {
		_mkdir((root + "HairModel/incidenceLayer_withSig").c_str());		//吸収係数ありの場合
		name = "HairModel/incidenceLayer_withSig/" + mField->getStringCellInfo();
	}

	_mkdir((root + name).c_str());	//ディレクトリの作成

	return name + "/";
}


//simple cuticle model
FazzyHair_cuti_Model::FazzyHair_cuti_Model(Field* f) :
	FazzyModel(f), ep1(1.6*1.6*EPSILON_0_S), ep2(EPSILON_0_S), alpha(3), cwidth(1), length(60), r(32), cmc(0.08)
	//alpha:キューティクルの角度(deg)  cwidth:キューティクルの厚さ(μm)  length:キューティクルの長さ(μm)	r:毛皮質範囲の半径(μm)  cmc:CMC範囲(μm)
{
	alphaR = alpha * PI / 180;

	int n = length * sin(alphaR) / (cwidth*cos(alpha));


	cout << "キューティクルの角度 : " + to_s(alpha) + "deg" << endl;
	cout << "キューティクル厚さ : " + to_s(cwidth) + "micro" << endl;
	cout << "キューティクル長さ : " + to_s(length) + "micro" << endl;
	cout << "ly : " + to_s(length * sin(alphaR)) + "micro" << endl;
	cout << "each length in the bottom : " + to_s((cwidth *cos(alpha)) / tan(alpha)) + "micro" << endl;
	//cout << "number of cuticle in the bottom : " + to_s(cn *cos(alpha)) + "micro" << endl;
	//cout << "キューティクル範囲幅 : " + to_s(length*sin(alphaR)) + "micro" << endl;
	cout << "キューティクル重なり枚数 : " + to_s(n) + "枚" << endl;
	//cout << "the height of cuticle inside : " + to_s(length * sin(alphaR) - cmc - cwidth) + "mu_m" << endl;
}
/*
struct Point
{
	float x;
	float y;
	Point(float x, float y)
	{
		this->x = x;
		this->y = y;
	}
};
// 计算 |p1 p2| X |p1 p|
float GetCross(Point& p1, Point& p2, Point& p)
{
	return (p2.x - p1.x) * (p.y - p1.y) - (p.x - p1.x) * (p2.y - p1.y);
}
//判断点是否在5X5 以原点为左下角的正方形内（便于测试）
bool IsPointInMatrix(Point& p, Point&p1, Point&p2, Point&p3, Point& p4)
{


	return GetCross(p1, p2, p) * GetCross(p3, p4, p) >= 0 && GetCross(p2, p3, p) * GetCross(p4, p1, p) >= 0;
	//return false;
}
*/
double FazzyHair_cuti_Model::calcEPS(const double& x, const double& y, enum INTEG f) {
	alphaR = alpha * PI / 180;
	ln = mField->nanoToCell(length * 100);
	lx = ln * cos(alphaR);
	ly = ln * sin(alphaR);
	rn = mField->nanoToCell(r * 100); //cortex
	cn = mField->nanoToCell(cwidth * 100);
	mn = mField->nanoToCell(cmc * 100);
	h_cn = abs(cn * cos(alpha));
	h_cnmn = abs((cn + mn)*cos(alpha));


	h_out = 0.1*ln*sin(alpha);
	l_out = abs(h_out / tan(alpha));
	//cout << "h_in範囲幅 : " + to_s(h_in) + "micro" << endl;
	//cout << "h_out範囲幅 : " + to_s(h_out) + "micro" << endl;

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2; //zhongdian
	//double cx = 0;
	double cy = mField->getNy() / 2;
	//cout << "i : " + to_s(floor(cx*2/l_out)) + "micro" << endl;

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML層



}

double FazzyHair_cuti_Model::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	rn = mField->nanoToCell(32 * 1000);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;	//PML層
	//if (my < cy)	my = 2 * cy - my;		//x軸に対して線対称


	if (my <= rn + cy) {
		int k = (int)(mField->cellToNano(my) - mField->cellToNano(cy) - 1500) % 4000;
		double l = (mField->cellToNano(my) - mField->cellToNano(cy) - 1500) / 4000;

		if (k > 0 && k <= 1000 && l < 8)	return 1.0;		//黒色の色素
		else	return 1.0;
	}
	else  return 1.0;


}

string FazzyHair_cuti_Model::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	string name;

	if (mField->sig == false) {
		_mkdir((root + "HairModel/incidenceLayer").c_str());				//吸収係数なしの場合
		name = "HairModel/incidenceLayer/" + mField->getStringCellInfo();
	}
	else if (mField->sig == true) {
		_mkdir((root + "HairModel/incidenceLayer_withSig").c_str());		//吸収係数ありの場合
		name = "HairModel/incidenceLayer_withSig/" + mField->getStringCellInfo();
	}

	_mkdir((root + name).c_str());	//ディレクトリの作成

	return name + "/";
}



/*--------------------------------------------------*/
/*-----------縦断面(多層膜キューティクル)-----------*/
/*--------------------------------------------------*/
FazzyHair_incidenceLayerModel::FazzyHair_incidenceLayerModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), alpha(5), cwidth(0.5), length(50), r(32), cmc(0.06)
	//alpha:キューティクルの角度(deg)  cwidth:キューティクルの厚さ(μm)  length:キューティクルの長さ(μm)	r:毛皮質範囲の半径(μm)  cmc:CMC範囲(μm)
{
	alphaR = alpha * PI / 180;
	int n = length * sin(alphaR) / (cmc + cwidth);

	cout << "キューティクルの角度 : " + to_s(alpha) + "deg" << endl;
	cout << "キューティクル厚さ : " + to_s(cwidth) + "micro" << endl;
	cout << "キューティクル長さ : " + to_s(length) + "micro" << endl;
	cout << "キューティクル1枚の露出幅 : " + to_s((cmc + cwidth) / sin(alphaR)) + "micro" << endl;//露出的斜边。
	cout << "キューティクル範囲幅 : " + to_s(length*sin(alphaR)) + "micro" << endl;
	cout << "キューティクル重なり枚数 : " + to_s(n) + "枚" << endl;
}
double FazzyHair_incidenceLayerModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	alphaR = alpha * PI / 180;
	ln = mField->nanoToCell(length * 1000);
	lx = ln * cos(alphaR);
	ly = ln * sin(alphaR);
	rn = mField->nanoToCell(r * 1000); //cortex
	cn = mField->nanoToCell(cwidth * 1000);
	mn = mField->nanoToCell(cmc * 1000);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2; //zhongdian
	double cy = mField->getNy() / 2;

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML層

	/**** 毛髪全体 ****/
	if (my < cy)	my = 2 * cy - my;		//x軸に対して線対称  //对称到上边
	if (my <= rn + cy)		return ep1;		//毛皮質部分
	my = my - rn - cy;		//軸移動 减去cortex？，只留cuti

	/****************** キューティクル部分のみ ******************/
	/* Field推奨サイズ                                          */
	/* Field(8000, 8000, 5, 10) Field(16000, 8000, 10, 10) など */
	/************************************************************/
//	if (my <= 100)		return ep1;		//毛皮質部分
//	my = my - 100;		//軸移動


	if (my > ly)	return ep2; //如果my在cuti上边：空气ep2
	for (int i = -(ly / (cn + mn)); (i*mn + i * cn) / tan(alphaR) < mField->getNx(); i++) {//这时候半径都斜着。还能除tan？（层数；pml界内；i++）
		if (my <= ly - cn) { //（cuti内cmc之上？）
			if (my >= tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) + 1 && my <= tan(alphaR) * (mx - (i*mn + i * cn) / tan(alphaR)) - 1)//斜着插秧。+1-1？
				return ep1;

			else if (my >= tan(alphaR) * (mx - ((i + 1)*mn + (i + 1)*cn) / tan(alphaR)) + 1 && my <= tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) - 1)
				return ep2;

			else if ((my > tan(alphaR) * (mx - (i*mn + i * cn) / tan(alphaR)) - 1 && my < tan(alphaR) * (mx - (i*mn + i * cn) / tan(alphaR)) + 1)
				|| (my > tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) - 1 && my < tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) + 1)) {

				double s = 0;
				for (double a = -16; a < 16; a += 1)
					for (double b = -16; b < 16; b += 1)
						if (my + b / 32.0 >= tan(alphaR) * (mx + a / 32.0 - (i*mn + (i + 1)*cn) / tan(alphaR)) && my + b / 32.0 <= tan(alphaR) * (mx + a / 32.0 - (i*mn + i * cn) / tan(alphaR)) && mx <= (ly + i * mn + i * cn) / tan(alphaR))
							s += 1;
				s /= 32.0*32.0;
				return ep1 * s + ep2 * (1 - s);

			}

		}
		else if (my > ly - cn) {
			if (mx >= (my + i * mn + i * cn + 1) / tan(alphaR) && mx <= (ly + i * mn + i * cn) / tan(alphaR))
				return ep1;

			else if (mx < (my + i * mn + i * cn + 1) / tan(alphaR) && mx >(my + i * mn + i * cn - 1) / tan(alphaR)) {

				double s = 0;
				for (double a = -16; a < 16; a += 1)
					for (double b = -16; b < 16; b += 1)
						if (mx + a / 32.0 >= (my + b / 32.0 + i * mn + i * cn) / tan(alphaR) && mx <= (ly + i * mn + i * cn) / tan(alphaR))
							s += 1;
				s /= 32.0*32.0;
				return ep1 * s + ep2 * (1 - s);
			}
		}
	}

	return ep2;
}

double FazzyHair_incidenceLayerModel::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	rn = mField->nanoToCell(32 * 1000);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;	//PML層
	if (my < cy)	my = 2 * cy - my;		//x軸に対して線対称

	if (my <= rn + cy) {
		int k = (int)(mField->cellToNano(my) - mField->cellToNano(cy) - 1500) % 4000;
		double l = (mField->cellToNano(my) - mField->cellToNano(cy) - 1500) / 4000;

		if (k > 0 && k <= 1000 && l < 8)	return 1.0;		//黒色の色素
		else	return 0;
	}
	else  return 0;
}

string FazzyHair_incidenceLayerModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	string name;

	if (mField->sig == false) {
		_mkdir((root + "HairModel/incidenceLayer").c_str());				//吸収係数なしの場合
		name = "HairModel/incidenceLayer/" + mField->getStringCellInfo();
	}
	else if (mField->sig == true) {
		_mkdir((root + "HairModel/incidenceLayer_withSig").c_str());		//吸収係数ありの場合
		name = "HairModel/incidenceLayer_withSig/" + mField->getStringCellInfo();
	}

	_mkdir((root + name).c_str());	//ディレクトリの作成

	return name + "/";
}

/*------------------------------------------------*/
/*-----------縦断面(キューティクルなし)-----------*/
/*------------------------------------------------*/
FazzyHair_NONcuticleModel::FazzyHair_NONcuticleModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), r(32)
	//r:毛の半径(μm)
{
	cout << "キューティクル : なし" << endl;
}

double FazzyHair_NONcuticleModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	rn = mField->nanoToCell(r * 1000);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML層
	if (my < cy)	my = 2 * cy - my;		//x軸に対して線対称

	if (my <= rn + cy)	return ep1;
	else  return ep2;
}

double FazzyHair_NONcuticleModel::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	rn = mField->nanoToCell(r * 1000);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;	//PML層
	if (my < cy)	my = 2 * cy - my;		//x軸に対して線対称

	if (my <= rn + cy) {
		int k = (int)(mField->cellToNano(my) - mField->cellToNano(cy) - 1500) % 4000;
		double l = (mField->cellToNano(my) - mField->cellToNano(cy) - 1500) / 4000;

		if (k > 0 && k <= 1000 && l < 8)	return 1.0;		//黒色の色素
		else	return 0;
	}
	else  return 0;
}

string FazzyHair_NONcuticleModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	string name;

	if (mField->sig == false) {
		_mkdir((root + "HairModel/NONcuticle").c_str());				//吸収係数なしの場合
		name = "HairModel/NONcuticle/" + mField->getStringCellInfo();
	}
	else if (mField->sig == true) {
		_mkdir((root + "HairModel/NONcuticle_withSig").c_str());		//吸収係数ありの場合
		name = "HairModel/NONcuticle_withSig/" + mField->getStringCellInfo();
	}

	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}

/*----------------------------*/
/*-----------横断面-----------*/
/*----------------------------*/
FazzyHair_normalModel::FazzyHair_normalModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), e(0.6), r(32)
	//a:離心率  r:毛の半径(μm)
{
	cout << "楕円の離心率 = " + to_s((double)e) << endl;
}

double FazzyHair_normalModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	rn = mField->nanoToCell(r * 1000);
	ax = rn;
	by = ax * sqrt(1 - e * e);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML層

	double _x = mx - 0.5*mField->getNx();//(N_X/2, N_Y/2)を原点にシフト
	double _y = my - 0.5*mField->getNy();

	//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
	double _ax = ax + 1, _by = by + 1;
	if ((_x*_x) / (_ax*_ax) + (_y*_y) / (_by*_by) >= 1)
		return ep2;

	_ax = ax - 1;
	_by = by - 1;
	//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外
	if ((_x*_x) / (_ax*_ax) + (_y*_y) / (_by*_by) <= 1)
		return ep1;

	double s = 0;

	double a = 1.0, b = 1.0;
	if (f == D_X) b = 0;
	if (f == D_Y) a = 0;
	for (double i = -16 + 0.5; i < 16; i += 1)
		for (double j = -16 + 0.5; j < 16; j += 1)
			if (pow(_x + a * i / 32.0, 2.0) / (ax*ax) + pow(_y + b * j / 32.0, 2.0) / (by*by) <= 1)
				s += 1;
	s /= 32.0*32.0;
	return ep1 * s + ep2 * (1 - s);
}

string FazzyHair_normalModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	_mkdir((root + "HairModel/normalplane").c_str());

	string name = "HairModel/normalplane/e=" + to_s((double)e);
	_mkdir((root + name).c_str());	//ディレクトリの作成

	name = "HairModel/normalplane/e=" + to_s((double)e) + "/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}

#include"Model.h"
#include"Field.h"


FazzyMieModel::FazzyMieModel(Field *f):
FazzyModel(f)
{
	ep = 1.6*1.6*EPSILON_0_S;		
		
}

string FazzyMieModel::mkdir(string root){
	_mkdir((root + "Mie").c_str());

	string name = "Mie/" + to_s((int)(mField->cellToNano(r))) +"nm,"+ mField->getStringCellInfo();
	_mkdir((root + name).c_str());	
	return name + "/";
}

double FazzyMieModel::calcEPS(const double& x, const double& y, enum INTEG f){

	double mx = x - mField->getNpml(); 
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;

	double _x = mx - 0.5*mField->getNx();
	double _y = my - 0.5*mField->getNy();


	if(_x*_x + _y*_y >= pow(r+1, 2.0))

		return EPSILON_0_S;


	if(_x*_x + _y*_y <= pow(r-1, 2.0))
		return ep;

	double s=0;

	double a = 1.0,b=1.0;
	if(f == D_X) b = 0;
	if(f == D_Y) a = 0;
	for(double i=-16+0.5; i<16; i+=1)
		for(double j=-16+0.5; j<16; j+=1)
			if(pow(_x+a*i/32.0, 2.0) + pow(_y+b*j/32.0, 2.0) <= r*r)
				s+=1; 
	s /= 32.0*32.0;
	return ep*s + EPSILON_0_S*(1-s);
}



Fazzy_amber::Fazzy_amber(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(1.4*1.4*EPSILON_0_S), width1(250), width2(50) {
}

double Fazzy_amber::calcEPS(const double& x, const double& y, enum INTEG f) {

	double mx = x - mField->getNpml(); 
	double my = y - mField->getNpml();

	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;


	double length = mField->nanoToCell(600);
	double Llength = mField->nanoToCell(100);
	double Hlength = mField->nanoToCell(500);
	double firstlayer = mField->nanoToCell(50);


	if (mx < mField->nanoToCell(0) || my < 0 || mx >= mField->getNx() - mField->nanoToCell(0) || my >= mField->getNy()) return EPSILON_0_S;


	double k = fmod(abs(my - cy), length);
	int l = abs(my - cy) / length; 





	if (my < cy && k > 0 && k <= Hlength && l < 6)
		return ep1;
	if (my<cy &&k > Hlength && l < 6)
		return ep2;
	else
		return EPSILON_0_S;



}

string Fazzy_amber::mkdir(string root) {
	_mkdir((root + "SlabModel").c_str());

	string name = "SlabModel/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	
	return name + "/";
}

double * getRandom(int n)
{
	static double  r[100];
	mt19937 gen((unsigned int)time(nullptr));
	uniform_real_distribution<double> dis(100, 400);


	for (int i = 0; i < 100; ++i)
	{

		r[i] = dis(gen);

	}

	return r;
}

double * getRandom2(int n)
{
	static double  r[100];

	srand((unsigned)time(NULL));
	for (int i = 0; i < 100; ++i)
	{
		r[i] = rand() % n + 1;

	}

	return r;
}


Fazzy_amber_noise::Fazzy_amber_noise(Field *f) :
	FazzyModel(f)
{

	ep1 = 1.73*1.73*EPSILON_0_S;	
	ep2 = 1.4*1.4*EPSILON_0_S;
	r1 = mField->nanoToCell(1000);
}

string Fazzy_amber_noise::mkdir(string root) {
	_mkdir((root + "Fazzy_amber_noise").c_str());

	string name = "Fazzy_amber_noise/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	
	return name + "/";
}


double Fazzy_amber_noise::calcEPS(const double& x, const double& y, enum INTEG f) {


	double *p;
	p = getRandom2(50);

	double a = mField->nanoToCell(125);
	double b = mField->nanoToCell(57.5);
	double La = mField->nanoToCell(50);
	double Lb = mField->nanoToCell(14.5);

	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;
	double Llength = mField->nanoToCell(29 );
	double Hlength = mField->nanoToCell(115 );

	double mx = x - mField->getNpml(); 
	double my = y - mField->getNpml(); 
	if (mx < mField->nanoToCell(0.0) || my < 0 || mx >= mField->getNx() - mField->nanoToCell(0.0) || my >= mField->getNy()) return EPSILON_0_S;

	double _x = mx - 0.5*mField->getNx();
	double _y = my - 0.5*mField->getNy();

	int numout = rand() % 100 + 1;
	if (my > 0.5 *mField->getNy() - (Hlength + Llength) * 4 && my < 0.5*mField->getNy()) {


		for (int j = 0; j < 5; j++) {
			for (int i = 0; i < 11; i++) {
				double index = 0.1 * i;
		
				int numx = rand() % 100 + 1;

				double _x = mx - index * mField->getNx();


				double _y = my - 0.5 *mField->getNy() + (Hlength + Llength)*j;

				if ((_x*_x) / pow(a + mField->nanoToCell(p[numx]), 2.0) + (_y * _y) / pow(b, 2.0) <= pow(1, 2.0)) { // 

					return ep1;
				}
			}
		}return ep2;

	}

	for (int i = 0; i < 11; i++) {
		double index = 0.1 * i;
		int numx = rand() % 100 + 1;
		double _x = mx - index * mField->getNx();
		double _y = my - 0.5 *mField->getNy();


		if ((_x*_x) / pow(a + mField->nanoToCell(p[numx]), 2.0) + (_y * _y) / pow(b, 2.0) <= pow(1, 2.0)) { // 

			return ep1;
		}
	}
	for (int i = 0; i < 11; i++) {
		double index = 0.1 * i;
	
		int numx = rand() % 100 + 1;
		double _x = mx - index * mField->getNx();
		double _y = my - 0.5 *mField->getNy() + (Hlength + Llength) * 4;


		if ((_x*_x) / pow(a + mField->nanoToCell(p[numx]), 2.0) + (_y * _y) / pow(b, 2.0) <= pow(1, 2.0)) {
			return ep1;
		}
	}
	for (int i = 0; i < 11; i++) {
		double index = 0.1 * i;
		
		int numx = rand() % 50 + 1;
		double _x = mx - index * mField->getNx();
		double _y = my - 0.5 *mField->getNy() + (Hlength + Llength) * 4 + Hlength * 0.5;

		if ((_x*_x) / pow(La + p[numx], 2.0) + (_y * _y) / pow(Lb, 2.0) <= pow(1, 2.0)) { 
			return ep2;
		}
	}
	return EPSILON_0_S;

}
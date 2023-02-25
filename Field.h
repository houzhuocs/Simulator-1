#ifndef _FIELD_INFO_H
#define _FIELD_INFO_H
#include<string>
#include"Object.h"

class Field{
private:
	const int WIDTH;	//����[nm]
	const int HEIGHT;	//�c��[nm]
	const double H_U;		//1�Z���̑傫��[nm]
	const int H_S;		//�Z���Ԋu(�킩��₷���悤��1�ɂ��Ă���)
	const int N_X;		//���̔z��
	const int N_Y;		//�c�̔z��
	const int N_PML;	//PML�ŗp���郌�C���̌���
	const int N_PX;		//PML���C����܂߂����̔z��
	const int N_PY;		//PML���C����܂߂��c�̔z��
	const int N_CELL;	//���Z����
public:
	Field(int width, int height, double h_u, int pml);

	bool sig;	//�z���W��(�F�f)�̗L��

	//�v�Z�̈�ɂ����A�N�Z�X�ł��Ȃ��C���f�b�N�X�֐�(���̐�����Ă�����΂�����)
	int index(const int& i, const int& j){
		//return (i+N_PML)*N_PY + (j+N_PML);
		return i*N_PY + j;
	}

	//pml�̈�ɂ�A�N�Z�X�ł���C���f�b�N�X�֐�
	int pmlIndex(const int &i, const int &j){
		return i*N_PY + j;
	}

	//�v�Z�̈�̏��,(�t�@�C�����p)
	string getStringCellInfo(){
		return "(" +to_s(H_U) + "nm,"+ to_s(N_X) + "cell" +  ")";
	}

	//nm�P�ʂ̕����ʂ�V�~�����[�V�����l�ɕϊ�
	double nanoToCell(const double &length){
		return length/H_U;
	}

	double cellToNano(const double &cell){
		return cell*H_U;
	}

	//�Q�b�^�[
	int getNx(){
		return N_X;
	}
	int getNy(){
		return N_Y;
	}
	int getNpx(){
		return N_PX;
	}
	int getNpy(){
		return N_PY;
	}
	int getNcel(){
		return N_CELL;
	}
	int getHu(){
		return H_U;
	}
	int getNpml(){
		return N_PML;
	}

	double sigmaX(const int &i, const int &j); //��x
	double sigmaY(const int &i, const int &j); //��y
};

#endif //_FIELD_INFO_H
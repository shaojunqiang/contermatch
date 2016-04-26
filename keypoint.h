


#pragma once



#include"CurveCSS.h"

struct Keypoint
{
	int octave; //�ؼ���������
	int interval;// �ؼ������ڲ�

	double offset_interval;//������Ĳ������

	int x; //x,y����,����octave��interval��ȡ�Ĳ���ͼ��
	int y;
	Point  position;
	Point  positionS;
	//scale = sigma0*pow(2.0, o+s/S)
	double scale; //�ռ�߶�����

	double dx; //���������꣬�����걻���ų�ԭͼ���С 
	double dy;

	double offset_x;
	double offset_y;

	//��˹���������ڸ���߶����꣬��ͬ�����ͬ���sigmaֵ��ͬ
	//�ؼ�������������ڳ߶�
	double octave_scale; //offset_i;

	double ori;//����

	int descr_length;
	vector<double> descriptor; //

	vector<double> globaldescriptor;
	double val;//��ֵ
};
double CalcCrossCorrelation(const vector<double>& x, const vector<double>& y);
vector<Keypoint>  GetKeypointDescriptor(vector<Point> & src, int N_neighbor, int nScale, int nBufferStep, int nScaleStart);

vector<Keypoint> GetKeypointGlobaldescriptorDescriptor(vector<Point> & src);

double MatchCurvesSmithWaterman(const vector<Keypoint >& a, const vector<Keypoint >& b,
	vector<Keypoint >& a_out, vector<Keypoint >& b_out,
	vector<Point>& traceback);


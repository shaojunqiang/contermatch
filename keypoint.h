


#pragma once



#include"CurveCSS.h"

struct Keypoint
{
	int octave; //关键点所在组
	int interval;// 关键点所在层

	double offset_interval;//调整后的层的增量

	int x; //x,y坐标,根据octave和interval可取的层内图像
	int y;
	Point  position;
	Point  positionS;
	//scale = sigma0*pow(2.0, o+s/S)
	double scale; //空间尺度坐标

	double dx; //特征点坐标，该坐标被缩放成原图像大小 
	double dy;

	double offset_x;
	double offset_y;

	//高斯金字塔组内各层尺度坐标，不同组的相同层的sigma值相同
	//关键点所在组的组内尺度
	double octave_scale; //offset_i;

	double ori;//方向

	int descr_length;
	vector<double> descriptor; //

	vector<double> globaldescriptor;
	double val;//极值
};
double CalcCrossCorrelation(const vector<double>& x, const vector<double>& y);
vector<Keypoint>  GetKeypointDescriptor(vector<Point> & src, int N_neighbor, int nScale, int nBufferStep, int nScaleStart);

vector<Keypoint> GetKeypointGlobaldescriptorDescriptor(vector<Point> & src);

double MatchCurvesSmithWaterman(const vector<Keypoint >& a, const vector<Keypoint >& b,
	vector<Keypoint >& a_out, vector<Keypoint >& b_out,
	vector<Point>& traceback);


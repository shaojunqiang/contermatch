#pragma once
#include<vector>
#include <map>
#include <opencv2/opencv.hpp>
#include<math.h>
#include"stdax.h"

#include"keypoint.h"

#include"SmithWaterman.h"
using namespace cv;
using namespace std;

void   changeBPoint(/*const*/ vector<Point>& a, vector<Point>& std, int n);

template<typename T, typename V>
void PolyLineSplit(const vector<Point_<T>>& pl, vector<V>& contourx, vector<V>& contoury)
{
	contourx.resize(pl.size());
	contoury.resize(pl.size());
	int size = pl.size();
	for (int j = 0; j<size; j++)
	{
		contourx[j] = (V)(pl[j].x);
		contoury[j] = (V)(pl[j].y);
	}
}
template<typename T, typename V>
void PolyLineMerge(vector<Point_<T> >& pl, const vector<V>& contourx, const vector<V>& contoury) {
	assert(contourx.size() == contoury.size());
	pl.resize(contourx.size());
	for (int j = 0; j<(int)contourx.size(); j++) {
		pl[j].x = (T)(contourx[j]);
		pl[j].y = (T)(contoury[j]);
	}
}
//转换曲线
template<typename T, typename V>
void ConvertCurve(const vector<Point_<T> >& curve, vector<Point_<V> >& output)
{
	output.clear();
	for (int j = 0; j<(int)curve.size(); j++)
	{
		output.push_back(Point_<V>(curve[j].x, curve[j].y));
	}
}

void ResampleCurve(const vector<double>& curvex, const vector<double>& curvey,
	vector<double>& resampleX, vector<double>& resampleY,
	int N, bool isOpen = false
	);

template<typename T, typename V>
void ResampleCurve(const vector<Point_<T> >& curve, vector<Point_<V> >& output,
	int N, bool isOpen = false
	)
{
	vector<double> curvex, curvey, resampleX, resampleY;
	PolyLineSplit(curve, curvex, curvey);
	ResampleCurve(curvex, curvey, resampleX, resampleY, N, isOpen);
	PolyLineMerge(output, resampleX, resampleY);
}
template<typename T>
void drawOpenCurve(Mat& img, const vector<Point_<T> >& curve, Scalar color, int thickness)
{
	if (curve.size() <= 0)
	{
		return;
	}
	vector<cv::Point> curve2i;
	ConvertCurve(curve, curve2i);
	for (int i = 0; i<(int)curve2i.size() - 1; i++)
	{
		line(img, curve2i[i], curve2i[i + 1], color, thickness);
	}
	/*line(img, curve2i[curve2i.size()-1], curve2i[0], color, thickness);*/
}


template<typename T>
void drawPoint(Mat& img, const vector<Point_<T> >& curve, vector<int> curvePoint, Scalar color, int thickness)
{
	if (curve.size() <= 0)
	{
		return;
	}
	vector<cv::Point> curve2i;
	ConvertCurve(curve, curve2i);
	for (int i = 0; i<(int)curvePoint.size() ; i++)
	{
		circle(img, curve2i[curvePoint[i]], 3, color, thickness);
	}
	/*line(img, curve2i[curve2i.size()-1], curve2i[0], color, thickness);*/
}
template<typename T>
void drawPoint(Mat& img, const vector<Point_<T> >& curve,  Scalar color, int thickness)
{
	if (curve.size() <= 0)
	{
		return;
	}
	vector<cv::Point> curve2i;
	ConvertCurve(curve, curve2i);
	for (int i = 0; i<(int)curve2i.size(); i++)
	{
		circle(img, curve2i[i], 3, color, thickness);
	}
	/*line(img, curve2i[curve2i.size()-1], curve2i[0], color, thickness);*/
}



void SmoothCurve(const vector<double>& curvex,
	const vector<double>& curvey,
	vector<double>& smoothX,
	vector<double>& smoothY,
	vector<double>& X,
	vector<double>& XX,
	vector<double>& Y,
	vector<double>& YY,
	double sigma,
	bool isOpen = false);


template<typename T, typename V>
void SimpleSmoothCurve(const vector<Point_<T> >& curve,
	vector<Point_<V> >& smooth,
	double sigma,
	bool isOpen = false) {
	vector<double> contourx(curve.size()), contoury(curve.size());
	PolyLineSplit(curve, contourx, contoury);

	vector<double> smoothx, smoothy, X, XX, Y, YY;
	SmoothCurve(contourx, contoury, smoothx, smoothy, X, XX, Y, YY, sigma, isOpen);

	PolyLineMerge(smooth, smoothx, smoothy);
}



//#pragma mark CSS Image

void ComputeCurveCSS(const vector<double>& curvex,
	const vector<double>& curvey,
	vector<double>& kappa,
	vector<double>& smoothX, vector<double>& smoothY,
	double sigma = 1.0,
	bool isOpen = false);

vector<int> FindCSSInterestPointsZero(const vector<double>& kappa);
vector<int> FindCSSInterestPointsPeak(const vector<double>& kappa);
vector<int> FindCSSInterestPointsPeak(const vector<double>& kappa, int n);
template<typename T>
vector<int> ComputeCSSImageMaximas(const vector<Point_<T> >& curve
	/*vector<Point_<T>>& smooth,bool isClosedCurve = true*/)
{
	vector<double> contourx_(curve.size()), contoury_(curve.size());
	PolyLineSplit(curve, contourx_, contoury_);
	vector<double> smoothx, smoothy;
	vector<int> ID;
	bool isClosedCurve = true;
	ID=ComputeCSSImageMaximas(contourx_,  contoury_, smoothx,  smoothy , isClosedCurve);
	/*PolyLineMerge(smooth, smoothx, smoothy)*/;

	return ID;
}
vector<int> ComputeCSSImageMaximas(const vector<double>& contourx_, const vector<double>& contoury_,
	vector<double>& contourx, vector<double>& contoury,
	bool isClosedCurve
	);

template<typename T>
void ComputeCurveCSS(const vector<Point_<T> >& curve,
	vector<double>& kappa,
	vector<Point_<T> >& smooth,
	double sigma,
	bool isOpen = false
	)
{
	vector<double> contourx(curve.size()), contoury(curve.size());
	PolyLineSplit(curve, contourx, contoury);

	vector<double> smoothx, smoothy;
	ComputeCurveCSS(contourx, contoury, kappa, smoothx, smoothy, sigma, isOpen);

	PolyLineMerge(smooth, smoothx, smoothy);
}
int getcurveSegments(const vector<Point> & in, vector<vector<Point2d> >& out, double sigma);

template<typename T>
double MatchCurvesSmithWaterman(const vector<T>& a, const vector<T >& b,
	vector<T>& a_out, vector<T>& b_out,
	vector<Point>& traceback);

int getcurveSegmentsmean(const vector<Point> & in, vector<vector<Point2d> >& out, double sigma, int n);

int getcurveSegmentsmean(const vector<Point> & in, vector<vector<double> >& out, double sigma, int n);

void getcurveSignature(const vector<Point> & in, vector<vector<double> >& out, int n, vector<vector<int>>&  KeyPoint);

void getcurveSignature(const vector<Point2d> & in, vector<vector<double> >& out, int n);




/* 根据4对坐标点计算最小二乘平面单应性变换矩阵
参数：
pts：坐标点数组
mpts：对应点数组，pts[i]与mpts[i]一一对应
n：pts和mpts数组中点的个数，pts和mpts中点的个数必须相同，一般是4
返回值：一个3*3的变换矩阵，将pts中的每一个点转换为mpts中的对应点，返回值为空表示失败
*/
/*
Calculates a least-squares planar homography from point correspondeces.
@param pts array of points
@param mpts array of corresponding points; each pts[i], i=0..n-1, corresponds to mpts[i]
@param n number of points in both pts and mpts; must be at least 4
@return Returns the 3 x 3 least-squares planar homography matrix that
transforms points in pts to their corresponding points in mpts or NULL if
fewer than 4 correspondences were provided
*/
CvMat* lsq_homog(vector<Point2d>& pts, vector<Point2d>& mpts, int n);

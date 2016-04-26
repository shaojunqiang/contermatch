


#pragma once
#include"stdax.h"

//double MatchCurvesSmithWaterman(const vector<Point2d> & a, const vector<Point2d> & b, vector<Point>& traceback);


template<typename T>
void   changeBPoint(/*const*/ vector<T>& a, int n)
{
	vector<T> out;
	out.clear();
	int j = 0;
	out.insert(out.begin(), a.begin() + n, a.end());
	for (int i = 0; i < n; i++)
	{
		out.push_back(a[i]);
	}
	a.clear();
	a.insert(a.begin(), out.begin(), out.end());
}
void   findMatchPointArray(vector<Point> & a, vector<Point> & b, int a_size, vector<Point> &out, int n);
Mat ImageMerge(Mat& imgMat, Mat& imgMat1);
Mat ScanImageAndReduceC(Mat& imgMat, Mat & imgMat1);
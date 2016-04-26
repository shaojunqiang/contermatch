#pragma once

#include"stdax.h"
#include"CurveCSS.h"

 template<typename T>
int closestPointOnCurveToPoint(const vector<cv::Point_<T> >& _tmp, const cv::Point& checkPt, const T cutoff) {
	vector<cv::Point_<T> > tmp = _tmp;
	Mat(tmp) -= Scalar(checkPt.x, checkPt.y);
	vector<float> tmpx, tmpy, tmpmag;
	PolyLineSplit(tmp, tmpx, tmpy);
	magnitude(tmpx, tmpy, tmpmag);
	double minDist = -1;
	cv::Point minLoc; minMaxLoc(tmpmag, &minDist, 0, &minLoc);
	if (minDist<cutoff)
		return minLoc.x;
	else
		return -1;
}
template<typename T>
void saveCurveToFile(const vector<Point_<T> >& curve) {
	static int curve_id = 0;

	stringstream ss; ss << "curves/curves_" << (curve_id++) << ".txt";
	while (fileExists(ss.str())) {
		ss.str("");
		ss << "curves/curves_" << (curve_id++) << ".txt";
	}

	ofstream ofs(ss.str().c_str());
	ofs << curve.size() << "\n";
	for (int i = 0; i<curve.size(); i++) {
		ofs << curve[i].x << " " << curve[i].y << "\n";
	}
	cout << "saved " << ss.str() << "\n";
}

template<typename T>
vector<Point_<T> > loadCurveFromFile(const string& filename) {
	vector<Point_<T> > curve;
	ifstream ifs(filename.c_str());
	int curve_size; ifs >> skipws >> curve_size;
	while (!ifs.eof()) {
		T x, y;
		ifs >> x >> y;
		curve.push_back(Point_<T>(x, y));
	}
	return curve;
}


template<typename V>
Mat_<double> Find2DRigidTransform(const vector<Point_<V> >& a, const vector<Point_<V> >& b,
	Point_<V>* diff = 0, V* angle = 0, V* scale = 0) {
	//use PCA to find relational scale
	Mat_<V> P; Mat(a).reshape(1, a.size()).copyTo(P);
	Mat_<V> Q; Mat(b).reshape(1, b.size()).copyTo(Q);
	PCA a_pca(P, Mat(), CV_PCA_DATA_AS_ROW), b_pca(Q, Mat(), CV_PCA_DATA_AS_ROW);
	double s = sqrt(b_pca.eigenvalues.at<V>(0)) / sqrt(a_pca.eigenvalues.at<V>(0));
	//	cout << a_pca.eigenvectors << endl << a_pca.eigenvalues << endl << a_pca.mean << endl;
	//	cout << b_pca.eigenvectors << endl << b_pca.eigenvalues << endl << b_pca.mean << endl;

	//convert to matrices and subtract mean
	//	Mat_<double> P(a.size(),2),Q(b.size(),2);
	Scalar a_m = Scalar(a_pca.mean.at<V>(0), a_pca.mean.at<V>(1));
	Scalar b_m = Scalar(b_pca.mean.at<V>(0), b_pca.mean.at<V>(1));
	//	for (int i=0; i<a.size(); i++) { P(i,0) = a[i].x - a_m[0]; P(i,1) = a[i].y - a_m[1]; }
	//	for (int i=0; i<b.size(); i++) { Q(i,0) = b[i].x - b_m[0]; Q(i,1) = b[i].y - b_m[1]; }
	P -= repeat((Mat_<V>(1, 2) << a_m[0], a_m[1]), P.rows, 1);
	Q -= repeat((Mat_<V>(1, 2) << b_m[0], b_m[1]), Q.rows, 1);

	//    cout << "new mean for a " << mean(P) << "\n";

	//from http://en.wikipedia.org/wiki/Kabsch_algorithm
	Mat_<double> A = P.t() * Q;
	SVD svd(A);
	Mat_<double> C = svd.vt.t() * svd.u.t();
	double d = (determinant(C) > 0) ? 1 : -1;
	Mat_<double> R = svd.vt.t() * (Mat_<double>(2, 2) << 1, 0, 0, d) * svd.u.t();
	Mat_<double> T = (Mat_<double>(3, 3) << 1, 0, b_m[0] / s, 0, 1, b_m[1] / s, 0, 0, 1) *
		(Mat_<double>(3, 3) << s, 0, 0, 0, s, 0, 0, 0, s) *
		(Mat_<double>(3, 3) << R(0, 0), R(0, 1), 0, R(1, 0), R(1, 1), 0, 0, 0, 1) *
		(Mat_<double>(3, 3) << 1, 0, -a_m[0], 0, 1, -a_m[1], 0, 0, 1)
		;
	if (diff != NULL) {
		diff->x = b_m[0] - a_m[0];
		diff->y = b_m[1] - a_m[1];
	}
	if (angle != NULL) {
		*angle = atan2(R(1, 0), R(0, 0));
	}
	if (scale != NULL) {
		*scale = s;
	}
	return T(Range(0, 2), Range::all());
}

template<typename T, typename V>
Mat_<T> ConvertToMat(const vector<vector<V> >& mnt_DB) {
	Mat_<T> mnt_DB_m(mnt_DB.size(), mnt_DB[0].size());
	for (int i = 0; i<mnt_DB.size(); i++) {
		for (int j = 0; j<mnt_DB[i].size(); j++) {
			mnt_DB_m(i, j) = (T)(mnt_DB[i][j]);
		}
	}
	return mnt_DB_m;
}
template<typename T>
void drawCurvePoints(Mat& img, const vector<Point_<T> >& curve_, const Scalar& color, int thickness) {
	vector<cv::Point> curve;
	ConvertCurve(curve_, curve);
	for (int i = 0; i<curve.size(); i++) {
		circle(img, curve[i], 3, color, thickness);
	}
}



template<class _Ty> inline
void swaps(_Ty& _Left, _Ty& _Right)
_NOEXCEPT_OP(is_nothrow_move_constructible<_Ty>::value
&& is_nothrow_move_assignable<_Ty>::value)
{	// exchange values stored at _Left and _Right
	_Ty _Tmp = _Move(_Left);
	_Left = _Move(_Right);
	_Right = _Move(_Tmp);
}
void GetCurveForImage(const Mat& filename, vector<Point>& whole, vector<Point>& curve_upper, vector<Point>& curve_lower) {
	assert(!filename.empty());
	Mat tmp; filename.copyTo(tmp);
	Mat gray, canny_output;
	if (tmp.type() == CV_8UC3)
		cvtColor(tmp, gray, CV_BGR2GRAY);
	else if (tmp.type() == CV_8UC1)
		gray = tmp;
	else
		cvError(-1, "GetCurveForImage", "unsupported image format", __FILE__, __LINE__);

	//形式：void cvThreshold( const CvArr* src, CvArr* dst, double threshold, double max_value, int threshold_type );
	threshold(gray, gray, 128, 255, THRESH_BINARY);//threshold_type：阈值类型 threshold_type=CV_THRESH_BINARY:
	//如果 src(x, y)>threshold, dst(x, y) = max_value; 否则, dst（x, y） = 0;

	vector<vector<Point> > contours, _cons;
	/*blur(gray, gray, Size(3, 3));*/

	//Canny(gray, canny_output, 80, 255, 3);
	vector<Vec4i> hierarchy;
	findContours(gray, _cons, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_NONE);

	for (int j = 0; j < (int)_cons.size(); j++)
	{
		cv::Mat contourMat = cv::Mat(_cons[j]);
		double area = cv::contourArea(contourMat);
		if (area < /*areaThresh*/ 1000) continue;



		contours.push_back(_cons[j]);
	}
	

	//vector<vector<Point> > contours;
	///*blur(gray, gray, Size(3, 3));*/

	////Canny(gray, canny_output, 80, 255, 3);
	//vector<Vec4i> hierarchy;
	//findContours(gray, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_NONE);

	//if (contours.size() <= 0) return;
	//Mat results(gray.size(), CV_8U, Scalar(0));

	//for (int i = 0; i < (int)contours.size(); i++)
	//{
	//	drawOpenCurve(results, (vector<Point>)contours[i], Scalar(100, 0, 0), 1);
	//}

	//namedWindow("results", CV_WINDOW_AUTOSIZE);
	//imshow("results", results); //显示图像

	vector<Point> upperCurve = contours[0];
	if (upperCurve.size() <= 50) {
		return;
	}

	//find minimal and maximal X coord
	vector<double> x, y;
	PolyLineSplit(contours[0], x, y);//x,y坐标值
	Point minxp, maxxp;
	minMaxLoc(x, 0, 0, &minxp, &maxxp);
	int minx = minxp.x, maxx = maxxp.x;
	if (minx > maxx)
	{
		swaps(minx, maxx);
	}

	//take lower and upper halves of the curve
	vector<Point> upper, lower;
	upper.insert(upper.begin(), contours[0].begin() + minx, contours[0].begin() + maxx);//复制向量
	lower.insert(lower.begin(), contours[0].begin() + maxx, contours[0].end());
	lower.insert(lower.end(), contours[0].begin(), contours[0].begin() + minx);

	//test which is really the upper part, by looking at the y-coord of the mid point

	if (lower[lower.size() / 2].y <= upper[upper.size() / 2].y) {
		curve_upper = lower;
		curve_lower = upper;
	}
	else {
		curve_upper = upper;
		curve_lower = lower;
	}

	//make sure it goes left-to-right
	if (curve_upper.front().x > curve_upper.back().x) { //hmmm, need to flip
		reverse(curve_upper.begin(), curve_upper.end());
	}

	whole.clear();
	/*whole.insert(whole.begin(),curve_upper.rbegin(),curve_upper.rend());
	whole.insert(whole.begin(),curve_lower.begin(),curve_lower.end());*/
	whole.insert(whole.begin(), contours[0].begin(), contours[0].end());
}

void GetCurveForImage(const Mat& filename, vector<Point>& curve, bool onlyUpper, bool getLower)
{
	vector<Point> whole, upper, lower;
	GetCurveForImage(filename, whole, upper, lower);
	if (onlyUpper) {
		if (getLower)
			curve = lower;
		else
			curve = upper;
	}
	else {
		curve = whole;
	}
}


void PrepareSignatureDB(const vector<Point>& curve_, vector<vector<double> >& DB, vector<Point>& DB_params)
{
	vector<Point> curve;
	if (curve_.size() != 200) {
		ResampleCurve(curve_, curve, 200, true);
	}
	else {
		curve = curve_;
	}


	vector<double> kappa;
	vector<Point2d> smooth;
	SimpleSmoothCurve(curve, smooth, 5.0, true);
	vector<Point2d> smalls;

	DB.clear(); DB_params.clear();
	for (int len = 50; len <(int)smooth.size() - 2; len += 5)
	{
		//iterate different curve sizes, starting at 20 points
		//		cout << "len " << len <<  endl;

		for (int off = (smooth.size() - len); off >= 0; off -= 5)
		{
			//iterate segments on A curve
			vector<Point2d> small_smooth_input(smooth.begin() + off, smooth.begin() + off + len);

			//resample to N points
			ResampleCurve(small_smooth_input, smalls, 200, true);

			//compute curvature
			vector<Point2d> small_smooth;
			ComputeCurveCSS(smalls, kappa, small_smooth, 0.66667, true);
			vector<double> kappa_(kappa.begin() + 1, kappa.end() - 1);

			DB.push_back(kappa_);
			DB_params.push_back(Point(len, off));
		}
	}

	cout << "DB size " << DB.size() << endl;
}




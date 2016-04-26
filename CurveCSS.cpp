




/*
*  CurveCSS.cpp
*  CurveMatching
*
*  Created by Roy Shilkrot on 11/28/12.
*
*/




#include "CurveCSS.h"
#include"keypoint.h"

#include"stdax.h"


#define TwoPi 6.28318530718

//#pragma mark Gaussian Smoothing and Curvature 

void   changeBPoint(/*const*/ vector<Point>& a, vector<Point>& std, int n)
{
	vector<Point> out;
	out.clear();
	int j = 0;
	out.insert(out.begin(), a.begin() + n, a.end());
	for (int i = 0; i < n; i++)
	{
		out.push_back(a[i]);
	}
	std.clear();
	std.insert(std.begin(), out.begin(), out.end());
}



/* 1st and 2nd derivative of 1D gaussian
*/
void getGaussianDerivs(double sigma, int M, vector<double>& gaussian, vector<double>& dg, vector<double>& d2g)
{
	//	static double sqrt_two_pi = sqrt(TwoPi);
	int L;
	if (sigma < 0)
	{
		M = 1;
		L = 0;
		dg.resize(M); d2g.resize(M); gaussian.resize(M);
		gaussian[0] = dg[0] = d2g[0] = 1.0;
		return;
	}

	L = (M - 1) / 2;
	dg.resize(M); d2g.resize(M); gaussian.resize(M);
	getGaussianKernel(M, sigma, CV_64F).copyTo(Mat(gaussian));

	double sigma_sq = sigma * sigma;
	double sigma_quad = sigma_sq*sigma_sq;
	for (double i = -L; i < L + 1.0; i += 1.0)
	{
		int idx = (int)(i + L);
		// from http://www.cedar.buffalo.edu/~srihari/CSE555/Normal2.pdf
		dg[idx] = (-i / sigma_sq) * gaussian[idx];
		d2g[idx] = (-sigma_sq + i*i) / sigma_quad * gaussian[idx];
	}
}

/* 1st and 2nd derivative of smoothed curve point */
void getdX(vector<double> x,
	int n,
	double sigma,
	double& gx,
	double& dgx,
	double& d2gx,
	const vector<double>& g,
	const vector<double>& dg,
	const vector<double>& d2g,
	bool isOpen = false)
{
	int L = (g.size() - 1) / 2;

	gx = dgx = d2gx = 0.0;
	//	cout << "Point " << n << ": ";
	for (int k = -L; k < L + 1; k++)
	{
		double x_n_k;
		if (n - k < 0)
		{
			if (isOpen)
			{
				//open curve - 
				//mirror values on border
				//				x_n_k = x[-(n-k)]; 
				//stretch values on border
				x_n_k = x.front();
			}
			else
			{
				//closed curve - take values from end of curve
				x_n_k = x[x.size() + (n - k)];
			}
		}
		else if (n - k > (int)x.size() - 1)
		{
			if (isOpen) {
				//mirror value on border
				//				x_n_k = x[n+k]; 
				//stretch value on border
				x_n_k = x.back();
			}
			else
			{
				x_n_k = x[(n - k) - (x.size())];
			}
		}
		else
		{
			//			cout << n-k;
			x_n_k = x[n - k];
		}
		//		cout << "* g[" << g[k + L] << "], ";

		gx += x_n_k * g[k + L]; //gaussians go [0 -> M-1]
		dgx += x_n_k * dg[k + L];
		d2gx += x_n_k * d2g[k + L];
	}
	//	cout << endl;
}


/* 0th, 1st and 2nd derivatives of whole smoothed curve */
void getdXcurve(vector<double> x,
	double sigma,
	vector<double>& gx,
	vector<double>& dx,
	vector<double>& d2x,
	const vector<double>& g,
	const vector<double>& dg,
	const vector<double>& d2g,
	bool isOpen = false)
{
	gx.resize(x.size());
	dx.resize(x.size());
	d2x.resize(x.size());
	for (int i = 0; i<(int)x.size(); i++)
	{
		double gausx, dgx, d2gx;
		getdX(x, i, sigma, gausx, dgx, d2gx, g, dg, d2g, isOpen);
		gx[i] = gausx;
		dx[i] = dgx;
		d2x[i] = d2gx;
	}
}

void ResampleCurve(const vector<double>& curvex, const vector<double>& curvey,
	vector<double>& resampleX, vector<double>& resampleY,
	int N,
	bool isOpen
	)
{
	assert(curvex.size()>0 && curvey.size()>0 && curvex.size() == curvey.size());

	vector<Point2d> resamplepl(N); resamplepl[0].x = curvex[0]; resamplepl[0].y = curvey[0];
	vector<Point2i> pl; PolyLineMerge(pl, curvex, curvey);

	double pl_length = arcLength(pl, true);
	double resample_size = pl_length / (double)N;
	int curr = 0;
	double dist = 0.0;
	for (int i = 1; i<N;)
	{
		assert(curr <(int)pl.size() - 1);
		double last_dist = norm(pl[curr] - pl[curr + 1]);
		dist += last_dist;
		//		cout << curr << " and " << curr+1 << "\t\t" << last_dist << " ("<<dist<<")"<<endl;
		if (dist >= resample_size)
		{
			//put a point on line
			double _d = last_dist - (dist - resample_size);
			Point2d cp(pl[curr].x, pl[curr].y), cp1(pl[curr + 1].x, pl[curr + 1].y);
			Point2d dirv = cp1 - cp; dirv = dirv * (1.0 / norm(dirv));
			//			cout << "point " << i << " between " << curr << " and " << curr+1 << " remaining " << dist << endl;
			assert(i <(int)resamplepl.size());
			resamplepl[i] = cp + dirv * _d;
			i++;

			dist = last_dist - _d; //remaining dist			

			//if remaining dist to next point needs more sampling... (within some epsilon)
			while (dist - resample_size > 1e-3)
			{
				//				cout << "point " << i << " between " << curr << " and " << curr+1 << " remaining " << dist << endl;
				assert(i <(int)resamplepl.size());
				resamplepl[i] = resamplepl[i - 1] + dirv * resample_size;
				dist -= resample_size;
				i++;
			}
		}

		curr++;
	}

	PolyLineSplit(resamplepl, resampleX, resampleY);
}

//#pragma mark CSS image

void SmoothCurve(const vector<double>& curvex,
	const vector<double>& curvey,
	vector<double>& smoothX,
	vector<double>& smoothY,
	vector<double>& X,
	vector<double>& XX,
	vector<double>& Y,
	vector<double>& YY,
	double sigma,
	bool isOpen)
{
	int M = round((10.0*sigma + 1.0) / 2.0) * 2 - 1;
	//	assert(M % 2 == 1); //M is an odd number

	vector<double> g, dg, d2g;
	getGaussianDerivs(sigma, M, g, dg, d2g);


	getdXcurve(curvex, sigma, smoothX, X, XX, g, dg, d2g, isOpen);
	getdXcurve(curvey, sigma, smoothY, Y, YY, g, dg, d2g, isOpen);
}

/* compute curvature of curve after gaussian smoothing
from "Shape similarity retrieval under affine transforms", Mokhtarian & Abbasi 2002
curvex - x position of points
curvey - y position of points
kappa - curvature coeff for each point
sigma - gaussian sigma
*/
void ComputeCurveCSS(const vector<double>& curvex,
	const vector<double>& curvey,
	vector<double>& kappa,
	vector<double>& smoothX, vector<double>& smoothY,
	double sigma,
	bool isOpen
	)
{
	vector<double> X, XX, Y, YY;
	SmoothCurve(curvex, curvey, smoothX, smoothY, X, XX, Y, YY, sigma, isOpen);

	kappa.resize(curvex.size());
	for (int i = 0; i<(int)curvex.size(); i++)
	{
		// Mokhtarian 02' eqn (4)
		kappa[i] = (X[i] * YY[i] - XX[i] * Y[i]) / pow(X[i] * X[i] + Y[i] * Y[i], 1.5);
	}
}

/* find zero crossings on curvature */
vector<int> FindCSSInterestPointsZero(const vector<double>& kappa)
{
	vector<int> crossings;
	for (int i = 0; i<(int)kappa.size() - 1; i++)
	{
		if ((kappa[i] < 0 && kappa[i + 1] > 0) || (kappa[i] > 0 && kappa[i + 1] < 0))
		{
			crossings.push_back(i);
		}
	}
	return crossings;
}
/* find zero crossings on curvature */
vector<int> FindCSSInterestPointsPeak(const vector<double>& kappa)
{
	vector<int> crossings;
	for (int i = 1; i<(int)kappa.size() - 1; i++)
	{
		
		if ((kappa[i] < kappa[i + 1] && kappa[i] < kappa[i - 1]) || (kappa[i] > kappa[i + 1] && kappa[i] > kappa[i - 1]))
		{
			crossings.push_back(i);
		}
		
	}
	return crossings;
}
vector<int> FindCSSInterestPointsPeak(const vector<double>& kappa, int n)
{
	vector<int> crossings;
	int m = round((n - 1) / 2.0);

	for (int i = m; i < (int)kappa.size() - m; i++)
	{
		bool done = true;
		int l = i + m;
		int p = i - m;
		for (int j = p; j <l; j++)
		{
			if (fabs(kappa[i])<fabs(kappa[j]))
			{
				done = false;	
			}

		}
		if (done == true)
		{
			crossings.push_back(i);
		}

	}

	return crossings;
}







vector<int> EliminateCloseMaximas(const vector<int>& maximasv, map<int, double>& maximas) {
	//eliminate degenerate segments (of very small length)
	vector<int> maximasvv;
	for (int i = 0; i<(int)maximasv.size(); i++)
	{
		if (i < (int)maximasv.size() - 1 &&
			maximasv[i + 1] - maximasv[i] <= 4)
		{
			//segment of small length (1 - 4) - eliminate one point, take largest sigma 
			if (maximas[maximasv[i]] > maximas[maximasv[i + 1]]) {
				maximasvv.push_back(maximasv[i]);
			}
			else {
				maximasvv.push_back(maximasv[i + 1]);
			}
			i++; //skip next element as well
		}
		else {
			maximasvv.push_back(maximasv[i]);
		}
	}
	return maximasvv;
}

/* compute the CSS image */
vector<int> ComputeCSSImageMaximas(const vector<double>& contourx_, const vector<double>& contoury_,
	vector<double>& contourx, vector<double>& contoury,
	bool isClosedCurve
	)
{
	ResampleCurve(contourx_, contoury_, contourx, contoury, 200, !isClosedCurve);
	vector<Point2d> pl; PolyLineMerge(pl, contourx, contoury);

	map<int, double> maximas;

	Mat_<Vec3b> img(500, 200, Vec3b(0, 0, 0)), contourimg(350, 350, Vec3b(0, 0, 0));
	bool done = false;
	//#pragma omp parallel for
	for (int i = 0; i<490; i++)
	{
		if (!done)
		{
			double sigma = 1.0 + ((double)i)*0.1;
			vector<double> kappa, smoothx, smoothy;
			ComputeCurveCSS(contourx, contoury, kappa, smoothx, smoothy, sigma);

			//			vector<vector<Point> > contours(1);
			//			PolyLineMerge(contours[0], smoothx, smoothy);
			//			contourimg = Vec3b(0,0,0);
			//			drawContours(contourimg, contours, 0, Scalar(255,255,255), CV_FILLED);

			vector<int> crossings = FindCSSInterestPointsZero(kappa);
			if (crossings.size() > 0)
			{
				for (int c = 0; c<crossings.size(); c++)
				{
					img(i, crossings[c]) = Vec3b(0, 255, 0);
					//					circle(contourimg, contours[0][crossings[c]], 5, Scalar(0,0,255), CV_FILLED);

					if (c < crossings.size() - 1) {
						if (fabs((float)crossings[c] - crossings[c + 1]) < 5.0)//fabs计算绝对值
						{
							//this is a maxima
							int idx = (crossings[c] + crossings[c + 1]) / 2;
							//#pragma omp critical
							maximas[idx] = (maximas[idx] < sigma) ? sigma : maximas[idx];

							circle(img, Point(idx, i), 3, Scalar(0, 0, 255), CV_FILLED);
						}
					}
				}
				//				char buf[128]; sprintf(buf, "evolution_%05d.png", i);
				//				imwrite(buf, contourimg);
				//				imshow("evolution", contourimg);
				//				waitKey(30);
			}
			else
			{
				done = true;
			}

		}
	}

	//find largest sigma
	double max_sigma = 0.0;
	for (map<int, double>::iterator itr = maximas.begin(); itr != maximas.end(); ++itr)
	{
		if (max_sigma < (*itr).second)
		{
			max_sigma = (*itr).second;
		}
	}
	//get segments with largest sigma
	vector<int> maximasv;
	for (map<int, double>::iterator itr = maximas.begin(); itr != maximas.end(); ++itr)
	{
		if ((*itr).second > max_sigma / 8.0)
		{
			maximasv.push_back((*itr).first);
		}
	}
	//eliminate degenerate segments (of very small length)
	vector<int> maximasvv = EliminateCloseMaximas(maximasv, maximas);	//1st pass
	maximasvv = EliminateCloseMaximas(maximasvv, maximas);				//2nd pass
	maximasv = maximasvv;
	for (vector<int>::iterator itr = maximasv.begin(); itr != maximasv.end(); ++itr) {
		cout << *itr << " - " << maximas[*itr] << endl;
	}
	//	Mat zoom; resize(img,zoom,Size(img.rows*2,img.cols*2));
	imshow("css image", img);
	//waitKey();
	return maximasv;
}
//#pragma mark Curve Matching

/* calculate the "centroid distance" for the curve */
void GetCurveSignature(const vector<Point2d>& a, vector<double>& signature) {
	signature.resize(a.size());
	Scalar a_mean = mean(a); Point2d a_mpt(a_mean[0], a_mean[1]);

	//centroid distance
	for (int i = 0; i<a.size(); i++) {
		signature[i] = norm(a[i] - a_mpt);
	}
}



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

CvMat* lsq_homog(vector<Point2d>& pts, vector<Point2d>& mpts, int n)
{
	CvMat* H, *A, *B, X;
	double x[9];//数组x中的元素就是变换矩阵H中的值
	int i;

	//输入点对个数不够4
	if (n < 4)
	{
		fprintf(stderr, "Warning: too few points in lsq_homog(), %s line %d\n",
			__FILE__, __LINE__);
		return NULL;
	}

	//将变换矩阵H展开到一个8维列向量X中，使得AX=B，这样只需一次解线性方程组即可求出X，然后再根据X恢复H
	/* set up matrices so we can unstack homography into X; AX = B */
	A = cvCreateMat(2 * n, 8, CV_64FC1);//创建2n*8的矩阵，一般是8*8
	B = cvCreateMat(2 * n, 1, CV_64FC1);//创建2n*1的矩阵，一般是8*1
	X = cvMat(8, 1, CV_64FC1, x);//创建8*1的矩阵，指定数据为x
	H = cvCreateMat(3, 3, CV_64FC1);//创建3*3的矩阵
	cvZero(A);//将A清零

	//由于是展开计算，需要根据原来的矩阵计算法则重新分配矩阵A和B的值的排列
	for (i = 0; i < n; i++)
	{
		cvmSet(A, i, 0, pts[i].x);//设置矩阵A的i行0列的值为pts[i].x
		cvmSet(A, i + n, 3, pts[i].x);
		cvmSet(A, i, 1, pts[i].y);
		cvmSet(A, i + n, 4, pts[i].y);
		cvmSet(A, i, 2, 1.0);
		cvmSet(A, i + n, 5, 1.0);
		cvmSet(A, i, 6, -pts[i].x * mpts[i].x);
		cvmSet(A, i, 7, -pts[i].y * mpts[i].x);
		cvmSet(A, i + n, 6, -pts[i].x * mpts[i].y);
		cvmSet(A, i + n, 7, -pts[i].y * mpts[i].y);
		cvmSet(B, i, 0, mpts[i].x);
		cvmSet(B, i + n, 0, mpts[i].y);
	}

	//调用OpenCV函数，解线性方程组
	cvSolve(A, B, &X, CV_SVD);//求X，使得AX=B
	x[8] = 1.0;//变换矩阵的[3][3]位置的值为固定值1
	X = cvMat(3, 3, CV_64FC1, x);
	cvConvert(&X, H);//将数组转换为矩阵

	cvReleaseMat(&A);
	cvReleaseMat(&B);
	return H;
}


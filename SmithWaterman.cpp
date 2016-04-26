



#include"SmithWaterman.h"
#include"CurveCSS.h"




//void  findMatchPointArray(vector<Point> & a, vector<Point> & b, int a_size, vector<Point> &out, int n)
//{
//
//	
//	double H = 0;
//	vector<Point> traceback;
//	
//	
//		vector<vector<double> > out_a, out_b;
//		double sigma = 5.0;
//		changeBPoint(b, 20);
//		//getcurveSignature(a, out_a,  n);
//		//getcurveSignature(b, out_b,  n);
//	   double  maxval=	MatchCurvesSmithWaterman(out_a, out_b, traceback);
//	   
//	   //if (H < maxval)
//	   //{
//		  // H = maxval;
//		  // out.clear();
//		  // out.insert(out.begin(), b.begin(), b.end());
//	   //}
//
//
//}


Mat ImageMerge(Mat& imgMat, Mat& imgMat1)
{






	CvSize size;
	size.width = imgMat.cols + imgMat1.cols + 10;
	size.height = (imgMat.rows > imgMat1.rows) ? imgMat.rows : imgMat1.rows;
	Mat   imgMat2 (size, CV_8U, Scalar(255));
	
	CvRect rect = cvRect(0, 0, imgMat.cols, imgMat.rows);
	Mat imageROI = imgMat2(rect);
	imgMat.copyTo(imageROI);
	
	/*CvRect    rect1 = cvRect(imgMat.cols + 10, 0, imgMat1.cols, imgMat1.rows);
	
	imgMat1.copyTo(imgMat2(rect1));*/
	return imageROI;
}
Mat ScanImageAndReduceC(Mat& imgMat, Mat & imgMat1)
{


	
	int col = imgMat.cols + imgMat1.cols + 10;
	int row = (imgMat.rows >= imgMat1.rows) ? imgMat.rows : imgMat1.rows;
	Mat end(row, col, imgMat.type());


	int channels = imgMat.channels();
	int nRows = imgMat.rows;
	int nCols = imgMat.cols* channels;
	//if (imgMat.isContinuous())
	//{
	//	nCols *= nRows;
	//	nRows = 1;
	//}
	
	uchar* p,*q;
	for (int i = 0; i < nRows; ++i)
	{
		p = end.ptr<uchar>(i);
		q = imgMat.ptr<uchar>(i);
		for (int j = 0; j < nCols; ++j)
		{
			p[j] = q[j];
		}
	}
	for (int i = 0; i < imgMat1.rows; ++i)
	{
		p = end.ptr<uchar>(i);
		q = imgMat1.ptr<uchar>(i);
		for (int j = 0; j < imgMat1.cols* channels; ++j)
		{
			p[j + 10 + nCols] = q[j];
		}
	}
	return end;
}
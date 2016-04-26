



#include"stdax.h"
#include"CurveCSS.h"
#include"CurveSignature.h"
#include"SmithWaterman.h"
#include <time.h> 
#include<stdlib.h>
#include <fstream>
//image
//unit(read)
//{
//	//image = cv::read()l..
//}
/////image
//UNIT(segmentation)
//{
//   //image
//}

void warp(vector<Point>& src, vector<Point>& std, Mat &H)
{
	std.clear();
	Point AA;
	for (int i = 0; i < src.size(); i++)
	{
		AA.x = round(H.at<double>(0, 0)*src[i].x + H.at<double>(0, 1)*src[i].y + H.at<double>(0, 2));
		AA.y = round(H.at<double>(1, 0)*src[i].x + H.at<double>(1, 1)*src[i].y + H.at<double>(1, 2));

		std.push_back(AA);
	}

}
void warp(vector<Point2d>& src, vector<Point>& std, Mat &H)
{
	std.clear();
	Point AA;
	for (int i = 0; i < src.size(); i++)
	{
		AA.x = round(H.at<double>(0, 0)*src[i].x + H.at<double>(0, 1)*src[i].y + H.at<double>(0, 2));
		AA.y = round(H.at<double>(1, 0)*src[i].x + H.at<double>(1, 1)*src[i].y + H.at<double>(1, 2));

		std.push_back(AA);
	}

}
#ifndef MAX_PATH
#define MAX_PATH          260
#endif

int main(int argc, char *argv[])
{


	//reading sample db 
	//for each frame k  and frame k+1
	  ///1.segmentation->object
	  ///2.contour
	  ///3.resample->points
	  ///4.for each point -> description
	  ///5.SW match -> points 1 - points 2
	  ///6.Using 3d points to compute transformation matrix(PCL ICP)
	  ///7.object patch rigisteration 

	//1.读取图像
		std::stringstream ImagesName1, ImagesName2;
		ImagesName1 <<"C:\\Users\\Lenovo\\Desktop\\1-2.jpg";
		ImagesName2 << "C:\\Users\\Lenovo\\Desktop\\1-1.jpg";
		Mat src = cv::imread(ImagesName1.str());
		Mat src1 = cv::imread(ImagesName2.str());
		if (src.empty() )
	{
		cout << "can't read image" << endl;
		exit(0);
	}

  //2.提取图像轮廓重采样	
	vector<Point> a, b, b_out;                          //a为图像1轮廓，b为图像2轮廓
	GetCurveForImage(src, a, false, false);
	//a重采样
	/*ResampleCurve(a, a,1000, false);*/
	//b重采样
	GetCurveForImage(src1, b, false, false);
	//ResampleCurve(b, b,1000, false);




	/*******************显示图像1轮廓高斯平滑后的曲线*************************/
	double nScaleStart=9;                    //与建立轮廓曲率描述符的起始数值相等
	vector<double>  kappa_; vector<Point>  small_smooth_;
	ComputeCurveCSS(b, kappa_, small_smooth_, nScaleStart, false);
	Mat bb_Curve(src.size(), CV_8UC3, Scalar(255, 255, 255));
	drawOpenCurve(bb_Curve, small_smooth_, Scalar(0, 255, 0), 1);
	namedWindow("平滑后的曲线", CV_WINDOW_AUTOSIZE);
	imshow("平滑后的曲线", bb_Curve); //显示轮廓图像
	/* **********************************************************/

	vector<Point> traceback, traceback2, traceback3;
	vector<vector<double> > out_a, out_b;

	vector<Keypoint>  a_KeyPoint, a_KeyPoint1;
	vector<Keypoint>   b_KeyPoint, b_KeyPoint1;

	//3.提取特征点

		a_KeyPoint.clear();
		b_KeyPoint.clear();
		traceback.clear();
		traceback2.clear();
		traceback3.clear();
		/*cout << "sigma= " << nScaleStart_<<endl;*/
		int N_neighbor_ = 3;                                 //特征提取时的临域
		int	nScale_ = 6;                                     //描述符的维度	
		double nScaleStart_ = nScaleStart;                   //起始尺度sigma
		int nBufferStep_ = 9;                                //曲率尺度间隔

		a_KeyPoint = GetKeypointDescriptor(a, N_neighbor_, nScale_, nBufferStep_, nScaleStart_);  //计算描述符
		b_KeyPoint = GetKeypointDescriptor(b, N_neighbor_, nScale_, nBufferStep_, nScaleStart_);

		cout << a_KeyPoint.size() << endl;



		//特征点配准


      //4.s-w方法
		double  maxval = MatchCurvesSmithWaterman(a_KeyPoint, b_KeyPoint, a_KeyPoint1, b_KeyPoint1, traceback);

		vector<Point2d> a_arry, b_arry;

		//5.获取配准序列 a_arry, b_arry
		for (int k = 0; k < traceback.size(); k++)
		{
			if (CalcCrossCorrelation(a_KeyPoint[traceback[k].y - 1].descriptor, b_KeyPoint[traceback[k].x - 1].descriptor)>0.3);
			{
				a_arry.push_back(a_KeyPoint[traceback[k].y - 1].positionS);
				b_arry.push_back(b_KeyPoint[traceback[k].x - 1].positionS);
			}
		}

		//最小二乘法获得转置矩阵   estimateRigidTransform
		Mat warp_mat(3, 3, CV_32FC1);
		warp_mat = Find2DRigidTransform(a_arry, b_arry);
		vector<Point> abc, abd;
		warp(a, abc, warp_mat);
		/*warp(a_arry, abd, warp_mat);*/
		Mat bbb_Curve(src1.size(), CV_8UC3, Scalar(255, 255, 255));
		drawOpenCurve(bb_Curve, abc, Scalar(0, 0, 255), 1);
		namedWindow("bbb_Curve", CV_WINDOW_AUTOSIZE);
		imshow("bbb_Curve", bb_Curve); //显示配准轮廓图像


	cv::waitKey(0);
	return 0;
}




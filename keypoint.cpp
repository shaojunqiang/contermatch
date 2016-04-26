

#include"keypoint.h"

vector<Keypoint>  GetKeypointDescriptor(vector<Point> & src, int N_neighbor, int nScale, int nBufferStep, int nScaleStart)
{
	bool done = false;
	vector<Keypoint> std;


			vector<double> kappa;
			vector<Point>  small_smooth;
			ComputeCurveCSS(src, kappa, small_smooth, nScaleStart, false);




			vector<int> crossings, crossings1;
			crossings = FindCSSInterestPointsZero(kappa);
			crossings1 = FindCSSInterestPointsPeak(kappa, N_neighbor);

			if (crossings.size() > 0 /*&& crossings1.size()<30*/)
			{
				int h = src.size();
				vector<vector<double> > out1(h);//曲率描述
				vector<vector<double> > globldiscript;//全局描述


				for (int h = 0; h < nScale; h++)
				{
					double sigma = (nScaleStart + ((double)h)*nBufferStep);
					vector<Point> small_smooth1;
					vector<double> kappa1, signature;
					ComputeCurveCSS(src, kappa1, small_smooth1, sigma, false);
					for (int j = 0; j < kappa.size(); j++)
					{
						out1[j].push_back(kappa1[j]);
					}

				}
				////计算全局描述符
				//globldiscript.clear();
				//for (int kk = 0; kk < crossings1.size(); kk++)
				//{
				//	vector<Point>  std;
				//	changeBPoint(small_smooth, std, crossings1[kk]);
				//	ResampleCurve(std, std, nScale, false);
				//	vector<double> distance, distance1;
				//	distance.clear();
				//	distance1.clear();
				//	double maxs=0;
				//	for (int nn = 1; nn < std.size();nn++)
				//	{
				//		double dist = sqrt((std[0].x - std[nn].x)*(std[0].x - std[nn].x) + (std[0].y - std[nn].y)*(std[0].y - std[nn].y));
				//		if (maxs<dist)
				//		{
				//			maxs = dist;
				//		}
				//		distance.push_back(dist);
				//	}
				//	
				//	for (int dd = 0; dd < std.size()-1; dd++)
				//	{
				//		double  dis = distance[dd] / maxs;
				//		distance1.push_back(dis);
				//	}
				//	globldiscript.push_back(distance1);

				//}
				for (int yy = 0; yy < crossings1.size(); yy++)
				{
					Keypoint nat;
					nat.position = src[crossings1[yy]];
					nat.positionS = small_smooth[crossings1[yy]];
					nat.descriptor = out1[crossings1[yy]];
					/*nat.globaldescriptor = globldiscript[yy];*/
					std.push_back(nat);
				}			
			}
			//if (crossings.size() == 0)
			//{
			//	done = true;
			//}
	return std;
}


//bijiao
/* from http://paulbourke.net/miscellaneous/correlate/ */
double CalcCrossCorrelation(const vector<double>& x, const vector<double>& y)
{
	assert(x.size() == y.size());
	int i, j, n = x.size();
	double mx, my, sx, sy, sxy, denom, r;

	/* Calculate the mean of the two series x[], y[] */
	mx = 0;
	my = 0;




	for (i = 0; i<n; i++) {
		mx += x[i];
		my += y[i];
	}
	mx /= n;
	my /= n;

	/* Calculate the denominator */
	sx = 0;
	sy = 0;
	for (i = 0; i<n; i++) {
		sx += (x[i] - mx) * (x[i] - mx);
		sy += (y[i] - my) * (y[i] - my);
	}
	denom = sqrt(sx*sy);

	/* Calculate the correlation series */
	////	for (delay=-maxdelay;delay<maxdelay;delay++) 
	int delay = 0;
	//{
		sxy = 0;
		for (i = 0; i<n; i++)
		{
			j = i + delay;
			//if (j < 0 || j >= n)
			//	continue;
			/*else*/
				sxy += (x[i] - mx) * (y[j] - my);
			/* Or should it be (?)
			if (j < 0 || j >= n)
			sxy += (x[i] - mx) * (-my);
			else
			sxy += (x[i] - mx) * (y[j] - my);
			*/
		}

		

			r = sxy / denom;
	
	//	/* r is the correlation coefficient at "delay" */
	//}
	return r;
}

/* calculate the similarity score between two curve segments
Mai 2010, "Affine-invariant shape matching and recognition under partial occlusion", section 4.1
*/
double MatchTwoSegments(const vector<double>& a_, const vector<double>& b_)
{
	double cc;
	vector<double> a = a_, b = b_;

	cc = CalcCrossCorrelation(a, b);


	//#if 0
	//	ShowMathGLCompareCurves(a_canon, b_canon, a_sig, b_sig, cc);
	//#endif

	return cc > 0.5? cc : 0.0;
}

/* match the two curves using adapted Smith-Waterman aligning algorithm
Mai 2010, "Affine-invariant shape matching and recognition under partial occlusion", section 4.2 */

Mat_<double> GetSmithWatermanHMatrix(const vector<Keypoint >& a, const vector<Keypoint >& b)
{
	int M = a.size() + 1;
	int N = b.size() + 1;

	//Smith-Waterman
	Mat_<double> H(M, N, 0.0);
	for (int i = 1; i < M; i++)
	{
		for (int j = 1; j < N; j++)
		{
			vector<double> v(4, 0.0);
			v[1] = H(i - 1, j - 1) + MatchTwoSegments(a[i - 1].descriptor, b[j - 1].descriptor);
			v[2] = H(i - 1, j) - 1.0;
			v[3] = H(i, j - 1) - 1.0;
			H(i, j) = *(max_element(v.begin(), v.end()));


		}
	}
	/*cout << H << endl;*/
	return H;
}


/* original Smith Waterman algorithm */
double MatchCurvesSmithWaterman(const vector<Keypoint >& a, const vector<Keypoint >& b,
	vector<Keypoint >& a_out, vector<Keypoint >& b_out,
	vector<Point>& traceback)
{
	
	Mat_<double> H = GetSmithWatermanHMatrix(a, b);
	Point maxp; double maxval;
	minMaxLoc(H, NULL, &maxval, NULL, &maxp);
	traceback.clear();
	while (H(maxp.y, maxp.x) != 0)
	{
		/*cout << "H(maxp.y-1,maxp.x-1) > H(maxp.y,maxp.x-1)" << H(maxp.y - 1, maxp.x - 1) << " > " << H(maxp.y, maxp.x - 1) << endl;*/
		if (H(maxp.y - 1, maxp.x - 1) > H(maxp.y, maxp.x - 1) &&
			H(maxp.y - 1, maxp.x - 1) > H(maxp.y - 1, maxp.x))
		{
			maxp = maxp - Point(1, 1);
			traceback.push_back(maxp);
		}
		else if (H(maxp.y - 1, maxp.x) > H(maxp.y - 1, maxp.x - 1) &&
			H(maxp.y - 1, maxp.x) > H(maxp.y, maxp.x - 1))
		{
			maxp.y--;
				traceback.push_back(maxp);
		}
		else
		if (H(maxp.y, maxp.x - 1) > H(maxp.y - 1, maxp.x - 1) &&
			H(maxp.y, maxp.x - 1) > H(maxp.y - 1, maxp.x))
		{
			maxp.x--;
				traceback.push_back(maxp);
		}
		else {
			break;
		}
	}

	//a_out.clear();
	//b_out.clear();
	//a_out.insert(a_out.begin(), a.begin() + traceback[0].y,a.end());
	//a_out.insert(a_out.end(), a.begin(), a.begin() + traceback[0].y);
	//b_out.insert(b_out.begin(), b.begin() + traceback[0].x, b.end());
	//b_out.insert(b_out.end(), b.begin(), b.begin() + traceback[0].x);


/*	for (int k = 0; k < traceback.size(); k++)
	{
		cout << traceback[k] << " -> ";
	}*/

	return maxval;
}



#include <iostream>
#include <opencv2/opencv.hpp>
using namespace cv;
using namespace std;

void resampling(Mat &img1, Mat &img2, double rate) 
{
	int x=0, y=0;
	
	for (y=0;y<img2.rows;y++) {
		for (x=0;x<img2.cols;x++) {
			int px = (int)(x/rate);
			int py = (int)(y/rate);

			double fx1 = (double)x/(double)rate - (double)px;
			double fx2 = 1 - fx1;
			double fy1 = (double)y/(double)rate - (double)py;
                        double fy2 = 1 - fy1;

			double w1 = fx2 * fy2;
			double w2 = fx1 * fy2;
			double w3 = fx2 * fy1;
			double w4 = fx1 * fy1;
			
			if (img1.channels() == 1)
			{
				uchar P1 = img1.at<uchar>(py,px);
				uchar P2 = img1.at<uchar>(py,px+1);
				uchar P3 = img1.at<uchar>(py+1,px);
				uchar P4 = img1.at<uchar>(py+1,px+1);
				img2.at<uchar>(y,x) = w1*P1 + w2*P2 + w3*P3 + w4
*P4;
			}
			else if (img1.channels() == 3)
			{
				Vec3b P1 = img1.at<Vec3b>(py,px);
				Vec3b P2 = img1.at<Vec3b>(py,px+1);
				Vec3b P3 = img1.at<Vec3b>(py+1,px+1);
				Vec3b P4 = img1.at<Vec3b>(py+1,px+1);
				
				img2.at<Vec3b>(y,x) = w1*P1 + w2*P2 + w3*P3 + w4*P4;
			}
		}
	}			
}

int main(int argc, char** argv)
{
    Mat image = imread("sample.jpeg");

    if (image.empty()) {
        cout << "Image file not found." << endl;

        return -1;
    }

    namedWindow("Input Image", WINDOW_AUTOSIZE);
    imshow("Input Image", image);

    int h = image.rows;
    int w = image.cols;

    int hNew = h * 2; 
    int wNew = w * 2;

    int imgColorState = (image.channels()==1) ? CV_8UC1 : CV_8UC3;
    int rate = 2;
    Mat result(hNew,wNew,imgColorState,Scalar(0));
    resampling(image,result,rate);
    
    namedWindow("Output Image", WINDOW_AUTOSIZE);
    imshow("Output Image",result); 
    
    cout << "Enter any key to exit." << endl;

    waitKey(0);
    destroyAllWindows();

    return 0;
}

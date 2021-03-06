#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <iostream>
#include <cstdio>
#include <libgen.h>

using namespace std;
using namespace cv;


int threshold_value = 190;
int threshold_type = 0;
int const max_value = 255;
int const max_type = 4; //     4: Threshold to Zero Inverted
  /* 0: Binary
     1: Binary Inverted
     2: Threshold Truncated
     3: Threshold to Zero
     4: Threshold to Zero Inverted
   */
int const max_BINARY_value = 255;

Mat src, src_gray, dst;
const char* trackbar_value = "Value";

//Scanner parameters - obtained through calibration (scan an image with known area)
const double pixelArea=0.282485876*0.282485876;

void apply_threshold( int, void* );

char faux[450];
char proc[450];

int main (int argc, char *argv[]) {
	
	const char * file = argv[1]; //File to open
	strcpy(faux,file);
	char *base = basename(faux);			//select output folder
	char *dir = dirname(faux);				//select output folder
	sprintf(proc, "%s/proc/%s", dir, base);	//select output folder

    Mat src = imread(file, CV_LOAD_IMAGE_COLOR); //open colored image
    
    if (src.empty())  { 
		fprintf(stderr,"Erro ao carregar imagem\n\n");
		return -2;
    }

	cvtColor( src, src_gray, CV_BGR2GRAY ); //converting image to gray scale
	
	namedWindow("Imagem 2",0);	//create a window to show image

	createTrackbar( trackbar_value,
					"Imagem 2", &threshold_value,
					max_value, apply_threshold );	//create a window to change threshould
    
    
    apply_threshold( 0, 0 ); //apply selected threshold

	return 0;
}

void apply_threshold( int, void* )
{
  /* 0: Binary
     1: Binary Inverted
     2: Threshold Truncated
     3: Threshold to Zero
     4: Threshold to Zero Inverted
   */

  threshold( src_gray, dst, threshold_value, max_BINARY_value,threshold_type );	//apply threshold
  
  int count_black = countNonZero(dst);
  int count=0;

  Size imSize = dst.size();
  for (int i = 0; i < imSize.height; i++) {		//counting image pixels
	for (int j = 0; j < imSize.width; j++) {	//counting image pixels

	   if(dst.at<uchar>(i,j) == 0) count++; //counting image pixels

	}
  }
  
  double area = count*pixelArea;	//calculate image area

  printf("%.6f\n", area/100);	//show image area
  
  
  imwrite(proc,dst);	//record processed image
  imshow( "Imagem 2", dst );	//show processed image
}

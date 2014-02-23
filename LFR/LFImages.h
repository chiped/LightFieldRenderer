#pragma once

#include <math.h>
#include <string>
#include <vector>

using namespace cv;
using namespace std;

class LFImages
{
public:
	static IplImage ***Images;
	static char* fileNames[16][16];

	static void readImages()
	{
		
		Images = (IplImage***)malloc(16*sizeof(IplImage**));
		for(int k=0;k<16;k++)
		{
			Images[k] = (IplImage**)malloc(16*sizeof(IplImage*));
			for(int l=0;l<16;l++)
			{
				Images[k][l] = (IplImage*)malloc(sizeof(IplImage));					
			}
		}
	}

	static void setPath(const char *path, int s, int t)
	{
		Images[s][t] = cvLoadImage(path);//path;
	}
};
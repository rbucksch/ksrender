#include "AllConfig.h"
#include "BMP.h"
#include "Wall.h"

#include <algorithm>
#include <iostream>
#include <cmath>

static constexpr double Ebene[6][4] = { {  1.0,  0.0,  0.0, 0.5 },
										{  0.0,  0.0,  1.0, 0.5 },
										{  0.0, -1.0,  0.0, 0.5 },
										{  0.0,  0.0, -1.0, 0.5 },
										{ -1.0,  0.0,  0.0, 0.5 },
										{  0.0,  1.0,  0.0, 0.5 } };

static constexpr double Bild[6][3] = { {  0.5,  0.5, -0.5 },
									   {  0.5,  0.5,  0.5 },
									   {  0.5, -0.5,  0.5 },
									   {  0.5, -0.5, -0.5 },
									   { -0.5, -0.5, -0.5 },
									   { -0.5,  0.5, -0.5 } };

Wall::Wall()
{

}

Wall::Wall(int index)
{
	ID = index;

	OS = Ebene[index][3];

	for(int i=0; i<3; ++i) {
		N[i] = Ebene[index][i];

		Start[i] = Bild[index][i];

		Dir[i] = Bild[(index + 1) % 6][i] - Start[i];
	}

	hPixels = 8;
	vPixels = 8;
	padPixels = 24;

	switch (index) {
	case 0:
		color1[0] = 180;	color1[1] = 180;	color1[2] = 255;
		color2[0] = 0;		color2[1] = 0;		color2[2] = 255;
		break;
	case 1:
		color1[0] = 255;	color1[1] = 255;	color1[2] = 180;
		color2[0] = 255;	color2[1] = 255;	color2[2] = 0;
		break;
	case 2:
		color1[0] = 0;		color1[1] = 255;	color1[2] = 255;
		color2[0] = 180;	color2[1] = 255;	color2[2] = 255;
		break;
	case 3:
		color1[0] = 255;	color1[1] = 180;	color1[2] = 180;
		color2[0] = 255;	color2[1] = 0;		color2[2] = 0;
		break;
	case 4:
		color1[0] = 255;	color1[1] = 180;	color1[2] = 255;
		color2[0] = 255;	color2[1] = 0;		color2[2] = 255;
		break;
	case 5:
		color1[0] = 180;	color1[1] = 255;	color1[2] = 180;
		color2[0] = 0;		color2[1] = 255;	color2[2] = 0;
	}
}

Wall::~Wall()
{

}

double Wall::dotp(double* a, double*b)
{
	double d = 0;
	for(int i=0; i<3; ++i)
		d += a[i]*b[i];

	return d;
}

double Wall::length(double* a)
{
	return sqrt(dotp(a,a));
}

double Wall::distance(double* V)
{
	return N[0]*V[1] + N[1]*V[2] + N[2]*V[3] - OS;
}

double Wall::endpoint(double* now, double* last)
{
	return (OS - N[0]*last[1] - N[1]*last[2] - N[2]*last[3])
		/ (N[0]*(now[1] - last[1]) + N[1]*(now[2] - last[2]) + N[2]*(now[3] - last[3]));
}

void Wall::get_BGR(unsigned char* BGR, double* data)
{
	double hitVec[3] = { data[1] - Start[0],
						 data[2] - Start[1],
						 data[3] - Start[2] };

	double hitDist = length(hitVec);
	double xHit = dotp(hitVec, Dir);
	double yHit = hitDist * sqrt(1.0 - pow(xHit/hitDist,2));

	int xPix = (int)floor(hPixels*xHit);
	int yPix = (int)ceil(vPixels*yHit);

	if ((yPix + xPix) % 2)
	{
		std::copy_n(color1, 3, BGR);
	}
	else
	{
		std::copy_n(color2, 3, BGR);
	}
}

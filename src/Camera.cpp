#include "AllConfig.h"
#include "Camera.h"
#include "KerrSchild.h"

//comment

#include <algorithm>
#include <cmath>
#include <iostream>

Camera::Camera(CamData &c)
{
	std::copy_n(c.CamPos, 3, camPos);
	std::copy_n(c.CamDir, 3, camDir);

	hPixels = c.hPixels;
	vPixels = c.vPixels;
	padPixels = (c.hPixels * 3 + 3) & (~3);

	scrImage = new unsigned char[vPixels*padPixels];

	scrWidth = 1.0;
	scrHeight = (double)vPixels / (double)hPixels;
	scrDistance = 0.375;

	get_hVec();
	get_vVec();
	get_scrStart();

	scaleVec(hVec, 1.0/hPixels);
	scaleVec(vVec, 1.0/hPixels);
}

Camera::~Camera()
{
	delete[] scrImage;
	scrImage = NULL;
}

void Camera::scaleVec(double* Vec, double d)
{
	double L = d / sqrt(Vec[0]*Vec[0] + Vec[1]*Vec[1] + Vec[2]*Vec[2]);

	for (int i=0; i<3; ++i)
	{
		Vec[i] *= L;
	}
}

void Camera::get_hVec()
{
	double d = (0.5-0.5/hPixels)*scrWidth / sqrt(camDir[0]*camDir[0]+camDir[1]*camDir[1]);

	hVec[0] = camDir[1] * d;
	hVec[1] = -camDir[0] * d;
	hVec[2] = 0.0;
}

void Camera::get_vVec()
{
	double d = (0.5 - 0.5/vPixels)*scrHeight;

	vVec[0] = camDir[1]*hVec[2] - camDir[2]*hVec[1];
	vVec[1] = camDir[2]*hVec[0] - camDir[0]*hVec[2];
	vVec[2] = camDir[0]*hVec[1] - camDir[1]*hVec[0];

	scaleVec(vVec, d);
}

void Camera::get_scrStart()
{
	scaleVec(camDir, 1.0);

	for (int i=0; i<3; ++i) {
		scrStart[i] = scrDistance*camDir[i] - hVec[i] - vVec[i];
	}
}

void Camera::get_pixDir(double* Vec, int hPos, int vPos)
{
	for(int i=0; i<3; ++i) {
		Vec[i] = scrStart[i] + hPos*hVec[i] + vPos*vVec[i];
	}
}

void Camera::InitialConditions(Kerrschild &H, double* dataVec, int vPos, int hPos)
{
	double Vec[3];
	get_pixDir(Vec, hPos, vPos);

	dataVec[0] = 0.0;
	dataVec[1] = camPos[0];
	dataVec[2] = camPos[1];
	dataVec[3] = camPos[2];
	dataVec[4] = 1.0;
	dataVec[5] = Vec[0];
	dataVec[6] = Vec[1];
	dataVec[7] = Vec[2];

	H.tiltDown(dataVec);
	H.updateGamma(dataVec);
	double scaler = -1.0 / H.gammaScale(dataVec + 5);
	H.tiltUp(dataVec);

	for (int i=5; i<8; ++i) {
		dataVec[i] *= scaler;
	}
}

{
	
}

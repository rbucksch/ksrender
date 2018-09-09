#pragma once

#include "KerrSchild.h"

struct CamData
{
	double CamPos[3];
	double CamDir[3];
	int hPixels, vPixels;
};

class Camera
{
private:
	double camPos[3];
	double scrStart[3];
	double hVec[3];
	double vVec[3];
	double camDir[3];

	int hPixels;
	int vPixels;
	double scrWidth;
	double scrHeight;
	double scrDistance;
	int padPixels;
public:
	Camera(CamData &cam);
	~Camera();

	unsigned char* scrImage;

	int get_hPixels() { return hPixels; }
	int get_vPixels() { return vPixels; }
	int get_padPixels() { return padPixels; }

	void scaleVec(double* Vec, double d);

	void get_hVec();
	void get_vVec();

	void get_scrStart();

	void get_pixDir(double* Vec, int hPos, int vPos);

	void InitialConditions(Kerrschild &H, double* dataVec, int vPos, int hPos);
};

#include "AllConfig.h"
#include "BMP.h"
#include "Camera.h"
#include "DormandPrince.h"
#include "KerrSchild.h"
#include "Raytracing.h"
#include "Wall.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

void traceRays(Camera &Cam, Wall Bounds[])
{
	unsigned char* img = Cam.scrImage;

	int hPixels = Cam.get_hPixels();
	int vPixels = Cam.get_vPixels();
	int padPixels = Cam.get_padPixels();

	int colorIter = 0;
	unsigned char BGR[64];

	for (int k = 0; k<vPixels; ++k) {
		for (int i = 0; i<hPixels; ++i) {

			int ObstacleHit = -1;
			double T = 0.0;
			double dt = 0.25;

			double dydx[8];
			double daten[8]; //array for calculation

			Kerrschild Blackhole;

			Cam.InitialConditions(Blackhole, daten, k, i);

			double last[8]; //saves the last point calculated

			Blackhole.tiltDown(daten);
			dp_derivs(Blackhole, daten, dydx);
			Blackhole.tiltUp(daten);

			while(true) {

				std::copy_n(daten, 8, last);

				Blackhole.tiltDown(daten);
				dopri(Blackhole, daten, dydx, dt, T, dopri_error);
				Blackhole.tiltUp(daten);
				
				ObstacleHit = ObstacleImpact(daten, last, Bounds);
				if (ObstacleHit != -1) {
					Bounds[ObstacleHit].get_BGR(BGR + colorIter, daten);
					break;
				}
				else if (T > T_MAX) {
					BGR[colorIter] = BGR[colorIter + 1] = BGR[colorIter + 2] = 0;
					break;
				}
			}

			if (colorIter == 57) {
				colorIter = 0;
				for (int color = 0; color<60; ++color) {
					img[k*padPixels + i*3 - 57 + color] = BGR[color];
				}
			}
			else
				colorIter += 3;

		}
	}
}

int ObstacleImpact(double* now, double* last, Wall Bounds[])
{
	double lin;

	double lin_temp = 0;
	int curWall = -1;

	for (int wallCount = 0; wallCount<6; ++wallCount) {
		if (Bounds[wallCount].distance(now)*Bounds[wallCount].distance(last) <= 0) {
			lin_temp = Bounds[wallCount].endpoint(now, last);

			if (curWall == -1) {
				if (lin_temp > 0.0 && lin_temp <= 1.0) {
					lin = lin_temp;
					curWall = wallCount;
				}
			}
			else {
				if (lin_temp > 0.0 && lin_temp <= 1.0 && lin_temp < lin) {
					lin = lin_temp;
					curWall = wallCount;
				}
			}
		}
	}

	if (curWall > -1) {
		for (int h = 0; h<8; h++) {
			now[h] = last[h] + lin * (now[h] - last[h]);
		}
	}

	return curWall;
}

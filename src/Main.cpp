#include "AllConfig.h"
#include "BMP.h"
#include "Camera.h"
#include "Raytracing.h"
#include "Wall.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace chrono;

int main(void)
{
	CamData c;
	c.CamPos[0] = -0.5;	c.CamPos[1] = 0.0;	c.CamPos[2] = 0.0;
	c.CamDir[0] = 1.0;	c.CamDir[1] = 0.0;	c.CamDir[2] = 0.0;
	c.hPixels = PIXEL_WIDTH;
	c.vPixels = PIXEL_HEIGHT;
	Camera Cam(c);

	Wall Bounds[6] = { Wall(0), Wall(1), Wall(2), Wall(3), Wall(4), Wall(5) };

	steady_clock::time_point begin = steady_clock::now();
	traceRays(Cam, Bounds);
	steady_clock::time_point end = steady_clock::now();	

	writeIMG(Cam, "C:/dev/ksOutput/img.bmp");

	duration<double> time_span = duration_cast<duration<double>>(end - begin);
	cout << "Elapsed Time = " << time_span.count() << "s\n";

	std::cin.get();
}

#pragma once

class Wall
{
public:
	Wall();
	Wall(int index);
	~Wall();

	double dotp(double* a, double* b);
	double length(double* a);
	double distance(double* a);
	double endpoint(double* now, double* last);

	int getID() { return ID; }

	void get_BGR(unsigned char* BGR, double* data);

private:
	int ID;

	double OS;  // offset of plane
	double N[3];  // normal vector of plane
	double Start[3]; // starting corner of picture
	double Dir[3]; // direction of the picture's bottom edge

	int hPixels;
	int vPixels;
	int padPixels;

	unsigned char color1[3], color2[3];
};

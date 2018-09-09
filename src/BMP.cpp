#include "BMP.h"
#include "Camera.h"

#include <iostream>
#include <fstream>
#include <stdint.h>

void writeIMG(Camera &Cam, const char* filename)
{

	int width = Cam.get_hPixels();
	int height = Cam.get_vPixels();
	unsigned char* image = Cam.scrImage;

	int width_padded = (width * 3 + 3) & (~3);

	int DATASIZE = sizeof(unsigned char) * height*width_padded;

	struct BITMAPFILEHEADER
	{
		uint16_t   bfType;
		uint32_t   bfSize;
		uint32_t   bfReserved;
		uint32_t   bfOffBits;
	} bmfh;

	struct BITMAPINFOHEADER
	{
		uint32_t  biSize;
		int32_t   biWidth;
		int32_t   biHeight;
		uint16_t  biPlanes;
		uint16_t  biBitCount;
		uint32_t  biCompression;
		uint32_t  biSizeImage;
		int32_t   biXPelsPerMeter;
		int32_t   biYPelsPerMeter;
		uint32_t  biClrUsed;
		uint32_t  biClrImportant;
	} bmih;

	bmfh.bfType = 0x4d42;
	bmfh.bfSize = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + DATASIZE;
	bmfh.bfReserved = 0;
	bmfh.bfOffBits = 0;

	bmih.biSize = sizeof(BITMAPINFOHEADER);
	bmih.biWidth = width;
	bmih.biHeight = height;
	bmih.biPlanes = 1;
	bmih.biBitCount = 24;
	bmih.biCompression = 0;
	bmih.biSizeImage = 0;
	bmih.biXPelsPerMeter = 0;
	bmih.biYPelsPerMeter = 0;
	bmih.biClrUsed = 0;
	bmih.biClrImportant = 0;

	std::fstream myfile (filename, std::fstream::out | std::fstream::binary );

	myfile.write(reinterpret_cast<char*>(&bmfh), 14);
	myfile.write(reinterpret_cast<char*>(&bmih), 40);
	for (int j = height - 1; j >= 0; j--) {
		myfile.write(reinterpret_cast<char*>(image + j * width_padded), width_padded * sizeof(unsigned char));
	}

	myfile.close();
}

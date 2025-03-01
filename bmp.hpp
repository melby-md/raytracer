#ifndef _BMP_HPP
#define _BMP_HPP
#include <stdio.h>

#include "common.hpp"

void WriteBMP(const char *file_name, int w, int h, byte *data)
{
	int filesize = 54 + 3*w*h;

	byte bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
	byte bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
	byte bmppad[3] = {0,0,0};

	bmpfileheader[ 2] = (byte)(filesize    );
	bmpfileheader[ 3] = (byte)(filesize>> 8);
	bmpfileheader[ 4] = (byte)(filesize>>16);
	bmpfileheader[ 5] = (byte)(filesize>>24);

	bmpinfoheader[ 4] = (byte)(       w    );
	bmpinfoheader[ 5] = (byte)(       w>> 8);
	bmpinfoheader[ 6] = (byte)(       w>>16);
	bmpinfoheader[ 7] = (byte)(       w>>24);
	bmpinfoheader[ 8] = (byte)(       h    );
	bmpinfoheader[ 9] = (byte)(       h>> 8);
	bmpinfoheader[10] = (byte)(       h>>16);
	bmpinfoheader[11] = (byte)(       h>>24);

	FILE *f = fopen(file_name, "wb");

	fwrite(bmpfileheader,1,14,f);
	fwrite(bmpinfoheader,1,40,f);
	for(int i=0; i<h; i++)
	{
		fwrite(data+(w*(h-i-1)*3),3,w,f);
		fwrite(bmppad,1,(4-(w*3)%4)%4,f);
	}

	fclose(f);
}

#endif

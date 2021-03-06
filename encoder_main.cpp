#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <bitset>
#include <math.h>


typedef __int16_t MyDataSize;
using namespace std;
typedef unsigned char BYTE;
typedef struct Block
{
    BYTE data[8][8];
} Block;

typedef struct Block16x16
{
    BYTE data[16][16];
	int x;
	int y;
} Block16x16;

typedef struct MV
{
    MyDataSize x;
    MyDataSize y;
} MV;

// Get the size of a file
#define mnint(a) ((a) < 0 ? (int)(a - 0.5) : (int)(a + 0.5))

//QCIF frame resolution
#define Y_frame_width 176
#define Y_frame_height 144
#define SIZE_IN_MB 3

//function to get the size of a file
//input file: file pointer
//returns file size
long getFileSize(FILE *file)
{
    long lCurPos, lEndPos;
    lCurPos = ftell(file);
    fseek(file, 0, 2);
    lEndPos = ftell(file);
    fseek(file, lCurPos, 0);
    return lEndPos;
}

//global variables
long total_number_of_frames;
int number_of_macroblocks_per_frame;
int y_buffer_size_bytes;
int u_buffer_size_bytes;
int v_buffer_size_bytes;
int framesize;
int numOfPFrames = 300;

int zigzag_order[8][8] = {
    {0, 1, 5, 6, 14, 15, 27, 28},
    {2, 4, 7, 13, 16, 26, 29, 42},
    {3, 8, 12, 17, 25, 30, 41, 43},
    {9, 11, 18, 24, 31, 40, 44, 53},
    {10, 19, 23, 32, 39, 45, 52, 54},
    {20, 22, 33, 38, 46, 51, 55, 60},
    {21, 34, 37, 47, 50, 56, 59, 61},
    {35, 36, 48, 49, 57, 58, 62, 63},
};

//////////////////////////////// ENCODE //////////////////////////

void writeMotionVector(MV mv,FILE *output_file)
{
    MyDataSize output_value;
    
	output_value = (MyDataSize)mv.x;
	//writing the value of mv.x
	fwrite(&output_value, sizeof(MyDataSize), 1, output_file);

    output_value = (MyDataSize)mv.y;
	//writing the value of mv.y
	fwrite(&output_value, sizeof(MyDataSize), 1, output_file);
}
void writeQuantization(MyDataSize Quant_parameter,FILE *output_file)
{
    MyDataSize output_value;
    
	output_value = (MyDataSize)Quant_parameter;
	//writing the value of mv.x
	fwrite(&output_value, sizeof(MyDataSize), 1, output_file);
}

void writeEncodedMacroblock(vector<int> rle, FILE *output_file)
{
    MyDataSize output_value;
    output_value = (MyDataSize)rle.size();
    //Writing the Size of each sequence before writing to the file
    fwrite(&output_value, sizeof(MyDataSize), 1, output_file);

    for (int j = 0; j < rle.size(); j++)
    {
	output_value = (MyDataSize)rle[j];
	fwrite(&output_value, sizeof(MyDataSize), 1, output_file);
    }
}


vector<int> Compute_VLC(int block[8][8])
{
    vector<int> rle;
    int zigzag_scaned_values[64];

    //first order coefficients in a single dimension array using zigzag scan order
    //coefficients are originally in blocks zigzag_scaned_values array
    //ordered coefficients will be in zigzag_scaned_values
    int zigzag_ordered_index;
    for (int col_index = 0; col_index < 8; col_index++)
	for (int row_index = 0; row_index < 8; row_index++)
	{
	    zigzag_ordered_index = zigzag_order[col_index][row_index];
	    zigzag_scaned_values[zigzag_ordered_index] = block[col_index][row_index];
	}

    //now perform run length coding on ordered coefficients in zigzag_scaned_values
    int run_num; //number of zeros before a non-zero value

    run_num = 0;
    rle.clear();

    for (int coefficient_index = 0; coefficient_index < 64; coefficient_index++)
    {
	if (zigzag_scaned_values[coefficient_index] == 0)
	{
	    run_num++;
	}
	else
	{
	    //store run-level pair
	    rle.push_back(run_num);
	    rle.push_back(zigzag_scaned_values[coefficient_index]);
	    //reset run counter
	    run_num = 0;
	}
    }

    return rle;
}

//Compute_Quantization
// Quantize an input block
//Use a fixed quantization of 8
//Input   8x8 DCT-Transformed block
//Output  8x8 Quantized block
void Compute_quantization(int inMatrix[8][8], int outMatrix[8][8], int Quant_parameter)
{
    for (int i = 0; i < 8; i++)
    {
		for (int j = 0; j < 8; j++)
		{
			outMatrix[i][j] = inMatrix[i][j] / Quant_parameter;
		}
    }
}


//Compute_SAD
// computes sad value between two blocks
//Input   16x16 block1
//Output  16x16 block2
//return  int sadValue
int Compute_SAD(BYTE block1[16][16], BYTE block2[16][16])
{
    int sadValue = 0;
    for (int i = 0; i < 16; i++)
    {
		for (int j = 0; j < 16; j++)
		{
			BYTE value1 = block1[i][j];
			BYTE value2 = block2[i][j];
			sadValue += abs(value1 - value2);
		}
    }
    return sadValue;
}

Block get8x8Block(BYTE *frame, int frameWidth, int startI, int startJ)
{
    Block block;
    for (int j = 0; j < 8; j++)
    {
		for (int i = 0; i < 8; i++)
		{
			block.data[i][j] = *(frame + (startI + i) + (startJ + j) * frameWidth);
		}
    }
    return block;
}

void set8x8Block(int block[8][8], BYTE *frame, int frameWidth, int offset)
{
    for (int j = 0; j < 8; j++)
    {
		for (int i = 0; i < 8; i++)
		{
			BYTE tobeSetValue = (BYTE)block[i][j];
			int pixelOffset = offset + i + j * frameWidth;
			*(frame +pixelOffset) = tobeSetValue; 
		}
    }
}

//startI, startJ are assumed to be multiples of 16
Block16x16 get16x16Block(BYTE *frame, int frameWidth, int startI, int startJ)
{
    Block16x16 block;
    for (int j = 0; j < 16; j++)
    {
		for (int i = 0; i < 16; i++)
		{
			block.data[i][j] = *(frame + (startI + i) + (startJ + j) * frameWidth);
		}
    }
	block.x = startI;
	block.y = startJ;
    return block;
}

//compute diffences
//matrix1 current frame block
//matrix2 deoceded reference frame  block
//outMatrix output 
void Compute_Diffrences(BYTE matrix1[8][8], BYTE matrix2[8][8], BYTE outMatrix[8][8])
{
    for (int i = 0; i < 8; i++)
    {
		for (int j = 0; j < 8; j++)
		{
			int val1 = matrix1[i][j];
			int val2 = matrix2[i][j];
			// BYTE res = (val1 > val2) ? val1- val2: 127+ (val2- val1);
			// BYTE res = (val1 > val2) ? val1- val2: 127+ (val1- val2);
			int diff = (val1- val2);
			bitset<32> x(diff);
			// cout << x << endl;
			
			int res = 127+ ((diff>>1));
			
			bitset<8> y(res);
			// cout << y << endl;
		
			outMatrix[i][j] = (BYTE)res;
			
			
			// outMatrix[i][j] = matrix1[i][j]- matrix2[i][j];
		}
    }
}

void Compute_Additions(int matrix1[8][8], BYTE matrix2[8][8], int outMatrix[8][8])
{
    for (int i = 0; i < 8; i++)
    {
		for (int j = 0; j < 8; j++)
		{
			BYTE asd = matrix1[i][j];
			bitset<8> x(1);
			// cout << x << endl;
			
			bitset<8> y(asd);
			// cout << y << endl;
			
			int16_t btn = -127;
			int16_t asds = btn+asd;
			int res =  asds<<1;
			outMatrix[i][j] = res + matrix2[i][j];
			outMatrix[i][j] = outMatrix[i][j]> 255? 255: outMatrix[i][j] < 0 ? 0: outMatrix[i][j];
		}
    }
}


MV LDSP(Block16x16 block, BYTE *frame, int isLDSP, int lPoints[9][2], int sPoints[5][2])
{
	int size;
	int (*points)[2];
	if(isLDSP){
		points = lPoints;
		size = 9;
	}
	else{
		points = sPoints;
		size = 5;
	}
	Block16x16 frameBlock;
	int offsetX;
	int offsetY;
	int minValue = -1;
    int minI = 0;
	for (int8_t i = 0; i < size; i++)
	{
		offsetX = points[i][0]*16;
		offsetY = points[i][1]*16;
		if(block.x + offsetX >= 0 && block.x + offsetX + 16 < Y_frame_width && block.y + offsetY >= 0 && block.y + offsetY + 16 < Y_frame_height)
		{
			frameBlock = get16x16Block(frame, Y_frame_width, block.x+offsetX, block.y+offsetY);
			int sad = Compute_SAD(block.data, frameBlock.data);
			if (minValue < 0 ||sad < minValue)
			{
				minI = i;
				minValue = sad;
			}
		}
	}
	block.x += points[minI][0]*16;
	block.y += points[minI][1]*16;
			
	if (minI == 0)
		if (isLDSP)
			return LDSP(block, frame, 0, lPoints, sPoints);//resturn SDSP
		else{
			MV mv = {(MyDataSize)(block.x/16), (MyDataSize)(block.y/16)};
    		return mv;
		}
	else 
		if (isLDSP)
			return LDSP(block, frame, isLDSP, lPoints, sPoints);
		else{
			MV mv = {(MyDataSize)(block.x/16), (MyDataSize)(block.y/16)};
    		return mv;
		}
}
MV Compute_MV(Block16x16 block, BYTE *frame)
{
	int lPoints[9][2] = {{0,0},{2,0},{0,2},{-2,0},{0,-2},{1,1},{-1,-1},{-1,1},{1,-1}};
	int sPoints[5][2] = {{0,0},{1,0},{0,1},{-1,0},{0,-1}};
	return LDSP(block, frame, 1, lPoints, sPoints);
}

//DCT on 8x8 block
//Input   8x8 spatial block
//Output  8x8 DCT-Transformed block
void Compute_DCT(BYTE input_block[8][8], int output_block[8][8])
{
    BYTE block[64];
    double coeff[64];

    int m = 0;
    for (int i = 0; i < 8; i++)
    {
	for (int j = 0; j < 8; j++)
	{
	    block[m] = input_block[i][j];
	    m++;
	}
    }

    int j1, i, j, k;
    vector<double> b(8);
    vector<double> b1(8);
    vector<vector<double> > d(8, vector<double>(8));
    double f0 = .7071068;
    double f1 = .4903926;
    double f2 = .4619398;
    double f3 = .4157348;
    double f4 = .3535534;
    double f5 = .2777851;
    double f6 = .1913417;
    double f7 = .0975452;

    for (i = 0, k = 0; i < 8; i++, k += 8)
    {
	for (j = 0; j < 8; j++)
	{
	    b[j] = block[k + j];
	}
	/* Horizontal transform */
	for (j = 0; j < 4; j++)
	{
	    j1 = 7 - j;
	    b1[j] = b[j] + b[j1];
	    b1[j1] = b[j] - b[j1];
	}
	b[0] = b1[0] + b1[3];
	b[1] = b1[1] + b1[2];
	b[2] = b1[1] - b1[2];
	b[3] = b1[0] - b1[3];
	b[4] = b1[4];
	b[5] = (b1[6] - b1[5]) * f0;
	b[6] = (b1[6] + b1[5]) * f0;
	b[7] = b1[7];
	d[i][0] = (b[0] + b[1]) * f4;
	d[i][4] = (b[0] - b[1]) * f4;
	d[i][2] = b[2] * f6 + b[3] * f2;
	d[i][6] = b[3] * f6 - b[2] * f2;
	b1[4] = b[4] + b[5];
	b1[7] = b[7] + b[6];
	b1[5] = b[4] - b[5];
	b1[6] = b[7] - b[6];
	d[i][1] = b1[4] * f7 + b1[7] * f1;
	d[i][5] = b1[5] * f3 + b1[6] * f5;
	d[i][7] = b1[7] * f7 - b1[4] * f1;
	d[i][3] = b1[6] * f3 - b1[5] * f5;
    }
    /* Vertical transform */
    for (i = 0; i < 8; i++)
    {
	for (j = 0; j < 4; j++)
	{
	    j1 = 7 - j;
	    b1[j] = d[j][i] + d[j1][i];
	    b1[j1] = d[j][i] - d[j1][i];
	}
	b[0] = b1[0] + b1[3];
	b[1] = b1[1] + b1[2];
	b[2] = b1[1] - b1[2];
	b[3] = b1[0] - b1[3];
	b[4] = b1[4];
	b[5] = (b1[6] - b1[5]) * f0;
	b[6] = (b1[6] + b1[5]) * f0;
	b[7] = b1[7];
	d[0][i] = (b[0] + b[1]) * f4;
	d[4][i] = (b[0] - b[1]) * f4;
	d[2][i] = b[2] * f6 + b[3] * f2;
	d[6][i] = b[3] * f6 - b[2] * f2;
	b1[4] = b[4] + b[5];
	b1[7] = b[7] + b[6];
	b1[5] = b[4] - b[5];
	b1[6] = b[7] - b[6];
	d[1][i] = b1[4] * f7 + b1[7] * f1;
	d[5][i] = b1[5] * f3 + b1[6] * f5;
	d[7][i] = b1[7] * f7 - b1[4] * f1;
	d[3][i] = b1[6] * f3 - b1[5] * f5;
    }
    for (i = 0; i < 8; i++)
    {
	for (j = 0; j < 8; j++)
	{
	    coeff[j + i * 8] = (d[i][j]);
	}
    }
    m = 0;
    for (int n = 0; n < 8; n++)
    {
	for (int z = 0; z < 8; z++)
	{
	    output_block[n][z] = (int)coeff[m];
	    m++;
	}
    }
}

//Inverse DCT function
void Compute_idct(int inblock[8][8], int outblock[8][8])
{

    vector<int> coeff(64);
    int m = 0;
    for (int i = 0; i < 8; i++)
    {
	for (int j = 0; j < 8; j++)
	{
	    coeff[m] = inblock[i][j];
	    m++;
	}
    }

    int j1, i, j;
    double b[8], b1[8], d[8][8];
    double f0 = .7071068;
    double f1 = .4903926;
    double f2 = .4619398;
    double f3 = .4157348;
    double f4 = .3535534;
    double f5 = .2777851;
    double f6 = .1913417;
    double f7 = .0975452;
    double e, f, g, h;
    int block[64];
    /* Horizontal */
    //    fprintf(Debugfptr,"Before IDCT\n");
    //	PrintMB(coeff,64);
    for (i = 0; i < 8; i++)
    {
	for (j = 0; j < 8; j++)
	    b[j] = coeff[j + i * 8];

	e = b[1] * f7 - b[7] * f1;
	h = b[7] * f7 + b[1] * f1;
	f = b[5] * f3 - b[3] * f5;
	g = b[3] * f3 + b[5] * f5;

	b1[0] = (b[0] + b[4]) * f4;
	b1[1] = (b[0] - b[4]) * f4;
	b1[2] = b[2] * f6 - b[6] * f2;
	b1[3] = b[6] * f6 + b[2] * f2;
	b[4] = e + f;
	b1[5] = e - f;
	b1[6] = h - g;
	b[7] = h + g;

	b[5] = (b1[6] - b1[5]) * f0;
	b[6] = (b1[6] + b1[5]) * f0;
	b[0] = b1[0] + b1[3];
	b[1] = b1[1] + b1[2];
	b[2] = b1[1] - b1[2];
	b[3] = b1[0] - b1[3];

	for (j = 0; j < 4; j++)
	{
	    j1 = 7 - j;
	    d[i][j] = b[j] + b[j1];
	    d[i][j1] = b[j] - b[j1];
	}
    }

    /* Vertical */

    for (i = 0; i < 8; i++)
    {
	for (j = 0; j < 8; j++)
	{
	    b[j] = d[j][i];
	}
	e = b[1] * f7 - b[7] * f1;
	h = b[7] * f7 + b[1] * f1;
	f = b[5] * f3 - b[3] * f5;
	g = b[3] * f3 + b[5] * f5;

	b1[0] = (b[0] + b[4]) * f4;
	b1[1] = (b[0] - b[4]) * f4;
	b1[2] = b[2] * f6 - b[6] * f2;
	b1[3] = b[6] * f6 + b[2] * f2;
	b[4] = e + f;
	b1[5] = e - f;
	b1[6] = h - g;
	b[7] = h + g;

	b[5] = (b1[6] - b1[5]) * f0;
	b[6] = (b1[6] + b1[5]) * f0;
	b[0] = b1[0] + b1[3];
	b[1] = b1[1] + b1[2];
	b[2] = b1[1] - b1[2];
	b[3] = b1[0] - b1[3];

	for (j = 0; j < 4; j++)
	{
	    j1 = 7 - j;
	    d[j][i] = b[j] + b[j1];
	    d[j1][i] = b[j] - b[j1];
	}
    }

    for (i = 0; i < 8; i++)
    {
	for (j = 0; j < 8; j++)
	{
	    block[i * 8 + j] = mnint(d[i][j]);
	}
    }

    m = 0;
    for (int n = 0; n < 8; n++)
    {
	for (int z = 0; z < 8; z++)
	{
	    outblock[n][z] = block[m];
	    m++;
	}
    }

    //	allFrames.push_back(block1);
}

//Inverse quantization

void Compute_Inverse_quantization(int inMatrix[8][8], int outMatrix[8][8], int Quant_parameter)
{ 
    for (int i = 0; i < 8; i++)
    {
		for (int j = 0; j < 8; j++)
		{
			outMatrix[i][j] = inMatrix[i][j] * Quant_parameter;
		}
    }
}
void Encode_Video_File()
{
    const char *inputFileName = "coastguard_qcif.yuv";
    BYTE *frameBuffer; // Pointer to current frame buffer
    //BYTE *outputEncodedStream;  //pointer to output stream buffer
    FILE *inputFileptr = NULL; // File pointer
    // Open the file in binary mode
    if ((inputFileptr = fopen(inputFileName, "rb")) == NULL)
    {
	cout << "Could not open specified file" << endl;
	exit(-1);
    }
    else
	cout << "File opened successfully" << endl;

    FILE *OutputFile;

    if ((OutputFile = fopen("output/encoded_coastguard_qcif.264", "wb")) == NULL)
    {
	cout << "Could not open output file" << endl;
	exit(-1);
    }
    else
	cout << "output File opened successfully" << endl;

    
	FILE *testFile;

    if ((testFile = fopen("output/test_reference_coastguard_qcif.yuv", "wb")) == NULL)
    {
	cout << "Could not open output file" << endl;
	exit(-1);
    }
    else
	cout << "output File opened successfully" << endl;


    //Set Y,U, and V frame sizes
    y_buffer_size_bytes = Y_frame_width * Y_frame_height;
    u_buffer_size_bytes = (Y_frame_width * Y_frame_height) / 4;
    v_buffer_size_bytes = (Y_frame_width * Y_frame_height) / 4;
    //Total YUV frame size
    framesize = y_buffer_size_bytes + u_buffer_size_bytes + v_buffer_size_bytes;

    // Get the size of the file in bytes
    long fileSize = getFileSize(inputFileptr);
    total_number_of_frames = fileSize / framesize;
    cout << "input file" << inputFileName << " size in bytes is " << fileSize << " & number of frames is " << total_number_of_frames << endl;

    // Allocate space in the buffer for a single frame
    frameBuffer = new BYTE[framesize];

    BYTE *referenceFrameBuffer = new BYTE[framesize];
    BYTE *uFrameStart, *vFrameStart;
    uFrameStart = frameBuffer + y_buffer_size_bytes;
    vFrameStart = uFrameStart + u_buffer_size_bytes;

    //Allocate space for output encoded stream
    //outputEncodedStream = new BYTE[total_number_of_frames*4000]; //asumption: assume average frame size is 4000 bytes

    int macroblock_Xpos, macroblock_Ypos; //hold top-left corner of current macroblock "to be encoded"
    Block current_blocks[6];
    vector<int> run_length_table;


	//array of frame widths
    int frameWidths[6] = {Y_frame_width, Y_frame_width, Y_frame_width, Y_frame_width, Y_frame_width/2, Y_frame_width/2};
    //yyyyuv each in the format of x,y
	int blockOffsetsXY[12] = {0,0,8,0,0,8,8,8,0,0,0,0};
	int mvMultipliers[6] = {2, 2, 2, 2, 1, 1};
	
    BYTE* frameOffsets[6] = {frameBuffer, frameBuffer, frameBuffer, frameBuffer, uFrameStart, vFrameStart};
    BYTE* refFrameOffsets[6] = {referenceFrameBuffer, referenceFrameBuffer, referenceFrameBuffer, referenceFrameBuffer, 
								referenceFrameBuffer + y_buffer_size_bytes, referenceFrameBuffer + y_buffer_size_bytes+ u_buffer_size_bytes};
    

    int isIframe;
    //start encoding loop
	int Quant_parameter = 16;
	int frameSizeInBytes;
	int expectedFrameSizeInBytes = (int)((SIZE_IN_MB*1024*1024)/(300)); 
    for (int frame_num = 0; frame_num < total_number_of_frames; frame_num++)
	{
		isIframe = (frame_num%numOfPFrames== 0) ? 1: 0;
		frameSizeInBytes = 0;
		writeQuantization(Quant_parameter, OutputFile);

		cout << "Encoding frame number: " << frame_num << " , quantization parameter : " << Quant_parameter << endl;
		//read a frame from input file into the frameBuffer
		fread(frameBuffer, framesize, 1, inputFileptr);

		// std::copy(std::begin(frameBuffer), std::end(frameBuffer), std::begin(lastIFrameBuffer));
		//loop accross blocks
		for (macroblock_Ypos = 0; macroblock_Ypos < Y_frame_height; macroblock_Ypos += 16)
			for (macroblock_Xpos = 0; macroblock_Xpos < Y_frame_width; macroblock_Xpos += 16)
			{
			//load individual blocks into current_blocks
			for (int block_y = 0; block_y < 8; block_y++)
				for (int block_x = 0; block_x < 8; block_x++)
				{
					//Y block 0
					current_blocks[0].data[block_x][block_y] = *(frameBuffer + (macroblock_Xpos + block_x) + (macroblock_Ypos + block_y) * Y_frame_width);

					//Y block 1
					current_blocks[1].data[block_x][block_y] = *(frameBuffer + (macroblock_Xpos + block_x + 8) + (macroblock_Ypos + block_y) * Y_frame_width);

					//Y block 2
					current_blocks[2].data[block_x][block_y] = *(frameBuffer + (macroblock_Xpos + block_x) + (macroblock_Ypos + block_y + 8) * Y_frame_width);

					//Y block 3
					current_blocks[3].data[block_x][block_y] = *(frameBuffer + (macroblock_Xpos + block_x + 8) + (macroblock_Ypos + block_y + 8) * Y_frame_width);

					//u block
					current_blocks[4].data[block_x][block_y] = *(uFrameStart + (macroblock_Xpos / 2 + block_x) + (macroblock_Ypos / 2 + block_y) * Y_frame_width / 2);

					//v block
					current_blocks[5].data[block_x][block_y] = *(vFrameStart + (macroblock_Xpos / 2 + block_x) + (macroblock_Ypos / 2 + block_y) * Y_frame_width / 2);
					
				}
			//macroblock processing
			//loop accross all blocks in the macroblock
			MV mv;
			mv.x = 0;
			mv.y = 0;
			int blockOffsets[6] = {
				0 + (macroblock_Xpos) + (macroblock_Ypos) * Y_frame_width,
				0 + (macroblock_Xpos + 8) + (macroblock_Ypos) * Y_frame_width,
				0 + (macroblock_Xpos) + (macroblock_Ypos + 8) * Y_frame_width,
				0 + (macroblock_Xpos + 8) + (macroblock_Ypos + 8) * Y_frame_width,
				y_buffer_size_bytes + (macroblock_Xpos / 2) + (macroblock_Ypos / 2) * Y_frame_width / 2,
				y_buffer_size_bytes + u_buffer_size_bytes + (macroblock_Xpos / 2) + (macroblock_Ypos / 2) * Y_frame_width / 2,
			};

			if (!isIframe)//TODO calculate mvs
			{
				mv = Compute_MV(get16x16Block(frameBuffer,Y_frame_width,macroblock_Xpos, macroblock_Ypos), referenceFrameBuffer);
				writeMotionVector(mv,OutputFile);
			}

			for (int block_index = 0; block_index < 6; block_index++)
			{
				int outBlock[8][8];
				BYTE diffBlock[8][8];
				//to store the block after it's  decoding
				int codec[8][8];

				Block refFrameBlock;
				if (!isIframe) //to mv or not to mv
				{
					refFrameBlock = get8x8Block(refFrameOffsets[block_index], frameWidths[block_index] ,  ((mv.x*mvMultipliers[block_index])*8+ blockOffsetsXY[block_index*2+0]),  ((mv.y*mvMultipliers[block_index])*8+  + blockOffsetsXY[block_index*2+1]));
					Compute_Diffrences(current_blocks[block_index].data, refFrameBlock.data, current_blocks[block_index].data);				
				}

				//encode then decode then store this current block to be used next frame
				//DCT
				Compute_DCT(current_blocks[block_index].data, outBlock);
				//Quantization
				Compute_quantization(outBlock, outBlock,Quant_parameter);
				//Zigzag and run length
				run_length_table = Compute_VLC(outBlock);
				// cout << "frame : " << frame_num << " macrobBlock [ "<< macroblock_Xpos<< ","<< macroblock_Ypos<< "] , blockType : "<< block_index << " coeffecients count: " << run_length_table.size() << endl;
				//write coefficient into output file
				writeEncodedMacroblock(run_length_table, OutputFile);
				frameSizeInBytes += run_length_table.size()*2;
				// inverse Quantization
				Compute_Inverse_quantization(outBlock, codec,Quant_parameter);
				//IDCT 
				Compute_idct(codec, codec);
				
				if (!isIframe) //to mv or not to mv
					Compute_Additions(codec, refFrameBlock.data, codec);
				
				set8x8Block(codec, frameBuffer, frameWidths[block_index], blockOffsets[block_index]);   				

				// encode the diffrence if it's an Iframe				
				

				
			}

			} //end macroblock loop
			fwrite(referenceFrameBuffer, framesize, 1, testFile);
			if(!isIframe)
			{
				if(abs(frameSizeInBytes - expectedFrameSizeInBytes) > 64)
				Quant_parameter = 1+(Quant_parameter*frameSizeInBytes/expectedFrameSizeInBytes);
			}
			copy(frameBuffer, frameBuffer + framesize, referenceFrameBuffer);
	
    } //end frame loop

    //free allocated memory
    delete frameBuffer;
    fclose(OutputFile);
	fclose(testFile);
    //	free(outputEncodedStream);

    return;
}



///////////////////////////////////////////////////////// DECODE //////////////////////////

//Inverse ZigZag scan
//read number of coefficients from rle (number = length), store the coefficients in output_block
void Compute_inverseZigzag(vector<int> rle, int num_coefficients, int output_block[8][8])
{
    memset(output_block, 0, sizeof(int) * 8 * 8); //initialzie the output_block to zeros
    int zigzag_scaned_values[64];

    //loop for all number of coefficients
    int run_num; //number of zeros before a non-zero value
    int level;   // non-zero value
    //int zigzag_scaned_values_index;
    int zigzag_ordered_index;

    //initalization
    zigzag_ordered_index = 0;
    memset(zigzag_scaned_values, 0, sizeof(int) * 64);

    for (int coefficient_index = 0; coefficient_index < num_coefficients; coefficient_index += 2) //read two at a time
    {
	run_num = rle[coefficient_index];
	level = rle[coefficient_index + 1];

	//fill the zigzag_scaned_values array
	for (int zero_coefficient_count = 0; zero_coefficient_count < run_num; zero_coefficient_count++)
	{
	    zigzag_scaned_values[zigzag_ordered_index] = 0;
	    zigzag_ordered_index++;
	}
	//fill level
	zigzag_scaned_values[zigzag_ordered_index] = level;
	zigzag_ordered_index++;

	//countinue with next coefficient pair
    }

    //now  zigzag_scaned_values has the coefficients in zigzag order, need to convert them to normal order in 8x8 block

    for (int col_index = 0; col_index < 8; col_index++)
	for (int row_index = 0; row_index < 8; row_index++)
	{
	    zigzag_ordered_index = zigzag_order[col_index][row_index];
	    output_block[col_index][row_index] = zigzag_scaned_values[zigzag_ordered_index];
	}
}


//Assume encoded file is QCIF resolution (176 x 144)
void Decode_Video_File()
{
    FILE *file;
    if ((file = fopen("output/encoded_coastguard_qcif.264", "rb")) == NULL)
    {
		cout << "Could not open specified file" << endl;
		return;
    }
    else
    {
		cout << "File opened successfully" << endl;
    }

    //open output file
    FILE *output_file;
    if ((output_file = fopen("output/decoded_coastguard_qcif.yuv", "wb")) == NULL)
    {
		cout << "Could not open specified file" << endl;
		return;
    }
    else
    {
		cout << "File opened successfully" << endl;
    }

    long fileSize = getFileSize(file);
    int coefficient_size_in_bytes = sizeof(MyDataSize);
    int total_num_encoded_coefficients = fileSize / coefficient_size_in_bytes;
    cout << "Number of encoded coefficients: " << total_num_encoded_coefficients << endl;
    MyDataSize *fileBuffer = new MyDataSize[total_num_encoded_coefficients];
    fread(fileBuffer, fileSize, 1, file);
    fclose(file);

    //Set Y,U, and V frame sizes
    y_buffer_size_bytes = Y_frame_width * Y_frame_height;
    u_buffer_size_bytes = (Y_frame_width * Y_frame_height) / 4;
    v_buffer_size_bytes = (Y_frame_width * Y_frame_height) / 4;
    //Total YUV frame size
    framesize = y_buffer_size_bytes + u_buffer_size_bytes + v_buffer_size_bytes;

    BYTE *yuv_frameBuffer; // Pointer to current frame buffer
    yuv_frameBuffer = new BYTE[framesize];
    BYTE *uFrameStart, *vFrameStart;
    uFrameStart = yuv_frameBuffer + y_buffer_size_bytes;
    vFrameStart = uFrameStart + u_buffer_size_bytes;
	//frame to store last frame for mv compansations
    BYTE *referenceFrameBuffer = new BYTE[framesize];


    int current_blocks[6][8][8];
    int macroblock_Xpos, macroblock_Ypos;					 //hold top-left corner of current macroblock "to be encoded"
    int number_macroblocks_per_frame = Y_frame_width / 16 * Y_frame_height / 16; //assume YUV 420

    //int i;
    int block_num;
    int macroblock_number;
    int number_block_coefficients;
    int coefficient_index;
    vector<int> inputBlock;
    coefficient_index = 0;
    int frame_num = 0;
	int isIframe;
	int mvX, mvY;
	BYTE* frameOffsets[6] = {yuv_frameBuffer, yuv_frameBuffer, yuv_frameBuffer, yuv_frameBuffer, uFrameStart, vFrameStart};
	BYTE* refFrameOffsets[6] = {referenceFrameBuffer, referenceFrameBuffer, referenceFrameBuffer, referenceFrameBuffer, 
								referenceFrameBuffer + y_buffer_size_bytes, referenceFrameBuffer + y_buffer_size_bytes+ u_buffer_size_bytes};
    
	int frameWidths[6] = {Y_frame_width, Y_frame_width, Y_frame_width, Y_frame_width, Y_frame_width/2, Y_frame_width/2};
	//yyyyuv each in the format of x,y
	int blockOffsetsXY[12] = {0,0,8,0,0,8,8,8,0,0,0,0};
	int mvMultipliers[6] = {2, 2, 2, 2, 1, 1};
	

    int num_macroblock_per_row = Y_frame_width / 16;

    int Quant_parameter = 0;
	while (coefficient_index < total_num_encoded_coefficients) //loop accross the whole input file
    {
		isIframe = (frame_num%numOfPFrames== 0) ? 1: 0;
		Quant_parameter = (int)fileBuffer[coefficient_index++];
	cout << "decoding frame number " << frame_num << endl;
	//frame loop: loop accross macroblocks in the whole QCIF frame
	for (macroblock_number = 0; macroblock_number < number_macroblocks_per_frame; macroblock_number++)
	{
	    macroblock_Xpos = (macroblock_number % num_macroblock_per_row) * 16;
	    macroblock_Ypos = (macroblock_number / num_macroblock_per_row) * 16;
	    cout << "macroblock#: " << macroblock_number << " xPos,yPos: " << macroblock_Xpos << "," << macroblock_Ypos << endl;
	    
		//get mvX, mvY per loop if it's not an I frame. 
		if(!isIframe)
		{
			mvX = fileBuffer[coefficient_index++];
			mvY = fileBuffer[coefficient_index++];
		}
		
		//macroblock loop: loop accross all blocks for a single macroblock
	    for (int block_index = 0; block_index < 6; block_index++)
	    {
		//RLE sequence for each block
		inputBlock.clear();
		number_block_coefficients = fileBuffer[coefficient_index];
		coefficient_index++; //advance index to point to next coefficient in the fileBuffer
		// cout << "reading block # " << block_index << " number of coefficients is: " << number_block_coefficients << endl;
		for (int j = 0; j < number_block_coefficients; j++)
		{
			inputBlock.push_back(fileBuffer[coefficient_index + j]);
		}
		Compute_inverseZigzag(inputBlock, number_block_coefficients, current_blocks[block_index]);
		Compute_Inverse_quantization(current_blocks[block_index], current_blocks[block_index],Quant_parameter);
		Compute_idct(current_blocks[block_index], current_blocks[block_index]);

		// add residuls to the mv block 
		Block bestMatchBlock; 
		if(!isIframe)
		{
			// Block bestMatchBlock = get8x8Block(frameOffsets[block_index], frameWidths[block_index] ,  ((mvX*mvMultipliers[block_index])*16+macroblock_Xpos+ blockOffsetsXY[block_index*2+0]),  ((mvY*mvMultipliers[block_index])*16+macroblock_Ypos+  + blockOffsetsXY[block_index*2+1]));
			Block bestMatchBlock = get8x8Block(refFrameOffsets[block_index], frameWidths[block_index] ,  ((mvX*mvMultipliers[block_index])*8+ blockOffsetsXY[block_index*2+0]),  ((mvY*mvMultipliers[block_index])*8+  + blockOffsetsXY[block_index*2+1]));
			Compute_Additions(current_blocks[block_index], bestMatchBlock.data, current_blocks[block_index]);	
		}

		char blockOffsets[6] = {
			0 + (macroblock_Xpos) + (macroblock_Ypos) * Y_frame_width,
			0 + (macroblock_Xpos + 8) + (macroblock_Ypos) * Y_frame_width,
			0 + (macroblock_Xpos) + (macroblock_Ypos + 8) * Y_frame_width,
			0 + (macroblock_Xpos + 8) + (macroblock_Ypos + 8) * Y_frame_width,
			y_buffer_size_bytes + (macroblock_Xpos / 2) + (macroblock_Ypos / 2) * Y_frame_width / 2,
			y_buffer_size_bytes + u_buffer_size_bytes + (macroblock_Xpos / 2) + (macroblock_Ypos / 2) * Y_frame_width / 2,
		};
		
		coefficient_index += number_block_coefficients;
		} //macroblock loop

		//copy macroblock into frame buffer
		//load individual blocks into current_blocks
		//might need clipping here
		for (int block_y = 0; block_y < 8; block_y++)
		for (int block_x = 0; block_x < 8; block_x++)
		{
			//Y block 0
			*(yuv_frameBuffer + (macroblock_Xpos + block_x) + (macroblock_Ypos + block_y) * Y_frame_width) = current_blocks[0][block_x][block_y];
			//Y block 1
			*(yuv_frameBuffer + (macroblock_Xpos + block_x + 8) + (macroblock_Ypos + block_y) * Y_frame_width) = current_blocks[1][block_x][block_y];
			//Y block 2
			*(yuv_frameBuffer + (macroblock_Xpos + block_x) + (macroblock_Ypos + block_y + 8) * Y_frame_width) = current_blocks[2][block_x][block_y];
			//Y block 3
			*(yuv_frameBuffer + (macroblock_Xpos + block_x + 8) + (macroblock_Ypos + block_y + 8) * Y_frame_width) = current_blocks[3][block_x][block_y];
			//u block
			*(uFrameStart + (macroblock_Xpos / 2 + block_x) + (macroblock_Ypos / 2 + block_y) * Y_frame_width / 2) = current_blocks[4][block_x][block_y];
			//v block
			*(vFrameStart + (macroblock_Xpos / 2 + block_x) + (macroblock_Ypos / 2 + block_y) * Y_frame_width / 2) = current_blocks[5][block_x][block_y];
		}

	} //frame loop completed

	isIframe = 0;
	//store frame into output file
	fwrite(yuv_frameBuffer, framesize, 1, output_file);
	//store last decoded frame to be used for motion compansation in the next frame
	copy(yuv_frameBuffer, yuv_frameBuffer + framesize, referenceFrameBuffer);
	frame_num++;
    }

    cout << "Decoding completed: total number of decoded frames is: " << frame_num << endl;
    fclose(output_file);
    delete[] fileBuffer;
    delete yuv_frameBuffer;
    //		cout << "Decoding done";
}

//Main function
//Encode inout video file coastguard_qcif into encoded.264
//decode encoded.264 into decoded,yuv
int main(int argc, char *argv[])
{
    Encode_Video_File();
    cout << "File encoding completed, press a key to continue";
    

    Decode_Video_File();
    cout << "File decoding completed, press a key to exit";
    cin.get();
}
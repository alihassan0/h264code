# AVP Report

## How to run 
```css
/*on unix machines*/
$ clone the repo
$ cd h264code/
$ clang++ encoder_main.cpp -o output/codec
$ ./output/codec

/*to run the output video using vlc*/
$ vlc  --rawvid-fps 25 --rawvid-width 176 --rawvid-height 144 --rawvid-chroma I420 output/decoded_coastguard_qcif.yuv 
```
---

## bitstream
- all frames are stored sequentially with the first one as IFrame and the rest as Pframes 
    - I frames
        - each frame is preceeded with Quant_parameter
        - macroblocks are stored sequentially consisting of 6 blocks
            - each block is turned into codingn coeffecients and is encoded as [size, coeffecients]
    - Non Iframes same as Iframes except
        - all macroblocks are preceeded with two DataSize of motion vector

## Implementaion

### Structs

- `Block`: 8*8 byte block used for sorting blocks 
```cpp
typedef struct Block
{
    BYTE data[8][8];
} Block;
```

- `Block16x16`: 16*16 byte block used for sorting macroblocks that are used in motion estimation 
```cpp 
typedef struct Block16x16
{
    BYTE data[16][16];
	int x;
	int y;
} Block16x16;
```

- `MV`: MV Struct used to store motion vector normalized 
```cpp 
typedef struct MV
{
    BYTE x;
    BYTE y;
} MV;
```

## **Writing data**

```cpp
/**
    - purpose:      writes motion vector in the bit stream
    - mv            motion vectore used for the cucrrent macroblocks
    - output_file   pointer to the ouputfile

**/
void writeMotionVector(MV mv,FILE *output_file);

/**
    - purpose:          writes Quant_parameter in the bit stream
    - Quant_parameter   motion vectore used for the cucrrent macroblocks
    - output_file       pointer to the ouputfile

**/

void writeQuantization(MyDataSize Quant_parameter, FILE *output_file)

/**
    - purpose:      writes run length coding coeffecients in the bit stream
    - rle           run length coding coeffecients
    - output_file   pointer to the ouputfile

**/
void writeEncodedMacroblock(vector<int> rle, FILE *output_file)
```

---
## **utils**
```cpp
/**
    - purpose:      gets 8*8 block used for motion estimation and compansation
    - frame         pointer to the start of Y or U or V frame
    - frameWidth    width of Y or U or V frame
    - startI        x index of the wanted block **must be multiple of 8**
    - startJ        y index of the wanted block **must be multiple of 8**
**/
Block get8x8Block(BYTE *frame, int frameWidth, int startI, int startJ);

/**
    - purpose:      gets 16*16 block used for motion estimation
    - frame         pointer to the start of Y or U or V frame
    - frameWidth    width of Y or U or V frame
    - startI        x index of the wanted block **must be multiple of 16**
    - startJ        y index of the wanted block **must be multiple of 16**
**/
Block16x16 get16x16Block(BYTE *frame, int frameWidth, int startI, int startJ)

/**
    - purpose:      sets 8*8 block used for storing encoded->decoded blocks to be used for reference
    - block         the 8*8 block to be set      
    - frame         pointer to the start of Y or U or V reference frame
    - frameWidth    width of Y or U or V frame
    - offset        startJ * width + startI 
**/
void set8x8Block(int block[8][8], BYTE *frame, int frameWidth, int offset)

```

---
## **motion estimation**

```cpp
/**
    - purpose:      calculated diamond search in the given frame recursively
    - block         the 16*16 that we are trying to find the best match
    - frame         pointer to the start of Y or U or V reference frame
    - isLDSP        whether to use small diamond search or large diamond search
    - frameWidth    width of Y or U or V frame
    - lPoints       points used for checking in LDSP 
    - sPoints       points used for checking in SDSP 
    - output        MV - > **motion vectors are positions and not directions** 
**/

MV LDSP(Block16x16 block, BYTE *frame, int isLDSP, int lPoints[9][2], int sPoints[5][2])

/**
    - purpose:      calculated motion vectos for a given macroblock
    - block         the 16*16 that we are trying to find the best match
    - frame         pointer to the start of Y or U or V reference frame
    - output        MV - > **motion vectors are positions and not directions** 
**/

MV Compute_MV(Block16x16 block, BYTE *frame)
```

---
## **Motion Compansation**

```cpp
/**
    - purpose:      calculates diffrences between two 8*8 blocks and stores it in a third block
    - matrix1       used for the current block
    - matrix2       used for the reference block
    - outMatrix     stores matrix1- matrix2 
**/
void Compute_Diffrences(BYTE matrix1[8][8], BYTE matrix2[8][8], BYTE outMatrix[8][8])


/**
    - purpose:      calculates additions between two 8*8 blocks and stores it in a third block
    - matrix1       used for the current block
    - matrix2       used for the bestMatch block
    - outMatrix     stores matrix1+ matrix2 
**/
void Compute_Additions(int matrix1[8][8], BYTE matrix2[8][8], int outMatrix[8][8])
```


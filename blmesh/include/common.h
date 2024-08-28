#ifndef _COMMOM_H
#define _COMMOM_H
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "BLVector.h"
#define OC_SUBNODES 8

typedef enum
{
	UNE = 0,
		UNW = 1,
		USW = 2,
		USE = 3,
		LNE = 4,
		LNW = 5,
		LSW = 6,
		LSE = 7
}OcCubeEnum;
 
typedef BLVector CUBECOORD;

struct OCCUBE
{
	OCCUBE() {}
	OCCUBE(CUBECOORD min, CUBECOORD max):upper(max),lower(min) {}
	CUBECOORD upper;
	CUBECOORD lower;
};


//for test (need to be deleted later)
extern int printflag;
////
#endif

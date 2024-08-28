#include "LiCheck.h"
#include "ComputInfo.h"
#include "md5.h"
#include <spdlog/spdlog.h> 
 #include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<time.h>


LiCheck::LiCheck(void)
{
}


LiCheck::~LiCheck(void)
{
}

void LiCheck::TransInteger(int num, char numcode[32])
{
	//char numcode[32];
	int r;
	memset(numcode, 0, sizeof(numcode));

	srand((unsigned)time(NULL));
	r = rand()%2;

	switch(num)
	{
	case 0:
		if(r == 0)
			strcpy(numcode, "JB3F");
		else
			strcpy(numcode, "BD3J");
		break;

	case 1:
		if(r == 0)
			strcpy(numcode, "T9SV");
		else
			strcpy(numcode, "V9SJ");
		break;

	case 2:
		if(r == 0)
			strcpy(numcode, "E1IR");
		else
			strcpy(numcode, "1ETR");
		break;

	case 3:
		if(r == 0)
			strcpy(numcode, "8ARZ");
		else
			strcpy(numcode, "Z8RA");
		break;

	case 4:
		if(r == 0)
			strcpy(numcode, "HQ4A");
		else
			strcpy(numcode, "ADQB");
		break;

	case 5:
		if(r == 0)
			strcpy(numcode, "W0ZD");
		else
			strcpy(numcode, "ODZP");
		break;

	case 6:
		if(r == 0)
			strcpy(numcode, "9EON");
		else
			strcpy(numcode, "NE0E");
		break;

	case 7:
		if(r == 0)
			strcpy(numcode, "YB2L");
		else
			strcpy(numcode, "BL2T");
		break;

	case 8:
		if(r == 0)
			strcpy(numcode, "T3JE");
		else
			strcpy(numcode, "JE9T");
		break;

	case 9:
		if(r == 0)
			strcpy(numcode, "3X8I");
		else
			strcpy(numcode, "X8I3");
		break;

	default:
		break;
	}
}

void LiCheck::TransString(char numcode[32], int *num)
{
	if(!strcmp(numcode, "JB3F") || !strcmp(numcode, "BD3J"))
		*num = 0;
	else if(!strcmp(numcode, "T9SV") || !strcmp(numcode, "V9SJ"))
		*num = 1;
	else if(!strcmp(numcode, "E1IR") || !strcmp(numcode, "1ETR"))
		*num = 2;
	else if(!strcmp(numcode, "8ARZ") || !strcmp(numcode, "Z8RA"))
		*num = 3;
	else if(!strcmp(numcode, "HQ4A") || !strcmp(numcode, "ADQB"))
		*num = 4;
	else if(!strcmp(numcode, "W0ZD") || !strcmp(numcode, "ODZP"))
		*num = 5;
	else if(!strcmp(numcode, "9EON") || !strcmp(numcode, "NE0E"))
		*num = 6;
	else if(!strcmp(numcode, "YB2L") || !strcmp(numcode, "BL2T"))
		*num = 7;
	else if(!strcmp(numcode, "T3JE") || !strcmp(numcode, "JE9T"))
		*num = 8;
	else if(!strcmp(numcode, "3X8I") || !strcmp(numcode, "X8I3"))
		*num = 9;
}

void LiCheck::DecodeDate(char date[64])
{
	int i, j, num;
	char tmp[64], str[8], strnum[8];

	memset(tmp, 0, sizeof(tmp));
	strcpy(tmp, date);
	memset(date, 0, sizeof(date));
	for (i=0; i<strlen(tmp);)
	{
		memset(str, 0, sizeof(str));
		memset(strnum, 0, sizeof(strnum));
		for (j=0; j<4; j++)
			str[j] = tmp[i++];
		TransString(str, &num);
		sprintf(strnum, "%d", num);
		strcat(date, strnum);
	}
}

void LiCheck::SplitDate(char date[64], int &year, int &month, int &day)
{
	char strnum[16];
	int i, num;

	memset(strnum, 0, sizeof(strnum));
	for(i=0; i<4; i++)
		strnum[i] = date[i];
	year = atoi(strnum);

	memset(strnum, 0, sizeof(strnum));
	for(i=4; i<6; i++)
		strnum[i-4] = date[i];
	month = atoi(strnum);

	memset(strnum, 0, sizeof(strnum));
	for(i=6; i<8; i++)
		strnum[i-6] = date[i];
	day = atoi(strnum);
}

int LiCheck::CheckLic(char *filename)
{
	int i, j, r1, r2;
	char title[256],strtp1[256], strtp2[256], strnum[8], strt[512];
	char username[64], md51[64], md52[64], startdate[64], expiredate[64];
	int expyear, expmonth, expday, curyear, curmonth, curday;
	int syear, smonth, sday;

	FILE *fin = NULL;
	fin = fopen(filename, "r");
	if(!fin)
	{
		//spdlog::info("Error: can not open file %s\n", filename);
		return N0_LICENSE_FILE;
	}
	memset(title, 0, sizeof(title));
	fscanf(fin, "%s\n", title);
	if (!strcmp(title, "ANONYMOUS"))
	{
		memset(strtp1, 0, sizeof(strtp1));
		memset(strtp2, 0, sizeof(strtp2));
		memset(md52, 0, sizeof(md52));
		fscanf(fin, "%s\n", strtp1);
		for (i=0; i<130; i++)
		{
			if(i<96)
				strtp2[i] = strtp1[i];
			else if(i == 96)
			{
				memset(strnum, 0, sizeof(strnum));
				strnum[0] = strtp1[i];
				//sprintf(strnum, "%s", strtp1[i]);
				r1 = atoi(strnum);
			}
			else if (i == 97)
			{
				memset(strnum, 0, sizeof(strnum));
				strnum[0] = strtp1[i];
				//sprintf(strnum, "%s", strtp1[i]);
				r2 = atoi(strnum);
			}
			else
			{
				md52[i-98] = strtp1[i];
			}
		}

		memset(strt, 0, sizeof(strt));
		sprintf(strt, "%s%s%d%d", title, strtp2, r1, r2);
		MD5 md5c;
		md5c.update(strt);
		memset(md51, 0, sizeof(md51));
		strcpy(md51, md5c.toString().c_str());
		strupr(md51);
		if (strcmp(md51, md52))
		{
			return INCORRECT_LICENSE_FILE;
		}

		memset(md51, 0, sizeof(md51));
		memset(startdate, 0, sizeof(startdate));
		memset(expiredate, 0, sizeof(expiredate));

		int idx = 0;
		if (r1<r2)
		{
			for (i=0; i<strlen(strtp2);)
			{
				if (i == r1)
				{
					for (j=0; j<32; j++)
						startdate[j] = strtp1[i++];
				}
				else if (i == r2+32)
				{
					for (j=0; j<32; j++)
						expiredate[j] = strtp1[i++];
				}
				else
					md51[idx++] = strtp2[i++];
			}
		} 
		else
		{
			for (i=0; i<strlen(strtp2);)
			{
				if (i == r1+32)
				{
					for (j=0; j<32; j++)
						startdate[j] = strtp1[i++];
				}
				else if (i == r2)
				{
					for (j=0; j<32; j++)
						expiredate[j] = strtp1[i++];
				}
				else
					md51[idx++] = strtp2[i++];
			}
		}

		//get hostname and hostmac, check it out with md51
		ComputInfo cpInfo;
		char struser[256], umd5[256];
		cpInfo.UserInfo(struser);
		MD5 md5user;
		md5user.update(struser);
		memset(umd5, 0, sizeof(umd5));
		strcpy(umd5, md5user.toString().c_str());
		strupr(umd5);
		if (strcmp(md51, umd5))
		{
			return INCORRECT_LICENSE_FILE;
		}

		//decodedate
		DecodeDate(startdate);
		DecodeDate(expiredate);
		SplitDate(startdate, syear, smonth, sday);
		SplitDate(expiredate, expyear, expmonth, expday);
		struct tm *calender;
		int months_left = 0;
		time_t now;

		expyear -= 1900;
		syear -= 1900;

		time(&now);
		calender = localtime(&now);

		if(calender->tm_year < syear)
			return INCORRECT_SYSTEM_TIME;
		else if (calender->tm_year == syear && calender->tm_mon+1 < smonth)
			return INCORRECT_SYSTEM_TIME;
		else if (calender->tm_year == syear && calender->tm_mon+1 == smonth && calender->tm_mday < sday)
			return INCORRECT_SYSTEM_TIME;

		if(calender->tm_year < expyear - 1)
			return CORRECT_LICENSE_FILE;
		else if(calender->tm_year == expyear - 1)
			months_left = expmonth + (11 - calender->tm_mon);
		else if(calender->tm_year == expyear)
			months_left = (expmonth - 1 - calender->tm_mon);
		else
			return EXPIRED_LICENSE;

		if(months_left > 0)
			spdlog::info("License expires in {} months .....\n", months_left);
		else if(months_left == 0)
		{
			if(expday >= calender->tm_mday)
				spdlog::info("License expires in {} days .....\n", expday-calender->tm_mday);
			else
			{
				return EXPIRED_LICENSE;
			}
		}
		else
			return EXPIRED_LICENSE;
	} 
	else
	{
	}

	fclose(fin);
	fin = NULL;
	return CORRECT_LICENSE_FILE;
}
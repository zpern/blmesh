#pragma once

#define N0_LICENSE_FILE -1
#define INCORRECT_LICENSE_FILE 0
#define CORRECT_LICENSE_FILE 1
#define EXPIRED_LICENSE 2
#define INCORRECT_SYSTEM_TIME 3

class LiCheck
{
public:
	LiCheck(void);
	~LiCheck(void);

	int CheckLic(char *filename);

private:
	void EncodeDate(char date[64]);
	void DecodeDate(char date[64]);
	void SplitDate(char date[64], int &year, int &month, int &day);
	void TransInteger(int num, char numcode[32]);
	void TransString(char numcode[32], int *num);
};


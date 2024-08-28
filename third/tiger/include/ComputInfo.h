#pragma once
#ifndef _WIN32_WINNT  
#define _WIN32_WINNT 0x0501  
#endif  

#include <windows.h>  
#include <winioctl.h>  

class ComputInfo
{
public:
	ComputInfo(void);
	~ComputInfo(void);

	void GetHostName();
	void GetMacAddress();
	void GetDiskSerialNo();
	void UserInfo(char userinfo[256]);
	void OutputLog(char *filename);

	//helper
	BOOL GetPhyDriveSerial(LPTSTR pModelNo, LPTSTR pSerialNo);  
	void ToLittleEndian(PUSHORT pWords, int nFirstIndex, int nLastIndex, LPTSTR pBuf);  
	void TrimStart(LPTSTR pBuf);
	void TcharToChar (const TCHAR * tchar, char * _char);

private:
	char m_sMacAddress[256];
	char m_sMacAddressWireless[256];
	char m_sHostName[256];
	char m_sDiskSerialNo[256];
	bool bWireless;
	bool bNowireless;
};


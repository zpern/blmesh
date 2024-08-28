#include "ComputInfo.h"
#include <spdlog/spdlog.h> 
 #include <iostream>  
#include <string>  
#include <vector>  
#include <atlstr.h>   
#include<atlbase.h>
#include<atlconv.h>
#include<iphlpapi.h>
#pragma comment(lib,"Iphlpapi.lib")
#pragma comment(lib,"ws2_32.lib")
#include "HardDriveSerialNumer.h"

ComputInfo::ComputInfo(void)
{
	memset(m_sHostName, 0, sizeof(m_sHostName));
	memset(m_sDiskSerialNo, 0, sizeof(m_sDiskSerialNo));
	memset(m_sMacAddress, 0, sizeof(m_sMacAddress));
	memset(m_sMacAddressWireless, 0, sizeof(m_sMacAddressWireless));
	strcpy(m_sMacAddress, "10-BF-48-B8-97-34");
	strcpy(m_sMacAddressWireless, "10-BF-48-B8-97-34");
	bWireless = false;
	bNowireless = false;
}


ComputInfo::~ComputInfo(void)
{
}

void ComputInfo::GetHostName()
{
	CString m_IP; 
	char PCnameBuffer[128];  

	WSAData data;  
	if(WSAStartup(MAKEWORD(1,1),&data)!=0)  
	{  
		spdlog::info("Error: Failed to initialize network!\n");
		return;
	}  
	else
	{
		if(0==gethostname(PCnameBuffer,128))  
		{
			strcpy(m_sHostName, PCnameBuffer);
			struct hostent* pHost;  
			int i;  
			pHost=gethostbyname(PCnameBuffer); 
			for (i=0;pHost!=NULL&&pHost->h_addr_list[i]!=NULL;i++)  
			{  
				LPCSTR psz = inet_ntoa(*(struct in_addr *)pHost->h_addr_list[i]);
				m_IP += psz;
			}
		}
		else
		{  
			spdlog::info("Error: failed to get host information!\n");
			return;
		}
	}
}

/*
void ComputInfo::GetMacAddress()
{
	char address[256], tmpstr[8];
	PIP_ADAPTER_INFO pAdapterInfo;
	PIP_ADAPTER_INFO pAdapter=NULL;
	DWORD dwRetVal=0;

	memset(address, 0, sizeof(address));
	memset(tmpstr, 0, sizeof(tmpstr));

	pAdapterInfo=(IP_ADAPTER_INFO*)malloc(sizeof(IP_ADAPTER_INFO));
	ULONG ulOutBufLen=sizeof(IP_ADAPTER_INFO);
	if(GetAdaptersInfo(pAdapterInfo,&ulOutBufLen)!=ERROR_SUCCESS)
	{
		GlobalFree(pAdapterInfo);
		pAdapterInfo=(IP_ADAPTER_INFO*)malloc(ulOutBufLen);
	}
	if((dwRetVal=GetAdaptersInfo(pAdapterInfo,&ulOutBufLen))==NO_ERROR)
	{
		pAdapter=pAdapterInfo;
		while(pAdapter)
		{
			if( strstr(pAdapter->Description,"PCI")>0	//pAdapter->Description中包含"PCI"为：物理网卡
				||pAdapter->Type==71					//pAdapter->Type是71为：无线网卡
				)
			{
// 				printf("AdapterName:\t%s\n",pAdapter->AdapterName);
// 				printf("AdapterDesc:\t%s\n",pAdapter->Description);
// 				spdlog::info("AdapterAddr:\t");
				for(UINT i=0;i<pAdapter->AddressLength;i++)
				{
					if (i==pAdapter->AddressLength-1)
					{
						//spdlog::info("%02x\n", pAdapter->Address[i]);
						sprintf(tmpstr, "%02X", pAdapter->Address[i]);
					}
					else
					{
						//spdlog::info("%02x-", pAdapter->Address[i]);
						sprintf(tmpstr, "%02X-", pAdapter->Address[i]);
					}
					strcat(address, tmpstr);
				}
// 				printf("AdapterType:\t%d\n",pAdapter->Type);
// 				printf("IPAddress:\t%s\n",pAdapter->IpAddressList.IpAddress.String);
// 				printf("IPMask:\t%s\n",pAdapter->IpAddressList.IpMask.String);
			}
			pAdapter=pAdapter->Next;
		}
		strcpy(m_sMacAddress, address);
	}
	else
	{
		spdlog::info("Error: failed to get MAC address!\n");
	}
}*/

void ComputInfo::GetMacAddress()
{
	ULONG BufferLength = 0;
	BYTE* pBuffer = 0;
	char address[256], tmpstr[8];

	memset(address, 0, sizeof(address));
	memset(tmpstr, 0, sizeof(tmpstr));

	if( ERROR_BUFFER_OVERFLOW == GetAdaptersInfo( 0, &BufferLength ))
	{
		// Now the BufferLength contain the required buffer length.
		// Allocate necessary buffer.
		pBuffer = new BYTE[ BufferLength ];
	}
	else
	{
		// Error occurred. handle it accordingly.
	}

	// Get the Adapter Information.
	PIP_ADAPTER_INFO pAdapterInfo =
		reinterpret_cast<PIP_ADAPTER_INFO>(pBuffer);
	
	if((GetAdaptersInfo( pAdapterInfo, &BufferLength ))==NO_ERROR)
	{
		while(pAdapterInfo)
		{
			memset(address, 0, sizeof(address));
			if( strstr(pAdapterInfo->Description,"PCI")>0	//pAdapter->Description中包含"PCI"为：物理网卡
				||pAdapterInfo->Type==71					//pAdapter->Type是71为：无线网卡
				)
			{
// 				printf("AdapterName:\t%s\n",pAdapter->AdapterName);
// 				printf("AdapterDesc:\t%s\n",pAdapter->Description);
// 				spdlog::info("AdapterAddr:\t");
				for(UINT i=0;i<pAdapterInfo->AddressLength;i++)
				{
					memset(tmpstr, 0, sizeof(tmpstr));
					if (i==pAdapterInfo->AddressLength-1)
					{
						//spdlog::info("%02x\n", pAdapter->Address[i]);
						sprintf(tmpstr, "%02X", pAdapterInfo->Address[i]);
					}
					else
					{
						//spdlog::info("%02x-", pAdapter->Address[i]);
						sprintf(tmpstr, "%02X-", pAdapterInfo->Address[i]);
					}
					strcat(address, tmpstr);
				}
// 				printf("AdapterType:\t%d\n",pAdapter->Type);
// 				printf("IPAddress:\t%s\n",pAdapter->IpAddressList.IpAddress.String);
// 				printf("IPMask:\t%s\n",pAdapter->IpAddressList.IpMask.String);
			}
			if(pAdapterInfo->Type==71)
			{
				bWireless = true;
				strcpy(m_sMacAddressWireless, address);
			}
			else
			{
				bNowireless = true;
				strcpy(m_sMacAddress, address);
			}

			pAdapterInfo=pAdapterInfo->Next;
		}

	}
	else
	{
		spdlog::info("Error: failed to get MAC address!\n");
	}
}
/*
void ComputInfo::GetDiskSerialNo()
{
	TCHAR szModelNo[48], szSerialNo[24]; 
	char cModelNo[48], cSerialNo[24];
	if(GetPhyDriveSerial(szModelNo, szSerialNo))  
	{  
		TcharToChar(szModelNo, cModelNo);
		TrimStart(szSerialNo);  
		TcharToChar(szSerialNo, cSerialNo);
		//spdlog::info("Model No: %s\n", cModelNo);
		//spdlog::info("Serial No: %s\n", cSerialNo);
		strcpy(m_sDiskSerialNo, cSerialNo);
	}  
	else 
	{  
		spdlog::info("Error: failed to get disk serial number.\n");  
	}  
}*/

std::string GetWText(std::vector<char> &v){
	int w_num = v.size() + 1;
	wchar_t* w_str = new wchar_t[w_num];	w_str[w_num - 1] = 0;
	//MultiByteToWideChar(CP_UTF8, 0, (char *) &v[0], v.size(), w_str, w_num);
	//return w_str;
	std::string str(&v[0],v.size());
	return str;
}

void ComputInfo::GetDiskSerialNo()
{
	// Get Hard Drive Serial Number
	std::vector<char> hwId;
	MasterHardDiskSerial diskSerial;
	diskSerial.GetSerialNo(hwId);
	if (hwId.empty())
	{
		throw std::runtime_error("Can't retrieve hardware serial number");
	}

	//wchar_t *w_hwId = GetWText(hwId);
// 	CString cs = CString(&hwId[0], hwId.size());
// 	const char *str =  (char*)(LPCTSTR)cs;
// 	USES_CONVERSION;
// 	char * sz = T2CA(cs);
	std::string str = GetWText(hwId);
	strcpy(m_sDiskSerialNo, str.c_str());
}

void ComputInfo::UserInfo(char userinfo[256])
{
	memset(userinfo, 0, sizeof(userinfo));

	GetHostName();
//	GetMacAddress();
	GetDiskSerialNo();

	if(strcmp(m_sHostName, ""))
		strcat(userinfo, m_sHostName);

// 	if(strcmp(m_sMacAddress,""))
// 		strcat(userinfo, m_sMacAddress);

	if(strcmp(m_sDiskSerialNo, ""))
		strcat(userinfo, m_sDiskSerialNo);
}

void ComputInfo::OutputLog(char *filename)
{
	FILE *fout = NULL;
	fout = fopen(filename, "w");
	if (!fout)
	{
		spdlog::info("Error: Cann't open file %s\n", filename);
		return;
	}

	if(strcmp(m_sHostName, ""))
		fprintf(fout, "%s\n", m_sHostName);
	else
		fprintf(fout, "Failed to get host name.\n");

	if(bNowireless){
		if(strcmp(m_sMacAddress,""))
			fprintf(fout, "%s\n", m_sMacAddress);
		else
			fprintf(fout, "Failed to get Mac address.\n");
	}
	else if (bWireless)
	{
		if(strcmp(m_sMacAddressWireless,""))
			fprintf(fout, "%s\n", m_sMacAddressWireless);
		else
			fprintf(fout, "Failed to get Mac address.\n");
	}

	if(strcmp(m_sDiskSerialNo, ""))
		fprintf(fout, "%s\n", m_sDiskSerialNo);
	else
		fprintf(fout, "Failed to get disk serial number.\n");

	fclose(fout);
	fout = NULL;
}

BOOL ComputInfo::GetPhyDriveSerial(LPTSTR pModelNo,LPTSTR pSerialNo)  
{  
	BYTE IdentifyResult[sizeof(SENDCMDOUTPARAMS) + IDENTIFY_BUFFER_SIZE - 1];  
	DWORD dwBytesReturned;  
	GETVERSIONINPARAMS get_version;  
	SENDCMDINPARAMS send_cmd = { 0 };  

	HANDLE hFile = CreateFile(L"\\\\.\\PHYSICALDRIVE0", GENERIC_READ | GENERIC_WRITE,      
		FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, OPEN_EXISTING, 0, NULL);  
	if(hFile == INVALID_HANDLE_VALUE)  
		return FALSE;  

	//get version  
	DeviceIoControl(hFile, SMART_GET_VERSION, NULL, 0,  
		&get_version, sizeof(get_version), &dwBytesReturned, NULL);  

	//identify device  
	send_cmd.irDriveRegs.bCommandReg = (get_version.bIDEDeviceMap & 0x10)? ATAPI_ID_CMD : ID_CMD;  
	DeviceIoControl(hFile, SMART_RCV_DRIVE_DATA, &send_cmd, sizeof(SENDCMDINPARAMS) - 1,  
		IdentifyResult, sizeof(IdentifyResult), &dwBytesReturned, NULL);  
	CloseHandle(hFile);  

	//adjust the byte order  
	PUSHORT pWords = (USHORT*)(((SENDCMDOUTPARAMS*)IdentifyResult)->bBuffer);  
	ToLittleEndian(pWords, 27, 46, pModelNo);  
	ToLittleEndian(pWords, 10, 19, pSerialNo);  
	return TRUE;  
}  

//把WORD数组调整字节序为little-endian，并滤除字符串结尾的空格。  
void ComputInfo::ToLittleEndian(PUSHORT pWords, int nFirstIndex, int nLastIndex, LPTSTR pBuf)  
{  
	int index;  
	LPTSTR pDest = pBuf;  
	for(index = nFirstIndex; index <= nLastIndex; ++index)  
	{  
		pDest[0] = pWords[index] >> 8;  
		pDest[1] = pWords[index] & 0xFF;  
		pDest += 2;  
	}      
	*pDest = 0;  

	//trim space at the endof string; 0x20: _T(' ')  
	--pDest;  
	while(*pDest == 0x20)  
	{  
		*pDest = 0;  
		--pDest;  
	}  
}  

//滤除字符串起始位置的空格  
void ComputInfo::TrimStart(LPTSTR pBuf)  
{  
	if(*pBuf != 0x20)  
		return;  

	LPTSTR pDest = pBuf;  
	LPTSTR pSrc = pBuf + 1;  
	while(*pSrc == 0x20)  
		++pSrc;  

	while(*pSrc)  
	{  
		*pDest = *pSrc;  
		++pDest;  
		++pSrc;  
	}  
	*pDest = 0;  
}  

void ComputInfo::TcharToChar (const TCHAR * tchar, char * _char)  
{  
	int iLength ;
	iLength = WideCharToMultiByte(CP_ACP, 0, tchar, -1, NULL, 0, NULL, NULL);  
	WideCharToMultiByte(CP_ACP, 0, tchar, -1, _char, iLength, NULL, NULL);   
}  
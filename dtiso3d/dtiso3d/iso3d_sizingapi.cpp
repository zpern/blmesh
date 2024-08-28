#include <spdlog/spdlog.h> 
 #include "iso3d_sizingapi.h"



DTSizingFunc::DTSizingFunc()
{
}


DTSizingFunc::~DTSizingFunc()
{
	API_Delete_SF_Object(sfid);
}

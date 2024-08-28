#ifndef  CHAMFER_BEHAVIOR_H_
#define CHAMFER_BEHAVIOR_H_
#include <string>
struct ChamferBehavior
{
	std::string input_filename; //-i

	std::string getOutputName(std::string postfix=".pls") {
		std::string prefix= input_filename.substr(0,input_filename.find_last_of('.'));
		//std::string postfix = input_filename.substr(input_filename.find_last_of('.'), input_filename.size());
		return prefix + 'm' + postfix;
	}
};
#endif // ! MESH_SPLITTER_H_


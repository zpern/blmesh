#ifndef _SINGLETON_TERMINATE_H_
#define _SINGLETON_TERMINATE_H_
#include "singleton.h"
class SingletonTerminate :public Singleton<SingletonTerminate> {
public:
	bool terminate;
	SingletonTerminate() : Singleton<SingletonTerminate>() { terminate = false; }

};
static bool checkterminate() {
	bool ans= SingletonTerminate::GetSingletonPtr()->terminate;
	if (ans == false)
		return ans;
	else {
		SingletonTerminate::GetSingletonPtr()->terminate = false;
		return true;
	}
}
static void  resetterminate() {
	SingletonTerminate::GetSingletonPtr()->terminate = true;
}
#endif

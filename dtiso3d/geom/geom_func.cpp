#include "geom_func.h"
#include <spdlog/spdlog.h> 
 #include <mutex>
namespace GEOM_FUNC
{
	double exactinit_threadSafe() {
		static std::mutex mtx;
		mtx.lock();
		double ret = exactinit();
		mtx.unlock();
		return ret;
	}
}
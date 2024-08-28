/*
* @brife ： 支持动态插入和删除的数组，适应变化的网格数据结构，降低内存消耗
* @author ： hfye@zju.edu.cn
*/
#ifndef  _DYNAMIC_ARRAY_H_
#define  _DYNAMIC_ARRAY_H_

#include "BLMesh_define.h"
#include <array>
#include <vector>
#include <bitset>        
using INDEX_TYPE = unsigned int;

template <typename T>
class DynamicArray
{
public:
	DynamicArray() {
	}
	~DynamicArray() {}



	/*
	* @ brife: 下标操作
	*/
	T& operator [](const INDEX_TYPE& idx) {
		return buffer_[idx].first;
	}


	/*
	* @ brife: 新增一个元素
	*/
	INDEX_TYPE AddElem(const T& obj) {
#ifdef USE_DYNAMIC_ARRAY
		INDEX_TYPE idx;

		if (deleted_index_.empty()) {
			buffer_.push_back(pair<T, bool>(obj, true));
			idx = static_cast<INDEX_TYPE>(buffer_.size() - 1);
		}
		else {
			idx = deleted_index_.back();
			deleted_index_.pop_back();
			buffer_[idx] = pair<T, bool>(obj, true);
		}
		return idx;
#else
		buffer_.push_back(pair<T, bool>(obj,true));
		return static_cast<INDEX_TYPE>(buffer_.size() - 1);
#endif
		
		
	}
	/*
	* @ brife: 标识一个元素已经失效
	*/
	void DeleteElem(const INDEX_TYPE& idx) {

		if (buffer_[idx].second) {
			deleted_index_.push_back(idx);
			buffer_[idx].second = false;
		}
	}
	/*
	* @ brife: 预分配内存
	*/
	inline void Reserve(const unsigned int& size) { buffer_.reserve(size); deleted_index_.reserve(size/2); }

	/*
	* @ brife: 获取当前数组大小
	*/
	unsigned int GetSize() { return buffer_.size()- deleted_index_.size(); }

	/*
	* @ brife: 获取当前分配内存个数
	*/
	unsigned int GetCapacity() { return buffer_.capacity()+deleted_index_.capacity(); };

	/*
	* @brife: 删除buffer等数组
	*/
	void FreeMemory() {
		buffer_.clear();
		buffer_.shrink_to_fit();
		deleted_index_.clear();
	}
protected:
	std::vector<pair<T,bool>> buffer_;

	std::vector<INDEX_TYPE> deleted_index_;
};

#endif // ! _DYNAMIC_ARRAY_H_

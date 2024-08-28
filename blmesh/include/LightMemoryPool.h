#pragma once
#ifndef  LIGHT_WIGHT_MEMORY_POOL_H_
#define  LIGHT_WIGHT_MEMORY_POOL_H_
#include <vector>
#include <memory>
#define INIT(TYPE) std::vector<TYPE*> TYPE::memory_pools_ = std::vector<TYPE*>();\
std::vector<TYPE*> TYPE::extra_buffer_ = std::vector<TYPE*>();\
TYPE* TYPE::buffer_ptr_ = nullptr; 

/*
* @author yhf 2021/4/24
* @usage Inherit this class with public, for example, class A:public LightWightMomeryPool<A>, then initialize the pool array
*/
template <typename  T>
class LightWightMomeryPool {
public:
	//a simple memory allocater
	void* operator new(size_t t) {
		if (memory_pools_.empty()) {
			extra_buffer_.push_back((T*)malloc(sizeof(T)));
			memory_pools_.push_back(extra_buffer_.back());
		}
		T *buffer = memory_pools_.back();
		memory_pools_.pop_back();
		return buffer;
	}

	//a simple memory deletor
	void operator delete(void* buffer) {

		memory_pools_.push_back((T*)buffer);

	}
	

	static void preAllocate(size_t size) {
		freeMemoryInPool();
	    buffer_ptr_ = (T*)malloc(sizeof(T)*size);
		memory_pools_.resize(size);
		for (size_t i = 0; i < size; i++) {
			memory_pools_[i] = buffer_ptr_ + i;
		}
	}
	static void freeMemoryInPool() {
		if (buffer_ptr_) {
			free(buffer_ptr_);
			buffer_ptr_ = nullptr;
		}
		memory_pools_.clear();
		memory_pools_.shrink_to_fit();

	}
	static std::vector<T*> memory_pools_;

	static T* buffer_ptr_;
	static std::vector<T*> extra_buffer_;

};

#endif // ! LIGHT_WIGHT_MEMORY_POOL_H_

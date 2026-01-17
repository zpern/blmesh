#pragma once
#ifndef  LIGHT_WIGHT_MEMORY_POOL_H_
#define  LIGHT_WIGHT_MEMORY_POOL_H_

#include <vector>
#include <memory>
#include <unordered_set>
#include <new>        // ::operator new/delete
#include <cstdint>    // uintptr_t
#include <cstdlib>    // malloc/free

#define INIT(TYPE) \
std::vector<TYPE*> TYPE::memory_pools_ = std::vector<TYPE*>();\
std::vector<TYPE*> TYPE::extra_buffer_ = std::vector<TYPE*>();\
std::unordered_set<TYPE*> TYPE::extra_set_ = std::unordered_set<TYPE*>();\
TYPE* TYPE::buffer_ptr_ = nullptr; \
size_t TYPE::pool_size_ = 0; \
thread_local bool TYPE::pool_enabled_ = true;


/*
* @author yhf 2021/4/24
* @usage Inherit this class with public, for example, class A:public LightWightMomeryPool<A>, then initialize the pool array
*/
template <typename  T>
class LightWightMomeryPool {
public:
    // 池开关：禁用时 new/delete 走系统堆；preAllocate/freeMemoryInPool 不做任何事
    static void SetPoolEnabled(bool on) { pool_enabled_ = on; }
    static bool PoolEnabled() { return pool_enabled_; }

    // a simple memory allocater
    void* operator new(size_t t) {
        if (!pool_enabled_) {
            return ::operator new(t);
        }

        if (memory_pools_.empty()) {
            // 原逻辑：malloc 一个，再塞回池。这里保留，但我们记录到 extra_set_ 方便 delete 识别
            T* p = (T*)std::malloc(sizeof(T));
            if (!p) throw std::bad_alloc();
            extra_buffer_.push_back(p);
            extra_set_.insert(p);
            return p;
        }

        T* buffer = memory_pools_.back();
        memory_pools_.pop_back();
        return buffer;
    }

    // a simple memory deletor
    void operator delete(void* buffer) noexcept {
        if (!buffer) return;

        // 不依赖 pool_enabled_，只看指针是不是来自池；否则会把堆指针塞回池导致灾难
        if (isFromPool(buffer)) {
            memory_pools_.push_back((T*)buffer);
        } else {
            ::operator delete(buffer);
        }
    }

    static void preAllocate(size_t size) {
        // 关键：禁用池时，绝对不能 free 外层池
        if (!pool_enabled_) return;

        freeMemoryInPool();
        buffer_ptr_ = (T*)std::malloc(sizeof(T) * size);
        if (!buffer_ptr_) throw std::bad_alloc();

        pool_size_ = size;
        memory_pools_.resize(size);
        for (size_t i = 0; i < size; i++) {
            memory_pools_[i] = buffer_ptr_ + i;
        }
    }

    static void freeMemoryInPool() {
        // 关键：禁用池时，不许释放（避免误伤外层）
        if (!pool_enabled_) return;

        if (buffer_ptr_) {
            std::free(buffer_ptr_);
            buffer_ptr_ = nullptr;
        }
        pool_size_ = 0;

        // 你原来这里没有释放 extra_buffer_，会泄漏；顺手修掉
        for (auto p : extra_buffer_) {
            std::free(p);
        }
        extra_buffer_.clear();
        extra_set_.clear();

        memory_pools_.clear();
        memory_pools_.shrink_to_fit();
    }

    static std::vector<T*> memory_pools_;
    static T* buffer_ptr_;
    static std::vector<T*> extra_buffer_;
    static std::unordered_set<T*> extra_set_;
    static size_t pool_size_;
    static thread_local bool pool_enabled_;

private:
    static bool isFromPool(void* p) {
        // 1) 是否落在 preAllocate 的大块 buffer 里
        if (buffer_ptr_ && pool_size_ > 0) {
            std::uintptr_t begin = (std::uintptr_t)buffer_ptr_;
            std::uintptr_t end   = begin + sizeof(T) * pool_size_;
            std::uintptr_t addr  = (std::uintptr_t)p;
            if (addr >= begin && addr < end) return true;
        }
        // 2) 是否是 pool empty 时 malloc 的“额外块”
        return extra_set_.find((T*)p) != extra_set_.end();
    }
};

#endif // ! LIGHT_WIGHT_MEMORY_POOL_H_

#pragma once
#ifndef LIGHT_WIGHT_MEMORY_POOL_H_
#define LIGHT_WIGHT_MEMORY_POOL_H_

#include <vector>
#include <mutex>
#include <unordered_set>
#include <new>        // ::operator new/delete
#include <cstdint>    // uintptr_t
#include <cstdlib>    // malloc/free
#include <atomic>

// 初始化宏保持接口兼容
#define INIT(TYPE) \
template <> std::vector<TYPE*> LightWightMomeryPool<TYPE>::memory_pools_ = std::vector<TYPE*>();\
template <> std::vector<TYPE*> LightWightMomeryPool<TYPE>::extra_buffer_ = std::vector<TYPE*>();\
template <> std::unordered_set<TYPE*> LightWightMomeryPool<TYPE>::extra_set_ = std::unordered_set<TYPE*>();\
template <> TYPE* LightWightMomeryPool<TYPE>::buffer_ptr_ = nullptr; \
template <> size_t LightWightMomeryPool<TYPE>::pool_size_ = 0; \
template <> bool LightWightMomeryPool<TYPE>::pool_enabled_ = true;\
template <> std::mutex LightWightMomeryPool<TYPE>::pool_mutex_ = std::mutex();

/*
* @author yhf 2021/4/24 (modified)
* @usage Inherit this class with public, for example, class A:public LightWightMomeryPool<A>, then initialize the pool array
*/
template <typename T>
class LightWightMomeryPool {
public:
    // 控制是否启用内存池，保持原接口
    static void SetPoolEnabled(bool on) { pool_enabled_ = on; }
    static bool PoolEnabled() { return pool_enabled_; }

    // 重载new操作符，修复内存申请逻辑
    void* operator new(size_t t) {
        // 检查申请大小是否匹配
        if (t != sizeof(T)) {
            return ::operator new(t);
        }

        if (!pool_enabled_) {
            return ::operator new(t);
        }

        std::lock_guard<std::mutex> lock(pool_mutex_);
        
        // 优先从内存池获取
        if (!memory_pools_.empty()) {
            T* p = memory_pools_.back();
            memory_pools_.pop_back();
            return p;
        }

        // 内存池为空时，申请新内存并记录
        T* p = static_cast<T*>(::operator new(sizeof(T)));
        if (!p) throw std::bad_alloc();
        extra_buffer_.push_back(p);
        extra_set_.insert(p);
        return p;
    }

    // 重载delete操作符，修复内存释放逻辑
    void operator delete(void* buffer) noexcept {
        if (!buffer) return;

        std::lock_guard<std::mutex> lock(pool_mutex_);

        // 判断是否来自内存池
        if (isFromPool(buffer)) {
            // 归还到内存池
            memory_pools_.push_back(static_cast<T*>(buffer));
        } else {
            // 非池内内存，直接释放
            ::operator delete(buffer);
        }
    }

    // 预分配内存池，保持原接口
    static void preAllocate(size_t size) {
        if (!pool_enabled_) return;

        std::lock_guard<std::mutex> lock(pool_mutex_);
        
        // 先释放已有内存
        freeMemoryInPoolUnlocked();
        
        // 分配连续内存块，使用operator new保证内存对齐
        buffer_ptr_ = static_cast<T*>(::operator new(sizeof(T) * size));
        if (!buffer_ptr_) throw std::bad_alloc();

        pool_size_ = size;
        memory_pools_.reserve(size);
        // 初始化内存池列表
        for (size_t i = 0; i < size; i++) {
            memory_pools_.push_back(buffer_ptr_ + i);
        }
    }

    // 释放内存池所有内存，保持原接口
    static void freeMemoryInPool() {
        if (!pool_enabled_) return;

        std::lock_guard<std::mutex> lock(pool_mutex_);

        freeMemoryInPoolUnlocked();
    }

    // 静态成员变量声明
    static std::vector<T*> memory_pools_;
    static T* buffer_ptr_;
    static std::vector<T*> extra_buffer_;
    static std::unordered_set<T*> extra_set_;
    static size_t pool_size_;
    static bool pool_enabled_;
    static std::mutex pool_mutex_;  // 新增：线程安全锁

private:
    // 判断指针是否来自内存池
    static bool isFromPool(void* p) {
        // 1) 检查是否在预分配的buffer范围内
        if (buffer_ptr_ && pool_size_ > 0) {
            std::uintptr_t begin = reinterpret_cast<std::uintptr_t>(buffer_ptr_);
            std::uintptr_t end   = begin + sizeof(T) * pool_size_;
            std::uintptr_t addr  = reinterpret_cast<std::uintptr_t>(p);
            if (addr >= begin && addr < end) return true;
        }
        // 2) 检查是否是额外申请的内存
        return extra_set_.find(static_cast<T*>(p)) != extra_set_.end();
    }

    // 不加锁的释放内存池
    static void freeMemoryInPoolUnlocked()
    {
        // 释放预分配的连续内存块（与 ::operator new 配对）
        if (buffer_ptr_) {
            ::operator delete(buffer_ptr_);
            buffer_ptr_ = nullptr;
        }
        pool_size_ = 0;

        // 释放额外申请的内存
        for (auto p : extra_buffer_) {
            ::operator delete(p);
        }
        extra_buffer_.clear();
        extra_set_.clear();

        // 清空内存池
        memory_pools_.clear();
        memory_pools_.shrink_to_fit(); // 可选：频繁调用可考虑去掉以提升性能
    }
};

#endif // !LIGHT_WIGHT_MEMORY_POOL_H_

#ifndef SIMPLEOBJECTPOOL_H
#define SIMPLEOBJECTPOOL_H

#include "tbb/concurrent_queue.h"
#include "tbb/concurrent_vector.h"

namespace icy {
    template<class T>
    class ConcurrentPool;
}

template<class T>
class icy::ConcurrentPool
{

public:
    ConcurrentPool(int initialSize);
    ~ConcurrentPool();
    T* take();
    void release(T* obj);
    void releaseAll();
    void printout(); // for testing

private:
    tbb::concurrent_queue<T*> available;      // items that are free to use
    tbb::concurrent_vector<T*> registry;       // all items of the pool
};

template<class T>
icy::ConcurrentPool<T>::ConcurrentPool(int initialSize)
{
    registry.reserve(initialSize*2);
    for(int i=0;i<initialSize;i++)
    {
        T* obj = new T;
        available.push(obj);
        registry.push_back(obj);
    }
}

template<class T>
icy::ConcurrentPool<T>::~ConcurrentPool()
{
    for(auto &x : registry) delete x;
}

template <class T>
T* icy::ConcurrentPool<T>::take()
{
    T *obj;
    bool result = available.try_pop(obj);
    if(!result)
    {
        obj = new T;
        registry.push_back(obj);
    }
    return obj;
}

template<class T>
void icy::ConcurrentPool<T>::releaseAll()
{
    available.clear();
    for(unsigned i=0;i<registry.size();i++) available.push(registry[i]);
}

template<class T>
void icy::ConcurrentPool<T>::release(T* obj)
{
    available.push(obj);
}

#endif

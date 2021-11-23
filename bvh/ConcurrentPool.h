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
    ConcurrentPool& operator=(ConcurrentPool&) = delete;

    T* take();
    void release(T* p) {available.push(p);}
    void release(std::vector<T*> &vec);
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
    T *p;
    bool result = available.try_pop(p);
    if(!result)
    {
        p = new T;
        registry.push_back(p);
    }
    return p;
}

template<class T>
void icy::ConcurrentPool<T>::releaseAll()
{
    available.clear();
    for(unsigned i=0;i<registry.size();i++) available.push(registry[i]);
}

template<class T>
void icy::ConcurrentPool<T>::release(std::vector<T*> &vec)
{
    for(T* p : vec) available.push(p);
    vec.clear();
}


#endif

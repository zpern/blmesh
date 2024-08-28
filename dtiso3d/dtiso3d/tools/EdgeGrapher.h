#pragma once
#include <math.h>
#include <algorithm>
#include <map>
#include <cassert>
#include <utility>
#include <algorithm>
#include <cmath>
using std::max;
using std::min;
template<class T>
class EdgeGrapher
{
public:
	void addEdgeOverride(int n1, int n2, const T& val);
	bool getEdgeVal(int n1, int n2, T& ret);
	T getEdgeVal(int n1, int n2); //不进行安全检查
	bool hasEdge(int n1, int n2);
	void clear() {
		map_.clear();
	}
	std::map<std::pair<int, int>, T> map_;
private:
	std::pair<int, int> getKey(int n1, int n2) {
		return std::pair<int,int>(min(n1, n2), max(n1, n2));
	}
};

template<class T>
T EdgeGrapher<T>::getEdgeVal(int n1, int n2)
{
	return map_[getKey(n1, n2)];
}

template<class T>
bool EdgeGrapher<T>::getEdgeVal(int n1, int n2, T& ret)
{
	assert(n1 != n2);
	if (!hasEdge(n1, n2))
		return false;
	ret = map_[getKey(n1, n2)];
	return true;
}

template<class T>
bool EdgeGrapher<T>::hasEdge(int n1, int n2)
{
	return map_.find(getKey(n1, n2)) != map_.end();
}


template<class T>
void EdgeGrapher<T>::addEdgeOverride(int n1, int n2, const T& val)
{
	assert(n1 != n2);
	map_[getKey(n1, n2)] = val;
}


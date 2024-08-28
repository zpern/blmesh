/*
* File: IntIntMap.cpp
* ------------------
* This file implements the IntIntMap.h interface.
*/
#include <stdio.h>
#include <stdlib.h>
#include <spdlog/spdlog.h> 
 #include "IntIntMap.h"

/*
* Implementation notes: Map constructor
* -------------------------------------
* The constructor allocates the array of buckets and initializes
* each bucket to the empty list.
*/
IntIntMap::IntIntMap() 
{
	nBuckets = INITIAL_SIZE;
	buckets = new cellT *[nBuckets];
	if (!buckets)
	{
		spdlog::info("Not enough memory.\n");
		exit(1);
	}
	for (int i = 0; i < nBuckets; i++) 
	{
		buckets[i] = NULL;
	}
	nEntries = 0;
}

IntIntMap::IntIntMap(int initSize) 
{
	nBuckets = initSize > INITIAL_SIZE ? initSize : INITIAL_SIZE;
	buckets = new cellT *[nBuckets];
	if (!buckets)
	{
		spdlog::info("Not enough memory.\n");
		exit(1);
	}
	for (int i = 0; i < nBuckets; i++) 
	{
		buckets[i] = NULL;
	}
	nEntries = 0;
}

/*
* Implementation notes: Map destructor
* -------------------------------------
* The destructor must deallocate every cell (which it can do by
* calling clear) and then free the dynamic bucket array.
*/
IntIntMap::~IntIntMap() 
{
#if 0
#ifdef _DEBUG
	int nonzero = 0, entry_num;
	CellT *cell = NULL;
	spdlog::info("IntIntMap:#Entries:{}.\n", nEntries);
	for (int i = 0; i < nBuckets; i++) 
	{
		entry_num = 0;
		cell = buckets[i];
		while (cell)
		{
			entry_num++;
			cell = cell->link;
		}
		if (entry_num > 0)
		{
			nonzero++;
			spdlog::info("Bucket {}: #Entries = {}.\n", i+1, entry_num); 
		}
	}
	spdlog::info("#Nonempty buckets = {}.\n", nonzero);
#endif
#endif
	clear();
	if (buckets)
		delete[] buckets;
}

/*
* Implementation notes: size, isEmpty
* -----------------------------------
* These methods can each be implemented in a single line
* because the size is stored in the nEntries instance variable.
*/
int IntIntMap::size() 
{
	return nEntries;
}

bool IntIntMap::isEmpty() 
{
	return nEntries <= 0;
}

/*
* Implementation notes: clear
* ---------------------------
* This method calls the recursive deleteChain method for each
* bucket chain.
*/
void IntIntMap::clear() 
{
	for (int i = 0; i < nBuckets; i++) 
	{
		deleteChain(buckets[i]);
	}
	nEntries = 0;
}

/*
* Implementation notes: put
* -------------------------
* This method first looks to see whether the key already
* exists in the map by calling the findCell method. If one
* exists, this method simply changes the value; if not, the
* implementation adds a new cell to the beginning of the chain.
*/
void IntIntMap::put(int key, int value) 
{
	int index = hash(key) % nBuckets;
	cellT *cell = findCell(buckets[index], key);
	if (cell == NULL) 
	{
		cell = new cellT;
		if (!cell)
		{
			spdlog::info("Not enough memory.\n");
			exit(1);
		}
		cell->key = key;
		cell->link = buckets[index];
		buckets[index] = cell;
		nEntries++;
	}
	cell->value = value;
}

/*
 * Implementation notes: get
 * --------------------------------------
 * These methods uses findCell to find the key in the map, which is
 * where all the real work happens.
 */
bool IntIntMap::get(int key, int *value) 
{
	cellT *cell = findCell(buckets[hash(key) % nBuckets], key);
	if (cell == NULL) {
		return false; /* an invalid value */
	}
	*value = cell->value;
	return true;
}

bool IntIntMap::containsKey(int key) 
{
	return findCell(buckets[hash(key) % nBuckets], key) != NULL;
}

/*
* Implementation notes: remove
* ----------------------------
* The remove method cannot use the findCell method as it
* stands because it needs a pointer to the previous entry.
* Because that code is used only in this method, the loop
* through the cells in a chain is reimplemented here and
* therefore does not add any cost to the get/put operations.
*/
void IntIntMap::remove(int key) 
{
	int index = hash(key) % nBuckets;
	cellT *prev = NULL;
	cellT *cp = buckets[index];
	while (cp != NULL && cp->key != key) 
	{
		prev = cp;
		cp = cp->link;
	}
	if (cp != NULL) 
	{
		if (prev == NULL) 
		{
			buckets[index] = cp->link;
		} 
		else 
		{
			prev->link = cp->link;
		}
		delete cp;
		nEntries--;
	}
}

/* Private methods */
/*
* Implementation notes: hash
* Usage: bucket = hash(key);
* --------------------------
* This function takes the key and uses it to derive a hash code,
* which is a nonnegative integer. The hash code is computed
* using a method called linear congruence.
*/
int IntIntMap::hash(int key) 
{
	const long MULTIPLIER = -1664117991L;
	unsigned long hashcode = 0;
	int charLength = sizeof(int), i, bitMove;
	unsigned long bitValue;
	for (i = 0; i < charLength; i++) 
	{
		bitMove = 8*i;
		bitValue = ((0xff << bitMove) & key) >> bitMove;
		hashcode = hashcode * MULTIPLIER +  bitValue;
	}

	return hashcode & ((unsigned) -1 >> 1);
}

/*
* Implementation notes: findCell
* Usage: cell = findCell(chain, key);
* -----------------------------------
* This function finds a cell in the chain that matches key.
* If a match is found, findCell returns a pointer to that cell;
* if not, findCell returns NULL.
*/
IntIntMap::cellT *IntIntMap::findCell(cellT *chain, int key) 
{
	for (cellT *cp = chain; cp != NULL; cp = cp->link) 
	{
		if (cp->key == key) 
			return cp;
	}
	return NULL;
}

/*
* Private method: deleteChain
* ---------------------------
* This method deletes all of the cells in a bucket chain.
* It operates recursively by freeing the rest of the chain
* and the freeing the current cell.
*/
void IntIntMap::deleteChain(cellT *chain) 
{
	if (chain != NULL) 
	{
		deleteChain(chain->link);
		delete chain;
	}
}
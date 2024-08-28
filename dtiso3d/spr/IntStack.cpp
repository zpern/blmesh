#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <spdlog/spdlog.h> 
 #include "IntStack.h"

#define STACK_SIZE_ADD_RATIO 1.2
#define STACK_SIZE_MIN_NUM 8192 /* 2**13 */
IntStack::IntStack() : top_(-1), elems_(NULL), maxSize_(0)
{
}

IntStack::IntStack(int maxStackSize) : top_(-1), maxSize_(maxStackSize)
{
	elems_ = maxStackSize > 0 ? (int*)malloc(sizeof(int)*maxStackSize) : NULL;
	if (!elems_ && maxStackSize > 0)
	{
		spdlog::info("Not enough memory.\n");
		exit(1);
	}
}

IntStack::~IntStack()
{
	if (elems_)
		free(elems_);
}

void IntStack::push(int elem)
{
	int newSize, *newElems = NULL;
    if (top_ + 1 >= maxSize_)
	{
		newSize = maxSize_*STACK_SIZE_ADD_RATIO > STACK_SIZE_MIN_NUM ? maxSize_*STACK_SIZE_ADD_RATIO : STACK_SIZE_MIN_NUM;
		assert(newSize > top_ + 1);
		newElems = (int*)realloc(elems_, sizeof(int)*newSize);
		if (!newElems)
		{
			spdlog::info("Not enough memory.\n");
			exit(1);
		}
		elems_ = newElems;
		maxSize_ = newSize;
	}

    elems_[++top_] = elem;
}

void IntStack::pop()
{
    if (top_ + 1 == 0)
    {
		spdlog::info("Illegal access in IntStack.\n");
		exit(1);
	}

    --top_;
}

int IntStack::top()
{
   if (top_ + 1 == 0)
   {
		spdlog::info("Illegal access in IntStack.\n");
		exit(1);
	}

    return elems_[top_];
}

bool IntStack::empty()
{
    return top_ + 1 == 0;
}

void IntStack::reset()
{
	top_ = -1;
}
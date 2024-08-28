#pragma once

// algebra.h
#include "math.h"
//#include "symbols.h"

#include <complex>
using namespace std;

#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__

////////////////////////////////////////////////////
void matrixXvector(int nRow,int nCol, complex<double>** a, complex<double>* v, complex<double>* w);
void matrixXvector_complex(int nRow,int nCol, double* a, complex<double>* v, complex<double>* w);
void matrixXvector(long nRow,long nCol, complex<double>* a, complex<double>* v, complex<double>* w);

void matrixAddXvector(int nRow,int nCol, complex<double>** a, complex<double>* v, complex<double>* w);
void matrixAddXvector(int nRow,int nCol, complex<double>** a, double* v, complex<double>* w);

void matrixMinusXvector(int nRow,int nCol, complex<double>** a, complex<double>* v,complex<double>* w);
void matrixMinusXvector(int nRow,int nCol, complex<double>** a, double* v, complex<double>* w);

void matrixScaleAddXvector(int nRow,int nCol, double alpha, complex<double>** a,complex<double>* v, complex<double>* w);
void matrixScaleAddXvector(int nRow,int nCol, double alpha, complex<double>** a, double* v, complex<double>* w);

void transMatrixXvector(int nRow,int nCol, complex<double>** a, complex<double>* v, complex<double>* w);
void transMatrixXvector(int nRow,int nCol, complex<double>** a, double* v, complex<double>* w);

void transMatrixAddXvector(int nRow,int nCol, complex<double>** a, complex<double>* v,complex<double>* w);
void transMatrixAddXvector(int nRow,int nCol, complex<double>** a, double* v, complex<double>* w);

void transMatrixMinusXvector(int nRow,int nCol, complex<double>** a, complex<double>* v, complex<double>* w);
void transMatrixMinusXvector(int nRow,int nCol, complex<double>** a, double* v, complex<double>* w);

void initVector(long n, complex<double> *v, int incx, complex<double> alfa);
void initVector0(long n, complex<double> *v, int incx);
void initVector0(long n,complex<double> *v);
void initVector(long n,complex<double> *v,complex<double> t);
void initVector(long n,complex<double> *v,double t);
void initVector(long n,long m, complex<double> **v);


void VectorXVector(long n, complex<double> *v1, double *v2);
complex<double> VectorXVector_complex(int n, complex<double> *v1,double *v2);

int cinv(complex<double> **a,long n,double eps);
int cinv1(complex<double>** a,long n,double eps);
double norm0(complex<double>* a,long n);


complex<double> getScatterExactResult(double costh,double ka,double kr);
complex<double> Eix(double a, int n);//求幂函数a^n
double Wigner3jsymbol(int j1,int j2,int j3,int m1,int m2,int m3,int i);//求Wigner3jsymbol的结果
double Wigner3jsymbol0(int j1,int j2,int j3);
complex<double> TheResultofWigner3j(int i,int j,int n,int m,int l);
double Factorial(int n);//求阶乘






//////////////////////////////////////////////////
double dmax(double a, double b);
double dmin(double a, double b);
double dmax(double a, double b, double c);
double dmin(double a, double b, double c);
double dmax(double a, double b, double c, double d);
double dmin(double a, double b, double c, double d);
double sign(double a);
int imax(int a, int b);
int imin(int a, int b);
int imax(int a, int b, int c);
int imin(int a, int b, int c);
int imax(int a, int b, int c, int d);
int imin(int a, int b, int c, int d);

//w=a*v
void matrixXvector(int nRow,int nCol, double** a, double* v, double* w);
//w+=a*v
void matrixAddXvector(int nRow,int nCol, double** a, double* v, double* w);
void matrixMinusXvector(int nRow,int nCol, double** a, double* v, double* w);
//w+=alpha*a*v
void matrixScaleAddXvector(int nRow,int nCol, double alpha, double** a, double* v, double* w);
//w=a'*v
void transMatrixXvector(int nRow,int nCol, double** a, double* v, double* w);
//w+=a'*v
void transMatrixAddXvector(int nRow,int nCol, double** a, double* v, double* w);
void transMatrixMinusXvector(int nRow,int nCol, double** a, double* v, double* w);
//w+=alpha*a'*v
void transMatrixScaleAddXvector(int nRow,int nCol, double alpha, double** a, double* v, double* w);
//a*=alpha
void scaleXmatrix(int nRow,int nCol, double** a, double alpha);
// a += b
void matrixAddTomatrix(int nRow,int nCol, double** a, double** b);
// c = a + b
void matrixAddmatrix(int nRow,int nCol, double** c, double** a, double** b);
/* c = a * b */
void matrixXmatrix(int nRow, int nCol, int K, double** c, double** a, double** b);
void float_matrixXmatrix(int nRow, int nCol, int K, double** c, double** a, double** b);
// c += a * b
void matrixAddXmatrix(int nRow, int nCol, int K, double** c, double** a, double** b);
double dnorm2(int n, double* v);

void initVector(long n, double *v, int incx, double alfa);
void initVector0(long n, double *v, int incx);
void initVector0(long n, double *v);
//v1[i]*=v2[i]
void VectorXVector(long n, double *v1, double *v2);
////////////////
int mxInvert(double a[],int n);
int eqGauss(double a[],double b[],int n);
int eqGausJdn(double a[],double b[],int n, int m);

void dswap(double* x,double* y);
int inv0(double** a,long n,double eps);
int inv1(double** a,long n,double eps);
double norm0(double* a,long n);
double norm1(double* a,long n);
int gmres(long n, double* A, double* b, int restart, double tol, int maxit,
		  double* M, double* x0, double* x, int& flag, int* iter);

int gmres_solver(long n, double* A, double* b, double tolerance, int& flag, int* iter);
//the following functions have been tested in (int CSurface::setMatrixH())
int ludcmp(double**a,int n,int* indx, double eps);
void lubksb(double **a, int n, int *indx, double* b);//A*x=b
void lubksb_lower(double **a, int n, int *indx, double* b);//L*x=b
void lubksb_lower2(double **a, int n, int *indx, double* b);//x*L=b
void lubksb_upper(double **a, int n, double* b);//U*x=b
void lubksb_upper2(double **a, int n, double* b);//x*U=b

int dsvdcmp(double **a, int m, int n, double* w, double **v);
void dsvbksb(double **u, double* w, double **v, int m, int n, double* b, double* x);

int choleskydcmp(double **a, int n, double* p);//choldc(...)
void choleskybksb(double **a, int n, double* p, double* b, double* x);//cholsl(...)

int qrdcmp(double **a, int n, double *c, double *d);
void qrbksb(double **a, int n, double* c, double* d, double* b);//qrsolv
int XuShiliangQR(double** a,int m,double** q);
int XuShiliangQR2(int m, int n, double** a, double** q);//m>=n

double dpythag(double a, double b);

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define IMIN(a,b) ((a) < (b) ? (a) : (b))
#define SQR(a) ((a) == 0.0 ? 0.0 : (a)*(a))

/////////////
const static double gsINT_11_g5[5]={-0.906179845938664,-0.538469310105683,
						             0.0,0.538469310105683,0.906179845938664};
const static double gsINT_11_w5[5]={0.236926885056189,0.478628670499366,
						             0.568888888888889,0.478628670499366,0.236926885056189};

const static double gsINT_11_g4[4]={-0.861136311594053,-0.339981043584856,0.339981043584856,0.861136311594053};
const static double gsINT_11_w4[4]={0.347854845137454,0.652145154862546,0.652145154862546,0.347854845137454};

const static double gsINT_11_g3[3]={-0.774596669241483,0.0,0.774596669241483};
const static double gsINT_11_w3[3]={0.555555555555556,0.888888888888889,0.555555555555556};

const static double gsINT_11_g2[2]={-0.577350269189626,0.577350269189626};
const static double gsINT_11_w2[2]={1.0,1.0};

const static double gsINT_11_g6[6]={-0.932469514203152,-0.661209386466265,-0.238619186083197,
                            0.238619186083197, 0.661209386466265, 0.932469514203152};
												
const static double gsINT_11_w6[6]={ 0.171324492379170, 0.360761573048139, 0.467913934572691,
						    0.467913934572691, 0.360761573048139, 0.171324492379170};

const static double gsINT_11_g8[8]={-0.960289856497536,-0.796666477413627,-0.525532409916329,-0.183434642495650,
                            0.183434642495650, 0.525532409916329, 0.796666477413627, 0.960289856497536};
                           
const static double gsINT_11_w8[8]={ 0.101228536290376, 0.222381034453374, 0.313706645877887, 0.362683783378362,
                            0.362683783378362, 0.313706645877887, 0.222381034453374, 0.101228536290376 };
					     
const static double gsINT_11_g10[10]={-0.973906528517172,-0.865063366688985,-0.679409568299024, -0.433395394129247, -0.148874338981631,
                              0.148874338981631, 0.433395394129247, 0.679409568299024,  0.865063366688985,  0.973906528517172 };
                            
const static double gsINT_11_w10[10]={0.066671344308688,  0.149451349150580, 0.219086362515982, 0.269266719309996, 0.295524224714753,
                             0.295524224714753,  0.269266719309996, 0.219086362515982, 0.149451349150580, 0.066671344308688 };

const static double gsINT_11_g1[1]={0.0};
const static double gsINT_11_w1[1]={2.0};

/////////////
const static double gsINT_01_g5[5]={0.046910077030668, 0.2307653449471585,
						             0.5, 0.7692346550528415, 0.953089922969332};
const static double gsINT_01_w5[5]={0.1184634425280945,0.239314335249683,
						             0.2844444444444444,0.239314335249683,0.1184634425280945};

const static double gsINT_01_g4[4]={0.0694318442029735,0.330009478207572,0.669990521792428,0.9305681557970265};
const static double gsINT_01_w4[4]={0.173927422568727,0.326072577431273,0.326072577431273,0.173927422568727};

const static double gsINT_01_g3[3]={0.1127016653792585,0.5,0.8872983346207415};
const static double gsINT_01_w3[3]={0.277777777777778,0.4444444444444444,0.277777777777778};

const static double gsINT_01_g2[2]={0.211324865405187,0.788675134594813};
const static double gsINT_01_w2[2]={0.5,0.5};

							
//const static double gsINT_01_g1[1]={0.5};
//const static double gsINT_01_w1[1]={1.0};

/////////////////////////////////////////////////////////////////////
const static double gsTrgl1_g3[3]={0.5, 0.0, 0.5};
const static double gsTrgl2_g3[3]={0.5, 0.5, 0.0};
const static double gsTrgl3_g3[3]={0.0, 0.5, 0.5};
const static double gsTrigl_w3[3]={0.3333333333, 0.3333333333, 0.3333333333};

const static double gsTrgl1_g4[4]={0.33333333333333, 0.6, 0.2, 0.2};
const static double gsTrgl2_g4[4]={0.33333333333333, 0.2, 0.6, 0.2};
const static double gsTrgl3_g4[4]={0.33333333333333, 0.2, 0.2, 0.6};
const static double gsTrigl_w4[4]={-0.5625, 0.52083333333333,0.52083333333333,0.52083333333333};
							      
const static double gsTrgl1_g7[7]={0.3333333333,0.0597158717,0.4701420641,
							       0.4701420641,0.7974269853,0.1012865073,0.1012865073};
const static double gsTrgl2_g7[7]={0.3333333333,0.4701420641,0.0597158717,
								   0.4701420641,0.1012865073,0.7974269853,0.1012865073};
const static double gsTrgl3_g7[7]={0.3333333333,0.4701420641,0.4701420641,
							       0.0597158717,0.1012865073,0.1012865073,0.7974269853};
const static double gsTrigl_w7[7]={0.2255,0.1323941527,0.1323941527,
							       0.1323941527,0.1259391805,0.1259391805,0.1259391805};

///////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////
template <class T, int size>
class CStack
{
public:

	CStack()//constructor
	{
		depth=size;
		top=-1;
	} 
	~CStack()
	{
		top=-1;
	};

	void reset() {top=-1;}
	int push(T& item)
	{
        if(!isFull())
		{
			buffer[++top] = item;
			return 0;
		}
		return 1;
	}
	int pop(T& item)
	{
        if(!isEmpty())
		{
			item=buffer[top--];
			return 0;
		}
		return 1;
	}
	int isEmpty() const { return top == -1;}
	int isFull() const { return top == depth-1;}

private:
	int top;
	int depth;
	T  buffer[size];
};


// Template class for a single concatenated list
template<class T>
class CSingleList
{
public:
  // class for a single item in the list
	struct item
	{
		T data;
		item* next;
		item(item* n, const T& e) : data(e), next(n){};
	};

	item *curItem;   // currently selected item
protected:
	// first element of list
	item *first_;
	// size of the list aka number of stored items
	long size_;

public:
	// constructors and destructors
	CSingleList() : first_(0x0), size_(0) { curItem=0x0;};

	~CSingleList()
	{
		clear();
	}

	// add a new object to the beginning of list
	void push_front(const T& x)
	{
		item *tmp = new item(first_, x);

		first_ = tmp;
		++size_;
	}

	// pop and erase the first element from the list
	//in this case the data must be a pointer to some class or struct
	T pop_front()
	{
        if(size_ == 0) return 0x0;
		item* tmp = first_->next;
		T dt=first_->data; delete first_;
		first_ = tmp;
		--size_;
		return dt;
	}

	// empty the list
	void clear()
	{
		item *next = first_, *tmp;
		while ((tmp=next)!= 0x0)
		{
			next = tmp->next;
			delete tmp;
		}

		first_ = 0x0;
		size_ = 0;
	}
	// remove all elements from list and clear the list
	//in this case the data must be a pointer to some class or struct
	void remove()
	{
		item *next = first_, *tmp;
		while ((tmp=next)!= 0x0)
		{
			next = tmp->next;
			delete tmp->data;
			delete tmp;
		}

		first_ = 0x0;
		size_ = 0;
	}

	int getSize()
	{
		return size_;
	}

	void iterateBegin(){ curItem=first_; }

	item *iterateCurrent()
	{
		return curItem;
	}
	void iterateNext()
	{
		ASSERT(curItem);
		curItem=curItem->next;
	}
	void deleteCurrent()
	{
		ASSERT(curItem);
        if(curItem == first_)
		{
			pop_front();
			curItem=first_;
		}
		else
		{
			item* tmp=first_;
			while (tmp)
			{
		        if(tmp->next == curItem)
				{
					tmp->next = curItem->next;
					delete curItem;
					curItem=tmp->next;
					--size_;
					break;
				}
				tmp=tmp->next;
			}
		}
	}
	void deleteItem(T var)//the values are unique in the list
	{
        if(var == first_->data)
		{
			pop_front();
		}
		else
		{
			item* tmp=first_;
			while (tmp->next)
			{
		        if(tmp->next->data == var)
				{
					curItem = tmp->next;
					tmp->next = curItem->next;
					delete curItem;
					--size_;
					break;
				}
				tmp=tmp->next;
			}
		}
	}
	//add to an array 
	void addToArrayOrderly(T *Tarray, long iStart)
	{
		T dt=pop_front();
		long k=0;
		while(dt)
		{
			Tarray[iStart+k]=dt; k++;
			dt=pop_front();
		}
	}
	void addToArrayConversely(T *Tarray, long iStart)
	{
		T dt=pop_front();
		long k=0;
		while(dt)
		{
			Tarray[iStart-k]=dt; k++;
			dt=pop_front();
		}
	}

};


#endif


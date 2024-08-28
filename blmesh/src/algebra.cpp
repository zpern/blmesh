// Written by Jianming Zhang (zhangjm@homer.shinshu-u.ac.jp)
// Department of Mechanical Systems engineering, Shinshu University
// Copyright (c) November 20, 2002
//
///////////////////////////////////////////////////////////////////

#include <spdlog/spdlog.h> 
 #include "algebra.h"

//////////////////////////////////////////////globals
double norm0(complex<double>* a,long n)//length: n
{
	  double temp;
	  temp=0.0;
	  for(long i=0; i<n; i++)
	  {
		  temp+=norm(a[i]);
	  }
	  temp=sqrt(temp);
	  return temp;
}
void matrixXvector(int nRow,int nCol, double* a, complex<double>* v, complex<double>* w)
{
	for(int i=0; i<nRow; i++)
	{
		w[i]=0.0;
		for(int j=0; j<nCol; j++)
			w[i]+=a[i+nCol*j]*v[j];
	}
}
void matrixXvector(long nRow,long nCol, complex<double>* a, complex<double>* v, complex<double>* w)
{
	for(int i=0; i<nRow; i++)
	{
		w[i]=0.0;
		for(int j=0; j<nCol; j++)
			w[i]+=a[i+nCol*j]*v[j];
	}
}
void matrixXvector(int nRow,int nCol, double ** a, complex<double>* v, complex<double>* w)
{
	for(int i=0; i<nRow; i++)
	{
		w[i]=0.0;
		for(int j=0; j<nCol; j++)
			w[i]+=a[i][j]*v[j];
	}
}
void matrixAddXvector(int nRow,int nCol, complex<double>** a, complex<double>* v, complex<double>* w)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
			w[i]+=a[i][j]*v[j];
	}
}void matrixAddXvector(int nRow,int nCol, complex<double>** a, double* v,complex<double>* w)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
			w[i]+=a[i][j]*v[j];
	}
}
void matrixMinusXvector(int nRow,int nCol, complex<double>** a, complex<double>* v, complex<double>* w)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
			w[i]-=a[i][j]*v[j];
	}
}
void matrixMinusXvector(int nRow,int nCol, complex<double>** a, double* v, complex<double>* w)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
			w[i]-=a[i][j]*v[j];
	}
}
void matrixScaleAddXvector(int nRow,int nCol, double alpha, complex<double>** a, complex<double>* v, complex<double>* w)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
			w[i]+=a[i][j]*v[j]*alpha;
	}
}
void matrixScaleAddXvector(int nRow,int nCol, double alpha, complex<double>** a, double* v, complex<double>* w)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
			w[i]+=a[i][j]*v[j]*alpha;
	}
}
void transMatrixXvector(int nRow,int nCol, complex<double>** a, complex<double>* v, complex<double>* w)
{
	for(int i=0; i<nCol; i++)
	{
		w[i]=0.0;
		for(int j=0; j<nRow; j++)
			w[i]+=v[j]*a[j][i];
	}
}
void transMatrixXvector(int nRow,int nCol, complex<double>** a, double* v, complex<double>* w)
{
	for(int i=0; i<nCol; i++)
	{
		w[i]=0.0;
		for(int j=0; j<nRow; j++)
			w[i]+=v[j]*a[j][i];
	}
}
void transMatrixAddXvector(int nRow,int nCol, complex<double>** a, complex<double>* v, complex<double>* w)
{
	for(int i=0; i<nCol; i++)
	{
		for(int j=0; j<nRow; j++)
			w[i]+=v[j]*a[j][i];
	}
}
void transMatrixAddXvector(int nRow,int nCol, complex<double>** a, double* v, complex<double>* w)
{
	for(int i=0; i<nCol; i++)
	{
		for(int j=0; j<nRow; j++)
			w[i]+=v[j]*a[j][i];
	}
}
void transMatrixMinusXvector(int nRow,int nCol, complex<double>** a, complex<double>* v, complex<double>* w)
{
	for(int i=0; i<nCol; i++)
	{
		for(int j=0; j<nRow; j++)
			w[i]-=v[j]*a[j][i];
	}
}
void transMatrixMinusXvector(int nRow,int nCol, complex<double>** a, double* v, complex<double>* w)
{
	for(int i=0; i<nCol; i++)
	{
		for(int j=0; j<nRow; j++)
			w[i]-=v[j]*a[j][i];
	}
}
void initVector(long n, complex<double> *v, int incx, complex<double> alfa)
{
	complex<double> *p=v;
	for(long i=0; i<n; i++)
	{
		*p=alfa; p+=incx;
	}
}
void initVector(long n, complex<double> *v,complex<double> t)
{
	complex<double> *p=v;
	for(long i=0; i<n; i++)
	{
		*p=t; p++;
	}
}
void initVector(long n,long m, complex<double> **v)
{
	long i,j;
	complex<double> **p=v;
	for(j=0;j<m;j++)
	for(i=0; i<n; i++)
	{
		*p[j]=0.0; p++;
	}
}
void initVector(long n, complex<double> *v,double t)
{
	complex<double> *p=v;
	for(long i=0; i<n; i++)
	{
		*p=t; p++;
	}
}
void initVector0(long n, complex<double> *v, int incx)
{
	complex<double> *p=v;
	for(long i=0; i<n; i++)
	{
		*p=0.0; p+=incx;
	}
}
void initVector0(long n, complex<double> *v)
{
	complex<double> *p=v;
	for(long i=0; i<n; i++)
	{
		*p=0.0; p++;
	}
}
complex<double> VectorXVector_complex(int n, complex<double> *v1,double *v2)
{
	complex<double> *p1=v1;
	double *p2=v2;
	complex<double> tmp;

	for(long i=0; i<n; i++)
	{
		tmp+=(*p1)*(*p2);
		 p1++; p2++;
	}
	return(tmp);
}
void VectorXVector(long n, complex<double> *v1,double *v2)
{
	complex<double> *p1=v1;
	double *p2=v2;

	for(long i=0; i<n; i++)
	{
		(*p1)*=*p2; p1++; p2++;
	}
}
int cinv(complex<double> **a,long n,double eps)
{ 
	int *is,*js,i,j,k,l,u,v,w;
    double p,q,s,t,d,b;
    is=new int [n];
    js=new int [n];

	double **ar,**ai;
	ar=new double *[n];
	ai=new double *[n];

	for(i=0;i<n;i++)
	{
		ar[i]=new double[n];
		ai[i]=new double[n];
	}
	for(i=0;i<n;i++)
	for(j=0;j<n;j++)
	{
		//ar[i][j]=ai[i][j]=0.0;
		ar[i][j]=a[i][j].real();
		ai[i][j]=a[i][j].imag();
	}

    for (k=0; k<=n-1; k++)
	{ d=0.0;
        for (i=k; i<=n-1; i++)
        for (j=k; j<=n-1; j++)
          { 
            p=ar[i][j]*ar[i][j]+ai[i][j]*ai[i][j];
            if (p>d) { d=p; is[k]=i; js[k]=j;}
          }
        if (d<=eps)
		{
			delete[](is); delete[] (js); spdlog::info("err**not inv\n");
            return(0);
          }
        if (is[k]!=k)
          for (j=0; j<=n-1; j++)
            { 
              t=ar[k][j]; ar[k][j]=ar[is[k]][j]; ar[is[k]][j]=t;
              t=ai[k][j]; ai[k][j]=ai[is[k]][j]; ai[is[k]][j]=t;
            }
        if (js[k]!=k)
          for (i=0; i<=n-1; i++)
            { 
              t=ar[i][k]; ar[i][k]=ar[i][js[k]]; ar[i][js[k]]=t;
              t=ai[i][k]; ai[i][k]=ai[i][js[k]]; ai[i][js[k]]=t;
            }
        l=k*n+k;
        ar[k][k]=ar[k][k]/d; ai[k][k]=-ai[k][k]/d;
        for (j=0; j<=n-1; j++)
          if (j!=k)
            { 
              p=ar[k][j]*ar[k][k]; q=ai[k][j]*ai[k][k];
              s=(ar[k][j]+ai[k][j])*(ar[k][k]+ai[k][k]);
              ar[k][j]=p-q; ai[k][j]=s-p-q;
            }
        for (i=0; i<=n-1; i++)
          if (i!=k)
            { v=i*n+k;
              for (j=0; j<=n-1; j++)
                if (j!=k)
                  { 
                    p=ar[k][j]*ar[i][k]; q=ai[k][j]*ai[i][k];
                    s=(ar[k][j]+ai[k][j])*(ar[i][k]+ai[i][k]);
                    t=p-q; b=s-p-q;
                    ar[i][j]=ar[i][j]-t;
                    ai[i][j]=ai[i][j]-b;
                  }
            }
        for (i=0; i<=n-1; i++)
          if (i!=k)
            { 
              p=ar[i][k]*ar[k][k]; q=ai[i][k]*ai[k][k];
              s=(ar[i][k]+ai[i][k])*(ar[k][k]+ai[k][k]);
              ar[i][k]=q-p; ai[i][k]=p+q-s;
            }
      }
    for (k=n-1; k>=0; k--)
      { if (js[k]!=k)
          for (j=0; j<=n-1; j++)
            { 
              t=ar[k][j]; ar[k][j]=ar[js[k]][j]; ar[js[k]][j]=t;
              t=ai[k][j]; ai[k][j]=ai[js[k]][j]; ai[js[k]][j]=t;
            }
        if (is[k]!=k)
          for (i=0; i<=n-1; i++)
            { 
              t=ar[i][k]; ar[i][k]=ar[i][is[k]]; ar[i][is[k]]=t;
              t=ai[i][k]; ai[i][k]=ai[i][is[k]]; ai[i][is[k]]=t;
            }
      }
	for(i=0;i<n;i++)
	for(j=0;j<n;j++)
	{
		a[i][j].real(ar[i][j]);
		a[i][j].imag(ai[i][j]);
	}

	for(k=0;k<n;k++) delete[] ar[k];
	delete[] ar;
	for(k=0;k<n;k++) delete[] ai[k];
	delete[] ai;

    delete []is; delete [](js);
    return 0;
  }
  int cinv1(complex<double>** a,long n,double eps)
{
	double **ar,**ai;
	ar=new double *[n];
	ai=new double *[n];

	int i,j,k;
	for(i=0;i<n;i++)
	{
		ar[i]=new double[n];
		ai[i]=new double[n];
	}
	for(i=0;i<n;i++)
	for(j=0;j<n;j++)
	{
		ar[i][j]=ai[i][j]=0.0;
		ar[i][j]=a[i][j].real();
		ai[i][j]=a[i][j].imag();
	}
	/*for (i=0; i<n; i++)
    for (j=0; j<n; j++)
     { 
		 br[i][j]=ar[i][j]; 
		 bi[i][j]=ai[i][j];
	}*/
   // k=cinv(n,ar,ai,eps);

	return k;
}
  complex<double> getScatterExactResult(double costh,double ka,double kr)
{
	int n=20,k;
	double *delta,*P;
    complex<double> *h2,*dh2;
	double tol=0.000001;

	complex<double> i=complex<double> (0.0,1.0);

	delta=new double[n+1];
	h2=new complex<double>[n+2];
	dh2=new complex<double>[n+1];
	P=new double[n+1];

	h2[0]=-1.0*exp(-i*ka)/(ka*i);
	h2[1]=-1.0*exp(-i*ka)*(1.0+1.0/(i*ka))/(ka);

	for(k=0;k<=n;k++)
	{    
		if(k>=1)
		{
			h2[k+1]=(2.0*k+1.0)*h2[k]/(ka)-h2[k-1];
		}
		dh2[k]=k*1.0*h2[k]/ka-h2[k+1];
		delta[k]=atan(dh2[k].real()/dh2[k].imag());
		if(abs(delta[k])<tol)delta[k]=0.0;
	}	

	h2[0]=-1.0*exp(-i*kr)/(kr*i);
	h2[1]=-1.0*exp(-i*kr)*(1.0+1.0/(i*kr))/(kr);
	for(k=1;k<=n-1;k++)
		h2[k+1]=(2.0*k+1.0)*h2[k]/(kr)-h2[k-1];

	P[0]=1.0;
	P[1]=costh;
	for(k=2;k<=n;k++)
	    P[k]=(costh*(2.0*k-1.0)*P[k-1]-(1.0*k-1.0)*P[k-2])/k;

    complex<double> temp=0.0,Sum=0.0;
	for(k=0;k<=n;k++)
	{
		temp=1.0*pow(i,k+1);
		temp*=((2.0*k+1.0)*sin(delta[k])*exp(i*delta[k])*h2[k]*P[k]);

		Sum+=temp;
	}

	delete [] P;	    
	delete [] h2;
	delete [] delta;	delete [] dh2;

	return Sum;
}
complex<double> Eix(complex<double> a,int n)//pow(double x,double y)表示x^y
{
	int i;
    complex<double> sum = 1.0;
	if(n>=0)
    for(i=0; i <=n; i++)
    {
       if(i == 0)
           return 1;
       sum *= a;
    }
	else
	{
		for(i=1;i<=-n;i++)
		{
			sum*=a;
		}
		sum=1.0/sum;
	}
    return sum;
}
double Wigner3jsymbol(int j1,int j2,int j3,int m1,int m2,int m3,int i)
{
	int k,n;
	double tmp1,tmp2,tmp3,tmp;
	tmp1=tmp2=tmp3=tmp=0.0;

	if(j1==j2&&m1+m2==0&&j3==0&&m3==0)
	{
		tmp1=pow(-1.0,j1-m1);
		tmp=1.0/sqrt(2.0*j1+1);
		tmp*=tmp1;
	}
	else
	{
		tmp1=Factorial(j1+j2-j3);
	    tmp=Factorial(j1-j2+j3);
	    tmp1*=tmp;
	    tmp=Factorial(-j1+j2+j3);
	    tmp1*=tmp;
	    tmp=Factorial(j1+j2+j3+1);
	    tmp1=tmp1/tmp;
	    tmp1=sqrt(tmp1);

	    tmp2=Factorial(j1+m1);
	    tmp=Factorial(j1-m1);
		tmp2*=tmp;
	    tmp=Factorial(j2+m2);
		tmp2*=tmp;
	    tmp=Factorial(j2-m2);
	    tmp2*=tmp;	 
		tmp=Factorial(j3+m3);
	    tmp2*=tmp;
	    tmp=Factorial(j3-m3);
	    tmp2*=tmp;
	    tmp2=sqrt(tmp2);

        tmp3=1.0;
        double temp,temp1,temp2,temp3,temp4,temp5;
        temp1=temp2=temp3=temp4=temp5=temp=0.0;

	    int min=dmin(j1+j2-j3,j1-m1,j2+m2);
	    int max=dmax(j2-j3-m1,j1-j3+m2,0);

        for(n=max;n<=min;n++)
       {
		  tmp3=Factorial(n);
	      temp1=Factorial(j1+j2-j3-n);
	      temp2=Factorial(j1-m1-n);
	      temp3=Factorial(j2+m2-n);
	      temp4=Factorial(j3-j2+m1+n);
	      temp5=Factorial(j3-j1-m2+n);

	      tmp3=tmp3*temp1*temp2*temp3*temp4*temp5;
	      tmp=pow(-1.0,n+j1-j2-m3);
	      tmp3=tmp/tmp3;
	      temp+=tmp3;
       }
		//tmp=sqrt(2.0*j3+1);
       tmp=tmp1*tmp2*temp;
       i=min;
	}
    return tmp;
}
double Factorial(int n)//负数的阶乘是无穷大
{
	int i;
	double tmp;
    tmp=1.0;
	if(n==0) tmp=1.0;
	else if(n>0)
	{
		for(i=1;i<n+1;i++)
		         tmp=tmp*i;
	}
	else
	{
		tmp=10000000000000000;//1e16
	}

	return tmp;
}

complex<double> TheResultofWigner3j(int i,int n,int j,int m,int l)
{
	double tmp1,tmp2;
	int k,p;
    k=0;
	complex<double> w,tmp;
	complex<double> a=complex<double>(0.0,1.0);;
	

	if(abs(m+j)>l)tmp1=0.0;
	else
	{
		tmp1=Wigner3jsymbol(i,n,l,j,m,-j-m,k);
	}
	tmp2=Wigner3jsymbol0(i,n,l);
	tmp=pow(a,n-i+l);
	w=(2.0*l+1.0)*tmp*tmp1*tmp2;
	return w;
}
double Wigner3jsymbol0(int j1,int j2,int j3)
{
	int k,n;
	double tmp1,tmp2,tmp3,tmp;
	tmp1=tmp2=tmp3=tmp=0.0;

	k=j1+j2+j3;
	
	if(k%2)tmp=0.0;
	else
	{
		tmp1=Factorial(j1+j2-j3);
	    tmp=Factorial(j1-j2+j3);
	    tmp1*=tmp;
	    tmp=Factorial(-j1+j2+j3);
	    tmp1*=tmp;
	    tmp=Factorial(k+1);
	    tmp1=tmp1/tmp;
	    tmp1=sqrt(tmp1);
   
	    tmp2=Factorial(k/2-j1);
	    tmp=Factorial(k/2-j2);
	    tmp2*=tmp;
	    tmp=Factorial(k/2-j3);
	    tmp2*=tmp;
	    tmp=Factorial(k/2);
	    tmp2=tmp/tmp2;

	    tmp=pow(-1.0,k/2);
	    tmp*=tmp1*tmp2;
	}
	return tmp;
}


////////////////////////////////////////////////////
double dmax(double a, double b)
{
	return a>b ? a:b;
}
double dmin(double a, double b)
{
	return a<b ? a:b;
}
double dmax(double a, double b, double c)
{
	double max=a;
	if(b>max) max=b;	if(c>max) max=c;
	return max;
}
double dmin(double a, double b, double c)
{
	double min=a;
	if(b<min) min=b;	if(c<min) min=c;
	return min;
}
double dmax(double a, double b, double c, double d)
{
	double max=a;
	if(b>max) max=b;	if(c>max) max=c;	if(d>max) max=d;
	return max;
}
double dmin(double a, double b, double c, double d)
{
	double min=a;
	if(b<min) min=b;	if(c<min) min=c;	if(d<min) min=d;
	return min;
}
double sign(double a)
{
	return a>0.0 ? 1.0 : -1.0;
}
int imax(int a, int b)
{
	return a>b ? a:b;
}
int imin(int a, int b)
{
	return a<b ? a:b;
}
int imax(int a, int b, int c)
{
	int max=a;
	if(b>max) max=b;	if(c>max) max=c;
	return max;
}
int imin(int a, int b, int c)
{
	int min=a;
	if(b<min) min=b;	if(c<min) min=c;
	return min;
}
int imax(int a, int b, int c, int d)
{
	int max=a;
	if(b>max) max=b;	if(c>max) max=c;	if(d>max) max=d;
	return max;
}
int imin(int a, int b, int c, int d)
{
	int min=a;
	if(b<min) min=b;	if(c<min) min=c;	if(d<min) min=d;
	return min;
}

//////////////////////////////////////////////////////
void matrixXvector(int nRow,int nCol, double** a, double* v, double* w)
{
	for(int i=0; i<nRow; i++)
	{
		w[i]=0.0;
		for(int j=0; j<nCol; j++)
			w[i]+=a[i][j]*v[j];
	}
}
void matrixAddXvector(int nRow,int nCol, double** a, double* v, double* w)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
			w[i]+=a[i][j]*v[j];
	}
}
void matrixMinusXvector(int nRow,int nCol, double** a, double* v, double* w)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
			w[i]-=a[i][j]*v[j];
	}
}
void matrixScaleAddXvector(int nRow,int nCol, double alpha, double** a, double* v, double* w)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
			w[i]+=a[i][j]*v[j]*alpha;
	}
}
void transMatrixXvector(int nRow,int nCol, double** a, double* v, double* w)
{
	for(int i=0; i<nCol; i++)
	{
		w[i]=0.0;
		for(int j=0; j<nRow; j++)
			w[i]+=v[j]*a[j][i];
	}
}
void transMatrixAddXvector(int nRow,int nCol, double** a, double* v, double* w)
{
	for(int i=0; i<nCol; i++)
	{
		for(int j=0; j<nRow; j++)
			w[i]+=v[j]*a[j][i];
	}
}
void transMatrixMinusXvector(int nRow,int nCol, double** a, double* v, double* w)
{
	for(int i=0; i<nCol; i++)
	{
		for(int j=0; j<nRow; j++)
			w[i]-=v[j]*a[j][i];
	}
}
void transMatrixScaleAddXvector(int nRow,int nCol, double alpha, double** a, double* v, double* w)
{
	for(int i=0; i<nCol; i++)
	{
		for(int j=0; j<nRow; j++)
			w[i]+=v[j]*a[j][i]*alpha;
	}
}
void scaleXmatrix(int nRow,int nCol, double** a, double alpha)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
			a[i][j]*=alpha;
	}
}
// c = a + b
void matrixAddmatrix(int nRow,int nCol, double** c, double** a, double** b)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
			c[i][j]=a[i][j]+b[i][j];
	}
}
// a += b
void matrixAddTomatrix(int nRow,int nCol, double** a, double** b)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
			a[i][j]+=b[i][j];
	}
}
/* c = a * b */
void matrixXmatrix(int nRow, int nCol, int K, double** c, double** a, double** b)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
		{
			c[i][j]=0.0;
			for(int k=0; k<K; k++)
				c[i][j]+=a[i][k]*b[k][j];
		}
	}
}
void float_matrixXmatrix(int nRow, int nCol, int K, double** c, double** a, double** b)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
		{
			c[i][j]=0.0;
			for(int k=0; k<K; k++)
				c[i][j]+=a[i][k]*b[k][j];
		}
	}
}
void matrixAddXmatrix(int nRow, int nCol, int K, double** c, double** a, double** b)
{
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
		{
			for(int k=0; k<K; k++)
				c[i][j]+=a[i][k]*b[k][j];
		}
	}
}
double dnorm2(int n, double* v)
{
	double norm=0;
	for(int i=0; i<n; i++) norm+=v[i]*v[i];
	return sqrt(norm);
}

void initVector(long n, double *v, int incx, double alfa)
{
	double *p=v;
	for(long i=0; i<n; i++)
	{
		*p=alfa; p+=incx;
	}
}
void initVector0(long n, double *v, int incx)
{
	double *p=v;
	for(long i=0; i<n; i++)
	{
		*p=0.0; p+=incx;
	}
}
void initVector0(long n, double *v)
{
	double *p=v;
	for(long i=0; i<n; i++)
	{
		*p=0.0; p++;
	}
}
//v1[i]*=v2[i]
void VectorXVector(long n, double *v1, double *v2)
{
	double *p1=v1;
	double *p2=v2;

	for(long i=0; i<n; i++)
	{
		(*p1)*=*p2; p1++; p2++;
	}
}

///////////////////////////////////////////////
int mxInvert(double a[],int n)
{
	int *is,*js,i,j,k,l,u,v;
    double d,p;
    is=new int[n];
    js=new int[n];
    for (k=0; k<=n-1; k++)
    { 		
		d=0.0;
        for (i=k; i<=n-1; i++)
        for (j=k; j<=n-1; j++)
        {
			l=i*n+j; p=fabs(a[l]);
            if (p>d) { d=p; is[k]=i; js[k]=j;}
        }
        if (d+1.0==1.0)
        { 
			delete []is; delete []js;
            return(1);
        }
        if (is[k]!=k)
          for (j=0; j<=n-1; j++)
          { 
			  u=k*n+j; v=is[k]*n+j;
              p=a[u]; a[u]=a[v]; a[v]=p;
          }
        if (js[k]!=k)
          for (i=0; i<=n-1; i++)
          { 
			  u=i*n+k; v=i*n+js[k];
              p=a[u]; a[u]=a[v]; a[v]=p;
          }
        l=k*n+k;
        a[l]=1.0/a[l];
        for (j=0; j<=n-1; j++)
          if (j!=k) { u=k*n+j; a[u]=a[u]*a[l];}
        for (i=0; i<=n-1; i++)
          if (i!=k)
            for (j=0; j<=n-1; j++)
              if (j!=k)
              { 
				  u=i*n+j;
                  a[u]=a[u]-a[i*n+k]*a[k*n+j];
              }
        for (i=0; i<=n-1; i++)
          if (i!=k) { u=i*n+k; a[u]=-a[u]*a[l];}
    }
    for (k=n-1; k>=0; k--)
    { 
		if (js[k]!=k)
          for (j=0; j<=n-1; j++)
          { 
			  u=k*n+j; v=js[k]*n+j;
              p=a[u]; a[u]=a[v]; a[v]=p;
          }
        if (is[k]!=k)
          for (i=0; i<=n-1; i++)
          { 
			  u=i*n+k; v=i*n+is[k];
              p=a[u]; a[u]=a[v]; a[v]=p;
          }
    }
    delete []is; delete []js;
    return(0);
}
//////////////////////////////////////////
int eqGauss(double a[],double b[],int n)
{
	int *js,l,k,i,j,is=0,p,q;
    double d,t;

    js=new int[n];
    l=1;
    for (k=0;k<=n-2;k++)
    { 
		d=0.0;
        for (i=k;i<=n-1;i++)
          for (j=k;j<=n-1;j++)
          { 
			  t=fabs(a[i*n+j]);
              if (t>d) { d=t; js[k]=j; is=i;}
          }
        if (d+1.0==1.0) l=0;
        else
        { 
			if (js[k]!=k)
              for (i=0;i<=n-1;i++)
              { 
				  p=i*n+k; q=i*n+js[k];
                  t=a[p]; a[p]=a[q]; a[q]=t;
              }
            if (is!=k)
            { 
				for (j=k;j<=n-1;j++)
                { 
					p=k*n+j; q=is*n+j;
                    t=a[p]; a[p]=a[q]; a[q]=t;
                }
                t=b[k]; b[k]=b[is]; b[is]=t;
            }
        }
        if (l==0)
        { 
			delete []js; 
            return(1);
        }
        d=a[k*n+k];
        for (j=k+1;j<=n-1;j++) { p=k*n+j; a[p]=a[p]/d;}
        b[k]=b[k]/d;
        for (i=k+1;i<=n-1;i++)
        { 
			for (j=k+1;j<=n-1;j++)
            { 
				p=i*n+j;
                a[p]=a[p]-a[i*n+k]*a[k*n+j];
            }
            b[i]=b[i]-a[i*n+k]*b[k];
        }
    }
    d=a[(n-1)*n+n-1];
    if (fabs(d)+1.0==1.0)
    {  
		delete []js; 
        return(1);
    }
    b[n-1]=b[n-1]/d;
    for (i=n-2;i>=0;i--)
    { 
		t=0.0;
        for (j=i+1;j<=n-1;j++) t=t+a[i*n+j]*b[j];
        b[i]=b[i]-t;
    }
    js[n-1]=n-1;
    for (k=n-1;k>=0;k--)
      if (js[k]!=k) { t=b[k]; b[k]=b[js[k]]; b[js[k]]=t;}
    delete []js;
    return(0);
}
//////////////////////////////////////////
int eqGausJdn(double a[],double b[],int n, int m)
{
	int *js,l,k,i,j,is=0,p,q;
    double d,t;
    js=new int[n];;
    l=1;
    for (k=0;k<=n-1;k++)
    { 
		d=0.0;
        for (i=k;i<=n-1;i++)
        for (j=k;j<=n-1;j++)
        { 
			t=fabs(a[i*n+j]);
            if (t>d) { d=t;  js[k]=j;  is=i;}
        }
        if (d+1.0==1.0) l=0;
        else
        { 
			if (js[k]!=k)
			  for (i=0;i<=n-1;i++)
			  { 
				  p=i*n+k; q=i*n+js[k];
				  t=a[p]; a[p]=a[q]; a[q]=t;
			  }
            if (is!=k)
            { 
				for (j=k;j<=n-1;j++)
                { 
					p=k*n+j; q=is*n+j;
                    t=a[p]; a[p]=a[q]; a[q]=t;
                }
                for (j=0;j<=m-1;j++)
                { 
					p=k*m+j; q=is*m+j;
                    t=b[p]; b[p]=b[q]; b[q]=t;
                }
            }
        }
        if (l==0)
        { 
			delete []js; 
            return(1);
        }
        d=a[k*n+k];
        for (j=k+1;j<=n-1;j++) { p=k*n+j; a[p]=a[p]/d;}
        for (j=0;j<=m-1;j++) { p=k*m+j; b[p]=b[p]/d;}
        for (j=k+1;j<=n-1;j++)
        for (i=0;i<=n-1;i++)
        { 
			p=i*n+j;
            if (i!=k) a[p]=a[p]-a[i*n+k]*a[k*n+j];
        }
        for (j=0;j<=m-1;j++)
        for (i=0;i<=n-1;i++)
        { 
			p=i*m+j;
            if (i!=k) b[p]=b[p]-a[i*n+k]*b[k*m+j];
        }
    }
    for (k=n-1;k>=0;k--)
      if (js[k]!=k)
        for (j=0;j<=m-1;j++)
        { 
			p=k*m+j; q=js[k]*m+j;
            t=b[p]; b[p]=b[q]; b[q]=t;
        }
    delete []js;
    return(0);
}

/*
template <class T>
Stack<T>::Stack(int s)
{
	size=s;
	top=-1;
	stackPtr=new T[size];
}
*/
//////////////////////////////////////////////globals
//////////////////////////////////////////////globals
void dswap(double* x,double* y)
{
	double temp;
	temp=*x;
	*x=*y;
	*y=temp;
}

int inv0(double** a,long n,double eps)//length n
{
  long k,i,j,*is,*js;
  double t,d;
  is=new long[n];
  js=new long[n];
  
  if(is==0x0||js==0x0)
  {
	  delete[] is;
	  delete[] js;
	  return 1;
  }

  for(k=1;k<=n;k++)
  {
	  d=0.0;
	  for(i=k;i<=n;i++)
	  {
		  for(j=k;j<=n;j++)
		  {
			  t=sqrt(a[i-1][j-1]*a[i-1][j-1]);
			  if(t>d)
			  {
				  d=t;is[k-1]=i-1;js[k-1]=j-1;
			  }
		  }
	  }

	  if(d<=eps)
	  {
		  delete[] is;
		  delete[] js;
		  return 1;
	  }

	  for(j=1;j<=n;j++)
	  {
		  i=is[k-1];
		  dswap(&a[k-1][j-1],&a[i][j-1]);
	  }

	  for(i=1;i<=n;i++)
	  {
		  j=js[k-1];
		  dswap(&a[i-1][k-1],&a[i-1][j]);
	  }

	  a[k-1][k-1]=1.0/a[k-1][k-1];

	  for(j=1;j<=n;j++)
	  {
		  if(j!=k)
			  a[k-1][j-1]*=a[k-1][k-1];
	  }

	  for(i=1;i<=n;i++)
	  {
		  if(i!=k)
		  {
			  for(j=1;j<=n;j++)
			  {
				  if(j!=k)
					  a[i-1][j-1]-=a[i-1][k-1]*a[k-1][j-1];
			  }
		  }
	  }//end for i

	  for(i=1;i<=n;i++)
	  {
		  if(i!=k)
			  a[i-1][k-1]*=-1.0*a[k-1][k-1];
	  }//end for i

  }//end for k

  for(k=n;k>=1;k--)
  {
	  for(j=1;j<=n;j++)
	  {
		  i=js[k-1];
		  dswap(&a[k-1][j-1],&a[i][j-1]);
	  }//end for j

	  for(i=1;i<=n;i++)
	  {
		  j=is[k-1];
		  dswap(&a[i-1][k-1],&a[i-1][j]);
	  }//end for i
	  
  }//end for k

  delete[] is;
  delete[] js;
  return 0;
}

int inv1(double** a,long n,double eps)//length n+1
{
  long k,i,j,*is,*js;
  double t,d;
  is=new long[n+1];
  js=new long[n+1];
  
  if(is==0x0||js==0x0)
  {
	  delete[] is;
	  delete[] js;
	  return -1;
  }

  for(k=1;k<=n;k++)
  {
	  d=0.0;
	  for(i=k;i<=n;i++)
	  {
		  for(j=k;j<=n;j++)
		  {
			  t=sqrt(a[i][j]*a[i][j]);
			  if(t>d)
			  {
				  d=t;is[k]=i;js[k]=j;
			  }
		  }//end for j
	  }  //end for i

	  if(d<=eps)
	  {
		  delete[] is;
		  delete[] js;
		  return -1;
	  }

	  for(j=1;j<=n;j++)
	  {
		  i=is[k];
		  dswap(&a[k][j],&a[i][j]);
	  }

	  for(i=1;i<=n;i++)
	  {
		  j=js[k];
		  dswap(&a[i][k],&a[i][j]);
	  }

	  a[k][k]=1.0/a[k][k];

	  for(j=1;j<=n;j++)
	  {
		  if(j!=k)
			  a[k][j]*=a[k][k];
	  }

	  for(i=1;i<=n;i++)
	  {
		  if(i!=k)
		  {
			  for(j=1;j<=n;j++)
			  {
				  if(j!=k)
					  a[i][j]-=a[i][k]*a[k][j];
			  }
		  }
	  }//end for i

	  for(i=1;i<=n;i++)
	  {
		  if(i!=k)
			  a[i][k]*=-1.0*a[k][k];
	  }//end for i

  }//end for k

  for(k=n;k>=1;k--)
  {
	  for(j=1;j<=n;j++)
	  {
		  i=js[k];
		  dswap(&a[k][j],&a[i][j]);
	  }//end for j

	  for(i=1;i<=n;i++)
	  {
		  j=is[k];
		  dswap(&a[i][k],&a[i][j]);
	  }//end for i
	  
  }//end for k

  delete[] is;
  delete[] js;
  return 1;

}

//求2范数
double norm1(double* a,long n)//length: n+1
{
	  double temp;
	  temp=0;
	  for(long i=1; i<=n; i++)
	  {
		  temp+=a[i]*a[i];
	  }
	  temp=sqrt(temp);
	  return temp;
}

double norm0(double* a,long n)//length: n
{
	  double temp;
	  temp=0;
	  for(long i=0; i<n; i++)
	  {
		  temp+=a[i]*a[i];
	  }
	  temp=sqrt(temp);
	  return temp;
}

//precondictioned GMRES迭代法
int gmres(long n, double* A, double* b, int restart, double tol, int maxIter,
		  double* M, double* x0, double* x, int& flag, int* iter)
{
	//A, b, M, x0, x:   are all from 0; whereas the locals from 1 with length of n+1

	//n:             方程阶数
	//A:             系数矩阵
	//b:             右手向量
	//restart:       内循环迭代次数
	//tol:           迭代允许误差
	//maxIter:       外循环次数
	//M:             预处理矩阵
	//x0:            初始向量
	//x:             结果向量
	//flag:          结果是否正确0，1，3
	//iter:          内外迭代次数    iter[0]:inner      iter[1]:outer
	long int i,j,k,l;          //used for loop
	double tolb;                        //tol*norm(b)
	double n2b;                         //b的2范数
	double* r=new double[n+1];          //r=b-A*x
	double* vh=new double[n+1];         //相当于r0
    double* u=new double[n+1];
	double* u2=new double[n+1];
	double* q=new double[restart+1];
	double normr;                       //r的2范数

	int stag;                           //是否中断
	double phibar;
	double rt,c,s,temp;

	double** V=new double*[n+1];        //Arnoldi 向量
	for(i=0;i<=n;i++) V[i]=new double[restart+2];
	double* h=new double[restart+2];    //上Hessenberg矩阵的一列 st A*V = V*H
	double** QT=new double*[restart+2]; //旋转正交矩阵 st QT*H = R
	for(i=0;i<=restart+1;i++) QT[i]=new double[restart+2];
	double** R=new double*[restart+1];  //旋转后上三角矩阵 st H = Q*R
	for(i=0;i<=restart;i++) R[i]=new double[restart+1];
	double* f=new double[restart+1];    //y = R\f => x = x0 + V*y
	double** W=new double*[n+1];        //W = V*inv(R)
	for(i=0;i<=n;i++) W[i]=new double[restart+1];

	n2b = norm0(b,n);                     
	if(n2b==0.0)
	{
		for(i=1;i<=n;i++) x[i-1]=0.0;
		flag=0;
		iter[0]=0;
		iter[1]=0;
		delete[] r; delete[] vh; delete[] u; delete[] u2; delete[] q; delete[] h; delete[] f;
	    for(k=0;k<=n;k++) delete[] V[k];
	    delete[] V;
	    for(k=0;k<=restart+1;k++) delete[] QT[k];
	    delete[] QT;
	    for(k=0;k<=restart;k++) delete[] R[k];
	    delete[] R;
	    for(k=0;k<=n;k++) delete[] W[k];
	    delete[] W;
		return 1;
	}//end if

	for(i=0;i<n;i++) x[i]=x0[i];

	flag=1;
	tolb=tol*n2b;

	for(i=1;i<=n;i++)
	{
		r[i]=b[i-1];
		for(j=1;j<=n;j++) r[i]-=A[(i-1)*n+j-1]*x[j-1];
	}//end for ii

	normr=norm1(r,n);
	if(normr<=tolb)
	{
		flag=0;
		iter[0]=0;
		iter[1]=0;
		delete[] r; delete[] vh; delete[] u; delete[] u2; delete[] q; delete[] h; delete[] f;
	    for(k=0;k<=n;k++) delete[] V[k];
	    delete[] V;
	    for(k=0;k<=restart+1;k++) delete[] QT[k];
	    delete[] QT;
	    for(k=0;k<=restart;k++) delete[] R[k];
	    delete[] R;
	    for(k=0;k<=n;k++) delete[] W[k];
	    delete[] W;
		return 1;
	}//end if


	stag=0;
	for(i=1;i<=maxIter;i++)  //外循环
	{
		for(j=1;j<=restart+1;j++)
			for(k=1;k<=restart+1;k++) QT[j][k]=0.0;

		for(j=1;j<=restart;j++)
			for(k=1;k<=restart;k++) R[j][k]=0.0;

		for(j=1;j<=n;j++)
		{
			vh[j]=0.0;
			for(k=1;k<=n;k++) vh[j]+=M[(j-1)*n+k-1]*r[k];
		}//vh = M \ r

		h[1]=norm1(vh,n);
		for(j=1;j<=n;j++) V[j][1]=vh[j]/h[1];  //V(:,1) = vh / h(1)
     
		QT[1][1]=1.0;
		phibar=h[1];

		for(j=1;j<=restart;j++)        //内循环
		{
			for(k=1;k<=n;k++)
			{
				u2[k]=0.0;
				for(l=1;l<=n;l++) u2[k]+=A[(k-1)*n+l-1]*V[l][j];
			}//u2=A*v(:,j)

			for(k=1;k<=n;k++)
			{
				u[k]=0.0;
				for(l=1;l<=n;l++) u[k]+=M[(k-1)*n+l-1]*u2[l];
			}//u=M*u2

			for(k=1;k<=j;k++)
			{
				h[k]=0.0;
				for(l=1;l<=n;l++) h[k]+=V[l][k]*u[l];//h(k) = V(:,k)' * u;
				for(l=1;l<=n;l++) u[l]-=h[k]*V[l][k];//u = u - h(k) * V(:,k)
			}

			h[j+1]=norm1(u,n);

			for(k=1;k<=n;k++) V[k][j+1]=u[k]/h[j+1];//V(:,j+1) = u / h(j+1);

			for(k=1;k<=j;k++)
			{
				R[k][j]=0.0;
				for(l=1;l<=j;l++) R[k][j]+=QT[k][l]*h[l];
			}//R(1:j,j) = QT(1:j,1:j) * h(1:j)

			rt=R[j][j];

	        if(h[j+1]==0.0)
			{
				c=1.0;
				s=0.0;
			}
			else if(fabs(h[j+1])>fabs(rt))
			{
				temp=rt/h[j+1];
                s=1.0/sqrt(1.0+fabs(temp)*fabs(temp));
                c=-temp*s;
			}
			else
			{
				temp=h[j+1]/ rt;
                c=1.0/sqrt(1.0+fabs(temp)*fabs(temp));
                s=-temp*c;
			}

			R[j][j]=c*rt-s*h[j+1];
			for(k=1;k<=j;k++) q[k]=QT[j][k];
			for(k=1;k<=j;k++) QT[j][k]=c*q[k];
			for(k=1;k<=j;k++) QT[j+1][k]=s*q[k];

			QT[j][j+1]=-s;
            QT[j+1][j+1]=c;
            f[j]=c*phibar;
            phibar=s*phibar;

	        if(j<restart-1)
			{
		        if(f[j]==0.0) stag=1;
		        if(j==1)
				{
					for(k=1;k<=n;k++) W[k][j]=V[k][j]/R[j][j];
				}
				else
				{
					for(k=1;k<=n;k++)
					{
						W[k][j]=V[k][j]/R[j][j];
						for(l=1;l<=j-1;l++) W[k][j]-=W[k][l]*R[l][j]/R[j][j];
					}
				}
				for(k=1;k<=n;k++) x[k-1]+=f[j]*W[k][j];
			}
			else
			{
				inv1(R,j,1e-15);
				for(k=1;k<=j;k++)
				{
					q[k]=0.0;
					for(l=1;l<=j;l++) q[k]+=R[k][l]*f[l];
				}
				for(k=1;k<=n;k++)
				{
					x[k-1]=0.0;
					for(l=1;l<=j;l++) x[k-1]+=V[k][l]*q[l];
					x[k-1]+=x0[k-1];
				}
			}

			for(k=1;k<=n;k++)
			{
				vh[k]=b[k-1];
				for(l=1;l<=n;l++) vh[k]-=A[(k-1)*n+l-1]*x[l-1];
			}

			normr=norm1(vh,n);
	        if(normr<=tolb)
			{
				flag = 0;
                iter[0]=i;
				iter[1]=j;
				delete[] r; delete[] vh; delete[] u; delete[] u2; delete[] q; delete[] h; delete[] f;
	            for(k=0;k<=n;k++) delete[] V[k];
	            delete[] V;
	            for(k=0;k<=restart+1;k++) delete[] QT[k];
	            delete[] QT;
	            for(k=0;k<=restart;k++) delete[] R[k];
	            delete[] R;
	            for(k=0;k<=n;k++) delete[] W[k];
	            delete[] W;
                return 1;
			}

	        if(stag==1)
			{
				flag=3;
				iter[0]=i;
				iter[1]=j;
				delete[] r; delete[] vh; delete[] u; delete[] u2; delete[] q; delete[] h; delete[] f;
	            for(k=0;k<=n;k++) delete[] V[k];
	            delete[] V;
	            for(k=0;k<=restart+1;k++) delete[] QT[k];
	            delete[] QT;
	            for(k=0;k<=restart;k++) delete[] R[k];
	            delete[] R;
	            for(k=0;k<=n;k++) delete[] W[k];
	            delete[] W;
				return -1;
			}

		}//end for j

        if(flag==1)
		{
			for(k=1;k<=n;k++) x0[k-1]=x[k-1];
			for(k=1;k<=n;k++)
			{
				r[k]=b[k-1];
				for(l=1;l<=n;l++) r[k]-=A[(k-1)*n+l-1]*x0[l-1];
			}
			iter[0]=i;
			iter[1]=j;
		}

	}//end for i

	delete[] r; delete[] vh; delete[] u; delete[] u2; delete[] q; delete[] h; delete[] f;
	for(i=0;i<=n;i++) delete[] V[i];
	delete[] V;
	for(i=0;i<=restart+1;i++) delete[] QT[i];
	delete[] QT;
	for(i=0;i<=restart;i++) delete[] R[i];
	delete[] R;
	for(i=0;i<=n;i++) delete[] W[i];
	delete[] W;

	return 1;
}

int gmres_solver(long n, double *A, double *b, double tolerance, int& flag, int* iter)
{
	long i, j;

	int restart=300;
	int maxIter=200;

	long nn=n*n;
	double *M=new double[nn];
	double *x0=new double[n];
	double *x=new double[n];

	for(i=0; i<nn; i++) M[i]=0.0;

	for(i=1; i<=n; i++)
	{
		x0[i-1]=1.0;
		for(j=1; j<=n; j++) M[(i-1)*n+i-1]+=A[(i-1)*n+j-1]*A[(i-1)*n+j-1];
		M[(i-1)*n+i-1]=1.0/sqrt(M[(i-1)*n+i-1]);
	}

	int ret=gmres(n, A, b, restart, tolerance, maxIter, M, x0, x, flag, iter);

	for(i=0; i<n; i++) b[i]=x[i];
	
	delete []M; delete []x0; delete []x;

	if(ret==-1) return 1;
	
	return 0;
}
///////////////////////////////////////////////////
int ludcmp(double**a,int n,int* indx, double eps)
{
	int i, imax=0, j, k;
	double big, dum, sum, temp;
	double *vv;

	vv=new double[n];

	for (i=0;i<n;i++)
	{
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big <= eps){ delete[] vv; return 1; }
		vv[i]=1.0/big;
	}
	for (j=0; j<n; j++)
	{
		for (i=0; i<j; i++)
		{
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j; i<n; i++)
		{
			sum=a[i][j];
			for (k=0; k<j; k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big)
			{
				big=dum;
				imax=i;
			}
		}
		if (j != imax)
		{
			for (k=0;k<n;k++)
			{
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (fabs(a[j][j]) <= eps) { delete[] vv; return 1; }///different from that in the book
		if (j != n-1)
		{
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
	delete[] vv;
	return 0;
}
/*this is decomposed to lower and upper cources; 
void lubksb(double **a, int n, int *indx, double* b)
{
	int i,ii=0,ip,j;
	double sum;

	for (i=0;i<n;i++)
	{
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
        if(ii)
			for (j=ii-1;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if(sum)
		{
			ii=i+1;
		}
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--)
	{
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
*/
void lubksb(double **a, int n, int *indx, double* b)//A*x=b
{
	lubksb_lower(a, n, indx, b);
	lubksb_upper(a, n, b);
}
void lubksb_lower(double **a, int n, int *indx, double* b)//L*x=b
{
	int i,ip,j;
	double sum;

	for (i=0;i<n;i++)
	{
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		for (j=0;j<=i-1;j++) sum -= a[i][j]*b[j];
		b[i]=sum;
	}
}
void lubksb_lower2(double **a, int n, int *indx, double* b)//x*L=b
{//to be tested
	int i, ip, j;
	double sum;

	for (i=n-1;i>=0;i--)
	{
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		for (j=i+1;j<n;j++) sum -= a[j][i]*b[j];
		b[i]=sum;
	}
}
void lubksb_upper(double **a, int n, double* b)//U*x=b
{
	int i,j;
	double sum;

	for (i=n-1;i>=0;i--)
	{
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
void lubksb_upper2(double **a, int n, double* b)//x*U=b
{
	int i,j;
	double sum;

	for (i=0;i<n;i++)
	{
		sum=b[i];
		for (j=0;j<i;j++) sum -= a[j][i]*b[j];
		b[i]=sum/a[i][i];
	}
}
////////////SVD/////////////////////////////////////
double dpythag(double a, double b)
{
	double absa,absb, ab;
	absa=fabs(a);  absb=fabs(b);
	if (absa > absb)
	{
		ab=absb/absa;
		return absa*sqrt(1.0+ab*ab);
	}
	else
	{
        if(absb == 0.0) return 0.0;
		else
		{
			ab=absa/absb;
			return absb*sqrt(1.0+ab*ab);
		}
	}
}
int dsvdcmp(double **a, int m, int n, double* w, double **v)
{
	int flag,i,its,j,jj,k,l=0,nm=0;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=new double[n];
	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i < m && i != n-1) {
			for (k=l;k<n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<m;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<n;k++) a[i][k] *= scale;
			}
		}
		anorm=dmax(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g) {
				for (j=l;j<n;j++) v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m-1,n-1);i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<m;j++) a[j][i] *= g;
		} else for (j=i;j<m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=dpythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30)
			{
				delete[] rv1;
				return 1;
			}//nrerror("no convergence in 30 dsvdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=dpythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=dpythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=dpythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	delete[] rv1;
	return 0;
}
void dsvbksb(double **u, double* w, double **v, int m, int n, double* b, double* x)
{
	int jj,j,i;
	double s,*tmp;

	tmp=new double[n];
	for (j=0;j<n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=0;i<m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=0;j<n;j++) {
		s=0.0;
		for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
	delete[] tmp;
}
int choleskydcmp(double **a, int n, double* p)//choldc(...)
{
	int i,j,k;
	double sum;

	for (i=0;i<n;i++) {
		for (j=i;j<n;j++) {
			for (sum=a[i][j],k=i-1;k>=0;k--) sum -= a[i][k]*a[j][k];
			if (i == j) {
				if (sum <= 0.0) return 1;
//					nrerror("choldc failed");
				p[i]=sqrt(sum);
			} else a[j][i]=sum/p[i];
		}
	}
	return 0;
}
void choleskybksb(double **a, int n, double* p, double* b, double* x)//cholsl(...)
{
	int i,k;
	double sum;

	for (i=0;i<n;i++) {
		for (sum=b[i],k=i-1;k>=0;k--) sum -= a[i][k]*x[k];
		x[i]=sum/p[i];
	}
	for (i=n-1;i>=0;i--) {
		for (sum=x[i],k=i+1;k<n;k++) sum -= a[k][i]*x[k];
		x[i]=sum/p[i];
	}
}

int qrdcmp(double **a, int n, double *c, double *d)
{
	int i,j,k;
	double scale,sigma,sum,tau;

	for (k=0;k<n-1;k++)
	{
		scale = 0.0;
		for (i=k;i<n;i++) scale=dmax(scale,fabs(a[i][k]));
		if (scale == 0.0) return 1;
		else
		{
			for (i=k;i<n;i++) a[i][k] /= scale;
			for (sum=0.0,i=k;i<n;i++) sum += SQR(a[i][k]);
			sigma=SIGN(sqrt(sum),a[k][k]);
			a[k][k] += sigma;
			c[k]=sigma*a[k][k];
			d[k] = -scale*sigma;
			for (j=k+1;j<n;j++)
			{
				for (sum=0.0,i=k;i<n;i++) sum += a[i][k]*a[i][j];
				tau=sum/c[k];
				for (i=k;i<n;i++) a[i][j] -= tau*a[i][k];
			}
		}
	}
	d[n-1]=a[n-1][n-1];
	if (d[n-1] == 0.0) return 1;
	else return 0;
}
void qrbksb(double **a, int n, double* c, double* d, double* b)//qrsolv
{
	int i,j;
	double sum,tau;

	for (j=0;j<n-1;j++)
	{//form Q^T*b
		for (sum=0.0,i=j;i<n;i++) sum += a[i][j]*b[i];
		tau=sum/c[j];
		for (i=j;i<n;i++) b[i] -= tau*a[i][j];
	}

	b[n-1] /= d[n-1];
	for (i=n-2;i>=0;i--)
	{
		for (sum=0.0,j=i+1;j<n;j++) sum += a[i][j]*b[j];
		b[i]=(b[i]-sum)/d[i];
	}
}

int XuShiliangQR(double** a,int m,double** q)
{
	int i,j,k;
	double u,alpha,w,t;

	for (i=0; i<m; i++)
	for (j=0; j<m; j++)
	{
		q[i][j]=0.0;
		if (i==j) q[i][j]=1.0;
	}

	for (k=0; k<m-1; k++)
	{ 
		u=0.0;
		for (i=k; i<m; i++)
		{
			w=fabs(a[i][k]);
			if (w>u) u=w;
		}
		alpha=0.0;
		for (i=k; i<m; i++)
		{
			t=a[i][k]/u;
			alpha=alpha+t*t;
		}
		if (a[k][k]>0.0) u=-u;

		alpha=u*sqrt(alpha);
		if (fabs(alpha)+1.0==1.0) return 1;

		u=sqrt(2.0*alpha*(alpha-a[k][k]));
		if ((u+1.0)!=1.0)
		{
			a[k][k]=(a[k][k]-alpha)/u;
			for (i=k+1; i<m; i++) a[i][k]=a[i][k]/u;

			for (j=0; j<m; j++)
			{
				t=0.0;
				for (i=k; i<m; i++) t+=a[i][k]*q[i][j];

				for (i=k; i<m; i++) q[i][j]-=2.0*t*a[i][k];
			}
			for (j=k+1; j<m; j++)
			{
				t=0.0;
				for (i=k; i<m; i++) t+=a[i][k]*a[i][j];

				for (i=k; i<m; i++)	a[i][j]-=2.0*t*a[i][k];
			}
			a[k][k]=alpha;
			for (i=k+1; i<m; i++) a[i][k]=0.0;
		}
	}

	for (i=0; i<m-1; i++)
	for (j=i+1; j<m;j++)
		dswap(&q[i][j], &q[j][i]);   

	return 0;
}

int XuShiliangQR2(int m, int n, double** a, double** q)//m>=n
{
	int i,j,k,l,nn, jj;
	double u,alpha,w,t;

	if (m<n) return 1;

	for (i=0; i<=m-1; i++)
	for (j=0; j<=m-1; j++)
	{
		q[i][j]=0.0;
		if (i==j) q[i][j]=1.0;
	}
	nn=n;
	if (m==n) nn=m-1;
	for (k=0; k<=nn-1; k++)
	{
		u=0.0; l=k*n+k;
		for (i=k; i<=m-1; i++)
		{
			w=fabs(a[i][k]);
			if (w>u) u=w;
		}
		alpha=0.0;
		for (i=k; i<=m-1; i++)
		{
			t=a[i][k]/u; alpha=alpha+t*t;
		}
		if (a[k][k]>0.0) u=-u;
		alpha=u*sqrt(alpha);
		if (fabs(alpha)+1.0==1.0) return 1;
		u=sqrt(2.0*alpha*(alpha-a[k][k]));
		if ((u+1.0)!=1.0)
		{
			a[k][k]=(a[k][k]-alpha)/u;
			for (i=k+1; i<=m-1; i++) a[i][k]=a[i][k]/u;
			for (j=0; j<=m-1; j++)
			{
				t=0.0;
				for (jj=k; jj<=m-1; jj++) t=t+a[jj][k]*q[jj][j];
				for (i=k; i<=m-1; i++) q[i][j]=q[i][j]-2.0*t*a[i][k];
			}
			for (j=k+1; j<=n-1; j++)
			{
				t=0.0;
				for (jj=k; jj<=m-1; jj++) t=t+a[jj][k]*a[jj][j];
				for (i=k; i<=m-1; i++) a[i][j]=a[i][j]-2.0*t*a[i][k];
			}
			a[k][k]=alpha;
			for (i=k+1; i<=m-1; i++) a[i][k]=0.0;
		}
	}
	for (i=0; i<=m-2; i++)
	for (j=i+1; j<=m-1;j++)
	{
		t=q[i][j]; q[i][j]=q[j][i]; q[j][i]=t;
	}

	return 0;
}

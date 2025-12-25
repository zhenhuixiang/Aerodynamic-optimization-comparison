/*
  CEC14 Test Function Suite 
  Jane Jing Liang (email: liangjing@zzu.edu.cn) 
  4th Dec. 2013
  1. Run the following command in Matlab window:
  mex cec14_func.cpp -DWINDOWS
  2. Then you can use the test functions as the following example:
  f = cec14_func(x,func_num); 
  Here x is a D*pop_size matrix.
*/
#include <WINDOWS.H>      
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <mex.h>

double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag,*SS;

#include <WINDOWS.H>      
#include <stdio.h>
#include <math.h>
#include <malloc.h>

#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

void sphere_func (double *, double *, int , double *,double *, int, int); /* Sphere */
void ellips_func(double *, double *, int , double *,double *, int, int); /* Ellipsoidal */
void step_func(double *, double *, int , double *,double *, int, int); /* Step */
void ackley_func (double *, double *, int , double *,double *, int, int); /* Ackley's */
void griewank_func (double *, double *, int , double *,double *, int, int); /* Griewank's  */
void rosenbrock_func (double *, double *, int , double *,double *, int, int); /* Rosenbrock's */
void rastrigin_func (double *, double *, int , double *,double *, int, int); /* Rastrigin's  */

void shiftfunc (double*,double*,int,double*);
void rotatefunc (double*,double*,int, double*);
void sr_func (double *, double *, int, double*, double*, double, int, int); /* shift and rotate */

void cec14_test_func(double *, double *,int,int,int);

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) 
{
	int  m, n,func_num;
	double  *f, *x;
	if ((nrhs < 2) || (nlhs < 1))
    {
		mexPrintf ("usage: f = cec14_func(x, func_num);\n");
		mexErrMsgTxt ("example: f= cec14_func([3.3253000e+000, -1.2835000e+000]', 1);");
    }
	n = mxGetM (prhs[0]);
	if (!(n==2||n==10||n==20||n==30))
    {
		mexPrintf ("usage: f = cec14_func(x, func_num);\n");
		mexErrMsgTxt ("Error: Test functions are only defined for D=2,10,20,30.");
    }
	m = mxGetN (prhs[0]);
	x = mxGetPr (prhs[0]);
	func_num= (int)*mxGetPr (prhs[1]);
	if (func_num>7)
    {
		mexPrintf ("usage: f = cec14_func(x, func_num);\n");
		mexErrMsgTxt ("Error: There are only 30 test functions in this test suite!");
    }

	plhs[0] = mxCreateDoubleMatrix (1, m, mxREAL);
	f = mxGetPr (plhs[0]);
	cec14_test_func(&x[0], &f[0], n,m,func_num);
}

/*
cec14_test_func

OShift
Matrix

*/

void cec14_test_func(double *x, double *f, int nx, int mx,int func_num)
{
	int cf_num=10,i;
	if (ini_flag==1)
	{
		if ((n_flag!=nx)||(func_flag!=func_num))
		{
			ini_flag=0;
		}
	}

	if (ini_flag==0)
	{
		FILE *fpt;
		char FileName[30];
		free(M);
		free(OShift);
		free(y);
		free(z);
		free(x_bound);
		y=(double *)malloc(sizeof(double)  *  nx);
		z=(double *)malloc(sizeof(double)  *  nx);
		x_bound=(double *)malloc(sizeof(double)  *  nx);
		for (i=0; i<nx; i++)
			x_bound[i]=100.0;

		if (!(nx==2||nx==10||nx==20||nx==30))
		{
			printf("\nError: Test functions are only defined for D=2,10,20,30.\n");
		}

		/* Load Matrix M*/
		sprintf(FileName, "input_data/M_%d_D%d.txt", func_num,nx);
		fpt = fopen(FileName,"r");
		if (fpt==NULL)
		{
		    printf("\n Error: Cannot open input file for reading \n");
		}
		
			M=(double*)malloc(nx*nx*sizeof(double));
			if (M==NULL)
				printf("\nError: there is insufficient memory available!\n");
			for (i=0; i<nx*nx; i++)
			{
				fscanf(fpt,"%Lf",&M[i]);
			}
		
		fclose(fpt);
		
		/* Load shift_data */
		sprintf(FileName, "input_data/shift_data_%d.txt", func_num);
		fpt = fopen(FileName,"r");
		if (fpt==NULL)
		{
			printf("\n Error: Cannot open input file for reading \n");
		}
		OShift=(double *)malloc(nx*cf_num*sizeof(double));
		if (OShift==NULL)
			printf("\nError: there is insufficient memory available!\n");
	
			for(i=0;i<nx;i++)
			{
				fscanf(fpt,"%Lf",&OShift[i]);
			}
	
		fclose(fpt);

		n_flag=nx;
		func_flag=func_num;
		ini_flag=1;
		//printf("Function has been initialized!\n");
	}


	for (i = 0; i < mx; i++)
	{
		switch(func_num)
		{
		case 1:	
			sphere_func(&x[i*nx],&f[i],nx,OShift,M,1,0);
			f[i]+=100.0;
			break;
		case 2:	
			ellips_func(&x[i*nx],&f[i],nx,OShift,M,1,0);
			f[i]+=200.0;
			break;
		case 3:	
			ellips_func(&x[i*nx],&f[i],nx,OShift,M,1,1);
			f[i]+=300.0;
			break;
		case 4:	
			step_func(&x[i*nx],&f[i],nx,OShift,M,1,0);
			f[i]+=400.0;
			break;
		case 5:	
			ackley_func(&x[i*nx],&f[i],nx,OShift,M,1,0);
			f[i]+=500.0;
			break;
		case 6:
			griewank_func(&x[i*nx],&f[i],nx,OShift,M,1,0);
			f[i]+=600.0;
			break;
		case 7:
			rosenbrock_func(&x[i*nx],&f[i],nx,OShift,M,1,0);
			f[i]+=700.0;
			break;
		case 8:	
			rastrigin_func(&x[i*nx],&f[i],nx,OShift,M,1,0);
			f[i]+=800.0;
			break;
		default:
			printf("\nError: There are only 28 test functions in this test suite!\n");
			f[i] = 0.0;
			break;
		}
		
	}

}

void sphere_func (double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag) /* Sphere */
{
	int i;
	f[0] = 0.0;
	sr_func (x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */	
	for (i=0; i<nx; i++)
	{					
		f[0] += z[i]*z[i];
	}

}



void ellips_func (double *x, double *f, int nx, double *Os,double *Mr, int s_flag, int r_flag) /* Ellipsoidal */
{
    int i;
	f[0] = 0.0;
	sr_func (x, z, nx, Os, Mr,1.0, s_flag, r_flag); /* shift and rotate */
	for (i=0; i<nx; i++)
	{
       f[0] += pow(10.0,6.0*i/(nx-1))*z[i]*z[i];
	}
}

void step_func (double *x, double *f, int nx, double *Os,double *Mr, int s_flag, int r_flag) /* Step */
{
    int i;
	
	sr_func (x, z, nx, Os, Mr,1.0, s_flag, r_flag); /* shift and rotate */

	f[0] = 0.0;
	for (i=0; i<nx; i++)
	{
       f[0] += pow(floor(z[i]+0.5),2.0);
	}
}

void rosenbrock_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Rosenbrock's */
{
    int i;
	double tmp1,tmp2;
	f[0] = 0.0;
	sr_func (x, z, nx, Os, Mr, 2.048/100.0, s_flag, r_flag); /* shift and rotate */
	z[0] += 1.0;//shift to orgin
	for (i=0; i<nx-1; i++)
	{
		z[i+1] += 1.0;//shift to orgin
		tmp1=z[i]*z[i]-z[i+1];
		tmp2=z[i]-1.0;
		f[0] += 100.0*tmp1*tmp1 +tmp2*tmp2;
	}
}

void ackley_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Ackley's  */
{
    int i;
    double sum1, sum2;
    sum1 = 0.0;
    sum2 = 0.0;

	sr_func (x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		sum1 += z[i]*z[i];
		sum2 += cos(2.0*PI*z[i]);
	}
	sum1 = -0.2*sqrt(sum1/nx);
	sum2 /= nx;
		f[0] =  E - 20.0*exp(sum1) - exp(sum2) +20.0;
}

void griewank_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Griewank's  */
{
    int i;
    double s, p;
    s = 0.0;
    p = 1.0;

	sr_func (x, z, nx, Os, Mr, 600.0/100.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		s += z[i]*z[i];
		p *= cos(z[i]/sqrt(1.0+i));
	}
	f[0] = 1.0 + s/4000.0 - p;
}

void rastrigin_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Rastrigin's  */
{
    int i;
	f[0] = 0.0;

	sr_func (x, z, nx, Os, Mr, 5.12/100.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		f[0] += (z[i]*z[i] - 10.0*cos(2.0*PI*z[i]) + 10.0);
	}
}


void shiftfunc (double *x, double *xshift, int nx,double *Os)
{
	int i;
    for (i=0; i<nx; i++)
    {
        xshift[i]=x[i]-Os[i];
    }
}

void rotatefunc (double *x, double *xrot, int nx,double *Mr)
{
	int i,j;
    for (i=0; i<nx; i++)
    {
        xrot[i]=0;
			for (j=0; j<nx; j++)
			{
				xrot[i]=xrot[i]+x[j]*Mr[i*nx+j];
			}
    }
}

void sr_func (double *x, double *sr_x, int nx, double *Os,double *Mr, double sh_rate, int s_flag,int r_flag) /* shift and rotate */
{
	int i;
	if (s_flag==1)
	{
		if (r_flag==1)
		{	
			shiftfunc(x, y, nx, Os);
			for (i=0; i<nx; i++)//shrink to the orginal search range
			{
				y[i]=y[i]*sh_rate;
			}
			rotatefunc(y, sr_x, nx, Mr);
		}
		else
		{
			shiftfunc(x, sr_x, nx, Os);
			for (i=0; i<nx; i++)//shrink to the orginal search range
			{
				sr_x[i]=sr_x[i]*sh_rate;
			}
		}
	}
	else
	{	

		if (r_flag==1)
		{	
			for (i=0; i<nx; i++)//shrink to the orginal search range
			{
				y[i]=x[i]*sh_rate;
			}
			rotatefunc(y, sr_x, nx, Mr);
		}
		else
		for (i=0; i<nx; i++)//shrink to the orginal search range
		{
			sr_x[i]=x[i]*sh_rate;
		}
	}
}


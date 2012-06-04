#include<stdio.h>
#include<math.h>
#include<R.h>
#include<Rmath.h>

//C code for kernel density estimation: the 3 routines refer to the use of 3 different kernel:
// gaussian, t_7

void c_kepdfopt(double *x, double *h, double *xeval, int *nx, int *ndim, int *neval, double *f, double *p)
{

 int i, j, k, _nx=*nx, _neval=*neval, _ndim=*ndim;
 double zero=0.0, s, p1=0.0, u=0.0, cons= pow( 0.3989423, _ndim )/_nx;

 for(i=0; i<_neval; i++)
 {
  s=zero;

   for(k=0; k<_nx; k++)
   {
    p1=zero;

    for(j=0; j<_ndim; j++)
    {
     //moltiplicare costa meno che dividere: in ingresso h=1/hm
     u  = ( x[k+j*_nx]-xeval[i+j*_neval] )*h[k+j*_nx]; 
     p1+=  u*u;
    }
    //il diviso 2 dentro exp() è stato portato fuori nel wrapper --> h = h/sqrt(2)
    //moltiplicare costa meno che dividere: in ingresso p=1/p
    s += exp(-p1)*p[k];
   }
    *f++=s*cons; 
 }
}




void c_kepdft7(double *x, double *h, double *xeval, int *nx, int
*ndim, int *neval, double *f, double *p)
{

 register int i, j, k;
 int _nx=*nx, _neval=*neval, _ndim=*ndim;
 double zero=0.0, uno=1.0, sette=7.0, s, prod1, prod, u, cons=pow(0.3849915, _ndim)/_nx, seven_raised;  seven_raised = pow(sette, (double)_ndim*4.0);

 for(i=0; i<_neval; i++)
 {
  s=zero;

   for(k=0; k<_nx; k++)
   {
	prod = uno;
	for(j=0; j<_ndim; j++)
{
  //moltiplicare costa meno che dividere: in ingresso h=1/hm
	u  = ( *(x + k + j*_nx) - *(xeval + i + j*_neval) )**(h + k + j*_nx);
    prod*=(sette+u*u);
}
	prod1 = prod*prod;
	s += (seven_raised/(prod1*prod1))**(p+k);
   }
    *f++=s*cons;
 }
}



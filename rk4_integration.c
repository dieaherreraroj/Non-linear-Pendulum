#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double w = 4.56;
double q = 0.98;
double u = 2.87;
double Fd = 0.95;

void copy_vect(int n, double *y, double *x);
void aux_interm(int n, int tag, double h, double *x, double *y, double *z);
void f(double t, double *y);
void update_kvect(int n, int tag, double h, double t, double *y, double *ki, double *kf);
void init_vector(int n, double *y, double x0, double v0);
void rk4_integ(int n, int NSTEP, double dt, double t0, double x0, double v0);

int main(int argc, char **argv){

  int n = 2;
  int NSTEP = 1000;

  double dt = 0.008;
  double t0 = 0.0;
  double x0 = 0.0;
  double v0 = 2.0;

  rk4_integ(n,NSTEP,dt,t0,x0,v0);
  
  return 0;
  
}

void copy_vect(int n, double *y, double *x){
  for(int ii = 0; ii<n; ii++)
    x[ii] = y[ii];
}

void aux_interm(int n, int tag, double h, double *x, double *y, double *z){
  if(tag == 1)
    copy_vect(n,y,z);
  if(tag == 2 || tag == 3)
    for(int ii = 0; ii<n; ii++)
      z[ii] = x[ii] + 0.5*h*y[ii];
  if(tag == 4)
    for(int ii = 0; ii<n; ii++)
      z[ii] = x[ii] + h*y[ii];
}

void f(double t, double *y){
  double a[2];
  a[0] = y[1];
  a[1] = -w*y[0];
  
  y[0] = a[0];
  y[1] = a[1];
}

void update_kvect(int n, int tag, double h, double t, double *y, double *ki, double *kf){
  double *aux = (double*) calloc(n,sizeof(double));
  aux_interm(n,tag,h,y,ki,aux);
  if(tag == 1)
    f(t,aux);
  if(tag == 2 || tag == 3)
    f(t+0.5*h,aux);
  if(tag == 4)
    f(t+h,aux);
  copy_vect(n,aux,kf);
  free(aux);
}

void init_vector(int n, double *y, double x0, double v0){
  y[0] = x0;
  y[1] = v0;
}

void rk4_integ(int n, int NSTEP, double dt, double t0, double x0, double v0){
  double *y = (double*) calloc(n,sizeof(double));
  double *k1 = (double*) calloc(n,sizeof(double));
  double *k2 = (double*) calloc(n,sizeof(double));
  double *k3 = (double*) calloc(n,sizeof(double));
  double *k4 = (double*) calloc(n,sizeof(double));
  
  double t = t0;
  init_vector(n,y,x0,v0);

  for(int ii = 0; ii<NSTEP; ii++){
    t = t0 +ii*dt;
    update_kvect(n,1,dt,t,y,y,k1);
    update_kvect(n,2,dt,t,y,k1,k2);
    update_kvect(n,3,dt,t,y,k2,k3);
    update_kvect(n,4,dt,t,y,k3,k4);
    for(int jj = 0; jj<n; jj++){
      y[jj] += (1.0/6.0)*dt*(k1[jj] + k4[jj] + 2.0*(k2[jj] + k3[jj]));
      printf("%d\t%4.7f\t%4.7f\t%4.7f\n",ii,t,y[0],y[1]);
    }
  }

  free(y);
  free(k1);
  free(k2);
  free(k3);
  free(k4);
}

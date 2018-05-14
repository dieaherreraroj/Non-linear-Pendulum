#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>

/*******************************GLOBAL VARIABLES OF MOTION*********************/

double w = 4.56;
double q = 0.98;
double u = 20.87;
double Fd = 0.95;

/***********************************AUXILIAR ROUTINES**************************/

void copy_vect(int n, double *y, double *x);
void aux_interm(int n, int tag, double h, double *x, double *y, double *z);
void f(double t, double *y);
void update_kvect(int n, int tag, double h, double t, double *y, double *ki, double *kf);
void init_vector(int n, double *y, double x0, double v0);

/***********************************ANALISIS ROUTINES**************************/

void rk4_integ(int n, int NSTEP, double dt, double t0, double x0, double v0, double *x);
void FFTW_Analize(int n, double *x, double *p);

/***********************************IMPLEMENTATIONS****************************/

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
  a[1] = -w*sin(y[0]) - q*y[1] + Fd*sin(u*t);

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

/*****************************MAIN IMPLEMENTATIONS*****************************/

void rk4_integ(int n, int NSTEP, double dt, double t0, double x0, double v0, double *x){
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
      //printf("%d\t%4.7f\t%4.7f\t%4.7f\n",ii,t,y[0],y[1]);
    }
    x[ii] = y[0];
  }

  free(y);
  free(k1);
  free(k2);
  free(k3);
  free(k4);
}

// We compute DFT using std library fftw3. Future feature: include mpi options
// This function is just a basic pilot in order to get started. By now it just
// computes DFT of a single array and returns discete power distribution in a
// different array.

 void FFTW_Analize(int n, double *x, double *p){
   fftw_complex *in, *out;
   fftw_plan plan;

   in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
   plan = fftw_plan_dft_1d(n,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

   for(int ii = 0; ii<n; ii++){
     in[ii][0] = x[ii];
     in[ii][1] = 0.0;
   }

   fftw_execute(plan);

   for(int ii = 0; ii<n; ii++)
     p[ii] = (out[ii][0]*out[ii][0]) + (in[ii][0]*in[ii][0]);

   fftw_destroy_plan(plan);
   fftw_free(in);
   fftw_free(out);
 }

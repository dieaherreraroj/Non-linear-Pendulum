#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>

/*******************************GLOBAL VARIABLES OF MOTION*********************/

// According to Baker and Gollub, a suitable range of parameters to produce a
// chaotic motion is q = 2.0, w = 2/3 and Fd in [0.5:1.5].

double w = 1.0;
double q = 0.5;
double wd = 2.0/3.0;
double Fd = 1.45;

/***********************************AUXILIAR ROUTINES**************************/

void copy_vect(int n, double *y, double *x);
void aux_interm(int n, int tag, double h, double *x, double *y, double *z);
double f(int comp, double t, double *y);
void update_kvect(int n, int tag, double h, double t, double *y, double *ki, double *kf);
void init_vector(int n, double *y, double x0, double v0);
void gen_animgif(int NSTEP, double t0, double dx, double *x);

double e_mec(double x, double v);
double power_mec(double a, double v);

/***********************************ANALISIS ROUTINES**************************/

void rk4_integ(int n, int NSTEP, double dt, double t0, double x0, double v0, double *x);
void FFTW_Analize(int n, double *x, double *p);

/******************************************************************************/

/***********************************IMPLEMENTATIONS****************************/

void copy_vect(int n, double *y, double *x){
  for(int ii = 0; ii<n; ii++)
    x[ii] = y[ii];
}

/******************************************************************************/

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

/******************************************************************************/

double f(int comp, double t, double *y){
  if(comp == 0)
    return y[1];
  if(comp == 1)
    return -w*w*sin(y[0]) - q*y[1] + Fd*cos(y[2]);
  if(comp == 2)
    return wd;
  else
    return 0.0;
}

/******************************************************************************/

void update_kvect(int n, int tag, double h, double t, double *y, double *ki, double *kf){
  double *aux = (double*) calloc(n,sizeof(double));
  double t_aux = 0.0;
  if(tag == 1) t_aux = t;
    else if(tag == 2 || tag == 3) t_aux = t + 0.5*h;
      else if(tag == 4) t_aux = t + h;
  else t_aux = 0.0;
  aux_interm(n,tag,h,y,ki,aux);
  for(int ii = 0; ii<n; ii++)
    kf[ii] = f(ii,t_aux,aux);
  free(aux);
}

/******************************************************************************/

void init_vector(int n, double *y, double x0, double v0){
  y[0] = x0;
  y[1] = v0;
  y[2] = 0.0;
}

/******************************************************************************/

void gen_animgif(int NSTEP, double t0, double dt, double *x){
  FILE *fptr;
  fptr = fopen("anim_gen.gp","w");
  double t;
  fprintf(fptr,"set terminal gif animate delay 2.5\n");
  fprintf(fptr,"set out 'nlp_Fd_1.45.gif'\n");
  fprintf(fptr,"set xrange [-1.5:1.5]\n");
  fprintf(fptr,"set yrange [-1.5:1.5]\n");
  for(int ii = 0; ii<NSTEP; ii+=20){
    t = t0 + ii*dt;
    if(q > 0.0){
      if(t > 10.0/q){
        fprintf(fptr,"plot '-' w l lw 4 lt 8\n");
        fprintf(fptr,"%4.7f\t  %4.7f\n",sin(x[ii]),-cos(x[ii]));
        fprintf(fptr,"%4.7f\t  %4.7f\n",0.0,0.0);
        fprintf(fptr,"e\n");
      }
    }
    else{
      fprintf(fptr,"plot '-' w l lw 4 lt 4\n");
      fprintf(fptr,"%4.7f\t  %4.7f\n",sin(x[ii]),-cos(x[ii]));
      fprintf(fptr,"%4.7f\t  %4.7f\n",0.0,0.0);
      fprintf(fptr,"e\n");
    }
  }
  fclose(fptr);
}

/******************************************************************************/

double e_mec(double x, double v){
  return 0.5*v*v + w*w*(1.0-cos(x));
}

/******************************************************************************/

double power_mec(double a, double v){
  return (Fd*sin(a) - q*v)*v;
}

/******************************************************************************/

/*****************************MAIN IMPLEMENTATIONS*****************************/

void rk4_integ(int n, int NSTEP, double dt, double t0, double x0, double v0, double *x){
  double *y = (double*) calloc(n,sizeof(double));
  double *k1 = (double*) calloc(n,sizeof(double));
  double *k2 = (double*) calloc(n,sizeof(double));
  double *k3 = (double*) calloc(n,sizeof(double));
  double *k4 = (double*) calloc(n,sizeof(double));
  FILE *fptr;
  fptr = fopen("ps_Fd_1.45.txt","w");
  double t = t0;
  init_vector(n,y,x0,v0);
  x[0] = x0;
  for(int ii = 1; ii<NSTEP; ii++){
    t = t0 +ii*dt;
    update_kvect(n,1,dt,t,y,y,k1);
    update_kvect(n,2,dt,t,y,k1,k2);
    update_kvect(n,3,dt,t,y,k2,k3);
    update_kvect(n,4,dt,t,y,k3,k4);
    for(int jj = 0; jj<n; jj++)
      y[jj] += (1.0/6.0)*dt*(k1[jj] + k4[jj] + 2.0*(k2[jj] + k3[jj]));
    double theta = atan2(sin(y[0]),cos(y[0]));
    if(q > 0.0){
      if(t > 10.0/q){
        fprintf(fptr,"%4.7f\t %4.7f\t %4.7f\t %4.7f\t %4.7f\n",t,theta,y[1],e_mec(theta,y[1])/(2*w*w),power_mec(y[2],y[1]));
        x[ii] = theta;
      }
      else
        x[ii] = 0.0;
    }
    else{
      fprintf(fptr,"%4.7f\t %4.7f\t %4.7f\t %4.7f\t %4.7f\n",t,theta,y[1],e_mec(theta,y[1])/(2*w*w),power_mec(y[2],y[1]));
      x[ii] = theta;
    }
  }
  fclose(fptr);
  free(y);
  free(k1);
  free(k2);
  free(k3);
  free(k4);
}

/******************************************************************************/

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
     p[ii] = (out[ii][0]*out[ii][0]) + (out[ii][1]*out[ii][1]);

   fftw_destroy_plan(plan);
   fftw_free(in);
   fftw_free(out);
 }

 /******************************************************************************/

#include "routines.h"

int main(int argc, char **argv){

  int n = 2;
  int NSTEP = 100000;

  double dt = 0.0008;
  double t0 = 0.0;
  double x0 = 3.1416*0.75;
  double v0 = 0.0;
  double f_samp = 1.0/dt;

  double *x = (double*) calloc(NSTEP,sizeof(double));
  double *p = (double*) calloc(NSTEP,sizeof(double));

  rk4_integ(n,NSTEP,dt,t0,x0,v0,x);
  FFTW_Analize(NSTEP,x,p);

  double f  = 0;
  for(int ii = 0; f<5.0; ii++){
    f = f_samp*((1.0*ii)/NSTEP);
    printf("%4.7f\t%4.7f\n",f,p[ii]);
  }

  return 0;

}

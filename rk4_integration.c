#include "routines.h"

int main(int argc, char **argv){

  int n = 3;
  int NSTEP = 0.0;

  if (q > 0.00) {
    NSTEP = 1000*((int) 10.0/q + 3*((int) (2*M_PI)/wd));
  } else {
    NSTEP = 1000*3*((int) (2*M_PI)/wd);
  }

  double dt = 0.001;
  double t0 = 0.0;
  double x0 = M_PI*0.80;
  double v_crit = 2.0*w*cos(x0/2.0);
  double v0 = 0.0*v_crit;
  double f_samp = 1.0/dt;

  double *x = (double*) calloc(NSTEP,sizeof(double));
  double *p = (double*) calloc(NSTEP,sizeof(double));

  rk4_integ(n,NSTEP,dt,t0,x0,v0,x);
  //gen_animgif(NSTEP,t0,dt,x);

  FFTW_Analize(NSTEP,x,p);

  FILE *fptr;
  fptr = fopen("fft_Fd_1.45.txt","w");
  double f  = 0;
  for(int ii = 0; f<6.0; ii++){
    f = f_samp*((1.0*ii)/NSTEP);
    fprintf(fptr,"%4.7f\t%4.7f\n",f,2.0*dt*p[ii]/NSTEP);
  }
  fclose(fptr);

  free(x);
  free(p);
  return 0;

}

//!########################################################################
//!                                                                       #
//! Copyright (C) University .                     #
//! This file is a modeling for 1d acoustic media  #
//!                                                                       #
//! 
//! Feb. 2014                                                      # 
//!                                                                       #
//!########################################################################
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "init.h"
#include "P0an.h"
//#include "Reflectiondata.h"

int main( int argc, char *argv[] )
{
  // grids dimensions 
  int nz    = 101; 
  int nx    = 401,nx2=(nx-1)/2;
  int ndims = nx*nz;      

  // physical grid size
  double dz = 5.0;     // modeling grid spacial step in z
  double dx = 5.0;    // modeling grid spacial step in x

 
  FILE* file1;
  char fname1[100];
  int ret;
  sprintf(fname1,"test_1D_2P_acoustic.out");
  file1 = fopen(fname1,"w");

  fprintf(file1, "--!\t                                     \n");
  fprintf(file1, "--!\t1D acoustic Volterra modeling research code\n");
  fprintf(file1, "--!\t                                     \n");

  ret = fflush(file1);

  const char * NameWeekDay[] = {"Sunday", "Monday", "Tuesday", "Wenesday", "Thursday", "Friday", "Saturday"};
  time_t timestamp;
  struct tm * tc;
    
  timestamp = time(NULL);
  tc = localtime(&timestamp);
  fprintf(file1,"--!\t                                     \n");
  fprintf(file1,"--!\tStart Time : %s, ", NameWeekDay[tc->tm_wday]);
  fprintf(file1,"%02u/%02u/%04u, ", tc->tm_mon+1, tc->tm_mday, 1900 + tc->tm_year);
  fprintf(file1,"Clock : %02uh %02umn %02usec.\n", tc->tm_hour, tc->tm_min, tc->tm_sec);

  ret = fflush(file1);

  //velocity & density initinalizaiton
  double c0=1500; 
  double rho0=1;
  double bulk0=c0*c0*rho0;

  double* vel  = malloc( nz*sizeof(double) );
  init_acou_layer_velocity(file1,  dx,dz,nx,nz, c0, vel); 
  
  double* den  = malloc( nz*sizeof(double) );
  init_acou_layer_density(file1, dx,dz,nx,nz,rho0,den);

  double* alpha  = malloc( nz*sizeof(double) );
  double* belta  = malloc( nz*sizeof(double) );
  double* gamma  = malloc( nz*sizeof(double) );
  double* beltad  = malloc( nz*sizeof(double) );
  FILE* file2;  
  char fname2[100];
  sprintf(fname2,"alpha.dat");
  file2 = fopen(fname2,"w");
  int ix,iz;
  double bulk;

  for (iz=0;iz<nz;iz++){
	bulk=vel[iz]*vel[iz]*den[iz];
	alpha[iz]=1-bulk0/bulk;
	belta[iz]=1-rho0/den[iz];
	gamma[iz]=1-c0*c0/vel[iz]/vel[iz];
        fprintf(file2," %d  %.12lf %.12lf %.12lf \n",iz, alpha[iz],belta[iz],gamma[iz]);                         	
  }
  fclose(file2);

  beltad[0]=0;
  for (iz=1;iz<nz;iz++){	
	beltad[iz]=(belta[iz]-belta[iz-1])/dz;	                          	
   }

	
  // time and frequency discretization parameter
  int   nt = 1000;        // number of steps
  double dt = 0.001;      // time step

  int nw=151;
  double dw=2;

  double* sourcet = malloc(nt*sizeof(double) );
  double* sourcefren  = malloc( (2*nw)*sizeof(double) );
  double* fren  = malloc( nw*sizeof(double) );
  init_source_ricker_fwps( file1, nt, dt,nw,dw, sourcet,fren, sourcefren); 
  
  double pi=3.1415926;
  int iw,ik,index,index2;
  double omega,k0,k02;
  int nkg=nx,nkg2=(nkg-1)/2;
  double dkg=2*pi/(nx-1)/dx/2;
  double kg,kg2,qg2,qg;
  double x,z;
  double* Pk1TRe = malloc(nw*nkg*nz*sizeof(double));
  double* Pk1TIm = malloc(nw*nkg*nz*sizeof(double));
  double* PkTRe = malloc(nw*nkg*nz*sizeof(double));
  double* PkTIm = malloc(nw*nkg*nz*sizeof(double));
  
  double* Pk1Re = malloc(nz*sizeof(double));
  double* Pk1Im = malloc(nz*sizeof(double));
  double* PkRe = malloc(nz*sizeof(double));
  double* PkIm = malloc(nz*sizeof(double));
 
 
  double* Rkr = malloc(nw*nkg*sizeof(double));
  double* Rki = malloc(nw*nkg*sizeof(double));
  
  double Rk[2];

  double* PsTRe = malloc(nw*nx*nz*sizeof(double));
  double* PsTIm = malloc(nw*nx*nz*sizeof(double));
  double sum1,sum2,arg;
  for (iw=0;iw<nw;iw++){
	omega=fren[iw];
        k0=omega/c0;
        k02=k0*k0;
	fprintf(file1, "--!\tFrequency_%f\n",omega);
        for(ik=0;ik<nkg;ik++){
		kg=(ik-nkg2)*dkg;
                kg2=kg*kg;
                qg2=k02-kg2;
                Rkr[iw*nkg+ik] = 0; 
                Rki[iw*nkg+ik] = 0; 
		for(iz=0;iz<nz;iz++){
			index=iw*(nkg*nz)+ik*nz+iz;
			Pk1TRe[index]=0;
                	Pk1TIm[index]=0;
			PkTRe[index]=0;
                	PkTIm[index]=0;
		}
		if (qg2>0){
			qg=sqrt(qg2);
			//ComputeRk( file1, nz, dz, qg,  k0, kg,alpha,belta,Rk,Pk1Re, Pk1Im, PkRe, PkIm );
			//ComputeRk2(file1, omega,nz,  dz, qg,  k0, kg,alpha,belta,beltad,vel, Rk, Pk1Re, Pk1Im, PkRe, PkIm );
                        ComputeRk3(file1, omega,nz,  dz, qg,  k0, kg,gamma,den,beltad,vel, Rk, Pk1Re, Pk1Im, PkRe, PkIm );
			Rkr[iw*nkg+ik]=Rk[0];
			Rki[iw*nkg+ik]=Rk[1];
			for(iz=0;iz<nz;iz++){
				index=iw*(nkg*nz)+ik*nz+iz;
				Pk1TRe[index]=Pk1Im[iz]*rho0/2/qg;
                		Pk1TIm[index]=+Pk1Re[iz]*rho0/2/qg;
				PkTRe[index]= PkIm[iz]*rho0/2/qg;
                		PkTIm[index]= +PkRe[iz]*rho0/2/qg;
				/*Pk1TRe[index]=Pk1Re[iz];
                		Pk1TIm[index]=Pk1Im[iz];
				PkTRe[index]= PkRe[iz];
                		PkTIm[index]= PkIm[iz];*/
			}
			
			
		}
	}
        if((iw%50)==0){
		sprintf(fname2,"P1_iw_%d.dat", iw);
		file2 = fopen(fname2,"w");
 		for (ik=0;ik<nkg;ik++){
			kg=(ik-nkg2)*dkg;
			for (iz=0;iz<nz;iz++){
				index=iw*(nkg*nz)+ik*nz+iz;
				z=iz*dz;		
                          	fprintf(file2," %f %f %.12lf %.12lf %.12lf %.12lf\n",kg, z, Pk1TRe[index],Pk1TIm[index], PkTRe[index],PkTIm[index]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);
	}
	
        for (iz=0;iz<nz;iz++){
		for (ix=0;ix<nx;ix++){
			x=(ix-nx2)*dx;
			index=iw*(nx*nz)+ix*nz+iz;
			sum1 = 0.0;
                     	sum2 = 0.0;
			//index2=iw*(nkg*nz)+nkg2*nz+iz;
			for (ik=0;ik<nkg;ik++){
				index2=iw*(nkg*nz)+ik*nz+iz;
				kg=(ik-nkg2)*dkg;
				arg=kg*x;
				sum1+=cos(arg)*PkTRe[index2]-sin(arg)*PkTIm[index2];
				sum2+=sin(arg)*PkTRe[index2]+cos(arg)*PkTIm[index2];
			}
			PsTRe[index]=(sum1*sourcefren[2*iw]-sum2*sourcefren[2*iw+1])*dkg;
			PsTIm[index]=(sum2*sourcefren[2*iw]+sum1*sourcefren[2*iw+1])*dkg;
		}
	}
	 if((iw%50)==0){
		sprintf(fname2,"Ps_iw_%d.dat", iw);
		file2 = fopen(fname2,"w");
 		for (ix=0;ix<nx;ix++){
			x=(ix-nx2)*dx;
			for (iz=0;iz<nz;iz++){
				index=iw*(nx*nz)+ix*nz+iz;
				z=iz*dz;		
                          	fprintf(file2," %f %f %.12lf %.12lf\n",x, z, PsTRe[index],PsTIm[index]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);
	}
							       	
	
  }
  sprintf(fname2,"Reflection.dat");
  file2 = fopen(fname2,"w");
   for (iw=0;iw<nw;iw++){
	omega=fren[iw];
        for(ik=0;ik<nkg;ik++){
		kg=(ik-nkg2)*dkg;
		fprintf(file2," %f %f %.12lf %.12lf\n",omega, kg, Rkr[iw*nkg+ik],Rki[iw*nkg+ik]);
	}
	 fprintf(file2," \n");
   }
    fclose(file2);

  nt=501;
  double* Pret = malloc(nz*nx*nt*sizeof(double) );  
  Pwtot(nw,nx, nz,  dw, dx, dz, dt, nt,fren,PsTRe, PsTIm, Pret);	

   fprintf(file1, "--!\t        Done                             \n");
   
}

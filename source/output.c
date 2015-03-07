#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



void Pwtot(int nw, int nx,  int nz,  double dw, double dx,double dz, double dt, int nt,  double* fren, double *P1r, double *P1i, double *Pret){
	int ix,iw,it,iz;
	double pi = 3.14159265358979323846; 
	double t,x,z, arg=0;
        for(iz=0;iz<nz;iz++){
		for(ix=0;ix<nx;ix++){
		x=(ix-(nx-1)/2)*dx;
		for(it=0;it<nt;it++){
			t=it*dt;
			Pret[iz*nt*nx+ix*nt+it]=0;
			for (iw=0;iw<nw;iw++){
				arg = fren[iw]*t;		
				Pret[iz*nt*nx+ix*nt+it]+=cos(arg)*P1r[iw*(nx*nz)+ix*nz+iz]-sin(arg)*P1i[iw*(nx*nz)+ix*nz+iz];
			}
			Pret[iz*nt*nx+ix*nt+it]*=dw*2/2/pi;
			
		}
		
		
	}
	}
	FILE* file2;  
  	char fname2[100];
	for(it=0;it<nt;it++){
		if((it%100)==0){
		sprintf(fname2,"snapshot_%d.dat", it);
		file2 = fopen(fname2,"w");
 		for (ix=0;ix<nx;ix++){
			x=(ix-(nx-1)/2)*dx;
			for (iz=0;iz<nz;iz++){
				z=iz*dz;		
                          	fprintf(file2," %f %f %.12lf \n", x, z, Pret[iz*nt*nx+ix*nt+it]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);
		}
	}
/*
	for(iz=0;iz<nz;iz++){
		if((iz%10)==0){
		sprintf(fname2,"trace_%d.dat", iz);
		file2 = fopen(fname2,"w");
 		for (ix=0;ix<nx;ix++){
			x=(ix-(nx-1)/2)*dx;
			for (it=0;it<nt;it++){
				t=it*dt;		
                          	fprintf(file2," %f %f %.12lf \n", x, t, Pret[iz*nt*nx+ix*nt+it]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);
		}
	}
*/	
}
void Pwtot2(int nw, int nx,  int nz,  double dw, double dx,double dz, double dt, int nt,  double* fren, double *P1r, double *P1i, double *Pret){
	int ix,iw,it,iz;
	double pi = 3.14159265358979323846; 
	double t,x,z, arg=0;
        for(iz=0;iz<nz;iz++){
		for(ix=0;ix<nx;ix++){
		x=(ix-(nx-1)/2)*dx;
		for(it=0;it<nt;it++){
			t=it*dt;
			Pret[iz*nt*nx+ix*nt+it]=0;
			for (iw=0;iw<nw;iw++){
				arg = fren[iw]*t;		
				Pret[iz*nt*nx+ix*nt+it]+=cos(arg)*P1r[iw*(nx*nz)+ix*nz+iz]-sin(arg)*P1i[iw*(nx*nz)+ix*nz+iz];
			}
			Pret[iz*nt*nx+ix*nt+it]*=dw*2/2/pi;
			
		}
		
		
	}
	}
	FILE* file2;  
  	char fname2[100];
	for(it=0;it<nt;it++){
		if((it%100)==0){
		sprintf(fname2,"snapshot_%d.dat", it);
		file2 = fopen(fname2,"w");
 		for (ix=0;ix<nx;ix++){
			x=(ix-(nx-1)/2)*dx;
			for (iz=0;iz<nz;iz++){
				z=iz*dz;		
                          	fprintf(file2," %f %f %.12lf \n", x, z, Pret[iz*nt*nx+ix*nt+it]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);
		}
	}
/*
	for(iz=0;iz<nz;iz++){
		if((iz%10)==0){
		sprintf(fname2,"trace_%d.dat", iz);
		file2 = fopen(fname2,"w");
 		for (ix=0;ix<nx;ix++){
			x=(ix-(nx-1)/2)*dx;
			for (it=0;it<nt;it++){
				t=it*dt;		
                          	fprintf(file2," %f %f %.12lf \n", x, t, Pret[iz*nt*nx+ix*nt+it]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);
		}
	}
*/	
}

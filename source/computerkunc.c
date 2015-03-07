#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <omp.h>


#include "operation_complex.h"
#include "computerkunc.h"


void ComputeRk( FILE *file1, int nz, double dz,double qg, double k0,double kg, double *alpha,double *belta, double *Rk,double* Pk1Re,double* Pk1Im,double* PkRe,double* PkIm ){
  int iz,iz2;
  double k02=k0*k0;
  double kg2=kg*kg;
  double qg2=qg*qg;
  double z=0.0 ;
  double I1temp[2], I2temp[2], Exptem1[2];
  double I1int[2], I2int[2];
  double Rkn[2], Rkd[2];

  double  Pk1Retemp, Pk1Imtemp, I1re,I1im,I2re,I2im;
  double* Pk1dzRe = malloc(nz*sizeof(double));
  double* Pk1dzIm = malloc(nz*sizeof(double));


     FILE* file2;
     char fname2[100];

     //sprintf(fname2,"rktk_.dat" );
     //file2 = fopen(fname2,"w");
        
     for(iz=0;iz<nz;iz++){
		z = iz*dz;
                if(iz<3){
                	Pk1Re[iz]=cos(qg*z);
                	Pk1Im[iz]=sin(qg*z);
			Pk1dzRe[iz]=-qg2*belta[iz]*cos(qg*z);
                	Pk1dzIm[iz]=-qg2*belta[iz]*sin(qg*z);
			//Pk1dzRe[iz]=0;
                	//Pk1dzIm[iz]=0;
          	}else{
			Pk1dzRe[iz-1]=(belta[iz-1]*(Pk1Re[iz-1]-Pk1Re[iz-2])-belta[iz-2]*(Pk1Re[iz-2]-Pk1Re[iz-3]))/dz/dz;
                	Pk1dzIm[iz-1]=(belta[iz-1]*(Pk1Im[iz-1]-Pk1Im[iz-2])-belta[iz-2]*(Pk1Im[iz-2]-Pk1Im[iz-3]))/dz/dz;
			//Pk1dzRe[iz-1]=0;
                	//Pk1dzIm[iz-1]=0;
                	Pk1Re[iz]=cos(qg*z);
               		Pk1Im[iz]=sin(qg*z);
                	for(iz2=1;iz2<iz;iz2++){
                		Pk1Retemp=dz*sin(qg*(iz-iz2)*dz)*((k02*alpha[iz2]-kg2*belta[iz2])*Pk1Re[iz2]+Pk1dzRe[iz2])/qg;
                		Pk1Imtemp=dz*sin(qg*(iz-iz2)*dz)*((k02*alpha[iz2]-kg2*belta[iz2])*Pk1Im[iz2]+Pk1dzIm[iz2])/qg;
                		Pk1Re[iz]+=Pk1Retemp;
                		Pk1Im[iz]+=Pk1Imtemp;
						
                	}
		}
	}

	I1re=0;
        I1im=0;
	I2re=0;
        I2im=0;
        for(iz=1;iz<nz-1;iz++){
		z=iz*dz;
        	I1temp[0]=(k02*alpha[iz]-kg2*belta[iz])*Pk1Re[iz]+Pk1dzRe[iz];
        	I1temp[1]=(k02*alpha[iz]-kg2*belta[iz])*Pk1Im[iz]+Pk1dzIm[iz];
		I2temp[0]=I1temp[0];
        	I2temp[1]=-I1temp[1];
        	Exptem1[0]=cos(qg*z);
        	Exptem1[1]=sin(qg*z);
        	prodcomp(I1temp,Exptem1,I1int);
		prodcomp(I2temp,Exptem1,I2int);
        	I1re+=I1int[0]*dz;
        	I1im+=I1int[1]*dz;
		I2re+=I2int[0]*dz;
        	I2im+=I2int[1]*dz;
                	
        }
          // compute Rk
          Rkn[0]=I1im/2/qg;
          Rkn[1]=-I1re/2/qg;
          Rkd[0]=1-I2im/2/qg;
          Rkd[1]=I2re/2/qg;
          divcomp(Rkn,Rkd,Rk);
	  for(iz=0;iz<nz;iz++){
		PkRe[iz]=Pk1Re[iz]+Rk[0]*Pk1Re[iz]+Rk[1]*Pk1Im[iz];
                PkIm[iz]=Pk1Im[iz]+Rk[1]*Pk1Re[iz]-Rk[0]*Pk1Im[iz];	
                	
           }         

}

void ComputeRk2( FILE *file1, double omega, int nz, double dz,double qg, double k0,double kg, double *alpha,double *belta,double *beltad,double *vel, double *Rk,double* Pk1Re,double* Pk1Im,double* PkRe,double* PkIm ){
  int iz,iz2;
  double kp2;
  double k02=k0*k0;
  double kg2=kg*kg;
  double qg2=qg*qg;
  double z=0.0 ;
  double I1temp[2], I2temp[2], Exptem1[2];
  double I1int[2], I2int[2];
  double Rkn[2], Rkd[2];

  double  Pk1Retemp, Pk1Imtemp, I1re,I1im,I2re,I2im;
  double* Pk1dzRe = malloc(nz*sizeof(double));
  double* Pk1dzIm = malloc(nz*sizeof(double));


     FILE* file2;
     char fname2[100];

     //sprintf(fname2,"rktk_.dat" );
     //file2 = fopen(fname2,"w");
        
     for(iz=0;iz<nz;iz++){
		z = iz*dz;
                if(iz<3){
                	Pk1Re[iz]=cos(qg*z);
                	Pk1Im[iz]=sin(qg*z);
			Pk1dzRe[iz]=-qg*beltad[iz]*sin(qg*z);
                	Pk1dzIm[iz]=-qg*beltad[iz]*cos(qg*z);
          	}else{
			Pk1dzRe[iz-1]=beltad[iz-1]*(Pk1Re[iz-1]-Pk1Re[iz-2])/dz;
                	Pk1dzIm[iz-1]=beltad[iz-1]*(Pk1Im[iz-1]-Pk1Im[iz-2])/dz;
                	Pk1Re[iz]=cos(qg*z);
               		Pk1Im[iz]=sin(qg*z);
                	for(iz2=1;iz2<iz;iz2++){
				kp2=omega*omega/vel[iz2]/vel[iz2];
                		Pk1Retemp=dz*sin(qg*(iz-iz2)*dz)*((k02*alpha[iz2]-kp2*belta[iz2])*Pk1Re[iz2]+Pk1dzRe[iz2])/qg;
                		Pk1Imtemp=dz*sin(qg*(iz-iz2)*dz)*((k02*alpha[iz2]-kp2*belta[iz2])*Pk1Im[iz2]+Pk1dzIm[iz2])/qg;
                		Pk1Re[iz]+=Pk1Retemp;
                		Pk1Im[iz]+=Pk1Imtemp;
						
                	}
		}
	}

	I1re=0;
        I1im=0;
	I2re=0;
        I2im=0;
        for(iz=1;iz<nz-1;iz++){
		z=iz*dz;
		kp2=omega*omega/vel[iz]/vel[iz];
        	I1temp[0]=(k02*alpha[iz]-kp2*belta[iz])*Pk1Re[iz]+Pk1dzRe[iz];
        	I1temp[1]=(k02*alpha[iz]-kp2*belta[iz])*Pk1Im[iz]+Pk1dzIm[iz];
		I2temp[0]=I1temp[0];
        	I2temp[1]=-I1temp[1];
        	Exptem1[0]=cos(qg*z);
        	Exptem1[1]=sin(qg*z);
        	prodcomp(I1temp,Exptem1,I1int);
		prodcomp(I2temp,Exptem1,I2int);
        	I1re+=I1int[0]*dz;
        	I1im+=I1int[1]*dz;
		I2re+=I2int[0]*dz;
        	I2im+=I2int[1]*dz;
                	
        }
          // compute Rk
          Rkn[0]=I1im/2/qg;
          Rkn[1]=-I1re/2/qg;
          Rkd[0]=1-I2im/2/qg;
          Rkd[1]=I2re/2/qg;
          divcomp(Rkn,Rkd,Rk);
	  for(iz=0;iz<nz;iz++){
		PkRe[iz]=Pk1Re[iz]+Rk[0]*Pk1Re[iz]+Rk[1]*Pk1Im[iz];
                PkIm[iz]=Pk1Im[iz]+Rk[1]*Pk1Re[iz]-Rk[0]*Pk1Im[iz];	
                	
           }         

}
// Compute Rk
// we suppose V is cst after zmax

/*
void ComputeRk_outgoing( FILE *file1, int nz,  double dz, double qg, double k0, double *alpha, double *belta, double *Rk,double* Pk1Re,double* Pk1Im,double* PkRe,double* PkIm ){

 int i=0, j=0, m=0, jj=0;
  double k=0.0, z=0.0, lz=0.0, kappa=0.0, kappa_scale=0.0;
  double z1[2], z2[2], zres[2];
  double Rkn[2], Rkd[2];
  double  Pk1Retemp, Pk1Imtemp;
 double matar[4], mataim[4];
 double rhsr[2], rhsim[2]; 


 // double arg,s,c,tkcoef;

  kappa_scale = sqrt(1.0-vexact[nz-1]); 

  int nzset = 5;
  int nzt = nz+2*nzset;

  double* Pk1Re = malloc(nzt*sizeof(double));
  double* Pk1Im = malloc(nzt*sizeof(double));


      fprintf(file1," rk outgoing \n");
      rkr[0] = (1-sqrt(1-vexact[nz-1]))/(1+sqrt(1-vexact[nz-1]));
      rki[0] = 0.0;
      fprintf(file2," %lf %lf %lf \n",  0.0, rkr[0], rki[0]);


      for(i=1;i<nk;i++){

      	k=i*dk+0.0;
        kappa = kappa_scale*k;
        //fprintf(stdout," k %lf \n", k);

        // computation of Pk1
       	for(j=0;j<nz;j++){

		z = j*dz+zs;
                if(j<2){
                	Pk1Re[j]=cos(k*z);
                	Pk1Im[j]=sin(k*z);
          	}else{
                	Pk1Re[j]=cos(k*z);
               		Pk1Im[j]=sin(k*z);
                	for(m=1;m<j;m++){
                		Pk1Retemp=k*dz*sin(k*(j-m)*dz)*vexact[m]*Pk1Re[m];
                		Pk1Imtemp=k*dz*sin(k*(j-m)*dz)*vexact[m]*Pk1Im[m];
                		Pk1Re[j]+=Pk1Retemp;
                		Pk1Im[j]+=Pk1Imtemp;
                	}
		}
                //fprintf(file2," %lf %lf %lf  \n",  z, Pk1Re[j], Pk1Im[j]);

      	 }

       rhsr[0]  = 0.0;
       rhsr[1]  = 0.0;
       rhsim[0] = 0.0;
       rhsim[1] = 0.0;

       matar[0]  = 0.0;
       matar[1]  = 0.0;
       matar[2]  = 0.0;
       matar[3]  = 0.0;

       mataim[0]  = 0.0;
       mataim[1]  = 0.0;
       mataim[2]  = 0.0;
       mataim[3]  = 0.0;

       // after zmax V(z)= cst
       // P1(z)+ rk P2(z) = Tk P3(z)
       // P2(z) = P1*(z) 
       // P3(z) = Tk exp(i kappa z)
       // kappa = k*c0/cj

       // compute first set
       // 1rst line of equation system
       	for(j=0;j<nzset;j++){

              jj = nz+j;
              z = jj*dz+zs;

                	Pk1Re[jj]=cos(k*z);
               		Pk1Im[jj]=sin(k*z);
                	for(m=1;m<nz;m++){
                		Pk1Retemp=k*dz*sin(k*(jj-m)*dz)*vexact[m]*Pk1Re[m];
                		Pk1Imtemp=k*dz*sin(k*(jj-m)*dz)*vexact[m]*Pk1Im[m];
                		Pk1Re[jj]+=Pk1Retemp;
                		Pk1Im[jj]+=Pk1Imtemp;
                	}
                	for(m=nz;m<jj;m++){
                		Pk1Retemp=k*dz*sin(k*(jj-m)*dz)*vexact[nz-1]*Pk1Re[m];
                		Pk1Imtemp=k*dz*sin(k*(jj-m)*dz)*vexact[nz-1]*Pk1Im[m];
                		Pk1Re[jj]+=Pk1Retemp;
                		Pk1Im[jj]+=Pk1Imtemp;
                        }


              // rhs = sum P1(zi)
              rhsr[0]  += Pk1Re[jj];
              rhsim[0] += Pk1Im[jj];

              // mata(i,0) = -sum P2(zi) = -sum P1*
              matar[0]  -= Pk1Re[jj];
              mataim[0] += Pk1Im[jj];

              // mata(i,1) = sum exp( i kappa z)
              matar[1]  += cos(kappa*z);
              mataim[1] += sin(kappa*z);

             //   fprintf(file2," %lf %lf %lf  \n",  z, Pk1Re[j], Pk1Im[j]);



        }        

       // compute second set
       // 2nd line of equation system 
       	for(j=0;j<nzset;j++){

              jj = nz+nzset+j;
              z = jj*dz+zs;

                	Pk1Re[jj]=cos(k*z);
               		Pk1Im[jj]=sin(k*z);
                	for(m=1;m<nz;m++){
                		Pk1Retemp=k*dz*sin(k*(jj-m)*dz)*vexact[m]*Pk1Re[m];
                		Pk1Imtemp=k*dz*sin(k*(jj-m)*dz)*vexact[m]*Pk1Im[m];
                		Pk1Re[jj]+=Pk1Retemp;
                		Pk1Im[jj]+=Pk1Imtemp;
                	}
                	for(m=nz;m<jj;m++){
                		Pk1Retemp=k*dz*sin(k*(jj-m)*dz)*vexact[nz-1]*Pk1Re[m];
                		Pk1Imtemp=k*dz*sin(k*(jj-m)*dz)*vexact[nz-1]*Pk1Im[m];
                		Pk1Re[jj]+=Pk1Retemp;
                		Pk1Im[jj]+=Pk1Imtemp;
                        }

              rhsr[1]  += Pk1Re[jj];
              rhsim[1] += Pk1Im[jj];

              // mata(i,0) = -sum P2(zi) = -sum P1*
              matar[2]  -= Pk1Re[jj];
              mataim[2] += Pk1Im[jj];

              // mata(i,1) = sum exp( i kappa z)
              matar[3]  += cos(kappa*z);
              mataim[3] += sin(kappa*z);

              //  fprintf(file2," %lf %lf %lf  \n",  z, Pk1Re[j], Pk1Im[j]);


        }    

    

         //         | rhs[0] a(0,1) |
         //         | rhs[1] a(1,1) |    
         // rk =    ---------------------
         //             det(A)
         //

         z1[0] = rhsr[0];
         z1[1] = rhsim[0];
         z2[0] = matar[3];
         z2[1] = mataim[3];

         prodcomp( z1, z2, zres);

         Rkn[0]=zres[0];
         Rkn[1]=zres[1];


         z1[0] = rhsr[1];
         z1[1] = rhsim[1];
         z2[0] = matar[1];
         z2[1] = mataim[1];

         prodcomp( z1, z2, zres);
         Rkn[0] -= zres[0];
         Rkn[1] -= zres[1];


         z1[0] = matar[0];
         z1[1] = mataim[0];
         z2[0] = matar[3];
         z2[1] = mataim[3];

         prodcomp( z1, z2, zres);

         Rkd[0]=zres[0];
         Rkd[1]=zres[1];


         z1[0] = matar[2];
         z1[1] = mataim[2];
         z2[0] = matar[1];
         z2[1] = mataim[1];

         prodcomp( z1, z2, zres);
         Rkd[0] -= zres[0];
         Rkd[1] -= zres[1];

         divcomp(Rkn,Rkd,Rk);
 	

    fprintf(file2," %lf %lf %lf\n",  k, Rk[0], Rk[1]);
                rkr[i] = Rk[0]; 
                rki[i] = Rk[1]; 
		
         

  }

  fclose(file2);

  free(Pk1Re);
  free(Pk1Im);

}
*/


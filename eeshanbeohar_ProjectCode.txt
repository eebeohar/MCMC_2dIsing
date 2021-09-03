/*2-d Ising Model Project Code Eeshan Beohar 19036762014 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//Functions
float rnd();
float rnda2b(float a, float b);
float E_total (int L,int spin[L+1][L+1]);
float M_total (int L,int spin[L+1][L+1]);
void RandomInitial(int N,int L,float T,float E[N],float M[N],int spin[L+1][L+1]);
void ColdInitial(int N,int L,float T,float E[N],float M[N],int spin[L+1][L+1]);
float DeltaE(int x,int y,int L,int spin[L+1][L+1]);
void data(float T,int Nequil,int mcs,int Nskip,int L,float E[mcs+Nequil],float M[mcs+Nequil]);
void Metropolis(int N,int L,int spin[L+1][L+1],float E[N],float T,float M[N]);
void save_config(int N,int L,float T,int spin[L+1][L+1],FILE*fp);

float rnd(){
	return ((float)rand()/(float)RAND_MAX);
}

float rnda2b(float a, float b){
		return (a + (b-a)*rnd());
}
void save_config(int N,int L,float T,int spin[L+1][L+1],FILE*fp)
   {
      int x,y,i;
      for (y=1; y <= L; y++){
         for (x=1; x <= L; x++){
            fprintf(fp,"%d\t",spin[x][y]);}
            fprintf(fp,"\n");}
   }

float E_total (int L,int spin[L+1][L+1]){
	
	int x,y,up,right,sum;
	float Eavg;
      for (y = 1; y <= L; y++)
        {
/*       periodic boundary conditions */
         if (y==(L)) 
           up = 1;
         else 
           up = y + 1;
         for (x = 1; x <= L; x++)
           {
/*           conditional expression */
             if (x==(L)) 
               right = 1; 
             else 
               right = x + 1;
             sum = spin[x][up] + spin[right][y];
             Eavg += (- spin[x][y]*sum);
            }
         }
         return Eavg;
}

float M_total (int L,int spin[L+1][L+1]){
	
	int x,y;
	float Mavg;
      for (y = 1; y <= L; y++)
        {
         for (x = 1; x <= L; x++)
           {
             Mavg += spin[x][y];
            }
         }
         return Mavg;
}





void ColdInitial(int N,int L,float T,float E[N],float M[N],int spin[L+1][L+1])
   {
      int x,y,up,right;
      float sum;  
    //random initial configuration (Hot Start) 
     /* for (y = 1; y <= *L; y++)
        for (x = 1; x <= *L; x++)
          {
            if (rnd() < 0.5) 
              spin[x][y] = 1;
            else  
              spin[x][y] = -1;
          }*/
          
        //Cold Start  
            for (y = 1; y <= L; y++)
			{
        		for (x = 1; x <= L; x++)
	          	{
	              spin[x][y] = 1;
	       		}
			 }
          E[0] = E_total(L,spin); M[0] = M_total(L,spin);
          //printf("Initial E = %f\tInitial M = %f\n",E,M);
}

void RandomInitial(int N,int L,float T,float E[N],float M[N],int spin[L+1][L+1])
   {
      int x,y,up,right;
      float sum;  
    //random initial configuration (Hot Start) 
     for (y = 1; y <= L; y++)
        for (x = 1; x <= L; x++)
          {
            if (rnd() < 0.5) 
              spin[x][y] = 1;
            else  
              spin[x][y] = -1;
          }
          E[0] = E_total(L,spin); M[0] = M_total(L,spin);
          //printf("Initial E = %f\tInitial M = %f\n",E,M);
}





float DeltaE(int x,int y,int L,int spin[L+1][L+1])
   {
/*    periodic boundary conditions */
      int dE,left,right,up,down;
      if (x == 1) 
        {
         left = spin[L][y];
         right = spin[2][y];
         }
      else if (x == L) 
         {
         left = spin[L-1][y];
         right = spin[1][y];
         }
      else
         {
         left = spin[x-1][y];
         right = spin[x+1][y];
         }
      if (y == 1) 
         {
         up = spin[x][2];
         down = spin[x][L];
         }
      else if (y == L) 
         {
         up = spin[x][1];
         down = spin[x][L-1];
         }
      else
         {
         up = spin[x][y+1];
         down = spin[x][y-1];
         }
      dE = 2*spin[x][y]*(left + right + up + down);
      return dE;
   }
             
void Metropolis(int mcs,int L,int spin[L+1][L+1],float E[mcs],float T,float M[mcs])
   {
      int x,y,i;
	  float dE,r;
	 
	  //printf("Met E = %f\tMet M = %f\n",E,M);
	  
	  FILE *fp;
	  fp = fopen("MvsE.txt","w");
	  
   	for (i=1;i<=mcs;i++){
     
	/*  random x and y coordinates for trial spin  */
	     x = (int)(rnda2b(1,L));
	     y = (int)(rnda2b(1,L));
	     
	     spin[x][y] = -spin[x][y];
	     
	     dE = DeltaE(x,y,L,spin);
	     
	     if(dE<0){
	     	M[i]= M[i-1]+2*spin[x][y];
			E[i] = E[i-1]+dE;		   		
		 }
		 else{
		 
		 r = rnd();
		 if (r <= exp(-dE/T))
	       {
	       	M[i]= M[i-1]+2*spin[x][y];
			E[i] = E[i-1]+dE;			   	
			}
			else{
				spin[x][y] = -spin[x][y];
				M[i]= M[i-1];
				E[i] = E[i-1];			     			
			}
		}
	//fprintf(fp,"%d\t%f\t%f\n",i,E[i]/(L*L),M[i]/(L*L));   /* for E and M MC history */
	}
}
                





main()
  {
  	//Initialize
	int i,pass,Nskip =10,mcs=10000,nequil=5000,L=23,N =L*L,spin[L+1][L+1];
	float E[mcs+nequil],M[mcs+nequil],T=2;
	
	FILE *fp_init;
	FILE *fp_equil;
	fp_init= fopen("Init_Config.txt","w");
	fp_equil= fopen("Equil_Config.txt","w");
	  
	ColdInitial(mcs+nequil,L,T,E,M,spin);
	//RandomInitial(N,L,T,E,M,spin); /* For Hot Start */
	save_config(N,L,T,spin,fp_init); // Initial Config
	//Thermalization 
	Metropolis(mcs+nequil,L,spin,E,T,M);
	save_config(N,L,T,spin,fp_equil); // Equilibrium Config
	// Discarding inital nequil cycles and Sampling mcs cycles
	data(T,nequil,mcs,Nskip,L,E,M);
}




void data(float T,int Nequil,int mcs,int Nskip,int L,float E[mcs+Nequil],float M[mcs+Nequil])
   {	
		int i,K;
		float Eavg=0.0,Mavg=0.0,E2avg=0.0,M2avg=0.0,HeatCap=0.0;
		FILE *fp;
		fp = fopen("DatvsE.txt","w");
		for(i=Nequil;i<mcs;i++){
			if (i%Nskip==0){
			Eavg+= E[i];
			Mavg+= M[i];
			M2avg+= M[i]*M[i];
			E2avg+= E[i]*E[i];
			K++;
			}
			//printf("%d\t%f\t%f\n",i,E[i]/(L*L),M[i]/(L*L)); 
		}
		
		float norm = (float)1/(L*L*K*Nskip);
		//printf("%d\n",K*Nskip); 
		Eavg = Eavg*norm;
		Mavg = Mavg*norm;
		E2avg = E2avg*norm;
		M2avg = M2avg*norm;
		
		//HeatCap = (1/L*L*T*T)*(E2avg-pow(Eavg,2));
			
	  printf(" \n");
	  printf("Effective Temperature (J/Kb) = %f\n",T);
      printf("Mean energy per spin = %f\n", Eavg);
      printf("Mean squared energy per spin = %1f\n", E2avg);
      printf("Mean magnetization per spin = %f\n",Mavg);
      printf("Mean squared mag. per spin = %f\n",M2avg);
      //printf("Heat Capacity = %.8f\n",HeatCap);
      //fprintf(fp,"%f\t%f\n",)  
		
   }

/* SAMPLE OUTPUT:
Effective Temperature (J/Kb) = 2.000000
Mean energy per spin = -0.513355
Mean squared energy per spin = 1394.213989
Mean magnetization per spin = 0.006511
Mean squared mag. per spin = 0.238087
 */



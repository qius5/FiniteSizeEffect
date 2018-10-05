// Dynamics of Neuro
// Dynamics of Neuro
// For theta model. This code is designed to get the stable solution independent of time.
#include <iostream>
using namespace std;
#include <math.h>
#include <pthread.h>
#include<time.h>
//#include <random>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <stack>
#define PI 3.1415926
//#define Ntot 5000 //assume total number of neuron is 1000
//#define Nsample 10000 // sample number
#define dtheta 2*PI/Ntheta //angle spacing is 2pi/1000
#define Iboot 1 
//#define Ilow -0.3
#define Ilast 10 
#define totime 100.
#define L 1
#define D 0
#define TimeBinSize 0.01
//#define NP 2 
//
double drand()   /* uniform distribution, (0..1] */
{
  return (rand()+1.0)/(RAND_MAX+1.0);
}

typedef struct SynapticDrive_data {
double value0; //the value as a function of z at time t0
double w; //transfer function w
double w0;
double theta;
double count;
double rand;
double I;
double trash;
int Ens;
int Id;
int ID;
} SynapticDrive_data;

typedef struct angle_data{
double theta;
int Id;
double cum;
double pdf;
} angle_data;

typedef struct stack_data{
int value;
} stack_data;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct params 
{ 
double *u; 
double *phase; 
double *theta; 
double *Ilow;
double *fc; //fc is the firing count vector
int tid; 
double t0; 
//size is same as Nsample*Ntot here;
int size;
int nthreads; 
int Ntot;
int Nsample;
double dt;
double beta;
double Ttot;
} params; 

void *compute_parallel(void *pinput)
{ 
double eps=0.0001;
params *p = (params*)pinput;
int tid = p->tid; 
int chunk_size = (p->size / p->nthreads); 
int start = tid * chunk_size; 
int end = start + chunk_size;
for(int i = start; i < end; i++) 
{
   if(p->Ilow[i]+p->u[i]>0)
     {
      double nu;
      nu=p->Ilow[i]+p->u[i];
//      double y;
      //y=0.5;
//      y=PI;

//      double fd;
//      fd=1;
//      std::random_device rd;
//      std::default_random_engine generator(rd());
//      std::uniform_real_distribution<double> distribution(0,1);
       double p1;
//      p1=distribution(generator);
      p1=drand();
//      printf ("The random number is %e\n", p1);
//      while(fd>eps)
//       {
//       double rho;
//       rho = sqrt(nu)/(1-cos(y)+nu*(1+cos(y)));
//       fd = atan(sin(y)/(sqrt(nu)*(1+cos(y))))/PI+.5-p1;
//       y -= fd/rho;
//       }
      p->phase[i]=2*atan(tan((p1-0.5)*PI)*sqrt(nu));
  //  cout << "p->phase[" << i << "] is " << p->phase[i] << endl;
     }
    else
    {
     double nu;
     nu=p->Ilow[i]+p->u[i];
     p->phase[i]=-2*atan(sqrt(nu));
    }
  //p->phase=&temp_phase;

}
//  p->phase=temp_phase;
 
} 

void initialize_phase(int Ntot, int Nsample, int NP, double *u, double *phase, double *Ilow)
{ 
int nthreads = NP; 
int size = Ntot*Nsample; 
pthread_t threads[nthreads]; 
//cout << "here? Inside the function?" << endl;
//array to hold thread information 
params *thread_params;
thread_params=(params*)malloc(nthreads*sizeof(params)); 
for(int i = 0; i < nthreads; i++)
{
//cout << "initializing" << endl; 
thread_params[i].phase = phase; 
thread_params[i].u = u; 
thread_params[i].Ilow = Ilow; 
thread_params[i].tid = i; 
thread_params[i].size = size; 
thread_params[i].nthreads = nthreads;
thread_params[i].Ntot=Ntot;
thread_params[i].Nsample=Nsample;
pthread_create(&threads[i], NULL, compute_parallel, (void*) &thread_params[i]); 
} 
for(int i = 0; i < nthreads; i++)
{
cout << "joining thread " << i << endl;  
pthread_join(threads[i], NULL); 
} 
free(thread_params); 
} 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This part is updating phase and synaptic drive on each tread with some set of samples, using same parameter structure
void *compute_parallel_updatepu(void *pinput)
{
double A1=30;
double A2=20;
double A=1.5;
double b=1;
params *p = (params*)pinput;
double beta=p->beta;
int Ntot=p->Ntot;
int Nsample=p->Nsample;
double t0=p->t0;
double dt=p->dt;
double dz=(double)L/(double)Ntot;
double local_w[Ntot];
double Ttot=p->Ttot;
//printf("\nDebug t0 is %e\n", t0);
for(int i=0; i<Ntot; i++)
{
 local_w[i]=100./(double)L*(A*exp(-A1*(i)*dz*i*dz)-b*exp(-A2*(i)*dz*i*dz)+A*exp(-A1*(Ntot-i)*dz*(Ntot-i)*dz)-b*exp(-A2*(Ntot-i)*dz*(Ntot-i)*dz));
 //local_w[i]=0.2;
 //local_w[i]=1;
 //local_w[i]=0.9+2.4*cos((i-1)*dz*2*PI/(double)L);
// local_w[i]=0.5+0.1*cos((i-1)*dz*2*PI/(double)L);
 //local_w[i]=100./(double)L*(A*exp(-10*A1*(i)*dz)-b*exp(-10*(i)*dz)+A*exp(-10*A1*(Ntot-i)*dz)-b*exp(-10*(Ntot-i)*dz))*dz;
 //local_w[i]=100./(double)L*(A*exp(-10*A1*(i)*dz)-exp(-10*(i)*dz))*dz;
}

double eps=0.0001;
int tid = p->tid;
int chunk_size = (p->size / p->nthreads);
int start = tid * chunk_size;
int end = start + chunk_size;
//define "sn" as the sample id, "spt" as sample per thread. 
int spt;
spt=(int)((double)Nsample/(double)p->nthreads);
for(int sn = tid*spt; sn < (tid+1)*spt; sn++)
{
std::stack<int> s;
int sum;
double sum1;
sum1=0;
int sweep;
sweep=0;
while(sweep*dt<Ttot)
{
//if(tid==0&&sweep%100==1)
//{
//printf("dealing with sample %d, time %e\n", sn, sweep*dt);
//}
sweep+=1;
for(int j=0; j<Ntot; j++)
{
int i;
i=sn*Ntot+j;
//cout << "index is " << i << endl;
//temp_phase[i]=p->phase[i];
// std::random_device rd;
// std::default_random_engine generator(rd());
// std::normal_distribution<double> distribution1(0,1);
double temp1;
  //temp1=dt*(1-cos(p->phase[i])+(1+cos(p->phase[i]))*(Ilow+p->u[i]-D/2.*sin(p->phase[i])))+(1+cos(p->phase[i]))*sqrt(D*dt)*distribution1(generator);
  temp1=dt*(1-cos(p->phase[i])+(1+cos(p->phase[i]))*(p->Ilow[i]+p->u[i]));
//if(sn==10&&j==250)
//{
//printf("\nDebug phase at time %e is %e with temp1 as %e\n", t0+sweep*dt, p->phase[i], temp1);
//}
  p->phase[i]+=temp1;
//if(t0<50&&p->phase[i]>=PI)
//{
//printf("\nDebug Why??\n");
//}
 if((PI-p->phase[i])<=0)
{
 s.push(j);
 p->fc[i]+=1;// adding up when neuron j of sample sn is firing.
 sum1+=1;
//printf("\nFiring, time is %e\n", t0+sweep*dt);
//printf("\nFiring %d %e %e %d\n", sn, sweep*dt, t0, j);
}
while(p->phase[i]>PI)
{
p->phase[i]-=2*PI;
}
while(p->phase[i]<-PI)
{
p->phase[i]+=2*PI;
}
//cout << "phase of " << i << "Neuron is " << p->phase[i] << "on "  << tid << "thread" << endl; 
}
//now we begin update u for specific sample
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
stack_data *data_stack;
data_stack=(stack_data*)malloc(sizeof(stack_data)*(int)sum1);
sum=0;
while(!s.empty())
{
data_stack[(int)sum1-1].value=s.top();
//if(t0<50)
//{
//printf("Debug empty or not?\n");
//}
//printf("\nFiring %d %e %e %d\n", sn, sweep*dt, t0, data_stack[(int)sum1-1].value);
//cout << "now, sum1 is " << sum1 << endl;
//cout << "now, data_stack value is " << s.top() << endl;
s.pop();
sum1--;
sum++;
}
//printf("sum is %d\n", sum);
//printf("place 2\n");

for(int j=0; j<Ntot; j++)
 {
int i=sn*Ntot+j;
p->u[i]-=dt*beta*p->u[i];
for(int k=0; k<sum; k++)
  {
   //cout << "The index is " << data_stack[k].value << endl; 
   if(fabs(data_stack[k].value-j)<Ntot)
   {
   p->u[i]+=beta*(double)L/(double)Ntot*local_w[abs(j-data_stack[k].value)]*exp(-beta*(PI+p->phase[sn*Ntot+abs(data_stack[k].value)])/2.);
   }
  }
 }

free(data_stack);
sum1=0;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
}
}

void update_phase_u(int Ntot, int Nsample, int NP, double *u, double *fc, double *phase, double dt, double Ttot, double beta, double *Ilow, double t0)
{
void *status;
int rc;
int nthreads = NP;
int size = Ntot*Nsample;
pthread_t threads[nthreads];
//array to hold thread information
params *thread_params;
thread_params=(params*)malloc(nthreads*sizeof(params));
for(int i = 0; i < nthreads; i++)
{
//cout << "updating phase and u in thread " << i << endl;
thread_params[i].phase = phase;
thread_params[i].beta = beta;
thread_params[i].u = u;
thread_params[i].fc = fc;
thread_params[i].t0 = t0;
thread_params[i].Ilow = Ilow;
thread_params[i].tid = i;
thread_params[i].size = size;
thread_params[i].nthreads = nthreads;
thread_params[i].Ntot=Ntot;
thread_params[i].Ttot=Ttot;
thread_params[i].dt=dt;
thread_params[i].Nsample=Nsample;
pthread_create(&threads[i], NULL, compute_parallel_updatepu, (void*)&thread_params[i]);
}

//cout << "finishe this step?" << endl;

for(int i = 0; i < nthreads; i++)
{
//cout << "joing thread " << i << endl;
//pthread_join(threads[i], NULL);
rc = pthread_join(threads[i], &status);
      if (rc){
         cout << "Error:unable to join," << rc << endl;
         exit(-1);
      }
}
free(thread_params);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char *argv[]){
if(argc<7)
{
cout << "Usage: " << argv[0] << " -NT <number of thread> -input <input file name> -Ntot <total number of Neurons> -Nsample <sample number> -dt <time step> -beta <beta value> -Iinput <I.dat>" << endl;
exit(0);
}
int Ntot=atoi(argv[3]);
int Nsample=atoi(argv[4]);
int NP=atoi(argv[1]);
double dt=atof(argv[5]);
double beta=atof(argv[6]);
double *Ilow;
Ilow=(double*)malloc(sizeof(double)*Ntot*Nsample);
//for(int j=0; j<Nsample; j++)
//{
//for(int i=0; i<Ntot; i++)
//{
 //Ilow[i+3*Ntot*j]=0.5+(double)i/(double)Ntot*10;
// Ilow[i+Ntot*j]=-0.3;
//}
//}
srand (time(NULL));
//std::random_device rd;
//std::default_random_engine generator(rd());
//std::uniform_real_distribution<double> distribution(0,1);
SynapticDrive_data *data_SD;
SynapticDrive_data *bumpu;
//angle_data *data_angle;
params *thread_params;
data_SD=(SynapticDrive_data*)malloc(sizeof(SynapticDrive_data)*Ntot*Nsample);
bumpu=(SynapticDrive_data*)malloc(sizeof(SynapticDrive_data)*Ntot);
//data_angle=(angle_data*)malloc(sizeof(angle_data)*Ntheta*Ntot*Nsample);
double corr[Ntot];
double corrfc[Ntot];
double corr1[Ntot];
double corrfc1[Ntot];
double corr2[Ntot];
double corr3[Ntot];
double fourmo[Ntot];
double avg[Ntot];
double Ifire=Ilast/dt;
double sweeptot=totime/dt;
double dz=(double)L/(double)Ntot;
double A1=2;
double A=2.5;
//double I=Iboot;
//double I=Ilow;
double gamma=1;
double c0, y0;
double ave[Ntot];
double avefc[Ntot];
double ave1[Ntot];
double avefc1[Ntot];
double mean[Ntot];
int sweep, k, l;
double temp, t1, c;
FILE *fp,*fp5, *fpbump, *fpcorr, *fpvarvar,*fpmean, *fpCorrTimeEvolve, *fpdata, *fpraster, *corrsnap, *corrfcsnap, *avesnap, *avefcsnap, *fpI;
//SynapticDrive_data *data_SD;
stack_data *data_stack;
fp =  fopen("OutFile", "w");
fp5 = fopen("OutFilep2", "w");
corrsnap=fopen("corrsnap","w");
corrfcsnap=fopen("corrfcsnap","w");
avesnap=fopen("avesnap","w");
avefcsnap=fopen("avefcsnap","w");
fpcorr = fopen("corr", "w");
fpvarvar = fopen("varvar", "w");
fpCorrTimeEvolve=fopen("CorrTimeEvolve", "w");
fpbump = fopen(argv[2], "r");
fpI = fopen(argv[7], "r");
fpmean = fopen("mean", "w");
fpdata = fopen("data", "w");
fpraster = fopen("raster", "w");
if (fp == NULL ){
     cout << "Error opening file!!" << endl;
}
int i, j;
double itep;
int timenum;

//here, we initialize the angle
double temp1;

for (i=0; i<Ntot; i++)
    {
    fscanf (fpbump, "%le%*[^\n]%*c", &bumpu[i].value0);
    }
    fclose(fpbump);
for (i=0; i<Ntot; i++)
    {
    fscanf (fpI, "%le%*[^\n]%*c", &Ilow[i]);
    }
for(j=1; j<Nsample; j++)
{
for(i=0; i<Ntot; i++)
 {
 Ilow[j*Ntot+i]=Ilow[i];
 }
}

    fclose(fpI);
//just make sure thesea re defined outside the loop. 
double *tempu;
double *tempphase;
double *tempfc;
tempu=(double*)malloc(sizeof(double)*Ntot*Nsample);
tempfc=(double*)malloc(sizeof(double)*Ntot*Nsample);
tempphase=(double*)malloc(sizeof(double)*Ntot*Nsample);
for(j=0; j<Nsample; j++)
{
for (i=0; i<Ntot; i++)
 {
// data_SD[j*Ntot+i].value0=bumpu[i].value0;
 tempu[j*Ntot+i]=bumpu[i].value0;
 tempphase[j*Ntot+i]=0;
 tempfc[j*Ntot+i]=0;
//   data_SD[j*Ntot+i].w0=100*(A*exp(-10*A1*(i)*dz)-exp(-10*(i)*dz)+A*exp(-10*A1*(Ntot-i)*dz)-exp(-10*(Ntot-i)*dz));
 }
}

cout << "here?" << endl;
initialize_phase(Ntot, Nsample, NP, tempu, tempphase, Ilow);


printf("Now the data at time 0 is set up.\n");
sweep=0;
int sa;
double count=0.0;
double t0;
int rad=(int)((double)Ntot/100./2.);
double tempvalue;
for(int snap=0; snap<(int)(10./TimeBinSize); snap+=1)
{
t0=(double)snap*TimeBinSize;
//printf("Debug main t0=%e\n", t0);

update_phase_u(Ntot, Nsample, NP, tempu, tempfc, tempphase, dt, TimeBinSize, beta, Ilow, t0);


cout << "Finished multithread updating and doing snapshot " << snap  << endl;

for(i=0; i<Ntot; i++)
{
 ave1[i]=0;
 corr1[i]=0;
 avefc1[i]=0;
 corrfc1[i]=0;
}
for(i=0; i<Ntot; i++)
{
ave[i]=0;
avefc[i]=0;
for(j=0; j<Nsample; j++)
{
ave[i]+=tempu[j*Ntot+i];
for(int in=0; in<rad; in++)
 {
 if(i-in>0&&i+in<Ntot)
 avefc[i]+=tempfc[j*Ntot+i-in]/(double)rad/2.+tempfc[j*Ntot+i+in]/(double)rad/2.;
 if(i-in<0)
 avefc[i]+=tempfc[j*Ntot+i-in+Ntot]/(double)rad/2.+tempfc[j*Ntot+i+in]/(double)rad/2.;
 if(i+in>Ntot)
 avefc[i]+=tempfc[j*Ntot+i-in]/(double)rad/2.+tempfc[j*Ntot+i+in-Ntot]/(double)rad/2.;
 }
}
ave[i]/=(double)Nsample;
avefc[i]/=(double)Nsample;
ave1[i]+=ave[i];
avefc1[i]+=avefc[i];
}


for(i=0; i<Ntot; i++)
{
corr[i]=0;
corrfc[i]=0;
for(j=0; j<Nsample; j++)
{
corr[i]+=(tempu[j*Ntot+i]-ave[i])*(tempu[j*Ntot+i]-ave[i]);
//reset the temporal variabel to be zero
tempvalue=0;
for(int in=0; in<rad; in++)
 {
 if(i-in>0&&i+in<Ntot)
 tempvalue+=tempfc[j*Ntot+i-in]/(double)rad/2.+tempfc[j*Ntot+i+in]/(double)rad/2.;
 if(i-in<0)
 tempvalue+=tempfc[j*Ntot+i-in+Ntot]/(double)rad/2.+tempfc[j*Ntot+i+in]/(double)rad/2.;
 if(i+in>Ntot)
 tempvalue+=tempfc[j*Ntot+i-in]/(double)rad/2.+tempfc[j*Ntot+i+in-Ntot]/(double)rad/2.;
 }
tempvalue-=avefc[i];
tempvalue*=tempvalue;
//corrfc[i]+=(tempfc[j*Ntot+i]-avefc[i])*(tempfc[j*Ntot+i]-avefc[i]);
 corrfc[i]+=tempvalue;
}
corr[i]/=((double)Nsample-1);
corrfc[i]/=((double)Nsample-1);
corr1[i]+=corr[i];
corrfc1[i]+=corrfc[i];
}

for(i=0; i<Ntot; i++)
{
fprintf(corrsnap, "%e %e\n", (double)i/(double)Ntot, corr[i]);
fprintf(avesnap, "%e %e\n", (double)i/(double)Ntot, ave[i]);
fprintf(corrfcsnap, "%e %e\n", (double)i/(double)Ntot, corrfc[i]);
fprintf(avefcsnap, "%e %e\n", (double)i/(double)Ntot, avefc[i]);
}

for(i=0; i<Ntot; i++)
{
 for(j=0; j<Nsample; j++)
  {
  tempfc[j*Ntot+i]=0; // need to return zero and get ready for measurement of next time interval.  
  }
}

}

count+=1.0;
//for(int cont=0; cont<100; cont++)
//{
//t0=10+cont*10./100;
//update_phase_u(Ntot, Nsample, NP, tempu, tempfc, tempphase, dt, 10./100., beta, Ilow, t0);



/*for(i=0; i<Ntot; i++)
{
ave[i]=0;
for(j=0; j<Nsample; j++)
{
ave[i]+=tempu[j*Ntot+i];
}
ave[i]/=(double)Nsample;
ave1[i]+=ave[i];
}
*/
/*
for(i=0; i<Ntot; i++)
{
corr[i]=0;
fourmo[i]=0;
for(j=0; j<Nsample; j++)
{
corr[i]+=(tempu[j*Ntot+i]-ave[i])*(tempu[j*Ntot+i]-ave[i]);
fourmo[i]+=(tempu[j*Ntot+i]-ave[i])*(tempu[j*Ntot+i]-ave[i])*(tempu[j*Ntot+i]-ave[i])*(tempu[j*Ntot+i]-ave[i]);
}
corr[i]/=((double)Nsample-1);
fourmo[i]/=(double)Nsample;
corr1[i]+=corr[i];
corr2[i]+=fourmo[i];
}
*/
//count+=1;

//}

/*for(i=0; i<Ntot; i++)
{
printf("%e %e\n", corr1[i], count);
ave1[i]/=count;
corr1[i]/=count;
corr2[i]/=count;
double tempNsample;
tempNsample=Nsample;
corr3[i]=(corr2[i]-(tempNsample-3)/(tempNsample-1)*corr1[i]*corr1[i])/tempNsample;
}
*/

// fprintf(fpCorrTimeEvolve, "%e %e\n", sweep*dt, corr[(int)((double)Ntot/2.)]);


//This is the end of time loop

//Here we print the corr and mean at the last time step, which is steady state.
/*for(i=0; i<Ntot; i++)
{
fprintf(fpcorr, "%e %e\n", i*dz, Ntot/(double)L*corr1[i]);
fprintf(fpmean, "%e %e\n", i*dz, ave1[i]);
fprintf(fpvarvar, "%e %e\n", i*dz, Ntot/(double)L*Ntot/(double)L*corr3[i]);
}
*/
for(sa=0; sa<Nsample; sa++)
{
{
 fprintf(fp, "%d %e %e\n", sa, sweep*dt, tempphase[sa*Ntot+(int)((double)Ntot/2.)]);
}
{
 fprintf(fp5, "%d %e %e\n", sa, sweep*dt, tempphase[sa*Ntot+3]);
}
}

for(sa=0; sa<Nsample; sa++) {
 for(i=0; i<Ntot; i++)
  {
  fprintf(fpdata, "%e ",tempu[sa*Ntot+i+Ntot]);
  }
  fprintf(fpdata, "\n");
}


free(data_SD);
fclose(fp);
fclose(fp5);
fclose(fpmean);
fclose(fpcorr);
fclose(fpvarvar);
fclose(corrsnap);
fclose(avesnap);
fclose(fpdata);
return 0;
}

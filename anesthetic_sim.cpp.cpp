//============================================================================
// Name        : atontransfer.cpp
// Author      : Bolaji
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

// so the only thing imay need to change is to decrease the inh firing rate/increase connectivity
// Also I need to increase the total input so there are no 1.5 gks neurons that dont fire
//
#include <iostream>
using namespace std;
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <time.h>
#include <string>
#include <cstdio>
#include <sstream>
//#include <array>
#include <iomanip>
#include <vector>
#include <stdlib.h>
#include <algorithm>
#include <stdio.h>
   // #include "matrixsize.h"
#include <random>
#define pi 3.14159
#include "network_dyn.hpp"


struct sim1 {
	double Vv;
	//double spiketimesdyn;
	double mv;
	double nv;
	double sv;
	//double zv;
	double hv;
	//double isyn;
};



//////////////////////gating
double hinf(double V){
	double x;
	x=pow((1+exp((V+53.0)/(7.0) )),(-1));
	return x;
};

double tauh(double V){
	double x;
		x= 0.37 + 2.78*pow((1 + exp( (V+40.5)/(6.0) ) ),(-1));
		return x;
};

double ninf(double V) {
	double x;
	x=pow((1+exp( - (V+30.0)/10.0 )),(-1));
	return x;
};
double minf(double V ){
	double x;
	x=pow( (1+exp( ((-V-30.0)/9.5 ) )),(-1));
	return x;
}

double taun(double V) {
	double x;
	x=0.37 + 1.85*pow((1 + exp( (V+27.0)/(15.0) ) ),(-1));
	return x;
};

double sinf(double V) {
	double x;
	x=pow((1+exp( - (V+39.0)/5.0 )),(-1));
	return x;
};

double taus(double V ){
	return 75.0;
}

double an (double V) {
	double x;
	x=0.01*(V+10.0)/(exp( (V+10.0)/10.0 ) -1);
	return x;
};

double bn (double V){

	double x;
	x= 0.125*exp(V/80.0);
	return x;
};

double am (double V){
	double x;
	x=0.1*(V+25.0)/(exp( (V+25.0)/10.0 ) -1);
	return x;
};

double bm (double V){
	double x;
	x=4*exp(V/18.0);
	return x;
};

double ah(double V){
	double x;
	x= 0.07*exp(V/20.0);
	return x;
};
double bh(double V){
	double x;
	x= 1/(exp((V+30.0)/10.0) +1);
	return x;
};

/////////////////////////////////////////////////////////////////////////////////////////////////


	double hd(double V,double hv){
		double x;
		x=(1.0/tauh(V))*(hinf(V)-hv);
		return x;
	}
	double nd(double V,double nv){
		double x;
		x=(1.0/taun(V))*(ninf(V)-nv);
		return x; }
	double sd(double V,double sv){
		double x;
		x=(1.0/taus(V))*(sinf(V)-sv);
		return x;
	}

	double md(double V,double mv){
		double x;
		x=0;
		return x;
	}
	double Vd(double V, double gna,double minf,double hv,double sv,double nv,double Vna,double Vk, double Vl, double gkdr,double gks,double gl,double Iext,double Isyn,double Id, double C){
		double Vout;
		Vout= (1.0/C)*(-gna*pow(minf,3)*hv*(V-Vna)-gkdr*pow(nv,4)*(V-Vk)-gks*sv*(V-Vk)-gl*(V-Vl)+Iext-Isyn +Id);
		return Vout;
	}

//gatingpm1 gatecomp;


sim1 V_integrate(double V1,double mv1, double nv1, double sv1, double hv1, double isyn1,double Iext,double Id,double dt,double gks){

	sim1 calcout;

	double C=1.0; // uF/cm^2
	double Is=0; //
	//double Id=2;

	double In=0;
	double Vna=55.0; //% mV
	double Vk= -90.0;// % mV
	double Vl=-60.0; //% mv
	double gkdr=3.0; //%mS/cm^2
	double gna= 24.0;// % mS/cm^2
	double gl=0.02; //% mS/cm^2


	double Vk1, hk1, nk1, sk1, halfdt, sixthdt;
	double Vk2, hk2, nk2, sk2;
	double Vk3, hk3, nk3, sk3;
	double Vk4, hk4, nk4, sk4;
	halfdt=(dt/2.0);
	sixthdt=(dt/6.0);


	//double dt=.05;


	Vk1=Vd(V1,gna,mv1,hv1,sv1,nv1,Vna,Vk,Vl,gkdr,gks,gl,Iext,isyn1,Id,C);
	hk1=hd(V1,hv1);
	nk1=nd(V1,nv1);
	sk1=sd(V1,sv1);


	Vk2=Vd(V1+halfdt*Vk1,gna,mv1,hv1+halfdt*hk1,sv1+halfdt*sk1,nv1+halfdt*nk1,Vna,Vk,Vl,gkdr,gks,gl,Iext,isyn1,Id,C);
	hk2=hd(V1+halfdt*Vk1,hv1+halfdt*hk1);
	nk2=nd(V1+halfdt*Vk1,nv1+halfdt*nk1);
	sk2=sd(V1+halfdt*Vk1,sv1+halfdt*sk1);

	Vk3=Vd(V1+halfdt*Vk2,gna,mv1,hv1+halfdt*hk2,sv1+halfdt*sk2,nv1+halfdt*nk2,Vna,Vk,Vl,gkdr,gks,gl,Iext,isyn1,Id,C);
	hk3=hd(V1+halfdt*Vk2,hv1+halfdt*hk2);
	nk3=nd(V1+halfdt*Vk2,nv1+halfdt*nk2);
	sk3=sd(V1+halfdt*Vk2,sv1+halfdt*sk2);

	Vk4=Vd(V1+dt*Vk3,gna,mv1,hv1+dt*hk3,sv1+dt*sk3,nv1+dt*nk3,Vna,Vk,Vl,gkdr,gks,gl,Iext,isyn1,Id,C);
	hk4=hd(V1+dt*Vk3,hv1+dt*hk3);
	nk4=nd(V1+dt*Vk3,nv1+dt*nk3);
	sk4=sd(V1+dt*Vk3,sv1+dt*sk3);


    //calcout = V_integrate(V[i],m[i],n[i],s[i],h[i],isyn[i],Iex );


	calcout.Vv= V1+ sixthdt*(Vk1 + 2.0*Vk2+2.0*Vk3+Vk4);
    calcout.hv=hv1+ sixthdt*(hk1 + 2.0*hk2+2.0*hk3+hk4);
    calcout.mv=minf(V1);//md(V1,mv1)*dt+mv1;
    calcout.nv=nv1+ sixthdt*(nk1 + 2.0*nk2+2.0*nk3+nk4);
    calcout.sv=sv1+ sixthdt*(sk1 + 2.0*sk2+2.0*sk3+sk4);
	return calcout;
}







///////////////////////////////////////////////////////////////////////////////////////////////
struct gatingpm2{
	double dt;
	double h(double V,double hv){
		double x;
		x=(ah(V)*(1-hv) -bh(V)*hv  )*dt + hv;
					return x;
		};
		double n(double V,double nv){
			double x;
			x= (an(V)*(1-nv) -bn(V)*nv  )*dt + nv;
			return x; };
		double s(double V,double sv){
			double x;
			x=0;
			return x; };

		double m(double V,double mv){
			double x;
			x=(am(V)*(1-mv) -bm(V)*mv  )*dt + mv;
			return x; };
};



void resetarray(int array12[],int n){

	for(int i=0;i<n;i++){
		array12[i]=0;
	}
}





int main(int argc, char* argv[]) {
        
        
        
	    
	    double pconee=atof(argv[1]);//.3;//.06;//.3
	    double pconii= atof(argv[2]);//.3;		
	    double pconie= atof(argv[3]);//.5;	
	    double pconei=atof(argv[4]);//.5;//.06;//.5;	
	
			
	    double wee_mult    =atof(argv[5]); ;// atof(argv[1]);
	    double wii_mult  =atof(argv[6]);;//atof(argv[2]);
	    double wie_mult  =atof(argv[7]); ; //atof(argv[3])x;
	    double wei_mult   =atof(argv[8]);; //atof(argv[4]);
	    double iE   =atof(argv[9]); ;//atof(argv[5]);
	    double iI    = atof(argv[10]);;//atof(argv[6]);
	    double gksexc = atof(argv[11]);	
	    double stim_val = atof(argv[12]);	
	    double w_mult = atof(argv[13]);
	    double wee_nmda_mult = atof(argv[14]);
	    double wei_nmda_mult = atof(argv[15]);
	    double log_norm_mu= atof(argv[16]);
	    double log_norm_std= atof(argv[17]);
	    double noise_freq=atof(argv[18]);
	    double noise_strength=atof(argv[19]);	    
	    double noise_dur=atof(argv[20]);	
	    double my_net=atof(argv[21]);
	    double my_runtime=atof(argv[22]);		
		double randee;
		
		double randii;
	    
		double randie;
		
		double randei;
		ofstream stimtext("stim.txt");
	    
        
	      std::default_random_engine generator;
	    //  std::lognormal_distribution<double> distribution(0.0,.9);
	      std::lognormal_distribution<double> distribution(log_norm_mu,log_norm_std);
	      

	     

	    int xl=1*time(NULL);

	    
	   /// cout<< argc<<endl;
	    if (argc != 23)
	    {
	    	cout << " wrong number of input arguments"<< endl;
	    	return -1;
	    }
	    
	    ofstream pfile("params.txt");
	    //pfile<< gsynee <<endl;
	    //pfile<< gsynii <<endl;
	    //pfile<< gsynie <<endl;
	    //pfile<< gsynei <<endl;
	    pfile<< iE <<endl; // super thresh=1.18	
	    pfile<< iI <<endl;	
	    pfile.close();
	    
	    
	    //double wie    = atof(argv[7]);
	    //double tauCaSlow = atof(argv[8]);
	    //double rmax = atof(argv[9]);

	    
	  /*  double gks    = atof(argv[1]);
	      double probE  = atof(argv[2]);
	     i double probI  = atof(argv[3]);
	      double rEEm   = atof(argv[4]);
	      double rIIm   = atof(argv[5]);
	      double wei    = atof(argv[6]);
	      double wie    = atof(argv[7]);
	      double tauCaSlow = atof(argv[8]);
	      double rmax = atof(argv[9]);*/
	    
	srand (my_net);
	Ran ICRand(xl);
	double dt=0.05; // ms
//	double runtime=20000;//500;// 2000.0; // ms
	//double runtime=70000;
	double noise_prob=(noise_freq)/(1000.0/dt);
	cout<<" noise_prob" <<noise_prob<< endl;
    double runtime=my_runtime;
	int ninh=200;
	int nexc=800;
	int ntot=ninh+nexc;
	double tauS_nmda= 200.0;
	double tauF_nmda= 0.3;
//double stim_val = 3.0;
//double t_stim= 0.75*runtime;
	double t_stim=2000.0;
	int t_stim_it= int(t_stim/dt);
	cout <<"t_stim_it" <<t_stim_it<< endl;


int stim_vect[nexc];
//int high_vec_weight[nexc];
    double stim_ratio=0.5;
	//double gsynee=0.0000625;
	//double gsynii=0.00025;	
	//double gsynie=0.00175;//mS/cm^2
	//double gsynei=0.00175;//mS/cm^2
		double fracei=1;
	double fracii=1;
	double fracie=2;
	double fracee=1;
	

	
	double inhoffset;
	double excoffset;
	
	//gsynee=gsynee/fracee;
	//gsynii=gsynii/fracii;
	//gsynie=gsynie/fracie;
	//gsynei=gsynei/fracei;
	//inhoffset=0;//.5;//.1
	//excoffset=0;//-4;//
	
	double simlength=runtime/dt;
	int iterj=0;
	double t;
	double rando2;
	double gks;
	double gksinh;
//	double gksexc;
	
	int inhtype=1;
	int exctype=2;
	double dcmaxexc;
	double dcminexc;
	double ia;
	double idcexc_mod=0.0;

	
	
	bool lognorm=false;
	if (inhtype==1){ 
			gksinh=0;
			ia=iI;
			//dcmaxexc=.68;
			//dcminexc=.46;
			//return -1;
			}
			else{ //gks=1.5;
			ia=iI;//ia=1.0//+inhoffset;
			gksinh=1.5;
			 //dcmaxexc=1.0*iE//8.68;//11.1+excoffset;
			 //dcminexc=0.9*iE;//8.7;//8.68+excoffset;
			}//1.5;
	
	/*
	if (exctype==1){
				gksexc=0;
				//gks=0;
				//ia=-0.2;
				//dcmaxexc=.68;
				//dcminexc=.46;
				//return -1;
				}
				else{ //gks=1.5;
				ia=iI;//ia=1.0//+inhoffset;
				gksexc=1.5;
				// dcmaxexc=1.0*iE;//8.68;//11.1+excoffset;
				 //dcminexc=0.9*iE;//8.7;//8.68+excoffset;
				}//1.5;
	*/
	
	
	/*if (type==1){ 
		//gks=0;
		//ia=-0.2;
		//dcmaxexc=.68;
		//dcminexc=.46;
		return -1;
		}
		else{ gks=1.5;
		ia=iI;//ia=1.0//+inhoffset;

		 dcmaxexc=1.0*iE;//8.68;//11.1+excoffset;
		 dcminexc=0.9*iE;//8.7;//8.68+excoffset;
		}//1.5;*/
	
	
	//gks=1.5;
	//ia=1.0;
	double idcexc[nexc];
	 //dcmaxexc=.68;
	// dcminexc=.46;
	 dcmaxexc=1*iE;//1.0*iE;//8.68;//11.1+excoffset;
		 dcminexc=.3*iE;//0.9*iE;
	double dcincrementexc=(dcmaxexc-dcminexc)/(nexc*1.0);
	 //dcmaxexc=1.2*iE;//1.0*iE;//8.68;//11.1+excoffset;
	 //dcminexc=.0*iE;//0.9*iE;
	double idcinh[ninh];
	double dcmaxinh=1.0*iI;//1.1*iI;
	double dcmininh=.9*iI;//.9*iI;
	double dcincrementinh=(dcmaxinh-dcmininh)/(ninh*1.0);

    

	
	double Esyni=-75.0;
	double Esyne=0.0;



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////end original w and strogatz

ofstream connectivityeetext("connectivityee.txt");
ofstream connectivityiitext("connectivityii.txt");
ofstream connectivityeitext("connectivityei.txt");
ofstream connectivityietext("connectivityie.txt");	
ofstream wietext("wie.txt");
ofstream weetext("wee.txt");
ofstream wiitext("wii.txt");
ofstream weitext("wei.txt");
ofstream weenmdatext("weenmda.txt");
ofstream weinmdatext("weinmda.txt");
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////begin ee strogatz


//double wee[nexc][nexc];
//	for (int i=0;i<nexc;i++ ){ for(int j=0;j<nexc;j++){ wee[i][j]=wmax/2;}; };// sets these to wmax/2
//prestrogatz
	



for(int i=0; i<(nexc); i++)
 //   for(int i=0; i<(num_ex); i++)
{
	if( i<(nexc/2.0))
{
    stim_vect[i]=1;
}
else{
	stim_vect[i]=0;
}
    
}





randomize(stim_vect,nexc,ICRand);

for(int i=0; i<nexc; i++)
{
    
stimtext <<  i <<" "<< stim_vect[i] << endl;

}

//--------------------------------------------------------------------

double **wie = new double*[nexc];

for(int i=0; i<nexc; i++)
{
wie[i]= new double[ninh];
}

double **wee = new double*[nexc];

for(int i=0; i<nexc; i++)
{
wee[i]= new double[nexc];
}


double **wei = new double*[ninh];

for(int i=0; i<ninh; i++)
{
wei[i]= new double[nexc];
}

double **wii = new double*[ninh];

for(int i=0; i<ninh; i++)
{
wii[i]= new double[ninh];
}


double **wee_nmda = new double*[nexc];

for(int i=0; i<nexc; i++)
{
wee_nmda[i]= new double[nexc];
}


double **wei_nmda = new double*[ninh];

for(int i=0; i<ninh; i++)
{
wei_nmda[i]= new double[nexc];
}

//--------------------------------------------------------------------






bool rewire;
int preconnectivityee[nexc][nexc];
	for(int i=0;i<nexc;i++)
	{
		for(int j=0;j<nexc;j++)
		{
			//cout<< randee<<endl;
			randee=((double)rand()/(RAND_MAX));
			rewire=randee<pconee;
			if ((rewire))
			{
				
			
				preconnectivityee[i][j]=1;
			
				
			}else
			{
				preconnectivityee[i][j]=0;
			}
			
		}
	}
	

	

int connectivityee[nexc][nexc];
	for(int i=0; i<nexc;i++)
	{
		for(int j=0; j<nexc;j++)
		{
			double log_norm_val = distribution(generator);
			if(log_norm_val>0.9)
                        {
                                log_norm_val=0.0;
                        }

			//cout<<preconnectivityee[j][i]<<"hi"<<endl;
			connectivityee[i][j]= preconnectivityee[j][i];
			connectivityeetext<<connectivityee[i][j];
			connectivityeetext<< " ";
			wee[i][j]=wee_mult*log_norm_val*connectivityee[i][j];
			weetext << wee[i][j];
			weetext << " ";


			wee_nmda[i][j]=wee_nmda_mult*log_norm_val*connectivityee[i][j];

			weenmdatext<< wee_nmda[i][j] << endl;
			weenmdatext<< " ";

		}
		connectivityeetext<<endl;
		weetext<<endl;
		weenmdatext << endl;
	}







/////////////////////////////////////////////////////////////////////////////////////////////////////////////end ee strogatz
////////////////////////////////////////////////////////////////////////////////////////////////////////////begin ii strogatz

//double wii[ninh][ninh];
	//for (int i=0;i<ninh;i++ ){ for(int j=0;j<ninh;j++){ wii[i][j]=wmax/2;}; };// sets these to wmax/2
//prestrogatz
	
//cout<< "geee 536"<<endl;
int preconnectivityii[ninh][ninh];
	for(int i=0;i<ninh;i++)
	{
		for(int j=0;j<ninh;j++)
		{
			randii=((double)rand()/(RAND_MAX));
			if ((randii<pconii))
			{
				preconnectivityii[i][j]=1;
			}else
			{
				preconnectivityii[i][j]=0;
			}

		}
	}

	
	
	int connectivityii[ninh][ninh];
	for(int i=0; i<ninh;i++)
	{
		for(int j=0; j<ninh;j++)
		{	
			 double log_norm_val = distribution(generator);
			if(log_norm_val>0.9)
			{
				log_norm_val=0.0;
			}	

			connectivityii[i][j]= preconnectivityii[j][i];
			connectivityiitext<<connectivityii[i][j];
			connectivityiitext<< " ";
			wii[i][j]=wii_mult*log_norm_val*connectivityii[i][j];
			wiitext<< wii[i][j];
			wiitext<< " ";
		}
		connectivityiitext<<endl;
		wiitext<<endl;
	}
	


	

//////////////////////////////////////////////////////////////////////////////////////////////////////////////end ii strogtaz
//////////////////////////////////////////////////////////////////////////////////////////////////////////////begin ei strogatz

// weight
//double wei[ninh][nexc];
//	for (int i=0;i<ninh;i++ ){ for(int j=0;j<nexc;j++){ wei[i][j]=wmax/2;}; };// sets these to wmax/2
//prestrogatz
	

int preconnectivityei[nexc][ninh];
	for(int i=0;i<nexc;i++)
	{
		for(int j=0;j<ninh;j++)
		{
			randei=((double)rand()/(RAND_MAX));
			if ((randei<pconei))
			{
				preconnectivityei[i][j]=1;
			}else
			{
				preconnectivityei[i][j]=0;
			}

		}
	}

int connectivityei[ninh][nexc];
	for(int i=0; i<ninh;i++)
	{
		for(int j=0; j<nexc;j++)
		{
			double log_norm_val= distribution(generator);
			if(log_norm_val>0.9)
                        {
                                log_norm_val=0.0;
                        }

			connectivityei[i][j]= preconnectivityei[j][i];
			connectivityeitext<<connectivityei[i][j];
			connectivityeitext<< " ";

			wei[i][j]= wei_mult*log_norm_val*connectivityei[i][j];
			weitext<< wei[i][j];
			weitext<< " ";

			wei_nmda[i][j]=wei_nmda_mult*log_norm_val*connectivityei[i][j];
			weinmdatext<< wei_nmda[i][j];
			weinmdatext<< " ";
		}
		connectivityeitext<<endl;
		weitext<<endl;
		weinmdatext<<endl;
	}




//cout<< "geee 649"<<endl;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////end ei strogatz
//////////////////////////////////////////////////////////////////////////////////////////////////////////////begin ie strogatz

// weight
//double wie[nexc][ninh];
	//for (int i=0;i<nexc;i++ ){ for(int j=0;j<ninh;j++){ wie[i][j]=wmax/2;}; };// sets these to wmax/2
//prestrogatz
	

int preconnectivityie[ninh][nexc];
	for(int i=0;i<ninh;i++)
	{
		for(int j=0;j<nexc;j++)
		{
			randie=((double)rand()/(RAND_MAX));
			if ((randie<pconie))
			{
				preconnectivityie[i][j]=1;
			}else
			{
				preconnectivityie[i][j]=0;
			}

		}
	}

int connectivityie[nexc][ninh];

	for(int i=0; i<nexc;i++)
	{
		for(int j=0; j<ninh;j++)
		{
			double log_norm_val =distribution(generator);
			if(log_norm_val>0.9)
                        {
                                log_norm_val=0.0;
                        }

			
			connectivityie[i][j]= preconnectivityie[j][i];
			connectivityietext<<connectivityie[i][j];
			connectivityietext<< " ";
			wie[i][j]=log_norm_val*wie_mult*connectivityie[i][j];
			wietext<< wie[i][j];
			wietext<< " ";
		}
		connectivityietext<<endl;
		wietext<<endl;
	}

//cout<< "geee 683"<<endl;








//
//
//








/////////////////////////////////////////////////////////////////////////////////////////////////////////////end ie strogatz
// make prestrogatz and then transpose to make this the implementation strogatz



/////////////////////////////////////////////network/voltage


double Vexc[nexc];//double V[ntot];
double hexc[nexc];//double h[ntot];
double nnexc[nexc];//double n[ntot];
double sexc[nexc];//double s[ntot];
double mexc[nexc];//double m[ntot];
double isynexc[nexc];//double isyn[ntot];
double Vinh[ninh];
double hinh[ninh];
double nninh[ninh];
double sinh[ninh];
double minh[ninh];
double isyninh[ninh];
ofstream vtextexc("vexc.txt");
ofstream rastertextexc("raster_datexc.txt");
ofstream mtextexc("mgateexc.txt");
ofstream ntextexc("ngateexc.txt");
ofstream stextexc("sgateexc.txt");
ofstream htextexc("hgateexc.txt");
ofstream isyntextexc("syncurrexc.txt");

ofstream vtextinh("vinh.txt");
ofstream rastertextinh("raster_datinh.txt");;
ofstream mtextinh("mgateinh.txt");
ofstream ntextinh("ngateinh.txt");
ofstream stextinh("sgateinh.txt");
ofstream htextinh("hgateinh.txt");
ofstream isyntextinh("syncurrinh.txt");
ofstream testtext("test.txt");



//cout<< "geee 731"<<endl;



int spikenumexc=0;
int spikenuminh=0;
sim1 calcoutexc;
sim1 calcoutinh;
for(int i=0;i<ninh;i++){

	idcinh[i]=dcmininh+(dcincrementinh)*i;
	//cout<<idcinh[i]<<endl;
}


if(lognorm==true)
{
for (int i=0; i<nexc; ++i) {
   idcexc[i] = distribution(generator);
  //cout<< idcexc[i]<<endl;
} 


//sort(idcexc, idcexc+nexc);
for (int i=0; i<nexc; ++i) {
   
  cout<< idcexc[i]<<endl;
}  
}
else{
	
for (int i=0; i<nexc; i++)
{
	//idcexc[i]=dcminexc+(dcincrementexc)*i;
	idcexc[i]=dcminexc+ ((double)rand()/(RAND_MAX))*(dcmaxexc) ; //rand(dcincrementexc)*i;
		//cout<<idcexc[i]<<endl;
}

}


//double spiketimes[ntot];
	//for(int i=0;i<ntot;i++)
	//{spiketimes[i]=0; }

double spiketimesexc[nexc];
double spiketimesinh[ninh];
for(int i=0;i<nexc;i++)
	{ spiketimesexc[i]=0;}
for(int i=0;i<ninh;i++)
	{ spiketimesinh[i]=0;}




double recentconnectedspiketimes[ntot];

double randexcV;
double randexcn;
double randexcs;
double randexcm;
double randexch;

double randinhV;
double randinhn;
double randinhs;
double randinhm;
double randinhh;
double isyne;
double isyni;
double gsynee_nmda_mod;
double B_NMDA;
double V_T;

for(int i=0;i<nexc;i++)
{

	randexcV=((double)rand()/((double)(RAND_MAX)));
	randexcn=((double)rand()/((double)(RAND_MAX)));
	randexcs=((double)rand()/((double)(RAND_MAX)));
	randexcm=((double)rand()/((double)(RAND_MAX)));
	randexch=((double)rand()/((double)(RAND_MAX)));
 Vexc[i]=-72+randexcV*40;
 nnexc[i]=.2 +randexcn*.4;
 sexc[i]=0.2+randexcs*.1;
 mexc[i]=0;
 hexc[i]=.2+randexch*.4;
}

for(int i=0;i<ninh;i++)
{
	randinhV=((double)rand()/((double)(RAND_MAX)));
	randinhn=((double)rand()/((double)(RAND_MAX)));
	randinhs=((double)rand()/((double)(RAND_MAX)));
	randinhm=((double)rand()/((double)(RAND_MAX)));
	randinhh=((double)rand()/((double)(RAND_MAX)));

 Vinh[i]=-72+randinhV*20;
 nninh[i]=.2+ randinhn*.1;
 sinh[i]=.15+randinhs*.1;
 minh[i]=0;
 hinh[i]=.2+randinhh*.1;
}

double gee;
double gie;
double gei;
double gii;
double idcinh_mod=0.0;
double gksinh_mod=0.0;
//Iex=Iext( t, idc, ampdc,Aosc, fosc);
double noise_time_last[ntot];
double noise_rand;
double noise_val;
double spike_rand_pre;
//double spike_pre_last[ntot][ntot];

double **spike_pre_last = new double*[ntot];

for(int i=0; i<ntot; i++)
{
spike_pre_last[i]= new double[ntot];
}

//time for loop j////////////////////////////////////////////////////////////////
for(int j=0;j<simlength;j++){
	t=(j/simlength)*runtime;
		if (j%1000==0){
			//cout<<(j/simlength)*runtime;
			//cout<<iterj;
			cout<<t<<endl;
			//cout<<spikenumexc<<endl;
		iterj++;
		}
//////////////////////////////////////////////////////////////////////////////////////////begin i increment of loop
	for(int i=0;i<ntot;i++)
	{
///////////////////////////////////////so bolaji add a conditional for this at different values
 noise_rand=((double)rand()/((double)(RAND_MAX)));
 	if(noise_rand<noise_prob)
	{

	noise_time_last[i]=t;


	}
	
//cout<<"got 1005"<<endl;
/*cout<<"got 1005"<endl;
 *  *if(noise_rand<prob_noise)
 *          { noise_time[i]=t;}
 *
 *           if ((t-noise_time[i])>noise_dur)
 *           {noise_val=noise_strength}
 *           else{noise_val=0}
 *
 *
 *           */

	for(int k=0;k<ntot;k++){
		
		
//	spike_rand_pre=((double)rand()/((double)(RAND_MAX)))*0.5 +0.5;	
//cout<<"got 1021"<<endl;
	//if for ee
	if( (i<nexc) && (k<nexc) )
	{
		if(connectivityee[i][k]){
		//isynexc[i]=isynexc[i]+gsynee*(  (exp(- (t-spiketimesexc[k])/(3.0) ) - exp(- (t-spiketimesexc[k])/(0.2) ) ) )*(Vexc[i]-Esyne);
			
	/*	if ( Vexc[i]>(-40))
		{
		gsynee_nmda_mod=gsynee_nmda;	
		}
		else
		{
		gsynee_nmda_mod=0;	
		}*/
		
		B_NMDA= (1.0)/( 1.0+exp(( -Vexc[i]-10.0  )/(16.13))  );

		/*
			if(k>(nexc-10))
			{	
			gee=gsynee*w_mult;
			wee[i][k]=gee;
			}
			else
			{

				gee=gsynee;
				wee[i][k]=gee;
			}*/

		//double tauS_nmda= 200.0;
		//double tauF_nmda= 20.0;
		isynexc[i]=isynexc[i]+spike_pre_last[i][k]*wee[i][k]*(  (exp(- (t-spiketimesexc[k])/(3.0) ) - exp(- (t-spiketimesexc[k])/(0.2) ) ) )*(Vexc[i]-Esyne) 
		+  wee_nmda[i][k]*B_NMDA*(  (exp(- (t-spiketimesexc[k])/(tauS_nmda) ) - exp(- (t-spiketimesexc[k])/(tauF_nmda) ) ) )*(Vexc[i]-Esyne)  ;
		}
	}
	//endif for ee
// spike_rand_pre=((double)rand()/((double)(RAND_MAX)))*0.5 +0.5;
//	cout<<"got 1060"<<endl;

//spike_rand_pre=((double)rand()/((double)(RAND_MAX)))*0.5 +0.5;
//if for ii
	if( (i<ninh) && (k<ninh) )
	{
		
		if((connectivityii[i][k]) ){
		
		/*
		if(i>(ninh+1))
			{	
			gii=gsynie;
			//wie[i][k]=gie;
			}
			else
			{

				gii=gsynii;
				//wie[i][k]=gie;
			}*/

//cout<< "799+k"<< 799+k<< endl;
//cout<< "799+i"<< 799+i<< endl;

		isyninh[i]=isyninh[i]+spike_pre_last[799+i][799+k]*wii[i][k]*(  (exp(- (t-spiketimesinh[k])/(5.5) ) - exp(- (t-spiketimesinh[k])/(0.2) ) ) )*(Vinh[i]-Esyni);
		}
		}
//	cout<<"got 1086"<<endl;	// if for ie
	if( (i<nexc) && (k<ninh) )
		{
		
		 if(connectivityie[i][k]) {

/*int connectivityie[nexc][ninh];

	for(int i=0; i<nexc;i++)
	{
		for(int j=0; j<ninh;j++)
		{
			connectivityie[i][j]= preconnectivityie[j][i];
			connectivityietext<<connectivityie[i][j];
			connectivityietext<< " ";
		}
		connectivityietext<<endl;
	}*/
			/*
			if(k>(ninh+1))
			{	
			gie=gsynie;//gsynee*w_mult;
			wie[i][k]=gie;
			}
			else
			{

				gie=gsynie;
				wie[i][k]=gie;
			}*/




			isynexc[i]=isynexc[i]+spike_pre_last[i][799+k]*wie[i][k]*(  (exp(- (t-spiketimesinh[k])/(5.5) ) - exp(- (t-spiketimesinh[k])/(0.2) ) ) )*(Vexc[i]-Esyni);
		}
		}
		// if for ei
//	cout<<"got 1124"<<endl;
 //spike_rand_pre=((double)rand()/((double)(RAND_MAX)))*0.5 +0.5;
if( (i<ninh) && (k<nexc)  )
			{
		
		 if(connectivityei[i][k]){

				//isyninh[i]=isyninh[i]+gsynei*(  (exp(- (t-spiketimesexc[k])/(3.0) ) - exp(- (t-spiketimesexc[k])/(0.2) ) ) )*(Vinh[i]-Esyne);

			/*
			if(i>(ninh-10))
			{	
			gei=gsynee;
			//wie[i][k]=gie;
			}
			else
			{

				gei=gsynei;
				//wie[i][k]=gie;
			}*/

				B_NMDA= (1.0)/( 1.0+exp(( -Vinh[i]-10.0  )/(16.13))  );

				isyninh[i]=isyninh[i]+spike_pre_last[799+i][k]*wei[i][k]*(  (exp(- (t-spiketimesexc[k])/(3.0) ) - exp(- (t-spiketimesexc[k])/(0.2) ) ) )*(Vinh[i]-Esyne)
			  + wei_nmda[i][k]*B_NMDA*(  (exp(- (t-spiketimesexc[k])/(tauS_nmda) ) - exp(- (t-spiketimesexc[k])/(tauF_nmda) ) ) )*(Vinh[i]-Esyne) ;
			}
			}


		/*if(strogatzmattrans[i][k])
		{
			//recentconnectedspiketimes[k]=spiketimes[k];
			isyn[i]=isyn[i]+ w[i][k]*(  (exp(- (t-spiketimes[k]-tauD)/tauS ) - exp(- (t-spiketimes[k]-tauD)/tauF ) ) )*(V[i]-Esyne);
		}*/








		



	}//end for k
	}// end for i

/* 
 *if(noise_rand<prob_noise) 
	{ noise_time[i]=t;}

 if ((t-noise_time[i])>noise_dur)
{noise_val=noise_strength} 
else{noise_val=0}
((d
   
   
   
   
   ((d
   
   
   
   cout<<"got 1005"<endl;
   

*/ //spike_rand_pre=((double)rand()/((double)(RAND_MAX)))*0.5 +0.5;
//cout<<"got 1194"<<endl;
for(int i=0;i<ntot;i++)
	{

	
	if ((j%t_stim_it)==0)
	{
		t_stim=t;

		//cout<< "t_stim"<< t_stim <<endl;
	}

	//Iex=Iext( t, idc[i], ampdc,Aosc, fosc);
	//V_integrate(double V1,double minf, double nv1, double sv1, double hv1, double isyn1,double Iext)
	//calcout = V_integrate(V[i],m[i],n[i],s[i],h[i],isyn[i], /*Iex*/,Id,dt,gks );

//cout<<"got 1210"<<endl;

	if(i<nexc){  
		if(t<300){isyne=0; idcexc_mod=idcexc[i]*1.0;}
		else{isyne=isynexc[i]; idcexc_mod=idcexc[i];}
        


				if( (t<(t_stim+4)) && (t >(t_stim)) && (stim_vect[i] ))//(curr_distance_prob && t<(200))//t<(t_stim*2) )
                {
                idcexc_mod=idcexc_mod+stim_val;//iE*3;
                //cout<< "x" <<radius_connect_ee[i][0] << "y" << radius_connect_ee[i][1] << endl;
        
                }

		
		if((t-noise_time_last[i])<noise_dur)
            	{noise_val=noise_strength;}
                else{noise_val=0;}
 
             
		mexc[i]=minf(Vexc[i]);
		calcoutexc = V_integrate(Vexc[i],mexc[i],nnexc[i],sexc[i],hexc[i],isyne,noise_val ,idcexc_mod/*Id*/,dt,gksexc ); 
	
		if( (calcoutexc.Vv>0) && (Vexc[i]<0) )
			{
				spiketimesexc[i]=t;
				rastertextexc<<i+1<<" "<<t<<endl;
				spikenumexc=spikenumexc+1;

 				for(int pp=0; pp<ntot; pp++)
                                {
                                        spike_pre_last[pp][i]=((double)rand()/((double)(RAND_MAX)))*0.5 +0.5;
                                }


			}
		
		
		Vexc[i]=calcoutexc.Vv;
			mexc[i]=calcoutexc.mv;
			nnexc[i]=calcoutexc.nv;
			sexc[i]=calcoutexc.sv;
			hexc[i]=calcoutexc.hv;

//	if(i==2){
//		vtextexc<<Vexc[i]<<" "<<t<<endl;
//	}	
		//if(j==120){cout<<isynexc[i]<<endl;}
	isynexc[i]=0;
	}
//	cout<<"got 1261"<<endl;
	if(i<ninh){
		minh[i]=minf(Vinh[i]);
		if(t<300){isyni=0;}
		else{isyni=isyninh[i];}
		

			if(i>(ninh+1))//if(i>(ninh-10))
			{	
			//idcinh_mod=idcinh[i]; //idcinh_mod=idcexc[i+599];
			//gksinh_mod=gks_inh; // gksexc;
			idcinh_mod=idcexc[i+599];
			gksinh_mod=gksexc;
			}
			else
			{

			//idcinh_mod=idcinh[i];
			//gksinh_mod=gksinh;	
			idcinh_mod=idcexc[i+599];
			gksinh_mod=gksexc;

			}

//cout<<"got 1285"<<endl;
    if((t-noise_time_last[799+i])<noise_dur)
                {noise_val=noise_strength;}
                else{noise_val=0;}	
	
		calcoutinh= V_integrate(Vinh[i],minh[i],nninh[i],sinh[i],hinh[i],isyni,/**/noise_val,idcinh_mod/* Id*/,dt,gksinh_mod);   
		





		if( (calcoutinh.Vv>0) && (Vinh[i]<0) )
		{
			spiketimesinh[i]=t;
			rastertextinh<<i+nexc+1<<" "<<t<<endl;
			spikenuminh=spikenuminh+1;


			for(int pp=0; pp<ntot; pp++)
		        	{
        	        		spike_pre_last[pp][799+i]=((double)rand()/((double)(RAND_MAX)))*0.5 +0.5;
        			}	

                           				
		}
		Vinh[i]=calcoutinh.Vv;
		minh[i]=calcoutinh.mv;
		nninh[i]=calcoutinh.nv;
		sinh[i]=calcoutinh.sv;
		hinh[i]=calcoutinh.hv;
		
		isyninh[i]=0;
		
		
	}
///xxxxxxxxxxxxxxxxxxxxxxxxxxooooooooooooooooooooooooooo// change these guys back form zero











///////////////////////////////////////////////////////////////////////////////// end of conditionals for different ei ie ee ii




	}//////////////////////////// end for i loop



////////////////////////////////////////////////////////////////////////////////////////////////// end i loop





}////////////////////////////////////////////////////////////// end for j time loop




 /*for(int i=0; i<ntot; i++)
	{
	 for(int j=0; j<ntot; j++){
	 wtext<<w[i][j];
	 wtext<< " ";
	 }
	 wtext<< endl;
	}*/




 /*
 
  for(int i=0; i<nexc; i++)
 	{
 	 for(int j=0; j<nexc; j++){
 	 strogatzobj<<strogatzmattrans[i][j];
 	 strogatzobj<< " ";
 	 }
 	 strogatzobj<< endl;
 	}
 strogatzobj.close();
 cout<<"done";

 }*/

/*






ofstream weetext("wee.txt");
ofstream wiitext("wii.txt");
ofstream wietext("wie.txt");
ofstream weitext("wei.txt");
*/



	/*for(int i=0; i<nexc;i++)
	{
		for(int j=0; j<ninh;j++)
		{
			//cout<<preconnectivityie[j][i]<<"hi"<<endl;

			
			wietext<<wie[i][j];
			wietext<< " ";
		}
		wietext<<endl;
	}

for(int i=0; i<nexc;i++)
	{
		for(int j=0; j<nexc;j++)
		{
			//cout<<preconnectivityie[j][i]<<"hi"<<endl;

			
			weetext<<wee[i][j];
			weetext<< " ";
		}
		weetext<<endl;
	}*/



 
 stimtext.close();
// wietext.close();
 //weetext.close();
 connectivityeetext.close();
 connectivityiitext.close();
 connectivityeitext.close();
 connectivityietext.close();
 
 
 vtextexc.close();
 rastertextexc.close();
 mtextexc.close();
 ntextexc.close();
 stextexc.close();
 htextexc.close();
 isyntextexc.close();
 
 vtextinh.close();
 rastertextinh.close();
 mtextinh.close();
 ntextinh.close();
 stextinh.close();
 htextinh.close();
 isyntextinh.close();
 
 testtext.close();

 
 weetext.close();
 wiitext.close();
 wietext.close();
 weitext.close();
 weinmdatext.close();
 weenmdatext.close();

 
 
 
 

cout<<" spikenumtot:"<<spikenumexc+spikenuminh<<endl;






}




///////////////////////////text files
/*ofstream connectobj;
connectobj.open("connectmat.txt");
 for(int i=0; i<ntot; i++)
	{
	 for(int j=0; j<ntot; j++){
	 connectobj<<connectmat[i][j];
	 connectobj<<" ";
	 }
	 connectobj<< endl;
	}
connectobj.close();*/





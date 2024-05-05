/*

PROGRAMMA DI DINAMICA MOLECOLARE: GAS 

*/
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <time.h> 
#include <TROOT.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAttLine.h>



using namespace std;

const double sigma = 3.4*pow(10, -10); // scala dell'argon in angstrom
const double eps = 120*1.38*pow(10,-23); // in Joule, -eps è il minimo dell'energia potenziale 
const double Kb = 1; //1.38*pow(10,-23); 
const double Nav = 6.022*pow(10,23);
//const double m = (n*39.88)/(N*1000); // massa in kg di un atomo di argon, 39.88 è la massa molare dell'argon 10^-26
//const double tau = sqrt((m*sigma*sigma)/eps); // unità di tempo, circa 10^-12 secondi
// definite queste costanti avremo a che fare con unità di energia, lunghezza e tempo in epsilon, sigma e tau che danno numeri ragionevolmente normali 
const double rCut = 2; // con rCut=6 il gas è ancora in regime di gas ideale
double Tc = 3; // 360 K in questo caso, temperatura caratteristica
double Pc = eps/pow(sigma,3);
const int nbin=5; // numero di bin in cui verrà scomposta la run per fare le misure
vector<double> T[nbin],Tbin; 
vector<double> P[nbin],Pbin;
vector<double> E[nbin];
vector<double> corrT1,corrT2; // funzione di autocorrelazione
vector<double> corrP1,corrP2;
vector<double> corrV;
double corr_timeT1,corr_timeT2;
double corr_timeP1,corr_timeP2;
double stdT,stdP,meanT,meanP;
double corr_timeV;
const double eta = 4*(1./pow(rCut,12)-1./pow(rCut,6));
double Ek_tot=0;
double Mom_tot; // momento totale
double uSum;
double T_true;
struct vec_3d {

	double x,y,z;

};

vec_3d L = { // in unità sigma
	.x = 30,
	.y = 30,
	.z = 30,
};

const int dmol = 7;
const int N = dmol*dmol*dmol;
double spacing = L.x / dmol;
double Vol = (L.x)*(L.y)*(L.z);
double rho = N/Vol; 



struct mol {

	vec_3d pos,vel,a;
};



vector<mol> *molecules = new vector<mol> (N);

double mean (vector<double> measure) {
	double mean = 0;
	for (int i=0;i<measure.size();i++) {
		mean += measure[i];
		
	}		
	mean = mean/measure.size();
	return mean;
}

double maxwell (double * v) {
	double vv = v[0]*v[0];
	double A = pow(1./2*M_PI*Kb*Tc,3/2);
	double e = exp(-vv/(2*Kb*Tc));
	return A*e;
}

void scaling_temp (double T) {
	double vvSum = 0;
	for (int i=0;i<molecules->size();i++) {
		vvSum += pow((*molecules)[i].vel.x,2)+pow((*molecules)[i].vel.y,2)+pow((*molecules)[i].vel.z,2);
	}
	double vvTrue = 3*T*N;
	double scale = sqrt(vvTrue/vvSum);
	for (int i=0;i<molecules->size();i++) {
		(*molecules)[i].vel.x *= scale;
		(*molecules)[i].vel.y *= scale;
		(*molecules)[i].vel.z *= scale;
	}
	
}
void ShowParameters() {

	cout << "L.x      L.y     L.z      V       T         " << endl;
	cout << L.x << "    " << "    " << L.y << "     " << L.z << "    " << Vol << "    " << T_true << endl; 
	
}


void ClearParametersOutStep() {
	
	for(int i=0;i<nbin;i++) {
		T[i].clear(); 
		P[i].clear();
	}
	
}






void init_coord () { // inizializzo posizioni e accellerazioni
	int n = 0;
	for (int ix=0;ix<dmol;ix++) { // POSIZIONI
		for(int iy=0;iy<dmol;iy++) {
			for(int iz=0;iz<dmol;iz++) {
				(*molecules)[n].pos.x = (L.x/dmol)*ix;
				(*molecules)[n].pos.y = (L.y/dmol)*iy;
				(*molecules)[n].pos.z = (L.z/dmol)*iz;
				n = n+1;
			}
		}
	}
	
	for(int i=0;i<molecules->size();i++) { //ACCELLERAZIONI
		(*molecules)[i].a.x = 0;
		(*molecules)[i].a.y = 0;
		(*molecules)[i].a.z = 0;
	}
}






void VRand()	{
	double phi,teta; // angoli che danno la direzione del versore
	double scale = 2*M_PI/RAND_MAX; // in questo modo genero numeri casuali tra 0 e 2PIgreco
	vec_3d * v;
	for(int i=0;i<molecules->size();i++) {
		v = &(*molecules)[i].vel;
		phi = (double)rand()*scale;
		teta = (double)rand()*scale;
		(*v).x = cos(teta)*sin(phi);
		(*v).y = sin(teta)*sin(phi);
		(*v).z = cos(phi);
	}

}





double vel_mag(double T) { 
	return sqrt(3*T);

}

vec_3d VSum = {
	.x = 0,
	.y = 0,
	.z = 0,
};

void init_vel (double T) {
	VRand();
	for (int i=0;i<molecules->size();i++) {
		(*molecules)[i].vel.x *= vel_mag(T);
		(*molecules)[i].vel.y *= vel_mag(T);
		(*molecules)[i].vel.z *= vel_mag(T);
		VSum.x += (*molecules)[i].vel.x;
		VSum.y += (*molecules)[i].vel.y;
		VSum.z += (*molecules)[i].vel.z;
		}
	for (int i=0;i<molecules->size();i++) {
		(*molecules)[i].vel.x -= (1./N)*VSum.x; // faccio in modo che il centro di massa abbia velocità nulla
		(*molecules)[i].vel.y -= (1./N)*VSum.y;
		(*molecules)[i].vel.z -= (1./N)*VSum.z;
	}
	VSum.x = 0;
	VSum.y = 0;
	VSum.z = 0;
	for (int i=0;i<molecules->size();i++) { 
		VSum.x += (*molecules)[i].vel.x;
		VSum.y += (*molecules)[i].vel.y;
		VSum.z += (*molecules)[i].vel.z;
		Ek_tot += 0.5*(pow((*molecules)[i].vel.x,2)+pow((*molecules)[i].vel.y,2)+pow((*molecules)[i].vel.z,2));
		// ENERGIA CINETICA TOTALE IN MD
	}
	
}



void SetParameters (double T) {
	init_coord (); // ho inizializzato posizione e acc
	init_vel(T); // ho inizializzato le velocità
} 

void acc_assignment (double fcVal, vec_3d dr, double j1, double j2) {
	vec_3d * a1;
	 a1 = &(*molecules)[j1].a;
	vec_3d * a2;
	a2 = &(*molecules)[j2].a;
	(*a1).x += +fcVal*dr.x;   // in MD units, vedi pag. 14 "art MD"
	(*a2).x += -fcVal*dr.x;
	(*a1).y += +fcVal*dr.y;
	(*a2).y += -fcVal*dr.y;
	(*a1).z += +fcVal*dr.z;
	(*a2).z += -fcVal*dr.z;
}

double Fvirial; // variabile utile per calcolare il viriale della forza 
void ComputeForces ()	{
	vec_3d dr; // distanza fra due molecole
	double rr,rrCut,rri,rri3,fcVal; // rispettivamente: distanza fra le due molecole al quadrato, cutoff al quadrato, inverso rr e (rri)^3
	int j1,j2,n;
	rrCut=rCut*rCut; 
	for (n=0;n<molecules->size();n++) { // questo ciclo è necessario perchè se dr>rCut allora l'accellerazione della molecola è nulla
		 (*molecules)[n].a.x =0;
		 (*molecules)[n].a.y =0;
		 (*molecules)[n].a.z =0; 
	}
	uSum = 0;
	Fvirial = 0;
	for (j1=0;j1<molecules->size();j1++) {
		for (j2=0;j2<j1;j2++) {
			dr.x= (*molecules)[j1].pos.x-(*molecules)[j2].pos.x; // sono già in unità di sigma perchè così inizializzate
			dr.y= (*molecules)[j1].pos.y-(*molecules)[j2].pos.y;
			dr.z= (*molecules)[j1].pos.z-(*molecules)[j2].pos.z;	
			
			if (dr.x >= 0.5*L.x) { // MINIMUM IMAGE CONVENTION X
				dr.x = dr.x - L.x;
				}
			if (dr.x <-0.5*L.x) {
				dr.x = dr.x + L.x;
			}
			if (dr.y >= 0.5*L.y) { // MINIMUM IMAGE CONVENTION X
				dr.y = dr.y - L.y;
			}
			if (dr.y <-0.5*L.y) {
				dr.y = dr.y + L.y;
			}
			
			if (dr.z >= 0.5*L.z) { // MINIMUM IMAGE CONVENTION X
				dr.z = dr.z - L.z;
			}
			if (dr.z <-0.5*L.z) {
				dr.z = dr.z + L.z;
			}
			
			rr = dr.x*dr.x+dr.y*dr.y+dr.z*dr.z;

			if (rr<rrCut) {
				rri=1./rr;
				rri3=pow(rri,3);
				fcVal = 48.*rri3*(rri3-0.5)*rri; // vedi lezione 4 pag.1, già scritta in MD units
				Fvirial += fcVal*rr;
				acc_assignment (fcVal,dr,j1,j2);
				uSum += 4.*rri3*(rri3-1.);//-eta; // già scritta in MD units
			}
		}
	}
	
}

void Boundary_Condition () { // una volta che i parametri della molecola sono stati aggiornati, applico le condizioni periodiche
	for (int i=0;i<molecules->size();i++) {
		if ((*molecules)[i].pos.x >= L.x) {
			(*molecules)[i].pos.x -= L.x;
		}
		
		if ((*molecules)[i].pos.x <= 0) {
			(*molecules)[i].pos.x += L.x;
		}
		
		if ((*molecules)[i].pos.y >= L.y) {
			(*molecules)[i].pos.y -= L.y;
		}
		
		if ((*molecules)[i].pos.y <= 0) {
			(*molecules)[i].pos.y += L.y;
		}
		
		if ((*molecules)[i].pos.z >= L.z) {
			(*molecules)[i].pos.z -= L.z;
		}
		
		if ((*molecules)[i].pos.z <= 0) {
			(*molecules)[i].pos.z += L.z;
		}
	
	
	}

}

double corr_func1 (vector<double> measure,double true_value,int t) { // t deve essere almeno 1
	int N = measure.size();
	double autocorr=0;

	for (int j=0;j<N-t;j++) {
		autocorr+=measure[j]*measure[j+t];
	}
	autocorr = autocorr/(N-t);
	autocorr = autocorr-true_value*true_value;
	return autocorr;
}

double corr_func2 (vector<double> measure,double true_value,int t) { // t deve essere almeno 1
	int N = measure.size();
	double autocorr=0;
	for (int j=0;j<N-t;j++) {
		autocorr+=(measure[j]-true_value)*(measure[j+t]-true_value);
	}
	autocorr = autocorr/N;
	return autocorr;
}

double corr_func3 (vector<double> measure,int t) { // t deve essere almeno 1
	int N = measure.size();
	double autocorr=0;
	double Xa =0;
	double Xb =0;
	for (int j=0;j<N-t;j++) {
		Xa += measure[j];
	}
	Xa = Xa/(N-t);
	for (int j=0;j<N-t;j++) {
		Xb += measure[j+t];
	}
	Xb = Xb/(N-t);
	for (int j=0;j<N-t;j++) {
		autocorr+=(measure[j]-Xa)*(measure[j+t]-Xb);
	}
	autocorr = autocorr/(N-t);
	return autocorr;
}

double autocorr_time (vector<double> corr_vec,vector<double> measure) {
	double sum = 0;
	for (int i=0;i<corr_vec.size();i++) {
		sum += abs(corr_vec[i]); // se metto abs nella somma totale viene sempre uno, giustamente, in quanto senza abs veniva sempre 0
	}
	double mean_square = 0;
	for (int i=0;i<measure.size();i++) {
		mean_square += measure[i]*measure[i];
	}
	mean_square = mean_square/measure.size();
	double variance = mean_square-mean(measure)*mean(measure);
	sum = 0.5*(1+2*sum/variance);
	return sum;
}


double Std (vector<double> measure) {
	double Mean = mean(measure);
	double mean_square = 0;
	for (int i=0;i<measure.size();i++) {
		mean_square += pow(measure[i],2);
		
	}	
	mean_square = mean_square/measure.size(); // media dei quadrati
	double std = sqrt((mean_square-Mean*Mean)/measure.size());
	return std;
}




void one_point (vec_3d * df, vec_3d * f, double dt) {
		
		(*f).x += (*df).x*dt;
		(*f).y += (*df).y*dt;
		(*f).z += (*df).z*dt;


}


void leapfrog (double dt) {
	ComputeForces();
	vec_3d * a;
	vec_3d * v;	
	vec_3d * x;
	for (int j=0;j<molecules->size();j++) {
		x = &(*molecules)[j].pos;
		v = &(*molecules)[j].vel;
		a= &(*molecules)[j].a;
		(*v).x += (*a).x*(dt/2);
		(*v).y += (*a).y*(dt/2);
		(*v).z += (*a).z*(dt/2);
		(*x).x += (*v).x*dt;
		(*x).y += (*v).y*dt;
		(*x).z += (*v).z*dt;
	}
	Boundary_Condition();
	ComputeForces();
	for (int j=0;j<molecules->size();j++) {
		x = &(*molecules)[j].pos;
		v = &(*molecules)[j].vel;
		a= &(*molecules)[j].a;
		(*v).x += (*a).x*(dt/2);
		(*v).y += (*a).y*(dt/2);
		(*v).z += (*a).z*(dt/2);
	}
}

TF1 maxw_distr (vector<double> measure, double width, int print,vector <double> T) { // measure sarà il vettore con le misure non correlate
	double mean=0;
	double min = measure[0];
	double max = measure[0]; // variabili che tengono traccia del valore max e min di M  
	for (int i=0; i<measure.size();i++) {
		if (measure[i]>max) {
			max = measure[i];
			}
		if (measure[i]<min) {
			min=measure[i];
		
		}
	}
	int n = int ((max-min)/width);
	TH1D * h = new TH1D ("v","maxw distr",n,min,max); // creo l'istogramma e lo passo al puntatore della funzione
	for (int i=0; i<measure.size();i++) {
		h->Fill(measure[i]);
	
	}
	h->Scale(1./h->Integral());
	TF1 maxw_f ("maxw-function","[0]*x*x*exp(-x*x/[1])",min,max); 
	gStyle->SetOptFit(1101);
	maxw_f.SetParameter(0,(4.*M_PI*rho)/pow(2*M_PI*T.back(),3/2));
	maxw_f.SetParameter(1,2*M_PI*T.back());
	h->Fit("maxw-function","RQ");
	if (print == 1) {
		TCanvas * myCanv = new TCanvas ("my canv","my canv");	
		cout << "p0= " << maxw_f.GetParameter(0) << " p0exp=  " << (4.*M_PI*rho)/pow(2*M_PI*T.back(),3/2) << endl;
		cout <<"p1= " << maxw_f.GetParameter(1) << " p1exp= " <<  2*T.back() << endl;
		h-> Draw();
		myCanv->Print("speed.gif","gif");
		delete gROOT->FindObject("my canv");
		}
	delete gROOT->FindObject("v");
	return maxw_f;
}
	
double Kvirial;
double virial; 
void integrator (double * t, double nStep, double dt, int k, int nbin) {
	*t =0;
	vector<double> V;
	ClearParametersOutStep();
	double Ek; //=0;
	int d = nStep/nbin;
	int n_equilibrium = 0;
	TF1 f;
	for(int j=0;j<nStep;j++) {
			virial = 0;
			Kvirial = 0;	
			Ek_tot = 0;
			VSum.x =0; VSum.y=0; VSum.z=0;
			V.clear();
			leapfrog(dt);
			for (int i = 0; i<molecules->size();i++) {
				VSum.x += (*molecules)[i].vel.x; 
				VSum.y += (*molecules)[i].vel.y; 
				VSum.z += (*molecules)[i].vel.z;
				V.push_back(sqrt(pow((*molecules)[i].vel.x,2)+pow((*molecules)[i].vel.y,2)+pow((*molecules)[i].vel.z,2)));
				Ek = 0.5*(pow((*molecules)[i].vel.x,2)+pow((*molecules)[i].vel.y,2)+pow((*molecules)[i].vel.z,2));// MDunits	
				Kvirial += 2*Ek;
				Ek_tot += Ek;
			}
			*t += dt;

			if (j%k == 0) { 
				for(int s=0;s<nbin;s++) {
					if (j/d==s) { // per questo ho questi due IF
						T[s].push_back ((2./3.)*(Ek_tot/(Kb*N))); // in MD units
						virial = Fvirial+Kvirial; 
						P[s].push_back ((1./3.)*virial/Vol); 
						E[s].push_back(Ek_tot+uSum);
						Mom_tot =  sqrt(pow(VSum.x,2)+pow(VSum.y,2)+pow(VSum.z,2));
						f=maxw_distr(V,0.15,0,T[s]);
						cout << "t        Ek         U         E          v          T           P" << endl;
						cout << *t << "  " << Ek_tot << "  " << uSum<< "  " << E[s].back() << "  "<< Mom_tot << "  ";
						cout << T[s].back() << "     " << P[s].back() << endl; 
						if (f.GetProb()>0.05 && n_equilibrium == 0 ) {
						n_equilibrium = j;
						}
					}
					
				}
			}
			if (j%d == 0) {
				scaling_temp(T_true);
			}
		}

		/*for (int t=1;t<T1.size();t++) {
			corrT1.push_back(corr_func2(T1,mean(T1),t));
			corrP1.push_back(corr_func2(P1,mean(P1),t));
		}*/
		double sigmaT[4];
		double sigmaP[4];
		double sigma_sq=0;
		double sigma_sqP=0;
		double sigma_sqT=0;
		for(int i=1;i<nbin;i++) { // escludo il primo bin
			Tbin.push_back(mean(T[i]));
			Pbin.push_back(mean(P[i]));
			sigmaT[i]=Std(T[i]);
			sigmaP[i]=Std(P[i]);
			sigma_sqT += sigmaT[i]*sigmaT[i];
			sigma_sqP += sigmaP[i]*sigmaP[i];
			}
		
		stdP = sqrt(sigma_sqP)/nbin;
		stdT = sqrt(sigma_sqT)/nbin;
		meanT = mean(Tbin);
		meanP = mean(Pbin);
		cout << "meanT        meanP        stdT        stdP" << endl; 
		cout << meanT << "   " << meanP;
		cout << "  " << stdT << "  " << stdP  << endl;	
		cout << "errT_rel    errP_rel" << endl;
		cout << 100*stdT/meanT<<"%" << "    " << 100*stdP/meanP<<"%" << endl;
		
		/*cout << "corr_timeT1 = " << autocorr_time (corrT1,T1) << endl;
		cout << "corr_timeP1 = " << autocorr_time (corrP1,P1) << endl;
		cout << "corr_timeT2 = " << autocorr_time (corrT2,T2) << endl;
		cout << "corr_timeP2 = " << autocorr_time (corrP2,P2) << endl;*/
		
		// COSTRUISCO L'ISTOGRAMMA DELLE VELOCITA'

		f=maxw_distr(V,0.15,1,T[nbin-1]);
		cout << "n_equilibrium_maxwell = " << n_equilibrium << endl;
		cout << f.GetProb() << endl;
}







int main () 	{

srand(time(NULL));
double t_in = 0;
int nStep = 20000; 
double dt = 0.001;
T_true = 3.8;
int nPoints= 15;

double * t = new double;
double measure_step = 50;
int k = measure_step;; // decido ogni quanti step fare la misura
*t = t_in;
double vol[25];
double Pmean[25];
double Tmean[4];			
TF1 fitfunz1("gas-curve1","343*[0]/x",3000,30000);
TF1 fitfunz2("gas-curve2","343*[0]/x",3000,30000);
TF1 fitfunz3("gas-curve3","343*[0]/x",3000,30000);
TF1 fitfunz4("gas-curve4","343*[0]/x",3000,30000);
TGraph pV[4],pV_copy[4];
TGraph *pv;
TGraph *pv_copy;
TMultiGraph * mg = new TMultiGraph();
TMultiGraph * mg_copy = new TMultiGraph();
for (int i=0;i<4;i++) {
	for (int j=0;j<nPoints;j++) {
			for(int l=0;l<2;l++) {	
				SetParameters(T_true);
				ShowParameters();
				integrator (t,nStep,dt,k,nbin);
				if(stdT/meanT < 0.01 && stdP/meanP <0.01) {
					cout << "EQUILIBRIUM REACHED BY ERROR" << endl;
					l=2;
				}
				else {
					if(l==1) {
						cout << "EQUILIBRIUM NOT REACHED BY ERROR" << endl;
					}
					else {
						k = k/2;
						Tbin.clear();
						Pbin.clear();
					}
				}
			}
	Tbin.clear();
	Pbin.clear();
	k = measure_step;
	pV[i].SetPoint(j,Vol,meanP);
	pV_copy[i].SetPoint(j,Vol,meanP);
	L.x *= 0.91; 
	L.y *= 0.91;
	L.z *= 0.91;
	Vol = L.x*L.y*L.z;
	rho = N/Vol;
	}
	Tmean[i] = T_true; // ho scelto T perchè che valore dovrei prendere per la media? Se calcolo il sistema per 15 diferrenti volumi
	pV[i].SetMarkerStyle(21+i);
	pV[i].SetMarkerColor(i+1);
	pV[i].SetMarkerSize(0.6);
	pV_copy[i].SetMarkerStyle(21+i);
	pV_copy[i].SetMarkerColor(i+1);
	pV_copy[i].SetMarkerSize(0.6);	
	pv = &pV[i];
	pv_copy = &pV_copy[i];
	mg->Add(pv);
	mg_copy->Add(pV_copy);
	L.x = 30;
	L.y = 30;
	L.z = 30;
	Vol = L.x*L.y*L.z;
	rho = N/Vol; 
	T_true -= 1;
	if (T_true<1) {
		nPoints==17;
	}
}
TCanvas pv_curves("canv","canv");
pv_curves.SetLogy(1);
mg -> SetTitle ("pV-curves;Volume;Pressure");
	fitfunz1.SetParameter(0,Tmean[0]);
	fitfunz2.SetParameter(0,Tmean[1]);
	fitfunz3.SetParameter(0,Tmean[2]);
	fitfunz4.SetParameter(0,Tmean[3]);
	fitfunz1.SetLineColor(kGreen);
	fitfunz2.SetLineColor(kRed);
	fitfunz3.SetLineColor(kBlue);
	fitfunz4.SetLineColor(kYellow);
	gStyle->SetOptFit(0000);
	pV[0].Fit("gas-curve1","R");
	pV[1].Fit("gas-curve2","R");
	pV[2].Fit("gas-curve3","R");
	pV[3].Fit("gas-curve4","R");

TLegend * legend1 = new TLegend (0.6,0.9,0.6,0.9);
legend1->AddEntry("gas-curve1","T=3.8","L");
legend1->AddEntry("gas-curve2","T=2.8","L");
legend1->AddEntry("gas-curve3","T=1.8","L");
legend1->AddEntry("gas-curve4","T=0.8","L");
mg->Draw("AP");
legend1->Draw();


pv_curves.Print("pVcurves.png","png");	



TCanvas pv_curves_real("canv2","canv2");
pv_curves_real.SetLogy(1);
mg -> SetTitle ("pV-curves-real;Volume;Pressure");
TF1 fitfunz_real1("gas-curve-real1","[0]/(x/343-[2])-[1]*pow(343/x,2)",3000,30000); 
TF1 fitfunz_real2("gas-curve-real2","[0]/(x/343-[2])-[1]*pow(343/x,2)",3000,30000); 
TF1 fitfunz_real3("gas-curve-real3","[0]/(x/343-[2])-[1]*pow(343/x,2)",3000,30000); 
TF1 fitfunz_real4("gas-curve-real4","[0]/(x/343-[2])-[1]*pow(343/x,2)",3000,30000); 
	
	fitfunz_real1.FixParameter(0,3.8);
	fitfunz_real2.FixParameter(0,2.8);
	fitfunz_real3.FixParameter(0,1.8);
	fitfunz_real4.FixParameter(0,0.8);
	/*fitfunz_real1.SetParameter(2,5);
	fitfunz_real2.SetParameter(2,5);
	fitfunz_real3.SetParameter(2,5);
	fitfunz_real4.SetParameter(2,5);
	fitfunz_real1.SetParameter(1,20);
	fitfunz_real2.SetParameter(1,20);
	fitfunz_real3.SetParameter(1,20);
	fitfunz_real4.SetParameter(1,20);*/
	fitfunz_real1.SetLineColor(kGreen);
	fitfunz_real2.SetLineColor(kRed);
	fitfunz_real3.SetLineColor(kBlue);
	fitfunz_real4.SetLineColor(kYellow);
	pV[0].Fit("gas-curve-real1","R");
	pV[1].Fit("gas-curve-real2","R");
	pV[2].Fit("gas-curve-real3","R");
	pV[3].Fit("gas-curve-real4","R");
	
	gStyle->SetOptFit(0000);
	


TLegend * legend2 = new TLegend (0.6,0.9,0.6,0.9);
legend2->AddEntry("gas-curve-real1","T=3.8","L");
legend2->AddEntry("gas-curve-real2","T=2.8","L");
legend2->AddEntry("gas-curve-real3","T=1.8","L");
legend2->AddEntry("gas-curve-real4","T=0.8","L");
mg->Draw("AP");	
legend2->Draw();


pv_curves_real.Print("pVcurves-real.png","png");	



return 0;
}





	

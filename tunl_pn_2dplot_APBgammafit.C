#include <TMath.h>
#include <stdio.h>		// Standard input and output
#include <string.h>		// String operations


//data carrying variables
TH1F *Data2 = new TH1F("tof_comb","Time of Flight Combined", 4096,0.,4095);
//TH1F *Data2;
TH1F *Data, *htemp, *BData;
TH1F *Data_rescale = new TH1F("tof","tof",24000,-463.72280,4.8851063); 
TFile *f;
TF1 *theory, *eback, *esig, *esig2, *esig3, *esig4;
Int_t npar = 8, i;

//fit function variables
Double_t params[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  
Double_t dparams[8];

//x (not parameter, is domain) ,AMPLitude, mew, sigma, sigmaGAMMA1, sigmaGAMMA2, sigmaGAMMA3, amplitudeGAMMA1, amplitudeGAMMA2, amplitudeGAMMA3, displacementGamma, localDis1, localDis2, 
//x (not parameter, is domain) ,AMPLitude, mew, sigma, sigmaGAMMA1, sigmaGAMMA2, sigmaGAMMA3, amplitudeGAMMA1, amplitudeGAMMA2, amplitudeGAMMA3, displacementGamma, localDis1, localDis2, 
Double_t Nparam[40] = {1.3, 0, 0.04, 0, 0, 0, 196, 388, 151, 461, 100, 2821, 2822, 2823, 2824, 2825, 3, 0, 0, 0, 0, 0.5, 1, 1, 0, 0, 3000, 3500, 800, 0, 3, -3, 0, 0, 0, 5, 10, 15, 20, 0};
Int_t bt = 3; //background type, 0 for linear, 1 for exp
// ***************SIGMA from GAMMA (bake in sigmas from gamma model)
Double_t Gsig1 = 0.8308; // peak 1
Double_t Gsig2 = 0.8958; // peak 2
Double_t Gsig3 = 2.265; //backgaus
// ***************AMPLITUDE from GAMMA
Double_t Ga1 = 2083.; // peak 1
Double_t Ga2 = 4003.; // peak 2
Double_t Ga3 =  1399.; // backgaus
// ***************DISPLACEMENTS of GAMMA Model (spacing of the gaussians within a single peak)(baked in)
Double_t Gdis1t2 = 3.;	//peak 1 to peak 2
Double_t Gdis1t3 = 0.0;	//peak 1 to backgaus
  




//left over variables from previous code
Double_t low, high; // this is the total ROI
Double_t low1, high1, low2, high2;
Double_t peak1, peak2, sum1, sum2;
const char *pname[] = {"Intercept", "Slope", "Constant 1", "Mean 1", "Sigma 1", "Constant 2", "Mean 2", "Sigma 2"}; 
 
//constants
Double_t mn = 939.7;//neutron mass (Mev/c^2)
Double_t Q = -0.8728;//136Xe(p,n)136Cs Q (MeV)
Double_t L = 13.0;//flight path (meters)
Double_t Eb = 7.0;//Beam energy (MeV)
Double_t Eloss = 0.0;//Energy loss in cell (MeV)
Double_t c = 0.2998;//speed of light (m/ns)
Double_t tofg = (L-0.1)/c;//gamma tof (ns)



//information from experinmental set up
//Table of Offsets
// channel#: offset#
// 0: 		+0
// 1: 		+4
// 2: 		+4
// 3: 		-7 
// 4: 		+9 
// 5: 		+7
// 6: 		-12
// 7: 		+1
// 8: 		-30 
// 9: 		+ test inconclusive: no peak in gamma, n or beam pulse; there does seem to be noise across all TOF
//10: 		+ test  inconclusive: no peak in gamma, n or beam pulse; data in relavent ranges is 0
//11: 		+ test no data in channel	
//12: 		+ test no data in channel	
//13: 		+ test no data in channel
//14: 		+ test no data in channel
//15: 		+ no active data
//16: 		-3 
//17: 		+2 
//18: 		+16 
//19: 		+4 
//20: 		+8
//21: 		+6

// loaddata takes data from file name passed to it and fills "Data2"
void loaddata(const char *filename)
{
//    char da[80];
    
//    sprintf(da, "%d", data_area);
    
    f = new TFile(filename);
    Data2 = (TH1F*)f->Get("tof_comb");
	//f->Close();
	/*
	BData = (TH1F*)f->Get("h_tof_bkg_0");
    Int_t numbins = Data->GetNbinsX();
   
    Int_t y=0;
    Double_t x = 0.0;

    for (Int_t i=0;i<numbins;i++){
       x = Data->GetBinCenter(i);
       y = Data->GetBinContent(i);
     printf("Bin  number %d, x: %f, y: %d \n",i,x,y); 
       for(Int_t j=0;j<y;j++){  
         Data_rescale->Fill(energy(x));
       }     
     }
*/
//Data = (TH1F*)f->Get(da);
    Data2->Draw();
   // Data_rescale->Draw();


}
	
//authored by previous researcher	
Double_t energy(Double_t tof) {
  Double_t ntof = tofg+tof;
  Double_t En = 0.5*mn*((L*L)/(ntof*ntof))/(c*c);//neutron energy (MeV)
  Double_t Ep = Eb - Eloss;//proton energy (MeV)
  Double_t answer = Ep+Q-En; //energy sum of proton and neutron (MeV)
return answer;
}

// -- The background function ----

Double_t Background(Double_t *x, Double_t *par) {
    Double_t val;
        
    if (bt == 0)
    	val = par[0] + par[1]*x[0];
    else if (bt == 1)
    	val = exp(par[0] + par[1]*x[0]);
	else if (bt == 2)
		val = par[0] + exp((x[0] - 2780.)/par[1]);
	else if (bt == 3)
		val = par[0] + par[1]*(x[0]-2780.) + par[2]*(x[0]-2780.)*(x[0]-2780.);
    else
    {
    	fprintf(stderr, "Huh?\n");
	exit(1);
    }
    
    return val;
}
Double_t Backgroundgamma(Double_t *x, Double_t *par) {
    Double_t val;
        
    if (bt == 0)
    	val = par[0] + par[1]*x[0];
    else if (bt == 1)
    	val = exp(par[0] + par[1]*x[0]);
	else if (bt == 2)
		val = par[0] + exp((x[0] - 2780.)/par[1]);
	else if (bt == 3)
		val = par[0] + par[1]*(x[0]-2780.) + par[2]*(x[0]-2780.)*(x[0]-2780.);
    else
    {
    	fprintf(stderr, "Huh?\n");
	exit(1);
    }
    
    return val;
}
// fit gaussian function for individual peaks



Double_t e1000(Double_t *x, Double_t *par){
   Double_t arg;
  if (par[2]){arg = (x[0] - par[1])/par[2];}
   return par[0]*TMath::Exp(-0.5*arg*arg);
}

Double_t e850(Double_t *x, Double_t *par){
  Double_t arg;
  if (par[2]){arg = (x[0] - par[1])/par[2];}
  return par[0]*TMath::Exp(-0.5*arg*arg);
}

Double_t e590(Double_t *x, Double_t *par){
    Double_t arg;
  if (par[2]){arg = (x[0] - par[1])/par[2];}
  return par[0]*TMath::Exp(-0.5*arg*arg);
}

Double_t e5839(Double_t *x, Double_t *par){
   Double_t arg;
  if (par[2]){arg = (x[0] - par[1])/par[2];}
   return par[0]*TMath::Exp(-0.5*arg*arg);
}

Double_t e5179(Double_t *x, Double_t *par){
  Double_t arg; 
 if (par[2]){arg = (x[0] - par[1])/par[2];}
   return par[0]*TMath::Exp(-0.5*arg*arg);
}

Double_t e1048(Double_t *x, Double_t *par){
   Double_t arg;
   if (par[2]){arg = (x[0] - par[1])/par[2];}
   return par[0]*TMath::Exp(-0.5*arg*arg);
}


Double_t e1200(Double_t *x, Double_t *par){
   Double_t arg;
  if (par[2]){arg = (x[0] - par[1])/par[2];}
   return par[0]*TMath::Exp(-0.5*arg*arg); 
}


Double_t e0260(Double_t *x, Double_t *par){
   Double_t arg;
   if(par[2]){arg = (x[0] - par[1])/par[2];}
   return par[0]*TMath::Exp(-0.5*arg*arg);
}

Double_t e1585(Double_t *x, Double_t *par){
    Double_t arg;
    if(par[2]){arg = (x[0] - par[1])/par[2];}
    return par[0]*TMath::Exp(-0.5*arg*arg);
}

// -- The signal function: a double gaussian -----

///Double_t Signal(Double_t *x, Double_t *par)

///{
//will fit peaks for following states in cs:
//1.000 MeV = 372.07 ns, 0.850 = 366.1 ns, 0.590 = 356.1 ns, 0.5839 = 355.8 ns, 0.5179 = 353.5 ns, 0.1048 = 339.6 ns
//dt1 = 6ns, dt2 = 15.7ns, dt3 = 16ns, dt4 = 18.3ns, dt5 = 32.3ns
    
   Double_t calctof_1000 = 372.07;
   Double_t calctof_0850 = 366.13;
   Double_t calctof_0590 = 356.41;
   Double_t calctof_05839 = 356.19; //not resolvable
   Double_t calctof_05179 = 353.83;
   Double_t calctof_01048 = 339.98;
   Double_t calctof_12000 = 380.41;  
   Double_t calctof_0260 = 345;  //found in spectrum
   

   Double_t dt1 = calctof_1000 - calctof_0850;
   Double_t dt2 = calctof_1000 - calctof_0590;
   Double_t dt3 = calctof_1000 - calctof_05839;
   Double_t dt4 = calctof_1000 - calctof_05179;
   Double_t dt5 = calctof_1000 - calctof_01048;
   Double_t dt6 = calctof_1000 - calctof_12000;
   Double_t dt7 = calctof_1000 - calctof_0260;
   
Double_t Signal(Double_t *x, Double_t *par){

    Double_t arg1 = 0;
    Double_t arg2 = 0;
    Double_t arg3 = 0;
    Double_t arg4 = 0;
    Double_t arg5 = 0;
    Double_t arg6 = 0;
    Double_t sig1 = 0;
    Double_t sig2 = 0;
    Double_t sig3 = 0;
    Double_t sig4 = 0;
    Double_t sig5 = 0;
    Double_t sig6 = 0;
    Double_t arg7 = 0;
    Double_t sig7 = 0;
    Double_t arg8 = 0;
    Double_t sig8 = 0;

   
    if (par[2])
    	arg1 = (x[0] - par[1])/par[2];
        sig1 = par[0]*TMath::Exp(-0.5*arg1*arg1);
   
    if (par[2])
    	arg2 = (x[0] - par[1] + dt1)/par[2];
        sig2 = par[3]*TMath::Exp(-0.5*arg2*arg2);

    if(par[2])
        arg3 = (x[0] - par[1] + dt2)/par[2];
        sig3 = par[4]*TMath::Exp(-0.5*arg3*arg3);

    if(par[2])
      arg4 = (x[0] - par[1] + dt3)/par[2];
      sig4 = par[5]*TMath::Exp(-0.5*arg4*arg4);

    if(par[2])
      arg5 = (x[0] - par[1] + dt4)/par[2];
      sig5 = par[6]*TMath::Exp(-0.5*arg5*arg5);

    if(par[2])
      arg6 = (x[0] - par[1] + dt5)/par[2];
      sig6 = par[7]*TMath::Exp(-0.5*arg6*arg6);

     if(par[2])
      arg7 = (x[0] - par[1] + dt6)/par[2]; 
      sig7 = par[8]*TMath::Exp(-0.5*arg7*arg7);
      
      if(par[2])
       arg8 = (x[0] - par[1] + dt7)/par[2];
       sig8 = par[9]*TMath::Exp(-0.5*arg8*arg8);


   
   return (sig1 + sig2 + sig3 + sig4 + sig5 + sig6 + sig7 + sig8);
}
Double_t SignalMAP(Double_t *x, Double_t *par){
	//par is Amplitude starting at 3
	//x is x
	


    Double_t arg1 = 0;
    Double_t arg2 = 0;
    Double_t arg3 = 0;
    Double_t arg4 = 0;
    
	Double_t sig1 = 0;
    Double_t sig2 = 0;
    Double_t sig3 = 0;
    Double_t sig4 = 0;
   
   // arg1 = (x - position)/width;
   //   sig1 = amplitude*TMath::Exp(-0.5*arg1*arg1);
    if (par[4])//peak 1
    	arg1 = (x[0] - par[8])/par[4];
        sig1 = par[0]*TMath::Exp(-0.5*arg1*arg1);
   
    if (par[5])//peak 2
    	arg2 = (x[0] - par[9])/par[5];
        sig2 = par[1]*TMath::Exp(-0.5*arg2*arg2);

    if(par[6])//peak 3
        arg3 = (x[0] - par[10])/par[6];
        sig3 = par[2]*TMath::Exp(-0.5*arg3*arg3);

    if(par[7])//peak 4
      arg4 = (x[0] - par[11])/par[7];
      sig4 = par[3]*TMath::Exp(-0.5*arg4*arg4);

    


   
   return (sig1 + sig2 + sig3 + sig4);
}

Double_t SignalMAPneutron(Double_t *x, Double_t *par){
	//par is Amplitude starting at 3
	//x is x
    Double_t arg1 = 0;
    Double_t arg2 = 0;
    Double_t arg3 = 0;
    Double_t arg4 = 0;
    
	Double_t sig1 = 0;
    Double_t sig2 = 0;
    Double_t sig3 = 0;
    Double_t sig4 = 0;
	Double_t sig5 = 0;
	Double_t sig6 = 0;
	
	//Neutron parameters mew, Nsig, Amp
	Double_t Nsig = par[2];
	Double_t mew = par[1];
	Double_t Nmew1 = mew;
	Double_t Nmew2 = mew - dt1;
	Double_t Nmew3 = mew - dt2;
	Double_t Nmew4 = mew - dt3;
	Double_t Nmew5 = mew - dt4;
	Double_t Nmew6 = mew - dt6;
	Double_t z1z1, z2z1, z3z1;
	Double_t z1z2, z2z2, z3z2;
	Double_t z1z3, z2z3, z3z3;
	Double_t z1z4, z2z4, z3z4;
	Double_t z1z5, z2z5, z3z5;
	Double_t z1z6, z2z6, z3z6;
	//Function put forward by Dr. Daniels for single neutron peak composed of 3 gaussians
	
	//peak  = Amp1( (TMath::Exp((-1/2)*((x-mew)/(TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig)))*(x-mew)/(TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig)))) + ((TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1)* (TMath::Exp((-1/2)*((x-(mew + Gdis1t2))/(TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig)))*((x-(mew + Gdis1t2))/(TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig)))))) + ((TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1)* (TMath::Exp((-1/2)*((x-(mew + Gdis1t3))/(TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig)))*((x-(mew + Gdis1t3))/(TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig)))))))
	
	
	
	//peak1= p1()+p2a2()*p3a2()*p4a2() + p2a3()*p3a3()*p4a3())
   
   // arg1 = (x - position)/width;
   //   sig1 = amplitude*TMath::Exp(-0.5*arg1*arg1);
    if (par[0])//peak 1
    	//arg1 = ((TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1)* (TMath::Exp((-1/2)*((x[0]-(Nmew1 + Gdis1t2))/(TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig)))*((x[0]-(Nmew1 + Gdis1t2))/(TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig)))))) + ((TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1)* (TMath::Exp((-1/2)*((x[0]-(Nmew1 + Gdis1t3))/(TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig)))*((x[0]-(Nmew1 + Gdis1t3))/(TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig))))));
		//sig1  = par[0]*( (TMath::Exp((-1/2)*((x[0]-Nmew1)/(TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig)))*(x[0]-Nmew1)/(TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig)))) + arg1);
		
		z1z1 = (x[0]-Nmew1 )/TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig);
		z2z1 = (x[0]-(Nmew1 + Gdis1t2))/TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig);
		z3z1 = (x[0]-(Nmew1 + Gdis1t3))/TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig);
		Double_t Zamp2z1 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1);
		Double_t Zamp3z1 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1);
		sig1  = par[0]*(  TMath::Exp(  (-0.5)*pow(z1z1, 2)) + Zamp2z1 * TMath::Exp((-0.5)*pow(z2z1, 2)) + Zamp3z1 * TMath::Exp((-0.5)*pow(z3z1, 2)));
		
		
    if (par[3])//peak 2
       //sig2  = par[3]*( (TMath::Exp((-1/2)*((x[0]-Nmew2)/(TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig)))*(x[0]-Nmew2)/(TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig)))) + ((TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1)* (TMath::Exp((-1/2)*((x[0]-(Nmew2 + Gdis1t2))/(TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig)))*((x[0]-(Nmew2 + Gdis1t2))/(TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig)))))) + ((TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1)* (TMath::Exp((-1/2)*((x[0]-(Nmew2 + Gdis1t3))/(TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig)))*((x[0]-(Nmew2 + Gdis1t3))/(TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig)))))));
		

		z1z2 = (x[0]-Nmew2)/TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig);
		z2z2 = (x[0]-(Nmew2 + Gdis1t2))/TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig);
		z3z2 = (x[0]-(Nmew2 + Gdis1t3))/TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig);
		Double_t Zamp2z2 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1);
		Double_t Zamp3z2 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1);
		sig2  = par[3]*(  TMath::Exp(  (-0.5)*pow(z1z2, 2)) + Zamp2z2 * TMath::Exp((-0.5)*pow(z2z2, 2)) + Zamp3z2 * TMath::Exp((-0.5)*pow(z3z2, 2)));
		
    if(par[4])//peak 3
       //sig3 = par[4]*( (TMath::Exp((-1/2)*((x[0]-Nmew3)/(TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig)))*(x[0]-Nmew3)/(TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig)))) + ((TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1)* (TMath::Exp((-1/2)*((x[0]-(Nmew3 + Gdis1t2))/(TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig)))*((x[0]-(Nmew3 + Gdis1t2))/(TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig)))))) + ((TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1)* (TMath::Exp((-1/2)*((x[0]-(Nmew3 + Gdis1t3))/(TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig)))*((x[0]-(Nmew3 + Gdis1t3))/(TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig)))))));
	

		z1z3 = (x[0]-Nmew3)/TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig);
		z2z3 = (x[0]-(Nmew3 + Gdis1t2))/TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig);
		z3z3 = (x[0]-(Nmew3 + Gdis1t3))/TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig);
		Double_t Zamp2z3 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1);
		Double_t Zamp3z3 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1);
		sig3  = par[4]*(  TMath::Exp(  (-0.5)*pow(z1z3, 2)) + Zamp2z3 * TMath::Exp((-0.5)*pow(z2z3, 2)) + Zamp3z3 * TMath::Exp((-0.5)*pow(z3z3, 2)));
		
    if(par[8])//peak 4
        //sig4  = par[8]*( (TMath::Exp((-1/2)*((x[0]-Nmew4)/(TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig)))*(x[0]-Nmew4)/(TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig)))) + ((TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1)* (TMath::Exp((-1/2)*((x[0]-(Nmew4 + Gdis1t2))/(TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig)))*((x[0]-(Nmew4 + Gdis1t2))/(TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig)))))) + ((TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1)* (TMath::Exp((-1/2)*((x[0]-(Nmew4 + Gdis1t3))/(TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig)))*((x[0]-(Nmew4 + Gdis1t3))/(TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig)))))));
	
		
		z1z4 = (x[0]-Nmew4)/TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig);
		z2z4 = (x[0]-(Nmew4 + Gdis1t2))/TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig);
		z3z4 = (x[0]-(Nmew4 + Gdis1t3))/TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig);
		Double_t Zamp2z4 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1);
		Double_t Zamp3z4 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1);
		sig4  = par[8]*(  TMath::Exp(  (-0.5)*pow(z1z4, 2)) + Zamp2z4 * TMath::Exp((-0.5)*pow(z2z4, 2)) + Zamp3z4 * TMath::Exp((-0.5)*pow(z3z4, 2)));
	
	if(par[9])//peak 5
        //sig5  = par[9]*( (TMath::Exp((-1/2)*((x[0]-Nmew5)/(TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig)))*(x[0]-Nmew5)/(TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig)))) + ((TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1)* (TMath::Exp((-1/2)*((x[0]-(Nmew5 + Gdis1t2))/(TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig)))*((x[0]-(Nmew5 + Gdis1t2))/(TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig)))))) + ((TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1)* (TMath::Exp((-1/2)*((x[0]-(Nmew5 + Gdis1t3))/(TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig)))*((x[0]-(Nmew5 + Gdis1t3))/(TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig)))))));
		z1z5 = (x[0]-Nmew5)/TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig);
		z2z5 = (x[0]-(Nmew5 + Gdis1t2))/TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig);
		z3z5 = (x[0]-(Nmew5 + Gdis1t3))/TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig);
		Double_t Zamp2z5 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1);
		Double_t Zamp3z5 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1);
		sig5  = par[9]*(  TMath::Exp(  (-0.5)*pow(z1z5, 2)) + Zamp2z5 * TMath::Exp((-0.5)*pow(z2z5, 2)) + Zamp3z5 * TMath::Exp((-0.5)*pow(z3z5, 2)));
		
	if(par[7])
		z1z6 = (x[0]-Nmew6)/TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig);
		z2z6 = (x[0]-(Nmew6 + Gdis1t2))/TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig);
		z3z6 = (x[0]-(Nmew6 + Gdis1t3))/TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig);
		Double_t Zamp2z6 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1);
		Double_t Zamp3z6 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1);
		sig6  = par[7]*(  TMath::Exp(  (-0.5)*pow(z1z6, 2)) + Zamp2z6 * TMath::Exp((-0.5)*pow(z2z6, 2)) + Zamp3z6 * TMath::Exp((-0.5)*pow(z3z6, 2)));
		
  return (sig1 + sig2 + sig3 + sig4 + sig5 + sig6);

	
	
}







//psudo code under construction

Double_t g1(Double_t *x, Double_t *par){
   Double_t arg;
  if (par[2]){arg = (x[0] - par[1])/par[2];}
   return par[0]*TMath::Exp(-0.5*arg*arg);
}
Double_t g2(Double_t *x, Double_t *par){
   Double_t arg;
  if (par[2]){arg = (x[0] - par[1])/par[2];}
   return par[0]*TMath::Exp(-0.5*arg*arg);
}
Double_t g3(Double_t *x, Double_t *par){
   Double_t arg;
  if (par[2]){arg = (x[0] - par[1])/par[2];}
   return par[0]*TMath::Exp(-0.5*arg*arg);
}
Double_t g4(Double_t *x, Double_t *par){
   Double_t arg;
  if (par[2]){arg = (x[0] - par[1])/par[2];}
   return par[0]*TMath::Exp(-0.5*arg*arg);
}



// peak 1 drawing function
Double_t n1(Double_t *x, Double_t *par){
	Double_t Nsig = par[2];
	Double_t Nmew1 = par[1];
	Double_t arg;
	Double_t z1, z2, z3;
	z1 = (x[0]-Nmew1)/TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig);
	z2 = (x[0]-(Nmew1+ Gdis1t2))/TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig);
	z3 = (x[0]-(Nmew1 + Gdis1t3))/TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig);
	Double_t Zamp2 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1);
	Double_t Zamp3 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1);
	arg  = par[0]*(  TMath::Exp(  (-0.5)*pow(z1, 2)) + Zamp2 * TMath::Exp((-0.5)*pow(z2, 2)) + Zamp3 * TMath::Exp((-0.5)*pow(z3, 2)));
   return arg;
}
// peak 2 drawing function
Double_t n2(Double_t *x, Double_t *par){
	Double_t Nsig = par[2];
	Double_t Nmew1 = par[1];
	Double_t arg;
	Double_t z1, z2, z3;
	z1 = (x[0]-Nmew1)/TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig);
	z2 = (x[0]-(Nmew1+ Gdis1t2))/TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig);
	z3 = (x[0]-(Nmew1 + Gdis1t3))/TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig);
	Double_t Zamp2 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1);
	Double_t Zamp3 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1);
	arg  = par[0]*(  TMath::Exp(  (-0.5)*pow(z1, 2)) + Zamp2 * TMath::Exp((-0.5)*pow(z2, 2)) + Zamp3 * TMath::Exp((-0.5)*pow(z3, 2)));
   return arg;
}
// peak 3 drawing function
Double_t n3(Double_t *x, Double_t *par){
	Double_t Nsig = par[2];
	Double_t Nmew1 = par[1];
	Double_t arg;
	Double_t z1, z2, z3;
	z1 = (x[0]-Nmew1)/TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig);
	z2 = (x[0]-(Nmew1+ Gdis1t2))/TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig);
	z3 = (x[0]-(Nmew1 + Gdis1t3))/TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig);
	Double_t Zamp2 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1);
	Double_t Zamp3 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1);
	arg  = par[0]*(  TMath::Exp(  (-0.5)*pow(z1, 2)) + Zamp2 * TMath::Exp((-0.5)*pow(z2, 2)) + Zamp3 * TMath::Exp((-0.5)*pow(z3, 2)));
   return arg;
}
// peak 4 drawing function
Double_t n4(Double_t *x, Double_t *par){
	Double_t Nsig = par[2];
	Double_t Nmew1 = par[1];
	Double_t arg;
	Double_t z1, z2, z3;
	z1 = (x[0]-Nmew1)/TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig);
	z2 = (x[0]-(Nmew1+ Gdis1t2))/TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig);
	z3 = (x[0]-(Nmew1 + Gdis1t3))/TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig);
	Double_t Zamp2 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1);
	Double_t Zamp3 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1);
	arg  = par[0]*(  TMath::Exp(  (-0.5)*pow(z1, 2)) + Zamp2 * TMath::Exp((-0.5)*pow(z2, 2)) + Zamp3 * TMath::Exp((-0.5)*pow(z3, 2)));
   return arg;
}
// peak 5 drawing function
Double_t n5(Double_t *x, Double_t *par){
	Double_t Nsig = par[2];
	Double_t Nmew1 = par[1];
	Double_t arg;
	Double_t z1, z2, z3;
	z1 = (x[0]-Nmew1)/TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig);
	z2 = (x[0]-(Nmew1+ Gdis1t2))/TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig);
	z3 = (x[0]-(Nmew1 + Gdis1t3))/TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig);
	Double_t Zamp2 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1);
	Double_t Zamp3 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1);
	arg  = par[0]*(  TMath::Exp(  (-0.5)*pow(z1, 2)) + Zamp2 * TMath::Exp((-0.5)*pow(z2, 2)) + Zamp3 * TMath::Exp((-0.5)*pow(z3, 2)));
   return arg;
   //previous code
   //arg  = par[0]*( (TMath::Exp((-0.5)*pow((x[0]-Nmew1)/(TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig))),2)))));// + ((TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1)* (TMath::Exp((-1/2)*((x[0]-(Nmew1 + Gdis1t2))/(TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig)))*((x[0]-(Nmew1 + Gdis1t2))/(TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig)))))) + ((TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1)* (TMath::Exp((-1/2)*((x[0]-(Nmew1 + Gdis1t3))/(TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig)))*((x[0]-(Nmew1 + Gdis1t3))/(TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig)))))));
	//arg  = par[0]*TMath::Exp((-0.5)*pow((x[0]-Nmew1)/TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig), 2));
}
Double_t n6(Double_t *x, Double_t *par){
	Double_t Nsig = par[2];
	Double_t Nmew1 = par[1];
	Double_t arg;
	Double_t z1, z2, z3;
	z1 = (x[0]-Nmew1)/TMath::Sqrt(Gsig1*Gsig1+Nsig*Nsig);
	z2 = (x[0]-(Nmew1+ Gdis1t2))/TMath::Sqrt(Gsig2*Gsig2+Nsig*Nsig);
	z3 = (x[0]-(Nmew1 + Gdis1t3))/TMath::Sqrt(Gsig3*Gsig3+Nsig*Nsig);
	Double_t Zamp2 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig2*Gsig2+Nsig*Nsig)))* (Ga2/Ga1)*(Gsig2/Gsig1);
	Double_t Zamp3 = (TMath::Sqrt((Gsig1*Gsig1+Nsig*Nsig)/(Gsig3*Gsig3+Nsig*Nsig)))* (Ga3/Ga1)*(Gsig3/Gsig1);
	arg  = par[0]*(  TMath::Exp(  (-0.5)*pow(z1, 2)) + Zamp2 * TMath::Exp((-0.5)*pow(z2, 2)) + Zamp3 * TMath::Exp((-0.5)*pow(z3, 2)));
   return arg;
}

// -- Combined background + signal ----
Double_t Total(Double_t *x, Double_t *par){
    Double_t tot = Background(x, par) + Signal(x, &par[3]);
    return tot;
}
Double_t TotalGAMMA(Double_t *x, Double_t *par){
    Double_t tot = Background(x, par) + SignalMAP(x, &par[3]);
    return tot;
}
Double_t TotalNEUTRON(Double_t *x, Double_t *par){
	Double_t tot = Background(x, par) + SignalMAPneutron(x, &par[3]);
    return tot;
}
	
void fit(Int_t l, Int_t h) { // parameters are start and end of domain

  TF1 *fitfunc = new TF1("fitfunc",Total,l,h,13);
	fitfunc->SetParameter(0,1.3);
	fitfunc->FixParameter(1,0);
	fitfunc->SetParameter(2,0.04);
  //fitfunc->FixParameter(0,0);
  //fitfunc->FixParameter(1,0);
  //fitfunc->SetParameter(3,371.2);

  fitfunc->SetParameter(3,169);
  fitfunc->SetParameter(4,2823.);
		fitfunc->SetParameter(5,4);
  fitfunc->SetParameter(6,388);
  fitfunc->SetParameter(7,151);
		fitfunc->FixParameter(8,0);
		fitfunc->FixParameter(9,0);
		fitfunc->FixParameter(10,0);
  fitfunc->SetParameter(11,461);
  fitfunc->SetParameter(12,59);

  fitfunc->SetParLimits(0,0,100);
  //fitfunc->SetParLimits(1,0,10);

  fitfunc->SetParLimits(3,0,1000);
  fitfunc->SetParLimits(6,0,1000);
  fitfunc->SetParLimits(7,0,1000);
  //fitfunc->SetParLimits(7,250,300);
  //fitfunc->SetParLimits(8,0,300);
  //fitfunc->SetParLimits(9,0,300);
  fitfunc->SetParLimits(11,0,1000);
  fitfunc->SetParLimits(12,0,1000);
  fitfunc->SetParLimits(5,0,8);
  
  printf("line 282 \n");
  
  
  
  
  Data2->Fit("fitfunc","R");
 
  fitfunc->GetParameters(params);



TF1 *e1200_obj = new TF1("e1200_obj","e1200",l,h,3);
TF1 *e1000_obj = new TF1("e1000_obj","e1000",l,h,3);
TF1 *e850_obj = new TF1("e850_obj","e850",l,h,3);
TF1 *e590_obj = new TF1("e590_obj","e590",l,h,3);
TF1 *e5839_obj = new TF1("e5839_obj","e5839",l,h,3);
TF1 *e5179_obj = new TF1("e5179_obj","e5179",l,h,3);
TF1 *e1048_obj = new TF1("e1048_obj","e1048",l,h,3);
TF1 *back = new TF1("back","Background",l,h,3);
TF1 *e0260_obj = new TF1("e0260_obj","e0260",l,h,3);

e1200_obj->SetParameter(0,params[11]);
e1200_obj->SetParameter(1,params[4]-dt6);
e1200_obj->SetParameter(2,params[5]);

e1000_obj->SetParameter(0,params[3]);  //amp
e1000_obj->SetParameter(1,params[4]);  //mean
e1000_obj->SetParameter(2,params[5]);  //sigma

e850_obj->SetParameter(0,params[6]);
e850_obj->SetParameter(1,params[4]-dt1);
e850_obj->SetParameter(2,params[5]);

e590_obj->SetParameter(0,params[7]);
e590_obj->SetParameter(1,params[4]-dt2);
e590_obj->SetParameter(2,params[5]);

e5839_obj->SetParameter(0,params[8]);
e5839_obj->SetParameter(1,params[4]-dt3);
e5839_obj->SetParameter(2,params[5]);

e5179_obj->SetParameter(0,params[9]);
e5179_obj->SetParameter(1,params[4]-dt4);
e5179_obj->SetParameter(2,params[5]);

e1048_obj->SetParameter(0,params[10]);
e1048_obj->SetParameter(1,params[6]-dt5);
e1048_obj->SetParameter(2,params[5]);

e0260_obj->SetParameter(0,params[12]);
e0260_obj->SetParameter(1,params[4]-dt7);
e0260_obj->SetParameter(2,params[5]);

back->SetParameter(0,params[0]);
back->SetParameter(1,params[1]);
back->SetParameter(2,params[2]);

e1000_obj->Draw("same");
e850_obj->Draw("same");
e590_obj->Draw("same");
e5839_obj->Draw("same");
e5179_obj->Draw("same");
e1048_obj->Draw("same");
e1200_obj->SetLineColor(3);
e1200_obj->Draw("same");
back->SetLineColor(4);
back->Draw("same");
e0260_obj->SetLineColor(3);
e0260_obj->Draw("same");
//BData->Draw("Same");

}
	
void twodplotalldetectors(Int_t det, Int_t det_2, Int_t det_3, Int_t det_4, Int_t det_5, Int_t det_6, Int_t det_7, Int_t det_8, Int_t det_9, Int_t det_10, Int_t det_11, Int_t det_12, Int_t det_13, Int_t det_14, Int_t det_15){// needs 0, 1, 2 as entries to select from portion of stack
	//TCanvas *canvas = new TCanvas("canvas", "detectors 'det', 'det+1', 'det+2'", 800, 600);
	Double_t ph,psd,ch,toff;
	Double_t det_off[22] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 000, 000, 000, 000, 000, 000, 000, 0, 0, 0, 0, 0, 0};

	//Double_t det_off[22] = {0, 4, 4, -7, 9, 7, -12, 1, -30, 000, 000, 000, 000, 000, 000, 000, -3, 2, 16, 4, 8, 6};
	
	TCanvas *c1 = new TCanvas("c1","detectors",1);
	TCanvas *c2 = new TCanvas("c2","detectors",1);
	TCanvas *c3 = new TCanvas("c3","detectors",1);
	
	TCanvas *c4 = new TCanvas("c4","detectors",1);
	TCanvas *c5 = new TCanvas("c5","detectors",1);
	TCanvas *c6 = new TCanvas("c6","detectors",1);
	
	TCanvas *c7 = new TCanvas("c7","detectors",1);
	TCanvas *c8 = new TCanvas("c8","detectors",1);
	TCanvas *c9 = new TCanvas("c9","detectors",1);
	
	TCanvas *c10 = new TCanvas("c10","detectors",1);
	TCanvas *c11 = new TCanvas("c11","detectors",800,600);
	TCanvas *c12 = new TCanvas("c12","detectors",800,600);
		
	TCanvas *c13 = new TCanvas("c13","detectors",800,600);
	TCanvas *c14 = new TCanvas("c14","detectors",800,600);
	TCanvas *c15 = new TCanvas("c15","detectors",800,600);
	
	//TH2F *twod = new TH2F("2D plot","PSD cut plot", 4096,0.,12.,4096,0.,16000.);
	
	TH1F *tof_top0 = new TH1F("tof_top","Time of Flight Top", 4096,0.,4095); //hist top
	TH1F *tof_middle0 = new TH1F("tof_middle","Time of Flight Middle", 4096,0.,4095); // hist mid
	TH1F *tof_bottom0 = new TH1F("tof_bottom","Time of Flight Bottom", 4096,0.,4095); //hist bot
	
	TH1F *tof_top3 = new TH1F("tof_top3","Time of Flight Top", 4096,0.,4095); //hist top
	TH1F *tof_middle3 = new TH1F("tof_middle3","Time of Flight Middle", 4096,0.,4095); // hist mid
	TH1F *tof_bottom3 = new TH1F("tof_bottom3","Time of Flight Bottom", 4096,0.,4095); //hist bot
	
	TH1F *tof_top6 = new TH1F("tof_top6","Time of Flight Top", 4096,0.,4095); //hist top
	TH1F *tof_middle6 = new TH1F("tof_middle6","Time of Flight Middle", 4096,0.,4095); // hist mid
	TH1F *tof_bottom6 = new TH1F("tof_bottom6","Time of Flight Bottom", 4096,0.,4095); //hist bot
	
	TH1F *tof_top9 = new TH1F("tof_top9","Time of Flight Top", 4096,0.,4095); //hist top
	TH1F *tof_middle9 = new TH1F("tof_middle9","Time of Flight Middle", 4096,0.,4095); // hist mid
	TH1F *tof_bottom9 = new TH1F("tof_bottom9","Time of Flight Bottom", 4096,0.,4095); //hist bot
	
	TH1F *tof_top12 = new TH1F("tof_top12","Time of Flight Top", 4096,0.,4095); //hist top
	TH1F *tof_middle12 = new TH1F("tof_middle12","Time of Flight Middle", 4096,0.,4095); // hist mid
	TH1F *tof_bottom12 = new TH1F("tof_bottom12","Time of Flight Bottom", 4096,0.,4095); //hist bot

	TChain cch1("data");
	cch1.Add("condor_run30175.root/ndet_events");
	cch1.Add("condor_run30176.root/ndet_events");
	cch1.Add("condor_run30177.root/ndet_events");
	cch1.Add("condor_run30178.root/ndet_events");
	cch1.Add("condor_run30179.root/ndet_events");
	cch1.Add("condor_run30180.root/ndet_events");

	cch1.SetBranchAddress("ndet_PH",&ph);
	cch1.SetBranchAddress("ndet_PSD",&psd);
	cch1.SetBranchAddress("ndet_ch",&ch);
	cch1.SetBranchAddress("ndet_tof",&toff);
	
	//establish data parameters
printf("detector %d at offset# %f , %f, %f \n", det, det_off[det], det_off[det+1], det_off[det+2]);
	
	for(Int_t i = 0; i <= cch1.GetEntries(); i++){
	//for(Int_t i = 0; i <= 10; i++){
		cch1.GetEvent(i);
		printf("\r%d %f ", i, ch);
		if(ch == det){
			//twod->Fill(psd,ph);
			tof_top0->Fill(toff+det_off[det]);
		}
		if(ch == (det_2)){
			//twod->Fill(psd,ph);
			tof_middle0->Fill(toff+det_off[det_2]);
		}
		if(ch == (det_3)){
			//twod->Fill(psd,ph);
			tof_bottom0->Fill(toff+det_off[det_3]);
		}
		//****** above fills stack 0 degree
		if(ch == det_4){
			//twod->Fill(psd,ph);
			tof_top3->Fill(toff+det_off[det_4]);
		}
		if(ch == (det_5)){
			//twod->Fill(psd,ph);
			tof_middle3->Fill(toff+det_off[det_5]);
		}
		if(ch == (det_6)){
			//twod->Fill(psd,ph);
			tof_bottom3->Fill(toff+det_off[det_6]);
		}
		//****** above fills stack 3 degree
		if(ch == det_7){
			//twod->Fill(psd,ph);
			tof_top6->Fill(toff+det_off[det_7]);
		}
		if(ch == (det_8)){
			//twod->Fill(psd,ph);
			tof_middle6->Fill(toff+det_off[det_8]);
		}
		if(ch == (det_9)){
			//twod->Fill(psd,ph);
			tof_bottom6->Fill(toff+det_off[det_9]);
		}
		//****** above fills stack 6 degree
		if(ch == det_10){
			//twod->Fill(psd,ph);
			tof_top9->Fill(toff+det_off[det_10]);
		}
		if(ch == (det_11)){
			//twod->Fill(psd,ph);
			tof_middle9->Fill(toff+det_off[det_11]);
		}
		if(ch == (det_12)){
			//twod->Fill(psd,ph);
			tof_bottom9->Fill(toff+det_off[det_12]);
		}
		//****** above fills stack 9 degree
		if(ch == det_13){
			//twod->Fill(psd,ph);
			tof_top12->Fill(toff+det_off[det_13]);
		}
		if(ch == (det_14)){
			//twod->Fill(psd,ph);
			tof_middle12->Fill(toff+det_off[det_14]);
		}
		if(ch == (det_15)){
			//twod->Fill(psd,ph);
			tof_bottom12->Fill(toff+det_off[det_15]);
		}
		
		//****** above fills stack 12 degree
	}
	
	//twod->Draw();
	c1->cd();
	tof_top0->Draw();
	
	c2->cd();
	tof_middle0->SetLineColor(kRed);
	tof_middle0->Draw();
	
	c3->cd();
	tof_bottom0->SetLineColor(kGreen);
	tof_bottom0->Draw();
	
	c4->cd();
	tof_top3->Draw();
	
	c5->cd();
	tof_middle3->SetLineColor(kRed);
	tof_middle3->Draw();
	
	c6->cd();
	tof_bottom3->SetLineColor(kGreen);
	tof_bottom3->Draw();
	
	c7->cd();
	tof_top6->Draw();
	
	c8->cd();
	tof_middle6->SetLineColor(kRed);
	tof_middle6->Draw();
	
	c9->cd();
	tof_bottom6->SetLineColor(kGreen);
	tof_bottom6->Draw();
	
	c10->cd();
	tof_top9->Draw();
	
	c11->cd();
	tof_middle9->SetLineColor(kRed);
	tof_middle9->Draw();
	
	c12->cd();
	tof_bottom9->SetLineColor(kGreen);
	tof_bottom9->Draw();
	
	c13->cd();
	tof_top12->Draw();
	
	c14->cd();
	tof_middle12->SetLineColor(kRed);
	tof_middle12->Draw();
	
	c15->cd();
	tof_bottom12->SetLineColor(kGreen);
	tof_bottom12->Draw();
	
	//***** summerizing a stack of detectors
	//TH1F *tof_combine = (TH1F*)tof_top->Clone();
	
	//tof_combine->SetName("Combined TOF");
	//tof_combine->Add(tof_middle);
	//tof_combine->Add(tof_bottom);
	//tof_combine->Draw("same");
	
	//Data2->Add(tof_top);
	//Data2->Add(tof_middle);
	//Data2->Add(tof_bottom);
	//Data2->Draw();
	//some kind of break in the code is happening in the 4 lines of code above for 3, 4, 5
	//seems to be some issue for 6, 7, 8, as well
	//seem to be fixed by implamenting the explicit initalization of global Data2
	//TFile* f = new TFile("rainbow.root","Recreate");
	//f->cd();
	//Data2->Write();
	//f->Close();
	
	//*******comparison file build
	// file number 1
	//TFile* fT = new TFile("dectector1T.root","Recreate");
	//f->cd();
	//tof_top->Write();
	//f->close();
	// file number 2
	//TFile* fM = new TFile("dectector1M.root","Recreate");
	//f->cd();
	//tof_middle->Write();
	//f->Draw();	
	// file number 3
	//TFile* fB = new TFile("dectector1B.root","Recreate");
	//f->cd();
	//tof_bottom->Write();
	//f->Draw();
}
// Look for largest 7 changes + and -; Then mark the distance between features
// looking for features in small fit domain by sorting derivitive of the 
/*
void peakFIND(Int_t l, Int_t h){
	Int_t rg = l - h;
	TH1F* DataDER = new TH1F("Derivative Short Period","Change in counts over small Time of Flight Combined", 4096,0.,4095);
	TH1F* DataDER1 = new TH1F("Derivative Short medium","Change in counts over medium Time of Flight Combined", 4096,0.,4095);
	TH1F* DataDER2 = new TH1F("Derivative Short long","Change in counts over long Time of Flight Combined", 4096,0.,4095);

	Double_t DER[7] = [0, 0, 0, 0, 0, 0, 0];
	Double_t NDER[7] = [0,0,0,0,0,0,0];
	for (int i, i > 100, i++)
	{
		for (int j, j> 4096, j++) {
			
			if (Data2[j] - Data2[j++] == i){
				for (int g, g > 7, g++){
					
				}
			}
			if (Data2[j++] - Data2[j] == i){
				NDER = i; 
			}
		}
		
		DataDER = Data2[i]-Data2[i++];
		
		
		DataDER = 
		
	}
}
*/
void twodplot(Int_t det, Int_t det_2, Int_t det_3){// needs 0, 1, 2 as entries to select from portion of stack
	//TCanvas *canvas = new TCanvas("canvas", "detectors 'det', 'det+1', 'det+2'", 800, 600);
	Double_t ph,psd,ch,toff;
	
	//Double_t det_off[22] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 000, 000, 000, 000, 000, 000, 000, 0, 0, 0, 0, 0, 0};
	Double_t det_off[22] = {0, 4, 4, -7, 9, 7, -12, 1, -30, 000, 000, 000, 000, 000, 000, 000, -3, 2, 16, 4, 8, 6};
	
	TCanvas *c1 = new TCanvas("c1","detectors",1);
	
	
	//TH2F *twod = new TH2F("2D plot","PSD cut plot", 4096,0.,12.,4096,0.,16000.);
	
	TH1F *tof_top0 = new TH1F("tof_top","Time of Flight Top", 4096,0.,4095); //hist top
	TH1F *tof_middle0 = new TH1F("tof_middle","Time of Flight Middle", 4096,0.,4095); // hist mid
	TH1F *tof_bottom0 = new TH1F("tof_bottom","Time of Flight Bottom", 4096,0.,4095); //hist bot
	

	TChain cch1("data");
	cch1.Add("condor_run30175.root/ndet_events");
	cch1.Add("condor_run30176.root/ndet_events");
	cch1.Add("condor_run30177.root/ndet_events");
	cch1.Add("condor_run30178.root/ndet_events");
	cch1.Add("condor_run30179.root/ndet_events");
	cch1.Add("condor_run30180.root/ndet_events");

	cch1.SetBranchAddress("ndet_PH",&ph);
	cch1.SetBranchAddress("ndet_PSD",&psd);
	cch1.SetBranchAddress("ndet_ch",&ch);
	cch1.SetBranchAddress("ndet_tof",&toff);
	
	//establish data parameters
printf("detector %d at offset# %f , %f, %f \n", det, det_off[det], det_off[det+1], det_off[det+2]);
	
	for(Int_t i = 0; i <= cch1.GetEntries(); i++){
	//for(Int_t i = 0; i <= 10; i++){
		cch1.GetEvent(i);
		printf("\r%d %f ", i, ch);
		if(ch == det){
			//twod->Fill(psd,ph);
			tof_top0->Fill(toff+det_off[det]);
		}
		if(ch == (det_2)){
			//twod->Fill(psd,ph);
			tof_middle0->Fill(toff+det_off[det_2]);
		}
		if(ch == (det_3)){
			//twod->Fill(psd,ph);
			tof_bottom0->Fill(toff+det_off[det_3]);
		}
		
		//****** above fills stack 12 degree
	}
	
	//twod->Draw();
	c1->cd();
	tof_top0->Draw();
	tof_middle0->SetLineColor(kRed);
	tof_middle0->Draw();
	tof_bottom0->SetLineColor(kGreen);
	tof_bottom0->Draw();
	
	
	//***** summerizing a stack of detectors
	TH1F *tof_combine = (TH1F*)tof_top0->Clone();
	
	tof_combine->SetName("Combined TOF");
	tof_combine->Add(tof_middle0);
	tof_combine->Add(tof_bottom0);
	tof_combine->SetLineColor(kRed);
	tof_combine->Draw("same");
	
	Data2->Add(tof_top0);
	Data2->Add(tof_middle0);
	//Data2->Add(tof_bottom);
	//Data2->Draw();
	//some kind of break in the code is happening in the 4 lines of code above for 3, 4, 5
	//seems to be some issue for 6, 7, 8, as well
	//seem to be fixed by implamenting the explicit initalization of global Data2
	//TFile* f = new TFile("rainbow.root","Recreate");
	//f->cd();
	//Data2->Write();
	//f->Close();
	
	//*******comparison file build
	// file number 1
	//TFile* fT = new TFile("dectector1T.root","Recreate");
	//f->cd();
	//tof_top->Write();
	//f->close();
	// file number 2
	//TFile* fM = new TFile("dectector1M.root","Recreate");
	//f->cd();
	//tof_middle->Write();
	//f->Draw();	
	// file number 3
	//TFile* fB = new TFile("dectector1B.root","Recreate");
	//f->cd();
	//tof_bottom->Write();
	//f->Draw();
}

void fitgamma(Int_t l, Int_t h) { // parameters are start and end of domain

  //set 13 parameters 
  TF1 *fitfunc = new TF1("fitfunc",TotalGAMMA,l,h,15);
  
  //***********background parameters
  fitfunc->SetParameter(0,1.3);
  fitfunc->FixParameter(1,0);
  fitfunc->SetParameter(2,0.04);
  //fitfunc->FixParameter(0,0);
  //fitfunc->FixParameter(1,0);
  //********************amplitudes
  fitfunc->SetParameter(3,3000.); //first peak
  fitfunc->SetParameter(4,3500.); //second peak
  fitfunc->SetParameter(5,800.); //backgaus 
  fitfunc->FixParameter(6, 0); //shoulder 2
  //********************width
  fitfunc->SetParameter(7,0.5);
  fitfunc->SetParameter(8,1);
  fitfunc->SetParameter(9,1);
  fitfunc->FixParameter(10,2);
  //********************position
  fitfunc->SetParameter(11,2451);
  fitfunc->SetParameter(12,2454);
  fitfunc->SetParameter(13,2448);
  fitfunc->FixParameter(14,2455);

  //define limits on background
  fitfunc->SetParLimits(0,0,100);
  //fitfunc->SetParLimits(1,0,10);
  //*************************Set limits on the fits
  fitfunc->SetParLimits(3,0,4000);
  fitfunc->SetParLimits(4,0,5000);
  fitfunc->SetParLimits(5,0, 1400);
  fitfunc->SetParLimits(9, 0, 10);
  //fitfunc->SetParLimits(6,0,1000);
  //fitfunc->SetParLimits(10, 0, 10);

  printf("line 282 \n");
  
  
  
  
  Data2->Fit("fitfunc","R");
 
  fitfunc->GetParameters(params);



TF1 *g1 = new TF1("g1","g1",l,h,3);
TF1 *g2 = new TF1("g2","g2",l,h,3);
TF1 *g3 = new TF1("g3","g3",l,h,3);
TF1 *g4 = new TF1("g4","g4",l,h,3);
TF1 *back = new TF1("back","Backgroundgamma",l,h,3);


g1->SetParameter(0,params[3]);
g1->SetParameter(1,params[11]);
g1->SetParameter(2,params[7]);

g2->SetParameter(0,params[4]);  //amp
g2->SetParameter(1,params[12]);  //mean
g2->SetParameter(2,params[8]);  //sigma

g3->SetParameter(0,params[5]);
g3->SetParameter(1,params[13]);
g3->SetParameter(2,params[9]);

g4->SetParameter(0,params[6]);
g4->SetParameter(1,params[14]);
g4->SetParameter(2,params[10]);

back->SetParameter(0,params[0]);
back->SetParameter(1,params[1]);
back->SetParameter(2,params[2]);

g1->SetLineColor(1);
g1->Draw("same");
g2->SetLineColor(2);
g2->Draw("same");
g3->SetLineColor(3);
g3->Draw("same");
g4->SetLineColor(3);
g4->Draw("same");
back->SetLineColor(4);
back->Draw("same");
//BData->Draw("Same");
}

//fit for neutron signal comprised of multiple gaussians from the gamma peak
void fitNEUTRON(Int_t l, Int_t h) { // parameters are start and end of domain
 TF1 *fitfunc = new TF1("fitfunc",TotalNEUTRON,l,h,13);

	//fitfunc->SetParameter(0,2.74);//background parameters
	//fitfunc->FixParameter(1,0);
	//fitfunc->SetParameter(2,0.04);
	
	fitfunc->FixParameter(0,10.);//background parameters
	fitfunc->FixParameter(1,0.);
	fitfunc->FixParameter(2,0.);
  //fitfunc->FixParameter(0,0);
  //fitfunc->FixParameter(1,0);
  //fitfunc->SetParameter(3,371.2);

  fitfunc->SetParameter(3,50.);//amp 1
  fitfunc->SetParameter(4,2830.);//mew
	fitfunc->SetParameter(5,1.);//sigma
  fitfunc->SetParameter(6,30.);//amp 2
  fitfunc->SetParameter(7,10.);//amp 3
		fitfunc->FixParameter(8,0);//empty
		fitfunc->FixParameter(9,0);//empty
  fitfunc->SetParameter(10,0.0);//amp 6 temp empty 
  fitfunc->SetParameter(11,10.);//amp 4
  fitfunc->SetParameter(12,10.);//amp 5

  //fitfunc->SetParLimits(0,0,100);
  fitfunc->SetParLimits(0,0,30);
  fitfunc->SetParLimits(4, 2810, 2840); //mew limits
  //fitfunc->SetParLimits(3,0,1000);
  //fitfunc->SetParLimits(6,0,1000);
  
  fitfunc->SetParLimits(3,0,500);// set parameter limits for all peaks for 0 to 1000
  fitfunc->SetParLimits(6,0,500);
  fitfunc->SetParLimits(7,0,500);
  fitfunc->SetParLimits(11,0,500);
  fitfunc->SetParLimits(12,0,500);
  
  //fitfunc->SetParLimits(7,250,300);
  //fitfunc->SetParLimits(11,0,1000);
  //fitfunc->SetParLimits(12,0,1000);
  fitfunc->SetParLimits(5,0,10);
   
  
  
  
  Data2->Fit("fitfunc","R");
 
  fitfunc->GetParameters(params);



TF1 *n1 = new TF1("n1","n1",l,h,3);
TF1 *n2 = new TF1("n2","n2",l,h,3);
TF1 *n3 = new TF1("n3","n3",l,h,3);
TF1 *n4 = new TF1("n4","n4",l,h,3);
TF1 *n5_obj = new TF1("n5_line","n5",l,h,3);
TF1 *n6_obj = new TF1("n6_line","n6",l,h,3);

//TF1 *e5179_obj = new TF1("e5179_obj","e5179",l,h,3);
//TF1 *e1048_obj = new TF1("e1048_obj","e1048",l,h,3);
TF1 *back = new TF1("back","Background",l,h,3);
//TF1 *e0260_obj = new TF1("e0260_obj","e0260",l,h,3);

n1->SetParameter(0,params[3]);
n1->SetParameter(1,params[4]);
n1->SetParameter(2,params[5]);

n2->SetParameter(0,params[6]);  //amp
n2->SetParameter(1,params[4]-dt1);  //mean
n2->SetParameter(2,params[5]);  //sigma

n3->SetParameter(0,params[7]);
n3->SetParameter(1,params[4]-dt2);
n3->SetParameter(2,params[5]);

n4->SetParameter(0,params[11]);
n4->SetParameter(1,params[4]-dt3);
n4->SetParameter(2,params[5]);

n5_obj->SetParameter(0,params[12]);
n5_obj->SetParameter(1,params[4]-dt4);
n5_obj->SetParameter(2,params[5]);

n6_obj->SetParameter(0,params[10]);
n6_obj->SetParameter(1,params[4]-dt6);
n6_obj->SetParameter(2,params[5]);

//n5_obj->SetParameter(0,461.);
//n5_obj->SetParameter(1,2800.);
//n5_obj->SetParameter(2,0.01);


back->SetParameter(0,params[0]);
back->SetParameter(1,params[1]);
back->SetParameter(2,params[2]);

n1->SetLineColor(1);
n1->Draw("same");
n2->SetLineColor(2);
n2->Draw("same");
n3->SetLineColor(3);
n3->Draw("same");
n4->SetLineColor(4);
n4->Draw("same");
n5_obj->SetLineColor(2);
n5_obj->Draw("same");
n6_obj->SetLineColor(2);
n6_obj->Draw("same");

back->SetLineColor(6);
back->Draw("same");
//BData->Draw("Same");
}
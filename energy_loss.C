#include <fstream>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include "TCutG.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSpectrum.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include "TVirtualFitter.h"
#include "TPolyMarker.h"
#include "eloss.cxx"
Int_t energy_loss()
{

ofstream sd1sim;
sd1sim.open("8Heexc.txt",ios::out);           //output

   // Stopping power table generated with LISE++
    std::string eLossFileP = "/home/gurmukh/simIris/simIris-master/LISE_dedx/lise_8He_in_P.txt";
    std::string eLossFileAl = "/home/gurmukh/simIris/simIris-master/LISE_dedx/lise_8He_in_Al.txt";
    std::string eLossFileAg = "/home/gurmukh/simIris/simIris-master/LISE_dedx/lise_8He_in_Ag.txt";
    std::string eLossFileSi = "/home/gurmukh/simIris/simIris-master/LISE_dedx/lise_8He_in_Si.txt";
    std::string eLossFileB = "//home/gurmukh/simIris/simIris-master/LISE_dedx/lise_8He_in_B.txt";
    std::string eLossFileSiO2 = "/home/gurmukh/simIris/simIris-master/LISE_dedx/lise_8He_in_SiO2.txt";
     std::string eLossFileD = "/home/gurmukh/simIris/simIris-master/LISE_dedx/lise_8He_in_D.txt";
     



    Double_t eD[100], dedxD[100];
    Double_t eAg[100], dedxAg[100];
    
    Double_t eAl[100], dedxAl[100];
    Double_t eB[100], dedxB[100];
    Double_t eSi[100], dedxSi[100];
    Double_t eP[100], dedxP[100];
    Double_t eSiO2[100], dedxSiO2[100];
    
    
    Double_t u=931.494061;
    Double_t mass1 = 8.033934390*u;   /// remeber to change the mass for different nucleus

    loadELoss(eLossFileAg,eAg,dedxAg,mass1);
    loadELoss(eLossFileAl,eAl,dedxAl,mass1);
    loadELoss(eLossFileP,eP,dedxP,mass1);
    loadELoss(eLossFileB,eB,dedxB,mass1);
    loadELoss(eLossFileSi,eSi,dedxSi,mass1);
    loadELoss(eLossFileD,eD,dedxD,mass1);
    
    loadELoss(eLossFileSiO2,eSiO2,dedxSiO2,mass1);


    double dE,dE1,E,E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11;
    
TString sd1s="cc.txt";           //input

FILE *sd1;
FILE *csiFile2;
    Int_t Chan=-1;
    double a,e,Chan1,a1,c1;
    Int_t r;
    Int_t ch=1001,ring[1001];
    double energy[1001],angle[1001];
    char buffer[1001];
       sd1 = fopen(sd1s.Data(), "r"); 
    printf("Reading input file '%s'\n",sd1s.Data());

     for (int i =0;i<1001;i++)
{
        fscanf(sd1,"%lf%lf",&a,&e);
        angle[i] = 1.*a;
        energy[i] = 1.*e;
      //  cout<<e<<endl;
        //printf("%d \t angle: %lf\t energy: %lf\n",i,angle[i],energy[i]);
        //printf(": %f\n",energy[i]);
 }



    double theta,M, S3Ring, hypot;
    TRandom rand;
    
    for(Int_t j =0;j<1001;j++)
    {
    
        E = energy[j];     
        theta = angle[j];
    
   
    //E = E-eloss(E,0.1*10.473*2.15/cos(theta*3.14/180),eAg,dedxAg); for silver runs
      //  cout<<E<<endl;
     E = E-eloss(E,0.1*0.201*27.5/cos(theta*3.14/180),eD,dedxD); // chnage for different thickness(half of thickness)
     E = E-eloss(E,0.1*2.702*1.5/cos(theta*3.14/180),eAl,dedxAl);
     E = E-eloss(E,0.1*2.32*3.5/cos(theta*3.14/180),eD,dedxSiO2);
     E = E-eloss(E,0.1*2.702*0.3/cos(theta*3.14/180),eAl,dedxAl);
  
     E = E-eloss(E,0.1*2.3502*0.5/cos(theta*3.14/180),eB,dedxB); //boron junction implant
    
  E1=E;
     E = E-eloss(E,0.1*2.3212*60/cos(theta*3.14/180),eSi,dedxSi);
    E2=E;
      
     E = E-eloss(E,0.1*1.8219*0.5/cos(theta*3.14/180),eP,dedxP);
    
     E = E-eloss(E,0.1*2.702*0.3/cos(theta*3.14/180),eAl,dedxAl);
     E = E-eloss(E,0.1*2.702*0.3/cos(theta*3.14/180),eAl,dedxAl);
     E = E-eloss(E,0.1*1.8219*0.5/cos(theta*3.14/180),eP,dedxP);
     E4=E;
      E = E-eloss(E,0.1*2.3212*500/cos(theta*3.14/180),eSi,dedxSi);
     
E3=E1-E2;

     

    //sd1sim<<theta<<"\t"<<E4<<"\t"<<E1<<endl;
        sd1sim<<theta<<"\t"<<E4<<endl;
  
   }

  sd1sim.close();
return 0;
}

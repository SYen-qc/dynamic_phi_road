#define xAna_new_cxx
#include "xAna_new.h"
#include <iostream>
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TProfile.h>
#include <TLine.h>
#include <TLorentzVector.h>

using namespace std;

void xAna_new:Loop()
{  
  //Gstyle->Setoptstat(111111);
  
  ifstream in;
  
  Char_t fname[200];
  Double_t es1=0,es2=0,es3=0,es4=0,res1=0,res2=0,res3=0,res4=0,res5=0,res6=0,res7=0,res8=0;
  Double_t Res1=0,Res2=0,Res3=0,Res4=0;
  //Double_t rhefP=0,rhefM=0,rherP=0,rherM=0;
  Double_t Ratio1=0,Ratio2=0,Ratio3=0,Ratio4=0;
  //Double_t sRatio1=0,sRatio2=0;
  Double_t dPhi=0,dEta=0,tempdPhi=0,tempdEta=0,tempresP1=0,tempresP2=0;
  Double_t deltaP=0,deltaE=0;
  Double_t stEta=0,stPhi=0,rhEta=0,rhPhi=0;
  Double_t dEold=0,dE024=0,dE020=0,dE021=0,dE022=0,dE023=0,dE025=0,dE026=0,dE027=0,dE028=0,dE029=0,dEnew=0,dEall=0,dE010=0,dE015=0,dE035=0,dE040=0,dE050=0,dEE=0,dEoldnew=0,dEalln=0;
  Int_t bin=0,bin1=0;
  Double_t saber=0,archer=0,lancer=0;
  Double_t suigintou=0,kanaria=0,suiseiseki=0;
  //Double_t r=0.12,s=0.15;
  Double_t MIP=0;
  Double_t Eold=0,Enew=0;
  Double_t f=0;
 
  Int_t Z, P, X, Y, S;
  Double_t t, p, x, y, z;
  Double_t eta[2][2][40][40][32];
  Double_t phi[2][2][40][40][32];

  sprintf(fname, "/data2/cmkuo/ES/preshower_strip_position_321.dat");
  in.open(fname);
  Int_t nlines = 0;
  while(1){
    in >> Z >> P >> X >> Y >> S >> t >> p >> x >> y >> z;
    if(Z == -1) Z = 0;
    if(Z == 1) Z = 1;
    if(!in.good()) break;
    eta[Z][P-1][X-1][Y-1][S-1] = t;
    phi[Z][P-1][X-1][Y-1][S-1] = p;
    nlines++;
  }

  if (fChain == 0) return;

  //sprintf(fname,"/afs/cern.ch/user/s/stseng/outputEX/reducedRH_pt%d_eta%d.root",PT_,ETA_);
  sprintf(fname,"/afs/cern.ch/user/s/stseng/outputEX/brem_pt%d_eta%d.root",PT_,ETA_);
  //sprintf(fname,"/home/Prisemriber/output/dphi_0.15_brem_pt%d_eta%d.root",PT_,ETA_);
  //sprintf(fname,"/home/Prisemriber/outputEX/pu_brem_pt%d_eta%.root",PT_,ETA_);
  //sprintf(fname,"/home/Prisemriber/output/MIP3_brem_pt%d_eta%d.root",PT_,ETA_);
  fout_ = new TFile(fname,"RECREATE");
  InitHists();

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;

  /*for (Long64_t jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    rhefP=0, rhefM=0, rherP=0, rherM=0;
    for (int k=0; k<nRH; ++k){
    if (rhZ[k]==1 && rhP[k]==1)
    rhefP += rhE[k];
    if (rhZ[k]==-1 && rhP[k]==1)
    rhefM += rhE[k];
    if (rhZ[k]==1 && rhP[k]==2)
    rherP += rhE[k];
    if (rhZ[k]==-1 && rhP[k]==2)
    rherM += rhE[k];
    }
    }*/
  
  ofstream fout;
  //sprintf(fname,"/home/Prisemriber/output/dEdPR1R2.txt");
  //fout.open(fname);
  TProfile *hProf1 = new TProfile("hProf1", "", 100, 0, 2900, 0, 800);
  TProfile *hProf2 = new TProfile("hProf2", "", 100, 0, 2900, 0, 800);
  TProfile *hProf3 = new TProfile("hProf3", "", 100, 0, 2900, 0, 800);
  
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    //for (Long64_t jentry=0; jentry<10;jentry++){
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    for (int i=0; i<nEle; ++i) {
      if (fabs(eleSCEta[i]) < 1.65) continue;
      if (fabs(eleSCEta[i]) > 2.5) continue;
      es1  = 0;
      es2  = 0;
      res1 = 0;
      res2 = 0;
      res3 = 0;
      res4 = 0;
      Res1 = 0;
      Res2 = 0;
      Res3 = 0;
      Res4 = 0;
      H1_->Reset();
      H2_->Reset();
      H3_->Reset();
      H4_->Reset();
      
      for (int j=0; j<nSC; ++j) {
        //cout << nBC[j] << endl;
        Int_t Matched =0;
        for (int n=0; n<nMC; n++) {
          Float_t dR = deltaR(scEta[j], scPhi[j], mcEta[n], mcPhi[n]);
          Float_t dPt = fabs(elePt[i] - mcPt[n])/mcPt[n];
	  if (dR<0.1 && dPt < 0.2) Matched=1;
        }
	if (Matched==0) continue; 
	
	if (st1x[j][0]==0 || st2x[j][0]==0) continue;
	if (fabs(eleSCEta[i] - scEta[j]) < 0.001 && fabs(eleSCPhi[i] - scPhi[j]) < 0.001) {
	  Float_t resmaxp1 = -9999;
	  Float_t resmaxn1 = -9999;
	  Float_t resmaxp2 = -9999;
	  Float_t resmaxn2 = -9999;
	  Float_t tempmaxP = -9999;
	  Float_t tempmaxE = -9999;
	  Float_t tempminP =  9999;
	  Float_t tempminE =  9999;
	  if (st1z[j][0]==-1) st1z[j][0]=0;
	  if (st2z[j][0]==-1) st2z[j][0]=0;

	  //ratio_old->Fill(scEn[j]/(PT_*cosh(ETA_/10.))); 
	  //f = scEn/(scRawEn[j]+scESEn[j]);
	  
	  for (int k=0; k<nBC[j]; ++k) {
	    tempdPhi = bcPhi[j][0] - bcPhi[j][k];
	    tempdEta = bcEta[j][0] - bcEta[j][k];
	    if (tempdPhi > TMath::Pi()) tempdPhi = tempdPhi - 2.*TMath::Pi();
	    if (tempdPhi < -TMath::Pi()) tempdPhi += 2.*TMath::Pi();
	    if (tempmaxP<tempdPhi) tempmaxP=tempdPhi;
	    if (tempmaxE<tempdEta) tempmaxE=tempdEta;
	    if (tempminP>tempdPhi) tempminP=tempdPhi;
	    if (tempminE>tempdEta) tempminE=tempdEta;
	    //if (fabs(tempdPhi)>0.28)
	    //cout << bcPhi[j][0] << " " << bcPhi[j][k] << " " << tempdPhi << endl;
	    for (int m=0; m<4; ++m) {
	      if (bcESEnX[j][k][m] > 0) es1 += bcESEnX[j][k][m];
	      if (bcESEnY[j][k][m] > 0) es2 += bcESEnY[j][k][m];
	    }
	  }  
	  
	  deltaP = tempmaxP-tempminP;
	  if (deltaP > TMath::Pi()) deltaP = -(deltaP - 2.*TMath::Pi());
	  if (deltaP < -TMath::Pi()) deltaP += 2.*TMath::Pi();
	  deltaE = tempmaxE-tempminE;
	  //fout << tempmaxE << " " << tempminE << " " << deltaE << endl;
	  //printf("max=%f, min=%f, dPhi=%f\n",tempmaxE,tempminE,deltaP);
          //cout<<tempmaxP<<" "<<tempminP<<" "<<tempmaxE<<" "<<tempminE<<endl;
	  if (nBC[j]==1) HH1_->Fill(deltaP+0.2);
	  if (nBC[j]==2) HH2_->Fill(deltaP+0.2);
	  if (nBC[j]==3) HH3_->Fill(deltaP+0.2);
	  if (nBC[j]==4) HH4_->Fill(deltaP+0.2);
	  if (nBC[j]==5) HH5_->Fill(deltaP+0.2);
	  if (nBC[j]==6) HH6_->Fill(deltaP+0.2);
	  if (nBC[j]==7) HH7_->Fill(deltaP+0.2);
	  if (nBC[j]==8) HH8_->Fill(deltaP+0.2);
	  //H_->Fill(deltaP,deltaE);
	  
	  for (int l=0; l<nRH; ++l) {
	    if (rhFlag[l]==1||rhFlag[l]==14||(rhFlag[l]<=10&&rhFlag[l]>=5)) continue;
	    if (rhZ[l]*scEta[j]>0 && rhP[l]==1){
              if (rhE[l]/(92.8*0.000001)>=MIP)  res3 += rhE[l];
	      if (rhZ[l]==1){
		NRHP1Zp_->Fill(nRH);
		stEta=eta[st1z[j][0]][st1p[j][0]-1][st1x[j][0]-1][st1y[j][0]-1][st1s[j][0]-1];
		rhEta=eta[rhZ[l]][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dEta = -stEta+rhEta;
		stPhi=phi[st1z[j][0]][st1p[j][0]-1][st1x[j][0]-1][st1y[j][0]-1][st1s[j][0]-1];
		rhPhi=phi[rhZ[l]][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dPhi = -stPhi+rhPhi;
		if (dPhi > TMath::Pi()) dPhi = dPhi - 2.*TMath::Pi();
		if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
	      }
	      if (rhZ[l]==-1){
		NRHP1Zn_->Fill(nRH);
		//rhZ[l]=0;
		stEta=eta[st1z[j][0]][st1p[j][0]-1][st1x[j][0]-1][st1y[j][0]-1][st1s[j][0]-1];
		rhEta=eta[0][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dEta = -stEta+rhEta;
		stPhi=phi[st1z[j][0]][st1p[j][0]-1][st1x[j][0]-1][st1y[j][0]-1][st1s[j][0]-1];
		rhPhi=phi[0][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dPhi = -stPhi+rhPhi;
		if (dPhi > TMath::Pi()) dPhi = dPhi - 2.*TMath::Pi();
		if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
	      }
	      if ((0<=dPhi && dPhi<=(tempmaxP+r_) && fabs(dEta)<=s_) 
		  || ((tempminP-r_)<=dPhi && dPhi<=0 && fabs(dEta)<=s_)){
		if (rhE[l]/(92.8*0.000001)>=MIP)  res1 += rhE[l];
		//cout <<rhZ[l]<<" "<<rhP[l]<<" "<<dEta<<" "<<dPhi<<" " <<"it is picked up"<<" "<<tempmaxE+s<<" "<<tempminE-s<<" "<<tempmaxP+r<<" "<<tempminP-r<<endl;
	      }
              //else cout <<rhZ[l]<<" "<<rhP[l]<<" "<<dEta<<" "<<dPhi<<" " <<"it is excluded"<<" "<<tempmaxE+s<<" "<<tempminE-s<<" "<<tempmaxP+r<<" "<<tempminP-r<<endl;
	      if (rhZ[l]==1) H1_->Fill(dEta,dPhi,rhE[l]/92.8e-6);
	      if (rhZ[l]==-1) H2_->Fill(dEta,dPhi,rhE[l]/92.8e-6);
	    }
	    if (rhZ[l]*scEta[j]>0 && rhP[l]==2){
              if (rhE[l]/(92.8*0.000001)>=MIP)  res4 += rhE[l];
	      if (rhZ[l]==1){
		NRHP2Zp_->Fill(nRH);
		stEta=eta[st2z[j][0]][st2p[j][0]-1][st2x[j][0]-1][st2y[j][0]-1][st2s[j][0]-1];
		rhEta=eta[rhZ[l]][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dEta = -stEta+rhEta;
		stPhi=phi[st2z[j][0]][st2p[j][0]-1][st2x[j][0]-1][st2y[j][0]-1][st2s[j][0]-1];
		rhPhi=phi[rhZ[l]][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dPhi = -stPhi+rhPhi;
		if (dPhi > TMath::Pi()) dPhi = dPhi - 2.*TMath::Pi();
		if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
	      }
	      if (rhZ[l]==-1){
		NRHP2Zn_->Fill(nRH);
		//rhZ[l]=0;
		stEta=eta[st2z[j][0]][st2p[j][0]-1][st2x[j][0]-1][st2y[j][0]-1][st2s[j][0]-1];
		rhEta=eta[0][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dEta = -stEta+rhEta;
		stPhi=phi[st2z[j][0]][st2p[j][0]-1][st2x[j][0]-1][st2y[j][0]-1][st2s[j][0]-1];
		rhPhi=phi[0][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dPhi = -stPhi+rhPhi;
		if (dPhi > TMath::Pi()) dPhi = dPhi - 2.*TMath::Pi();
		if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
	      } 
	      if ((0<=dPhi && dPhi<=(tempmaxP+r_) && fabs(dEta)<=s_) 
		  || ((tempminP-r_)<=dPhi && dPhi<=0 && fabs(dEta)<=s_)){
		if (rhE[l]/(92.8*0.000001)>=MIP)  res2 += rhE[l];
		//cout <<rhZ[l]<<" "<<rhP[l]<<" "<<" "<<dEta<<" "<<dPhi<<" " <<"it is picked up"<<" "<<tempmaxE+s<<" "<<tempminE-s<<" "<<tempmaxP+r<<" "<<tempminP-r<<endl;
	      }  
	      //else cout <<rhZ[l]<<" "<<rhP[l]<<" "<<" "<<dEta<<" "<<dPhi<<" " <<"it is excluded"<<" "<<tempmaxE+s<<" "<<tempminE-s<<" "<<tempmaxP+r<<" "<<tempminP-r<<endl;
	      if (rhZ[l]==1)  H3_->Fill(dEta,dPhi,rhE[l]/92.8e-6);
	      if (rhZ[l]==-1) H4_->Fill(dEta,dPhi,rhE[l]/92.8e-6);
	    }
	  }

          tempresP1=res3-res1;
          tempresP2=res4-res2;
	  if ((resmaxp1<tempresP1)&&scEta[j]/fabs(scEta[j])==1) resmaxp1=tempresP1;
	  //if (res3==0&&res1==0&&scEta[j]/fabs(scEta[j])==1) resmaxp1=0;
	  if ((resmaxn1<tempresP1)&&scEta[j]/fabs(scEta[j])==-1) resmaxn1=tempresP1;
	  //if (res3==0&&res1==0&&scEta[j]/fabs(scEta[j])==-1) resmaxn1=0;
	  if ((resmaxp2<tempresP2)&&scEta[j]/fabs(scEta[j])==1) resmaxp2=tempresP2;
	  //if (res4==0&&res2==0&&scEta[j]/fabs(scEta[j])==1) resmaxp2=0;
	  if ((resmaxn2<tempresP2)&&scEta[j]/fabs(scEta[j])==-1) resmaxn2=tempresP2;
	  //if (res4==0&&res2==0&&scEta[j]/fabs(scEta[j])==-1) resmaxn2=0;
	  
	  for (int a=0; a<nRH; ++a) {
            if (rhFlag[a]==1||rhFlag[a]==14||(rhFlag[a]<=10&&rhFlag[a]>=5)) continue;
	    if (rhZ[a]==1 && rhP[a]==1){ 
	      if (rhE[a]/(92.8*0.000001)>=MIP)  Res1 += rhE[a];
	    }
	    if (rhZ[a]==1 && rhP[a]==2){
	      if (rhE[a]/(92.8*0.000001)>=MIP)  Res3 += rhE[a];
	    }
	    if (rhZ[a]==-1 && rhP[a]==1){
	      if (rhE[a]/(92.8*0.000001)>=MIP)  Res2 += rhE[a];
	    }
	    if (rhZ[a]==-1 && rhP[a]==2){
	      if (rhE[a]/(92.8*0.000001)>=MIP)  Res4 += rhE[a];
	    }
	  }  
	  
	  if (eleSCEta[i]>0){
	    //Ratio1 = res1/shES_efP;
	    Ratio1 = es1/res1;
	    //Ratio2 = res2/shES_erP;
	    Ratio3 = es2/res2;
	  }
	  else if (eleSCEta[i]<0){
	    //Ratio1 = res1/shES_efM;
	    Ratio2 = es1/res1;
	    //Ratio2 = res2/shES_erM;
	    Ratio4 = es2/res2;
	  }
	  
	  if (eleBrem[i]<0) continue;
	  bin = eleBrem[i]/0.25;
	  dEold = scRawEn[j] + scESEn[j] - PT_*cosh(ETA_/10.);
	  h1_->Fill(dEold);
	  dE024 = scRawEn[j] + 0.020*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h3_->Fill(dE024);
	  dE020 = scRawEn[j] + 0.020*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h5_->Fill(dE020);
	  dE021 = scRawEn[j] + 0.021*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h6_->Fill(dE021);
	  dE022 = scRawEn[j] + 0.022*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h7_->Fill(dE022);
	  dE023 = scRawEn[j] + 0.023*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h8_->Fill(dE023);
	  dE025 = scRawEn[j] + 0.025*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h9_->Fill(dE025);
	  dE026 = scRawEn[j] + 0.026*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h10_->Fill(dE026);
	  dE027 = scRawEn[j] + 0.027*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h11_->Fill(dE027);
	  dE028 = scRawEn[j] + 0.028*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h12_->Fill(dE028);
	  dE029 = scRawEn[j] + 0.029*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h13_->Fill(dE029);
	  dE010 = scRawEn[j] + 0.01*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h14_->Fill(dE010);
	  dE015 = scRawEn[j] + 0.015*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h15_->Fill(dE015);
	  dE035 = scRawEn[j] + 0.035*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h16_->Fill(dE035);
	  dE040 = scRawEn[j] + 0.04*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h17_->Fill(dE040);
	  dE050 = scRawEn[j] + 0.05*((res1/(92.8e-6))+(0.7*(res2/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h19_->Fill(dE050);
	  dEE = scRawEn[j] - PT_*cosh(ETA_/10.);
	  h18_->Fill(dEE);
	  hfb1_[bin]->Fill(dEold);
	  hfb2_[bin]->Fill(dE024);
	  hfb3_[bin]->Fill(dE020);
	  hfb4_[bin]->Fill(dE021);
	  hfb5_[bin]->Fill(dE022);
	  hfb6_[bin]->Fill(dE023);
	  hfb7_[bin]->Fill(dE025);
	  hfb8_[bin]->Fill(dE026);
	  hfb9_[bin]->Fill(dE027);
	  hfb10_[bin]->Fill(dE028);
	  hfb11_[bin]->Fill(dE029);
	  hfb13_[bin]->Fill(dE010);
	  hfb14_[bin]->Fill(dE015);
	  hfb15_[bin]->Fill(dE035);
	  hfb16_[bin]->Fill(dE040);
	  hfb17_[bin]->Fill(dEE);
	  hfb18_[bin]->Fill(dE050);
	  hProf1->Fill((es1/(92.8e-6))+(0.7*(es2/(92.8e-6))), scEn[j]-scESEn[j]);
	  hProf2->Fill((res1/(92.8e-6))+(0.7*(res2/(92.8e-6))), scRawEn[j]);
	  hProf3->Fill((res3/(92.8e-6))+(0.7*(res4/(92.8e-6))), scRawEn[j]);
	  potato_old_->Fill((es1/(92.8e-6))+(0.7*(es2/(92.8e-6))), scEn[j]-scESEn[j]);
	  potato_new_->Fill((res1/(92.8e-6))+(0.7*(res2/(92.8e-6))), scRawEn[j]);
	  potato_all_->Fill((res3/(92.8e-6))+(0.7*(res4/(92.8e-6))), scRawEn[j]);
          Eh1_->Fill(Ratio1);
          Eh2_->Fill(Ratio2);
          Eh3_->Fill(Ratio3);
          Eh4_->Fill(Ratio4);
          
          /*TCanvas *c2 = new TCanvas("c2","c2",1300,650);
	  gStyle->SetOptStat(0);
	  gStyle->SetPadLeftMargin(0.14);
	  gStyle->SetPadRightMargin(0.15);
	  c2->Divide(2,1);
	  c2->cd(1);
	  H1_->SetTitle("");
	  H1_->GetXaxis()->SetTitle("#Delta#eta");
	  H1_->GetYaxis()->SetTitle("#Delta#phi");
	  H2_->SetTitle("");
	  H2_->GetXaxis()->SetTitle("#Delta#eta");
	  H2_->GetYaxis()->SetTitle("#Delta#phi");
	  if (scEta[j]/fabs(scEta[j])==1) H1_->Draw("colz");
	  if (scEta[j]/fabs(scEta[j])==-1) H2_->Draw("colz");
	  TLine *L1;
	  L1 = new TLine (s, tempminP-r, s, tempmaxP+r); L1->Draw();
	  L1 = new TLine (-s, tempminP-r, -s, tempmaxP+r); L1->Draw();
	  L1 = new TLine (s, tempminP-r, -s, tempminP-r); L1->Draw();
	  L1 = new TLine (s, tempmaxP+r, -s, tempmaxP+r); L1->Draw();
	  c2->cd(2);
	  H3_->SetTitle("");
	  H3_->GetXaxis()->SetTitle("#Delta#eta");
	  H3_->GetYaxis()->SetTitle("#Delta#phi");
	  H4_->SetTitle("");
	  H4_->GetXaxis()->SetTitle("#Delta#eta");
	  H4_->GetYaxis()->SetTitle("#Delta#phi");
	  if (scEta[j]/fabs(scEta[j])==1) H3_->Draw("colz");
	  if (scEta[j]/fabs(scEta[j])==-1) H4_->Draw("colz");
	  TLine *L2;
	  L2 = new TLine (s, tempminP-r, s, tempmaxP+r); L2->Draw();
	  L2 = new TLine (-s, tempminP-r, -s, tempmaxP+r); L2->Draw();
	  L2 = new TLine (s, tempminP-r, -s, tempminP-r); L2->Draw();
	  L2 = new TLine (s, tempmaxP+r, -s, tempmaxP+r); L2->Draw();
          c2->Print("event_display.pdf");
	  cin.get(); 
	  delete c2; 
	  delete L1;
	  delete L2;*/
          
	  /*if (es1/res1>1||es2/res2>1){
	  //cout<< deltaE <<" "<< deltaP <<endl;
	  TCanvas *c1 = new TCanvas("c1","c1", 1300,650);
	  gStyle->SetOptStat(0);
	  gStyle->SetPadLeftMargin(0.14);
	  gStyle->SetPadRightMargin(0.15);
	  c1->Divide(2,1);
	  c1->cd(1);
	  H1_->SetTitle("");
	  H1_->GetXaxis()->SetTitle("#Delta#eta");
	  H1_->GetYaxis()->SetTitle("#Delta#phi");
	  H2_->SetTitle("");
	  H2_->GetXaxis()->SetTitle("#Delta#eta");
	  H2_->GetYaxis()->SetTitle("#Delta#phi");
	  //if (scEta[j]/fabs(scEta[j])==1&&(res3-res1)/92.8e-6>35&&(res4-res2)/92.8e-6>35){
	  if (scEta[j]/fabs(scEta[j])==1){
	  cout<<scEta[j]/fabs(scEta[j])<<" "<<deltaP<<" "<<deltaE<<" "<<es1/92.8e-6<<" "<<res1/92.8e-6<< " "<<es1/res1<<" "<<res3/92.8e-6<<" "<<es2/92.8e-6<<" "<<res2/92.8e-6<<" "<<es2/res2<<" "<<res4/92.8e-6<<endl;
	  //fout<<scEta[j]/fabs(scEta[j])<<" "<<deltaP<<" "<<deltaE<<" "<<es1<<" "<<res1<< " "<<es1/res1<<" "<<res3<<" "<<es2<<" "<<res2<<" "<<es2/res2<<" "<<res4<<endl;
	  //H1_->Draw("colztext");
	  H1_->Draw("colz");
	  } 
	  //if (scEta[j]/fabs(scEta[j])==-1&&(res3-res1)/92.8e-6>35&&(res4-res2)/92.8e-6>35){
	  if (scEta[j]/fabs(scEta[j])==-1){
	  cout<<scEta[j]/fabs(scEta[j])<<" "<<deltaP<<" "<<deltaE<<" "<<es1/92.8e-6<<" "<<res1/92.8e-6<< " "<<es1/res1<<" "<<res3/92.8e-6<<" "<<es2/92.8e-6<<" "<<res2/92.8e-6<<" "<<es2/res2<<" "<<res4/92.8e-6<<endl;
	  //fout<<scEta[j]/fabs(scEta[j])<<" "<<deltaP<<" "<<deltaE<<" "<<es1<<" "<<res1<< " "<<es1/res1<<" "<<res3<<" "<<es2<<" "<<res2<<" "<<es2/res2<<" "<<res4<<endl;
	  //H2_->Draw("colztext");
	  H2_->Draw("colz");
	  }
	  TLine *L1;
	  L1 = new TLine (s, tempminP-r, s, tempmaxP+r); L1->Draw();
	  L1 = new TLine (-s, tempminP-r, -s, tempmaxP+r); L1->Draw();
	  L1 = new TLine (s, tempminP-r, -s, tempminP-r); L1->Draw();
	  L1 = new TLine (s, tempmaxP+r, -s, tempmaxP+r); L1->Draw();
	  c1->cd(2);
	  H3_->SetTitle("");
	  H3_->GetXaxis()->SetTitle("#Delta#eta");
	  H3_->GetYaxis()->SetTitle("#Delta#phi");
	  H4_->SetTitle("");
	  H4_->GetXaxis()->SetTitle("#Delta#eta");
	  H4_->GetYaxis()->SetTitle("#Delta#phi");
	  //if (scEta[j]/fabs(scEta[j])==1&&(res3-res1)/92.8e-6>35&&(res4-res2)/92.8e-6>35){
	  if (scEta[j]/fabs(scEta[j])==1){
	  //fout<<scEta[j]/fabs(scEta[j])<<" "<<deltaP<<" "<<deltaE<<" "<<es1<<" "<<res1<< " "<<es1/res1<<" "<<es2<<" "<<res2<<" "<<es2/res2<<endl;
	  //H3_->Draw("colztext");
	  H3_->Draw("colz");
	  }
	  //if (scEta[j]/fabs(scEta[j])==-1&&(res3-res1)/92.8e-6>35&&(res4-res2)/92.8e-6>35){
	  if (scEta[j]/fabs(scEta[j])==-1){
	  //fout<<scEta[j]/fabs(scEta[j])<<" "<<deltaP<<" "<<deltaE<<" "<<es1<<" "<<res1<< " "<<es1/res1<<" "<<es2<<" "<<res2<<" "<<es2/res2<<endl;
	  //H4_->Draw("colztext");
	  H4_->Draw("colz");
	  }
	  TLine *L2;
	  L2 = new TLine (s, tempminP-r, s, tempmaxP+r); L2->Draw();
	  L2 = new TLine (-s, tempminP-r, -s, tempmaxP+r); L2->Draw();
	  L2 = new TLine (s, tempminP-r, -s, tempminP-r); L2->Draw();
	  L2 = new TLine (s, tempmaxP+r, -s, tempmaxP+r); L2->Draw();
	  //c1->Update();
	  c1->Print("event_display.pdf");
	  cin.get();
	  delete c1;
	  delete L1;
	  delete L2;
	  }*/ 
	  //if (scEta[j]/fabs(scEta[j])==1) cout<<resmaxp1/92.8e-6<<" "<<resmaxp2/92.8e-6<<" "<<res3/92.8e-6<<" "<<res1/92.8e-6<<" "<<res4/92.8e-6<<" "<<res2/92.8e-6<<endl;
	  //if (scEta[j]/fabs(scEta[j])==-1) cout<<resmaxn1/92.8e-6<<" "<<resmaxn2/92.8e-6<<" "<<res3/92.8e-6<<" "<<res1/92.8e-6<<" "<<res4/92.8e-6<<" "<<res2/92.8e-6<<endl;
	}
      } 
    }
  }
  //fout.close();
  
  TF1 *f1 = new TF1("f1", "pol1", 0, 2900);
  TF1 *f2 = new TF1("f2", "pol1", 0, 2900);
  TF1 *f3 = new TF1("f3", "pol1", 0, 2900);
  f1->SetLineColor(1);
  f1->SetLineWidth(3);
  f2->SetLineColor(2);
  f2->SetLineWidth(3);
  f3->SetLineColor(3);
  f3->SetLineWidth(3);
  if (PT_==20 && ETA_==17){
    hProf1->Fit("f1", "MRW0", "", 0, 200);
    hProf2->Fit("f2", "MRW0", "", 0, 200);
    hProf3->Fit("f3", "MRW0", "", 0, 200);
  }
  if (PT_==20 && ETA_==19){
    hProf1->Fit("f1", "MRW0", "", 0, 400);
    hProf2->Fit("f2", "MRW0", "", 0, 400);
    hProf3->Fit("f3", "MRW0", "", 0, 400);
  }
  if (PT_==20 && ETA_==21){
    hProf1->Fit("f1", "MRW0", "", 100, 570);
    hProf2->Fit("f2", "MRW0", "", 100, 570);
    hProf3->Fit("f3", "MRW0", "", 100, 570);
  }
  if (PT_==20 && ETA_==23){
    hProf1->Fit("f1", "MRW0", "", 0, 570);
    hProf2->Fit("f2", "MRW0", "", 0, 570);
    hProf3->Fit("f3", "MRW0", "", 0, 570);
  }
  if (PT_==40 && ETA_==17){ 
    hProf1->Fit("f1", "MRW0", "", 0, 630);
    hProf2->Fit("f2", "MRW0", "", 0, 630);
    hProf3->Fit("f3", "MRW0", "", 0, 630);
  }
  if (PT_==40 && ETA_==19){
    hProf1->Fit("f1", "MRW0", "", 0, 900);
    hProf2->Fit("f2", "MRW0", "", 0, 900);
    hProf3->Fit("f3", "MRW0", "", 0, 900);
  }
  if (PT_==40 && ETA_==21){
    hProf1->Fit("f1", "MRW0", "", 200, 850);
    hProf2->Fit("f2", "MRW0", "", 200, 850);
    hProf3->Fit("f3", "MRW0", "", 200, 850);
  }
  if (PT_==40 && ETA_==23){ 
    hProf1->Fit("f1", "MRW0", "", 50, 800);
    hProf2->Fit("f2", "MRW0", "", 50, 800);
    hProf3->Fit("f3", "MRW0", "", 50, 800);
  }
  if (PT_==60 && ETA_==17){ 
    hProf1->Fit("f1", "MRW0", "", 0, 950);
    hProf2->Fit("f2", "MRW0", "", 0, 950);
    hProf3->Fit("f3", "MRW0", "", 0, 950);
  }
  if (PT_==60 && ETA_==19){ 
    hProf1->Fit("f1", "MRW0", "", 50, 1350);
    hProf2->Fit("f2", "MRW0", "", 50, 1350);
    hProf3->Fit("f3", "MRW0", "", 50, 1350);
  }
  if (PT_==60 && ETA_==21){ 
    hProf1->Fit("f1", "MRW0", "", 200, 1200);
    hProf2->Fit("f2", "MRW0", "", 200, 1200);
    hProf3->Fit("f3", "MRW0", "", 200, 1200);
  }
  if (PT_==60 && ETA_==23){ 
    hProf1->Fit("f1", "MRW0", "", 200, 1250);
    hProf2->Fit("f2", "MRW0", "", 200, 1250);
    hProf3->Fit("f3", "MRW0", "", 200, 1250);
  }
  if (PT_==80 && ETA_==17){ 
    hProf1->Fit("f1", "MRW0", "", 200, 1300);
    hProf2->Fit("f2", "MRW0", "", 200, 1300);
    hProf3->Fit("f3", "MRW0", "", 200, 1300);
  }
  if (PT_==80 && ETA_==19){ 
    hProf1->Fit("f1", "MRW0", "", 200, 1750);
    hProf2->Fit("f2", "MRW0", "", 200, 1750);
    hProf3->Fit("f3", "MRW0", "", 200, 1750);
  }
  if (PT_==80 && ETA_==21){ 
    hProf1->Fit("f1", "MRW0", "", 200, 1350);
    hProf2->Fit("f2", "MRW0", "", 200, 1350);
    hProf3->Fit("f3", "MRW0", "", 200, 1350);
  }
  if (PT_==80 && ETA_==23){ 
    hProf1->Fit("f1", "MRW0", "", 200, 1400);
    hProf2->Fit("f2", "MRW0", "", 200, 1400);
    hProf3->Fit("f3", "MRW0", "", 200, 1400);
  }
  //Double_t chi2 = f1->GetChisquare();
  Double_t gamma1 = f1->GetParameter(1);
  Double_t egamma1 = f1->GetParError(1);
  Double_t gamma2 = f2->GetParameter(1);
  Double_t egamma2 = f2->GetParError(1);
  Double_t gamma3 = f3->GetParameter(1);
  Double_t egamma3 = f3->GetParError(1);
  Char_t text[200];
  sprintf(text, "Fit slope #gamma = %.02g #pm %.02g (GeV/MIP)", gamma1, egamma1);
  sprintf(text, "Fit slope #gamma = %.02g #pm %.02g (GeV/MIP)", gamma2, egamma2);
  sprintf(text, "Fit slope #gamma = %.02g #pm %.02g (GeV/MIP)", gamma3, egamma3);
  
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    for (int i=0; i<nEle; ++i) {
      
      if (fabs(eleSCEta[i]) < 1.65) continue;
      if (fabs(eleSCEta[i]) > 2.5) continue;
      es3 = 0;
      es4 = 0;
      res5 = 0;
      res6 = 0;
      res7 = 0;
      res8 = 0;
      
      for (int j=0; j<nSC; ++j) {
	Int_t Matched =0;
	for (int n=0; n<nMC; n++) {
	  Float_t dR = deltaR(scEta[j], scPhi[j], mcEta[n], mcPhi[n]);
          Float_t dPt = fabs(elePt[i] - mcPt[n])/mcPt[n];
	  if (dR<0.1 && dPt < 0.2) Matched=1;
	}
        if (Matched==0) continue; 
        if (st1x[j][0]==0 || st2x[j][0]==0) continue;
        if (fabs(eleSCEta[i] - scEta[j]) < 0.001 && fabs(eleSCPhi[i] - scPhi[j]) < 0.001) {
          Float_t tempmaxP = -9999;
          Float_t tempmaxE = -9999;
          Float_t tempminP =  9999;
          Float_t tempminE =  9999;
          if (st1z[j][0]==-1) st1z[j][0]=0;
          if (st2z[j][0]==-1) st2z[j][0]=0;

	  //ratio_old->Fill(scEn[j]/(PT_*cosh(ETA_/10.)));
          //f = scEn[j]/(scRawEn[j]+scESEn[j]);
          
	  for (int k=0; k<nBC[j]; ++k) {
	    tempdPhi = bcPhi[j][0] - bcPhi[j][k];
	    tempdEta = bcEta[j][0] - bcEta[j][k];
	    if (tempdPhi > TMath::Pi()) tempdPhi = tempdPhi - 2.*TMath::Pi();
	    if (tempdPhi < -TMath::Pi()) tempdPhi += 2.*TMath::Pi();
	    if (tempmaxP<tempdPhi) tempmaxP=tempdPhi;
	    if (tempmaxE<tempdEta) tempmaxE=tempdEta;
	    if (tempminP>tempdPhi) tempminP=tempdPhi;
	    if (tempminE>tempdEta) tempminE=tempdEta;
            for (int m=0; m<4; ++m) {
	      if (bcESEnX[j][k][m] > 0) es3 += bcESEnX[j][k][m];
	      if (bcESEnY[j][k][m] > 0) es4 += bcESEnY[j][k][m];
            }
	  }
	  
	  for (int l=0; l<nRH; ++l) {
            if (rhFlag[l]==1||rhFlag[l]==14||(rhFlag[l]<=10&&rhFlag[l]>=5)) continue;
            if (rhZ[l]*scEta[j]>0 && rhP[l]==1){
              if (rhZ[l]==1){
		stEta=eta[st1z[j][0]][st1p[j][0]-1][st1x[j][0]-1][st1y[j][0]-1][st1s[j][0]-1];
		rhEta=eta[rhZ[l]][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dEta = stEta-rhEta;
		stPhi=phi[st1z[j][0]][st1p[j][0]-1][st1x[j][0]-1][st1y[j][0]-1][st1s[j][0]-1];
		rhPhi=phi[rhZ[l]][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dPhi = stPhi-rhPhi;
		if (dPhi > TMath::Pi()) dPhi = dPhi - 2.*TMath::Pi();
		if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
              }
              if (rhZ[l]==-1){
		//rhZ[l]=0;
		stEta=eta[st1z[j][0]][st1p[j][0]-1][st1x[j][0]-1][st1y[j][0]-1][st1s[j][0]-1];
		rhEta=eta[0][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dEta = stEta-rhEta;
		stPhi=phi[st1z[j][0]][st1p[j][0]-1][st1x[j][0]-1][st1y[j][0]-1][st1s[j][0]-1];
		rhPhi=phi[0][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dPhi = stPhi-rhPhi;
		if (dPhi > TMath::Pi()) dPhi = dPhi - 2.*TMath::Pi();
		if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
              }
              if ((0<=dPhi && dPhi<=(tempmaxP+r_) && fabs(dEta)<=s_) 
		  || ((tempminP-r_)<=dPhi && dPhi<0 && fabs(dEta)<=s_)){
		if (rhE[l]/(92.8*0.000001)>MIP)  res5 += rhE[l];
              }
	    }
            if (rhZ[l]*scEta[j]>0 && rhP[l]==2){
              if (rhZ[l]==1){
		stEta=eta[st2z[j][0]][st2p[j][0]-1][st2x[j][0]-1][st2y[j][0]-1][st2s[j][0]-1];
		rhEta=eta[rhZ[l]][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dEta = stEta-rhEta;
		stPhi=phi[st2z[j][0]][st2p[j][0]-1][st2x[j][0]-1][st2y[j][0]-1][st2s[j][0]-1];
		rhPhi=phi[rhZ[l]][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dPhi = stPhi-rhPhi;
		if (dPhi > TMath::Pi()) dPhi = dPhi - 2.*TMath::Pi();
		if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
              }
              if (rhZ[l]==-1){
		//rhZ[l]=0;
		stEta=eta[st2z[j][0]][st2p[j][0]-1][st2x[j][0]-1][st2y[j][0]-1][st2s[j][0]-1];
		rhEta=eta[0][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dEta = stEta-rhEta;
		stPhi=phi[st2z[j][0]][st2p[j][0]-1][st2x[j][0]-1][st2y[j][0]-1][st2s[j][0]-1];
		rhPhi=phi[0][rhP[l]-1][rhX[l]-1][rhY[l]-1][rhS[l]-1];
		dPhi = stPhi-rhPhi;
		if (dPhi > TMath::Pi()) dPhi = dPhi - 2.*TMath::Pi();
		if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
              }
              if ((0<=dPhi && dPhi<=(tempmaxP+r_) && fabs(dEta)<=s_) 
		  || ((tempminP-r_)<=dPhi && dPhi<0 && fabs(dEta)<=s_)){
		if (rhE[l]/(92.8*0.000001)>MIP)  res6 += rhE[l];
              }  
	    }
	  }
	  
          for (int a=0; a<nRH; ++a) {
            if (rhFlag[a]==1||rhFlag[a]==14||(rhFlag[a]<=10&&rhFlag[a]>=5)) continue;
	    if (rhZ[a]*scEta[j]>0 && rhP[a]==1){
	      if (rhE[a]/(92.8*0.000001)>MIP)  res7 += rhE[a];
	    }
	    if (rhZ[a]*scEta[j]>0 && rhP[a]==2){
	      if (rhE[a]/(92.8*0.000001)>MIP)  res8 += rhE[a];
	    }
	  }
	  
	  if (eleBrem[i]<0) continue;
	  bin1 = eleBrem[i]/0.25;
	  ratio_old->Fill(scEn[j]/(PT_*cosh(ETA_/10.)));
          f = scEn[j]/(scRawEn[j]+scESEn[j]);
	  dEoldnew = scRawEn[j] + 0.020*((es3/(92.8e-6))+(0.7*(es4/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h20_->Fill(dEoldnew);
	  hfb19_[bin1]->Fill(dEoldnew);
	  dEnew = scRawEn[j] + 0.020*((res5/(92.8e-6))+(0.7*(res6/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h2_->Fill(dEnew);
	  hfb12_[bin1]->Fill(dEnew);
	  dEalln = scRawEn[j] + 0.020*((res7/(92.8e-6))+(0.7*(res8/(92.8e-6)))) - PT_*cosh(ETA_/10.);
	  h21_->Fill(dEalln);
	  hfb20_[bin1]->Fill(dEalln);
	  
	  if (eleEta[i]>0){ 
	    h25_->Fill(res5/res7);
	    h26_->Fill(res6/res8);
	  }
	  if (eleEta[i]<0){
	    h27_->Fill(res5/res7);
	    h28_->Fill(res6/res8);
	  }
	  
	  ratio_new->Fill(f*(0.020*((res5/92.8e-6)+(0.7*(res6/92.8e-6)))+scRawEn[j])/(PT_*cosh(ETA_/10.)));
	}
      }
    }
  }
  fout_->Write();
  fout_->Close();
}

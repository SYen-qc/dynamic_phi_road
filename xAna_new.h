//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 19 09:58:30 2012 by ROOT version 5.27/06b
// from TTree tree/tree
// found on file: /data2/cmkuo/particleGun/job_electron_gun_pt60_eta23.root
//////////////////////////////////////////////////////////

#ifndef xAna_new_h
#define xAna_new_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>

class xAna_new {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           lumis;
   Int_t           nVtx;
   Int_t           IsVtxGood;
   Int_t           nRH;
   Int_t           rhZ[100000];   //[nRH]
   Int_t           rhP[100000];   //[nRH]
   Int_t           rhX[100000];   //[nRH]
   Int_t           rhY[100000];   //[nRH]
   Int_t           rhS[100000];   //[nRH]
   Float_t         rhT[100000];   //[nRH]
   Float_t         rhE[100000];   //[nRH]
   Int_t           rhFlag[100000];   //[nRH]   
   Int_t           nMC;
   Int_t           mcPID[2];   //[nMC]
   Int_t           mcMomPID[2];   //[nMC]
   Float_t         mcPt[2];   //[nMC]
   Float_t         mcEta[2];   //[nMC]
   Float_t         mcPhi[2];   //[nMC]
   Float_t         mcE[2];   //[nMC]
   Int_t           nSC;
   Int_t           nBC[100000];   //[nSC]
   Float_t         scEn[100000];   //[nSC]
   Float_t         scEta[100000];   //[nSC]
   Float_t         scPhi[100000];   //[nSC]
   Float_t         scRawEn[100000];   //[nSC]
   Float_t         scESEn[100000];   //[nSC]
   Float_t         bcEta[100000][10];   //[nSC]
   Float_t         bcPhi[100000][10];   //[nSC]
   Float_t         bcX[100000][10];   //[nSC]
   Float_t         bcY[100000][10];   //[nSC]
   Float_t         bcESEnX[100000][10][4];   //[nSC]
   Float_t         bcESEnY[100000][10][4];   //[nSC]
   Float_t         bcEn[100000][10];   //[nSC]
   Int_t           st1z[100000][10];   //[nSC]
   Int_t           st1p[100000][10];   //[nSC]
   Int_t           st1x[100000][10];   //[nSC]
   Int_t           st1y[100000][10];   //[nSC]
   Int_t           st1s[100000][10];   //[nSC]
   Int_t           st2z[100000][10];   //[nSC]
   Int_t           st2p[100000][10];   //[nSC]
   Int_t           st2x[100000][10];   //[nSC]
   Int_t           st2y[100000][10];   //[nSC]
   Int_t           st2s[100000][10];   //[nSC]
   Int_t           nBCRH[100000][10];   //[nSC]
   Int_t           bcRHz[100000][10][30];   //[nSC]
   Int_t           bcRHx[100000][10][30];   //[nSC]
   Int_t           bcRHy[100000][10][30];   //[nSC]
   Int_t           nEle;
   Int_t           eleClass[3];   //[nEle]
   Float_t         elePt[3];   //[nEle]
   Float_t         eleEn[3];   //[nEle]
   Float_t         eleESEn[3];   //[nEle]
   Float_t         eleEta[3];   //[nEle]
   Float_t         elePhi[3];   //[nEle]
   Float_t         eleSCRawEn[3];   //[nEle]
   Float_t         eleSCEn[3];   //[nEle]
   Float_t         eleSCEta[3];   //[nEle]
   Float_t         eleSCPhi[3];   //[nEle]
   Float_t         eleBrem[10];   //[nEle]
   Float_t         eleHoverE[3];   //[nEle]
   Float_t         eleEoverP[3];   //[nEle]
   Float_t         eledEtaAtVtx[3];   //[nEle]
   Float_t         eledPhiAtVtx[3];   //[nEle]
   Float_t         eleSigmaIEtaIEta[3];   //[nEle]
   Float_t         eleIsoTrk[3];   //[nEle]
   Float_t         eleIsoEcal[3];   //[nEle]
   Float_t         eleIsoHcal[3];   //[nEle]
   Float_t         elePFSCRawEn[3];   //[nEle]
   Float_t         elePFESEn[3];   //[nEle]
   Int_t           npfEle;
   Float_t         pfeleESEn[3];   //[npfEle]
   Float_t         pfeleEta[3];   //[npfEle]
   Float_t         pfelePhi[3];   //[npfEle]
   Float_t         pfeleSCRawEn[3];   //[npfEle]
   Float_t         pfeleSCEn[3];   //[npfEle]
   Int_t           nPFEle;
   Float_t         PFEle_energyES1P;
   Float_t         PFEle_energyES2P;
   Float_t         PFEle_energyES1M;
   Float_t         PFEle_energyES2M;
   Float_t         PFEle_ESenergy_fromcode[3];   //[npfEle]
   Int_t           nPho;
   Float_t         phoEn[1];   //[nPho]
   Float_t         phoESEn[1];   //[nPho]
   Float_t         phoEta[1];   //[nPho]
   Float_t         phoPhi[1];   //[nPho]
   Float_t         phoEnSC[1];   //[nPho]
   Float_t         phoEtaSC[1];   //[nPho]
   Float_t         phoPhiSC[1];   //[nPho]
   Int_t           phoPixelSeed[1];   //[nPho]
   Float_t         phoR9[1];   //[nPho]
   Float_t         phoHoverE[1];   //[nPho]
   Float_t         phoSigmaIEtaIEta[1];   //[nPho]
   Float_t         phoIsoTrkSolid[1];   //[nPho]
   Float_t         phoIsoTrkHollow[1];   //[nPho]
   Float_t         phoIsoEcal[1];   //[nPho]
   Float_t         phoIsoHcal[1];   //[nPho]
   Int_t           nSHES;
   Int_t           shES_z[100000];   //[nSHES]
   Int_t           shES_p[100000];   //[nSHES]
   Int_t           shES_x[100000];   //[nSHES]
   Int_t           shES_y[100000];   //[nSHES]
   Int_t           shES_s[100000];   //[nSHES]
   Float_t         shES_e[100000];   //[nSHES]
   Float_t         shES_efP;
   Float_t         shES_erP;
   Float_t         shES_efM;
   Float_t         shES_erM;
   Float_t         shES_Mside;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_IsVtxGood;   //!
   TBranch        *b_nRH;   //!
   TBranch        *b_rhZ;   //!
   TBranch        *b_rhP;   //!
   TBranch        *b_rhX;   //!
   TBranch        *b_rhY;   //!
   TBranch        *b_rhS;   //!
   TBranch        *b_rhT;   //!
   TBranch        *b_rhE;   //!
   TBranch        *b_rhFlag;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_nSC;   //!
   TBranch        *b_nBC;   //!
   TBranch        *b_scEn;   //!
   TBranch        *b_scEta;   //!
   TBranch        *b_scPhi;   //!
   TBranch        *b_scRawEn;   //!
   TBranch        *b_scESEn;   //!
   TBranch        *b_bcEta;   //!
   TBranch        *b_bcPhi;   //!
   TBranch        *b_bcX;   //!
   TBranch        *b_bcY;   //!
   TBranch        *b_bcESEnX;   //!
   TBranch        *b_bcESEnY;   //!
   TBranch        *b_bcEn;   //!
   TBranch        *b_st1z;   //!
   TBranch        *b_st1p;   //!
   TBranch        *b_st1x;   //!
   TBranch        *b_st1y;   //!
   TBranch        *b_st1s;   //!
   TBranch        *b_st2z;   //!
   TBranch        *b_st2p;   //!
   TBranch        *b_st2x;   //!
   TBranch        *b_st2y;   //!
   TBranch        *b_st2s;   //!
   TBranch        *b_nBCRH;   //!
   TBranch        *b_bcRHz;   //!
   TBranch        *b_bcRHx;   //!
   TBranch        *b_bcRHy;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleClass;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleESEn;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEta;   //!
   TBranch        *b_eleIsoTrk;   //!
   TBranch        *b_eleIsoEcal;   //!
   TBranch        *b_eleIsoHcal;   //!
   TBranch        *b_elePFSCRawEn;   //!
   TBranch        *b_elePFESEn;   //!
   TBranch        *b_npfEle;   //!
   TBranch        *b_pfeleESEn;   //!
   TBranch        *b_pfeleEta;   //!
   TBranch        *b_pfelePhi;   //!
   TBranch        *b_pfeleSCRawEn;   //!
   TBranch        *b_pfeleSCEn;   //!
   TBranch        *b_nPFEle;   //!
   TBranch        *b_PFEle_energyES1P;   //!
   TBranch        *b_PFEle_energyES2P;   //!
   TBranch        *b_PFEle_energyES1M;   //!
   TBranch        *b_PFEle_energyES2M;   //!
   TBranch        *b_PFEle_ESenergy_fromcode;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoEn;   //!
   TBranch        *b_phoESEn;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoEnSC;   //!
   TBranch        *b_phoEtaSC;   //!
   TBranch        *b_phoPhiSC;   //!
   TBranch        *b_phoPixelSeed;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoIsoTrkSolid;   //!
   TBranch        *b_phoIsoTrkHollow;   //!
   TBranch        *b_phoIsoEcal;   //!
   TBranch        *b_phoIsoHcal;   //!
   TBranch        *b_nSHES;   //!
   TBranch        *b_shES_z;   //!
   TBranch        *b_shES_p;   //!
   TBranch        *b_shES_x;   //!
   TBranch        *b_shES_y;   //!
   TBranch        *b_shES_s;   //!
   TBranch        *b_shES_e;   //!
   TBranch        *b_shES_efP;   //!
   TBranch        *b_shES_erP;   //!
   TBranch        *b_shES_efM;   //!
   TBranch        *b_shES_erM;   //!
   TBranch        *b_shES_Mside;   //!

   xAna_new(Int_t PT=0,Int_t ETA=0,Double_t r=0,Double_t s=0);
   virtual ~xAna_new();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Double_t deltaPhi(Double_t scPhi, Double_t mcPhi);
   virtual Double_t deltaR(Double_t scEta, Double_t scPhi, Double_t mcEta, Double_t mcPhi);
   virtual void     InitHists();
   Int_t PT_;
   Int_t ETA_;
   Double_t r_;
   Double_t s_; 
   TFile *fout_;
   TH2F *H1_;
   TH2F *H2_;
   TH2F *H3_;
   TH2F *H4_;
   TH2F *potato_old_;
   TH2F *potato_new_;
   TH2F *potato_all_; 
   TH1F *NRHP1Zp_;
   TH1F *NRHP1Zn_;
   TH1F *NRHP2Zp_;
   TH1F *NRHP2Zn_;
   TH1F *Massgen_;
   TH1F *Massold_;
   TH1F *Massnew_;
   TH1F *HH1_;
   TH1F *HH2_;
   TH1F *HH3_;
   TH1F *HH4_;
   TH1F *HH5_;
   TH1F *HH6_;
   TH1F *HH7_;
   TH1F *HH8_;
   TH1F *h1_;
   TH1F *h2_;
   TH1F *h3_;
   TH1F *h4_;
   TH1F *h5_;
   TH1F *h6_;
   TH1F *h7_;
   TH1F *h8_;
   TH1F *h9_;
   TH1F *h10_;
   TH1F *h11_;
   TH1F *h12_;
   TH1F *h13_;
   TH1F *h14_;
   TH1F *h15_;
   TH1F *h16_;
   TH1F *h17_;
   TH1F *h18_;
   TH1F *h19_;
   TH1F *h20_;
   TH1F *h21_;
   TH1F *h22_;
   TH1F *h23_;
   TH1F *h24_;
   TH1F *h25_;
   TH1F *h26_;
   TH1F *h27_;
   TH1F *h28_;
   TH1F *Eh1_;
   TH1F *Eh2_;
   TH1F *Eh3_;
   TH1F *Eh4_;
   TH1F *F1_;
   TH1F *F2_;
   TH1F *F3_;
   TH2F *F4_;
   TH1F *hfb1_[10];
   TH1F *hfb2_[10];
   TH1F *hfb3_[10];
   TH1F *hfb4_[10];
   TH1F *hfb5_[10];
   TH1F *hfb6_[10];
   TH1F *hfb7_[10];
   TH1F *hfb8_[10];
   TH1F *hfb9_[10];
   TH1F *hfb10_[10];
   TH1F *hfb11_[10];
   TH1F *hfb12_[10];
   TH1F *hfb13_[10];
   TH1F *hfb14_[10];
   TH1F *hfb15_[10];
   TH1F *hfb16_[10];
   TH1F *hfb17_[10];
   TH1F *hfb18_[10];
   TH1F *hfb19_[10];
   TH1F *hfb20_[10];
   
   TH1F *ratio_old;
   TH1F *ratio_new;

};

#endif

#ifdef xAna_new_cxx
xAna_new::xAna_new(Int_t PT,Int_t ETA, Double_t r, Double_t s)
{
   PT_  = PT;
   ETA_ = ETA;
   r_ = 0.12;
   s_ = 0.15;
   Char_t fname[200];
   //sprintf(fname, "/data2/cmkuo/particleGun/reducedRH/job_electron_gun_pt%d_eta%d.root",PT_,ETA_);
   sprintf(fname, "/data2/cmkuo/particleGun/job_electron_gun_pt%d_eta%d.root",PT_,ETA_);
   //sprintf(fname, "/data2/cmkuo/particleGun/pu/job_electron_gun_pt%d_eta%d.root",PT_,ETA_);
   //sprintf(fname, "/data2/cmkuo/Zee/collisionLite_mc.root");
   TChain *chain = new TChain("demo/tree");
   chain->Add(fname);
   Init(chain);
}
 
xAna_new::~xAna_new()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t xAna_new::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t xAna_new::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void xAna_new::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("IsVtxGood", &IsVtxGood, &b_IsVtxGood);
   fChain->SetBranchAddress("nRH", &nRH, &b_nRH);
   fChain->SetBranchAddress("rhZ", rhZ, &b_rhZ);
   fChain->SetBranchAddress("rhP", rhP, &b_rhP);
   fChain->SetBranchAddress("rhX", rhX, &b_rhX);
   fChain->SetBranchAddress("rhY", rhY, &b_rhY);
   fChain->SetBranchAddress("rhS", rhS, &b_rhS);
   fChain->SetBranchAddress("rhT", rhT, &b_rhT);
   fChain->SetBranchAddress("rhE", rhE, &b_rhE);
   fChain->SetBranchAddress("rhFlag", rhFlag, &b_rhFlag);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcMomPID", mcMomPID, &b_mcMomPID);
   fChain->SetBranchAddress("mcPt", mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcEta", mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", mcE, &b_mcE);
   fChain->SetBranchAddress("nSC", &nSC, &b_nSC);
   fChain->SetBranchAddress("nBC", nBC, &b_nBC);
   fChain->SetBranchAddress("scEn", scEn, &b_scEn);
   fChain->SetBranchAddress("scEta", scEta, &b_scEta);
   fChain->SetBranchAddress("scPhi", scPhi, &b_scPhi);
   fChain->SetBranchAddress("scRawEn", scRawEn, &b_scRawEn);
   fChain->SetBranchAddress("scESEn", scESEn, &b_scESEn);
   fChain->SetBranchAddress("bcEta", bcEta, &b_bcEta);
   fChain->SetBranchAddress("bcPhi", bcPhi, &b_bcPhi);
   fChain->SetBranchAddress("bcX", bcX, &b_bcX);
   fChain->SetBranchAddress("bcY", bcY, &b_bcY);
   fChain->SetBranchAddress("bcESEnX", bcESEnX, &b_bcESEnX);
   fChain->SetBranchAddress("bcESEnY", bcESEnY, &b_bcESEnY);
   fChain->SetBranchAddress("bcEn", bcEn, &b_bcEn);
   fChain->SetBranchAddress("st1z", st1z, &b_st1z);
   fChain->SetBranchAddress("st1p", st1p, &b_st1p);
   fChain->SetBranchAddress("st1x", st1x, &b_st1x);
   fChain->SetBranchAddress("st1y", st1y, &b_st1y);
   fChain->SetBranchAddress("st1s", st1s, &b_st1s);
   fChain->SetBranchAddress("st2z", st2z, &b_st2z);
   fChain->SetBranchAddress("st2p", st2p, &b_st2p);
   fChain->SetBranchAddress("st2x", st2x, &b_st2x);
   fChain->SetBranchAddress("st2y", st2y, &b_st2y);
   fChain->SetBranchAddress("st2s", st2s, &b_st2s);
   fChain->SetBranchAddress("nBCRH", nBCRH, &b_nBCRH);
   fChain->SetBranchAddress("bcRHz", bcRHz, &b_bcRHz);
   fChain->SetBranchAddress("bcRHx", bcRHx, &b_bcRHx);
   fChain->SetBranchAddress("bcRHy", bcRHy, &b_bcRHy);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleClass", eleClass, &b_eleClass);
   fChain->SetBranchAddress("elePt", elePt, &b_elePt);
   fChain->SetBranchAddress("eleEn", eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleESEn", eleESEn, &b_eleESEn);
   fChain->SetBranchAddress("eleEta", eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleSCRawEn", eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEn", eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("eleSCEta", eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleBrem", eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eleHoverE", eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleEoverP", eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("eledEtaAtVtx", eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eleSigmaIEtaIEta", eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   fChain->SetBranchAddress("eleIsoTrk", eleIsoTrk, &b_eleIsoTrk);
   fChain->SetBranchAddress("eleIsoEcal", eleIsoEcal, &b_eleIsoEcal);
   fChain->SetBranchAddress("eleIsoHcal", eleIsoHcal, &b_eleIsoHcal);
   fChain->SetBranchAddress("elePFSCRawEn", elePFSCRawEn, &b_elePFSCRawEn);
   fChain->SetBranchAddress("elePFESEn", elePFESEn, &b_elePFESEn);
   fChain->SetBranchAddress("npfEle", &npfEle, &b_npfEle);
   fChain->SetBranchAddress("pfeleESEn", pfeleESEn, &b_pfeleESEn);
   fChain->SetBranchAddress("pfeleEta", pfeleEta, &b_pfeleEta);
   fChain->SetBranchAddress("pfelePhi", pfelePhi, &b_pfelePhi);
   fChain->SetBranchAddress("pfeleSCRawEn", pfeleSCRawEn, &b_pfeleSCRawEn);
   fChain->SetBranchAddress("pfeleSCEn", pfeleSCEn, &b_pfeleSCEn);
   fChain->SetBranchAddress("nPFEle", &nPFEle, &b_nPFEle);
   fChain->SetBranchAddress("PFEle_energyES1P", &PFEle_energyES1P, &b_PFEle_energyES1P);
   fChain->SetBranchAddress("PFEle_energyES2P", &PFEle_energyES2P, &b_PFEle_energyES2P);
   fChain->SetBranchAddress("PFEle_energyES1M", &PFEle_energyES1M, &b_PFEle_energyES1M);
   fChain->SetBranchAddress("PFEle_energyES2M", &PFEle_energyES2M, &b_PFEle_energyES2M);
   fChain->SetBranchAddress("PFEle_ESenergy_fromcode", PFEle_ESenergy_fromcode, &b_PFEle_ESenergy_fromcode);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoEn", &phoEn, &b_phoEn);
   fChain->SetBranchAddress("phoESEn", &phoESEn, &b_phoESEn);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoEnSC", &phoEnSC, &b_phoEnSC);
   fChain->SetBranchAddress("phoEtaSC", &phoEtaSC, &b_phoEtaSC);
   fChain->SetBranchAddress("phoPhiSC", &phoPhiSC, &b_phoPhiSC);
   fChain->SetBranchAddress("phoPixelSeed", &phoPixelSeed, &b_phoPixelSeed);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoIsoTrkSolid", &phoIsoTrkSolid, &b_phoIsoTrkSolid);
   fChain->SetBranchAddress("phoIsoTrkHollow", &phoIsoTrkHollow, &b_phoIsoTrkHollow);
   fChain->SetBranchAddress("phoIsoEcal", &phoIsoEcal, &b_phoIsoEcal);
   fChain->SetBranchAddress("phoIsoHcal", &phoIsoHcal, &b_phoIsoHcal);
   fChain->SetBranchAddress("nSHES", &nSHES, &b_nSHES);
   fChain->SetBranchAddress("shES_z", shES_z, &b_shES_z);
   fChain->SetBranchAddress("shES_p", shES_p, &b_shES_p);
   fChain->SetBranchAddress("shES_x", shES_x, &b_shES_x);
   fChain->SetBranchAddress("shES_y", shES_y, &b_shES_y);
   fChain->SetBranchAddress("shES_s", shES_s, &b_shES_s);
   fChain->SetBranchAddress("shES_e", shES_e, &b_shES_e);
   fChain->SetBranchAddress("shES_efP", &shES_efP, &b_shES_efP);
   fChain->SetBranchAddress("shES_erP", &shES_erP, &b_shES_erP);
   fChain->SetBranchAddress("shES_efM", &shES_efM, &b_shES_efM);
   fChain->SetBranchAddress("shES_erM", &shES_erM, &b_shES_erM);
   fChain->SetBranchAddress("shES_Mside", &shES_Mside, &b_shES_Mside);
   Notify();
}

Bool_t xAna_new::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void xAna_new::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t xAna_new::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Double_t xAna_new::deltaPhi(Double_t Phi1, Double_t Phi2)
{
   Double_t delPhi = Phi1 - Phi2;
   if (delPhi > TMath::Pi()) delPhi -= 2.*TMath::Pi();
   if (delPhi < -TMath::Pi()) delPhi += 2.*TMath::Pi();

   return delPhi;
}

Double_t xAna_new::deltaR(Double_t Eta1, Double_t Phi1, Double_t Eta2, Double_t Phi2)
{
   Double_t delEta, delPhi;
   delEta = Eta1 - Eta2;
   delPhi = deltaPhi(Phi1,Phi2);

   return sqrt(delEta*delEta+delPhi*delPhi);
}

void xAna_new::InitHists()
{
  ratio_old = new TH1F("ratio_old","ratio_old",200,0,2);
  ratio_new = new TH1F("ratio_new","ratio_new",200,0,2);
  
  NRHP1Zp_ = new TH1F("NRHP1Zp","NRHP1Zp",1000,0,1000);
  NRHP1Zn_ = new TH1F("NRHP1Zn","NRHP1Zn",1000,0,1000);
  NRHP2Zp_ = new TH1F("NRHP2Zp","NRHP2Zp",1000,0,1000);
  NRHP2Zn_ = new TH1F("NRHP2Zn","NRHP2Zn",1000,0,1000);
  Massgen_ = new TH1F("Massgen","Massgen",1000,-200,200);
  Massold_ = new TH1F("Massold","Massold",1000,-200,200);
  Massnew_ = new TH1F("Massnew","Massnew",1000,-200,200);
  h1_ = new TH1F("h1","dEold",1000,-200,200);
  //h2_ = new TH1F("h2","dEnew",4000,-5,5);
  h2_ = new TH1F("h2","dEnew",1000,-200,200);
  h3_ = new TH1F("h3","dE024",1000,-200,200);
  h4_ = new TH1F("h4","dEES",1000,-200,200);
  //h20_ = new TH1F("h20","dEoldnew",4000,-5,5);
  h20_ = new TH1F("h20","dEoldnew",1000,-200,200);
  //h21_ = new TH1F("h21","dEalln",4000,-5,5);
  h21_ = new TH1F("h21","dEalln",1000,-200,200);
  h5_ = new TH1F("h5","dE020",800,-100,200);
  h6_ = new TH1F("h6","dE021",800,-100,200);
  h7_ = new TH1F("h7","dE022",800,-100,200);
  h8_ = new TH1F("h8","dE023",800,-100,200);
  h9_ = new TH1F("h9","dE025",800,-100,200);
  h10_ = new TH1F("h10","dE026",800,-100,200);
  h11_ = new TH1F("h11","dE027",800,-100,200);
  h12_ = new TH1F("h12","dE028",800,-100,200);
  h13_ = new TH1F("h13","dE029",800,-100,200);
  h14_ = new TH1F("h14","dE010",800,-100,200);
  h15_ = new TH1F("h15","dE015",800,-100,200);
  h16_ = new TH1F("h16","dE035",800,-100,200);
  h17_ = new TH1F("h17","dE040",800,-100,200);
  h18_ = new TH1F("h18","dEE",800,-100,200);
  h19_ = new TH1F("h19","dE050",800,-100,200);
  h22_ = new TH1F("h22","saber",1000,0,1);
  h23_ = new TH1F("h23","archer",1000,0,1);
  h24_ = new TH1F("h24","lancer",1000,0,1);
  h25_ = new TH1F("h25","Ratio1",1000,0,50);
  h26_ = new TH1F("h26","Ratio2",1000,0,50);
  h27_ = new TH1F("h27","Ratio3",1000,0,50);
  h28_ = new TH1F("h28","Ratio4",1000,0,50);
  Eh1_ = new TH1F("Eh1","ESEn+ZP1",1000,0,2);
  Eh2_ = new TH1F("Eh2","ESEn-ZP1",1000,0,2);
  Eh3_ = new TH1F("Eh3","ESEn+ZP2",1000,0,2);
  Eh4_ = new TH1F("Eh4","ESEn-ZP2",1000,0,2);
  F1_ = new TH1F("F1","F1",1000,-5,5);
  F2_ = new TH1F("F2","F2",1000,-5,5);
  F3_ = new TH1F("F3","F3",1000,-5,5);
  F4_ = new TH2F("F4","F4",250,0,2500,1000,-1,10);
  //H1_ = new TH2F("H1","dPhivsdEta1",100,-0.15,0.15,100,-0.15,0.15);
  //H2_ = new TH2F("H2","dPhivsdEta2",100,-0.15,0.15,100,-0.15,0.15);
  //H3_ = new TH2F("H3","dPhivsdEta3",100,-0.15,0.15,100,-0.15,0.15);
  //H4_ = new TH2F("H4","dPhivsdEta4",100,-0.15,0.15,100,-0.15,0.15);
  H1_ = new TH2F("H1","dPhivsdEta1",100,-1,1,100,-3.2,3.2);
  H2_ = new TH2F("H2","dPhivsdEta2",100,-1,1,100,-3.2,3.2);
  H3_ = new TH2F("H3","dPhivsdEta3",100,-1,1,100,-3.2,3.2);
  H4_ = new TH2F("H4","dPhivsdEta4",100,-1,1,100,-3.2,3.2);
  potato_old_ = new TH2F("old","old",1000,0,1000,1000,0,1000);
  potato_new_ = new TH2F("new","new",1000,0,1000,1000,0,1000);
  potato_all_ = new TH2F("all","all",1000,0,1000,1000,0,1000);
  HH1_ = new TH1F("dPhi1","dPhi1",100,0.15,0.8);
  HH2_ = new TH1F("dPhi2","dPhi2",100,0.15,0.8);
  HH3_ = new TH1F("dPhi3","dPhi3",100,0.15,0.8);
  HH4_ = new TH1F("dPhi4","dPhi4",100,0.15,0.8);
  HH5_ = new TH1F("dPhi5","dPhi5",100,0.15,0.8);
  HH6_ = new TH1F("dPhi6","dPhi6",100,0.15,0.8);
  HH7_ = new TH1F("dPhi7","dPhi7",100,0.15,0.8);
  HH8_ = new TH1F("dPhi8","dPhi8",100,0.15,0.8);
  Char_t hname[200];
  /*for (int i=0; i<10; ++i){
    sprintf(hname,"hfb%d",i+1);
    hfb1_[i] = new TH1F(hname,hname,12,0,1.2);
    sprintf(hname,"hfb%d",i+11);
    hfb2_[i] = new TH1F(hname,hname,12,0,1.2);
    }*/
  for (int i=0; i<4; i++){
    sprintf(hname,"hfb%d",i+1);
    hfb1_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+11);
    hfb2_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+21);
    hfb3_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+31);
    hfb4_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+41);
    hfb5_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+51);
    hfb6_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+61);
    hfb7_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+71);
    hfb8_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+81);
    hfb9_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+91);
    hfb10_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+101);
    hfb11_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+111);
    hfb12_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+121);
    hfb13_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+131);
    hfb14_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+141);
    hfb15_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+151);
    hfb16_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+161);
    hfb17_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+171);
    hfb18_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+181);
    hfb19_[i] = new TH1F(hname,hname,800,-100,200);
    sprintf(hname,"hfb%d",i+191);
    hfb20_[i] = new TH1F(hname,hname,800,-100,200);
  }
}  
#endif // #ifdef xAna_new_cxx

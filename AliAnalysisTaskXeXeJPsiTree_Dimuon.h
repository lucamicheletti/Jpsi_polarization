#ifndef AliAnalysisTaskXeXeJPsiTree_Dimuon_H
#define AliAnalysisTaskXeXeJPsiTree_Dimuon_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include "TTreeStream.h"//why?

class TObjArray;
class AliVParticle;
class AliAODEvent;
class TLorentzVector;
//class TList;

class AliAnalysisTaskXeXeJPsiTree_Dimuon : public AliAnalysisTaskSE {
  public:

  AliAnalysisTaskXeXeJPsiTree_Dimuon();
  AliAnalysisTaskXeXeJPsiTree_Dimuon(const char *name);
  virtual ~AliAnalysisTaskXeXeJPsiTree_Dimuon();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  void Terminate(Option_t *);
  //virtual void NotifyRun(); 
  
  void SetBeamEnergy(Double_t en) {fBeamEnergy=en;}
  void SetAnalysisType(const char* type) {fkAnalysisType=type;}
  void SetPeriod(TString period) {fPeriod=period;}
    
 private:
  AliAnalysisTaskXeXeJPsiTree_Dimuon(const AliAnalysisTaskXeXeJPsiTree_Dimuon&);
  AliAnalysisTaskXeXeJPsiTree_Dimuon& operator=(const AliAnalysisTaskXeXeJPsiTree_Dimuon&);
     
 //protected:
     
  Double_t fBeamEnergy;   // Energy of the beam (required for the CS angle)    
  const char* fkAnalysisType; //ESD or AOD based analysis
  TString fPeriod; //period
  Int_t fCountTotEv;   // counter
  Int_t fCountTrigger;   // counter
  Int_t fCountCINT7;   // counter
  Int_t fCountCMUL7; //counter
  Int_t fCountCMLL7; //counter
  Int_t fCountCMSL7; //counter
  Int_t fCountCMSH7; //counter
  
  Double_t      fNevt		;      // event counter
  
  TTree         *fOutputTree	;      //! tree output
  
  char   	fTrigClass[500];	// fired trigger classes
  UInt_t        finpmask;               // trigger input mask

  Int_t		fNMuons;		// muon tracks in the event
  Int_t		fNTracklets;		// spd tracklets
  Int_t		fNContributors;		// n contributors
  Double_t	fVertex[3];		// x,y,z vertex
  Float_t       fPercentV0M;            // percentile V0
  Float_t       fPercentCL0;            // percentile n. cluster layer 0
  Float_t       fPercentCL1;             // percentile n. cluster layer 1
  Float_t       fPercentV0A;            // percentile V0A
  Float_t       fPercentV0C;            // percentile V0C
  Float_t       fPercentZNA;             // percentile ZNA
  Float_t       fPercentZNC;             // percentile ZNC
  Double_t	fPt[300];		 // single mu pT
  Double_t	fE[300];			// single mu E
  Double_t	fPx[300];		// single mu px
  Double_t	fPy[300];		// single mu py
  Double_t	fPz[300];		// single mu pz
  Double_t	fY[300];		// single mu y
  Double_t	fEta[300];		// single mu eta
  Int_t		fMatchTrig[300];		// single mu match trigger
  Double_t	fTrackChi2[300];		// single mu chi2 track
  Double_t	fMatchTrigChi2[300];	// single mu chi2 of match trigger
  Double_t	fDCA[300];		// single mu DCA
  Int_t	        fCharge[300];		// single mu charge
  Double_t	fRAtAbsEnd[300];		// single mu distance from beam center at end abs

  Int_t		fNDimu;			// dimuons in the event
  Int_t		fDimuMu[3000][2];	// reference to single mus
  Double_t	fDimuPt[3000];			    // dimuon pT
  Double_t	fDimuPx[3000]; 		    // dimuon px
  Double_t	fDimuPy[3000]; 		    // dimuon py
  Double_t	fDimuPz[3000]; 		    // dimuon pz
  Double_t	fDimuY[3000];			// dimuon y
  Double_t	fDimuMass[3000];		// dimuon invariant mass
  Int_t	        fDimuCharge[3000];		// dimuon charge
  Int_t	        fDimuMatch[3000];		// dimuon match
  Bool_t        fIsPileupFromSPDInMultBins;   //check on pile-up
  Bool_t        fIsPhysSelected;   //check on phys selection 
  
  AliAODEvent* fAODEvent;      //! AOD event  //tolgo !
//  TList *fOutput;  //!< List of histograms for data
  
 ClassDef(AliAnalysisTaskXeXeJPsiTree_Dimuon,1);
};

#endif

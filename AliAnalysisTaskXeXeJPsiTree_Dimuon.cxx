/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliAnalysisTaskXeXeJPsiTree_Dimuon.cxx $ */

//-----------------------------------------------------------------------------
// Analysis task to compute muon/dimuon kinematic distributions
// The output is a list of histograms.
// The macro class can run on AOD or in the train with the ESD filter.
// R. Arnaldi
//
//-----------------------------------------------------------------------------

//#ifndef AliAnalysisTaskXeXeJPsiTree_Dimuon_CXX
//#define AliAnalysisTaskXeXeJPsiTree_Dimuon_CXX

// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TList.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"

#include "AliInputEventHandler.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"  //for Centrality

#include "AliAnalysisTaskXeXeJPsiTree_Dimuon.h"

ClassImp(AliAnalysisTaskXeXeJPsiTree_Dimuon)
//__________________________________________________________________________
AliAnalysisTaskXeXeJPsiTree_Dimuon::AliAnalysisTaskXeXeJPsiTree_Dimuon() :
  AliAnalysisTaskSE(),
  fOutputTree(0x0),
  fNevt(0x0),
  fBeamEnergy(0.),
  fkAnalysisType(0x0),
  fPeriod(0x0),
  fCountTotEv(0x0),
  fCountTrigger(0x0),
  fCountCINT7(0x0),
  fCountCMUL7(0x0),
  fCountCMLL7(0x0),
  fCountCMSL7(0x0),
  fCountCMSH7(0x0),
  fNMuons(0x0),
  fNTracklets(0x0),
  fNContributors(0x0),
  fNDimu(0x0),
  fPercentV0M(0x0),       
  fPercentCL0(0x0),       
  fPercentCL1(0x0),       
  fPercentV0A(0x0),       
  fPercentV0C(0x0),       
  fPercentZNA(0x0),       
  fPercentZNC(0x0),    
  fIsPileupFromSPDInMultBins(0x0),
  fIsPhysSelected(0x0),
  fAODEvent(0x0)
{
  //
  //Default ctor
  //  
  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<300;i++){
    fPt[i]=999.;
    fE[i]=999.;
    fPx[i]=999; 
    fPy[i]=999; 
    fPz[i]=999; 
    fY[i]=999.; 
    fEta[i]=999.; 
    fMatchTrig[i]=999.; 
    fTrackChi2[i]=999.; 
    fMatchTrigChi2[i]=999.;
    fDCA[i]=999.;
    fCharge[i]=999;
    fRAtAbsEnd[i]=999;
  }
  for(Int_t i=0; i<3000;i++){  
    fDimuPt[i]=999.; 
    fDimuPx[i]=999.; 
    fDimuPy[i]=999.; 
    fDimuPz[i]=999.; 
    fDimuY[i]=999.; 
    fDimuMass[i]=999.;
    fDimuCharge[i]=999;
    fDimuMatch[i]=999;
    for(Int_t k=0;k<2;k++) fDimuMu[i][k]=999;
  } 
  
//  DefineOutput(1,TTree::Class());  //added

}

//__________________________________________________________________________
AliAnalysisTaskXeXeJPsiTree_Dimuon::AliAnalysisTaskXeXeJPsiTree_Dimuon(const char *name) :
  AliAnalysisTaskSE(name),
  fOutputTree(0x0),
  fNevt(0x0),
  fBeamEnergy(0.),
  fkAnalysisType(0x0),
  fPeriod(0x0),
  fCountTotEv(0x0),
  fCountTrigger(0x0),
  fCountCINT7(0x0),
  fCountCMUL7(0x0),
  fCountCMLL7(0x0),
  fCountCMSL7(0x0),
  fCountCMSH7(0x0),
  fNMuons(0x0),
  fNTracklets(0x0),
  fNContributors(0x0),
  fNDimu(0x0),
  fPercentV0M(0x0),       
  fPercentCL0(0x0),       
  fPercentCL1(0x0),       
  fPercentV0A(0x0),       
  fPercentV0C(0x0),       
  fPercentZNA(0x0),       
  fPercentZNC(0x0),       
  fIsPileupFromSPDInMultBins(0x0),
  fIsPhysSelected(0x0),
  fAODEvent(0x0)
{
 //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskXeXeJPsiTree_Dimuon","Calling Constructor");
  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<300;i++){
    fPt[i]=999.;
    fE[i]=999.;
    fPx[i]=999; 
    fPy[i]=999; 
    fPz[i]=999; 
    fY[i]=999.; 
    fEta[i]=999.; 
    fMatchTrig[i]=999.; 
    fTrackChi2[i]=999.; 
    fMatchTrigChi2[i]=999.;
    fDCA[i]=999.;
    fCharge[i]=999;
    fRAtAbsEnd[i]=999;
  }
  for(Int_t i=0; i<3000;i++){  
    fDimuPt[i]=999.; 
    fDimuPx[i]=999.; 
    fDimuPy[i]=999.; 
    fDimuPz[i]=999.; 
    fDimuY[i]=999.; 
    fDimuMass[i]=999.;
    fDimuCharge[i]=999;
    fDimuMatch[i]=999;
    for(Int_t k=0;k<2;k++) fDimuMu[i][k]=999;
  } 
  
  DefineOutput(1,TTree::Class());
}

//___________________________________________________________________________
AliAnalysisTaskXeXeJPsiTree_Dimuon& AliAnalysisTaskXeXeJPsiTree_Dimuon::operator=(const AliAnalysisTaskXeXeJPsiTree_Dimuon& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
   // fNevt = c.fNevt ;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskXeXeJPsiTree_Dimuon::AliAnalysisTaskXeXeJPsiTree_Dimuon(const AliAnalysisTaskXeXeJPsiTree_Dimuon& c) :
  AliAnalysisTaskSE(c),  
  fOutputTree(c.fOutputTree),
  fNevt(c.fNevt),
  fBeamEnergy(c.fBeamEnergy),
  fkAnalysisType(c.fkAnalysisType),
  fPeriod(c.fPeriod),
  fCountTotEv(c.fCountTotEv),
  fCountTrigger(c.fCountTrigger),
  fCountCINT7(c.fCountCINT7),
  fCountCMUL7(c.fCountCMUL7),
  fCountCMLL7(c.fCountCMLL7),
  fCountCMSL7(c.fCountCMSL7),
  fCountCMSH7(c.fCountCMSH7),
  fNMuons(c.fNMuons),
  fNTracklets(c.fNTracklets),
  fNContributors(c.fNContributors),
  fNDimu(c.fNDimu),
  fPercentV0M(c.fPercentV0M),       
  fPercentCL0(c.fPercentCL0),       
  fPercentCL1(c.fPercentCL1),       
  fPercentV0A(c.fPercentV0A),       
  fPercentV0C(c.fPercentV0C),       
  fPercentZNA(c.fPercentZNA),       
  fPercentZNC(c.fPercentZNC),       
  fIsPileupFromSPDInMultBins(c.fIsPileupFromSPDInMultBins),
  fIsPhysSelected(c.fIsPhysSelected),
  fAODEvent(c.fAODEvent)
 {
  //
  // Copy Constructor									
  //
}

//___________________________________________________________________________
AliAnalysisTaskXeXeJPsiTree_Dimuon::~AliAnalysisTaskXeXeJPsiTree_Dimuon() {
  //
  //destructor
  //
  Info("~AliAnalysisTaskXeXeJPsiTree_Dimuon","Calling Destructor");
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis) delete fOutputTree;
}


//___________________________________________________________________________
void AliAnalysisTaskXeXeJPsiTree_Dimuon::UserCreateOutputObjects(){

  // devo mettere la phys selection nel tree???
    
   if (fOutputTree) return; //already initialised ADDED
   
  OpenFile(1,"RECREATE");
  fOutputTree = new TTree("XeXeTree","Data Tree");

  fOutputTree->Branch("FiredTriggerClasses",fTrigClass,"FiredTriggerClasses/C");
  fOutputTree->Branch("inpmask",&finpmask,"inpmask/i"); //unsigned integer

  fOutputTree->Branch("NMuons",&fNMuons,"NMuons/I");
  fOutputTree->Branch("NContributors",&fNContributors,"NContributors/I");
  fOutputTree->Branch("NTracklets",&fNTracklets,"NTracklets/I");
  fOutputTree->Branch("Vertex",fVertex,"Vertex[3]/D");
  fOutputTree->Branch("PercentV0M",&fPercentV0M,"PercentV0M/F");
  fOutputTree->Branch("PercentCL0",&fPercentCL0,"PercentCL0/F");
  fOutputTree->Branch("PercentCL1",&fPercentCL1,"PercentCL1/F");
  fOutputTree->Branch("PercentV0A",&fPercentV0A,"PercentV0A/F");
  fOutputTree->Branch("PercentV0C",&fPercentV0C,"PercentV0C/F");
  fOutputTree->Branch("PercentZNA",&fPercentZNA,"PercentZNA/F");
  fOutputTree->Branch("PercentZNC",&fPercentZNC,"PercentZNC/F");

  fOutputTree->Branch("Pt",fPt,"Pt[NMuons]/D");
  fOutputTree->Branch("E",fE,"E[NMuons]/D");
  fOutputTree->Branch("Px",fPx,"Px[NMuons]/D");
  fOutputTree->Branch("Py",fPy,"Py[NMuons]/D");
  fOutputTree->Branch("Pz",fPz,"Pz[NMuons]/D");
  fOutputTree->Branch("Y",fY,"Y[NMuons]/D");
  fOutputTree->Branch("Eta",fEta,"Eta[NMuons]/D");
  fOutputTree->Branch("MatchTrig",fMatchTrig,"MatchTrig[NMuons]/I");
  fOutputTree->Branch("TrackChi2",fTrackChi2,"TrackChi2[NMuons]/D");
  fOutputTree->Branch("MatchTrigChi2",fMatchTrigChi2,"MatchTrigChi2[NMuons]/D");
  fOutputTree->Branch("DCA",fDCA,"DCA[NMuons]/D");
  fOutputTree->Branch("Charge",fCharge,"Charge[NMuons]/I");
  fOutputTree->Branch("RAtAbsEnd",fRAtAbsEnd,"RAtAbsEnd[NMuons]/D");
 
  fOutputTree->Branch("NDimu",&fNDimu,"NDimu/I");
  fOutputTree->Branch("DimuMu",fDimuMu,"DimuMu[NDimu][2]/I");
  fOutputTree->Branch("DimuPt",fDimuPt,"DimuPt[NDimu]/D");
  fOutputTree->Branch("DimuPx",fDimuPx,"DimuPx[NDimu]/D");
  fOutputTree->Branch("DimuPy",fDimuPy,"DimuPy[NDimu]/D");
  fOutputTree->Branch("DimuPz",fDimuPz,"DimuPz[NDimu]/D");
  fOutputTree->Branch("DimuY",fDimuY,"DimuY[NDimu]/D");
  fOutputTree->Branch("DimuMass",fDimuMass,"DimuMass[NDimu]/D");
  fOutputTree->Branch("DimuCharge",fDimuCharge,"DimuCharge[NDimu]/I");
  fOutputTree->Branch("DimuMatch",fDimuMatch,"DimuMatch[NDimu]/I");
  fOutputTree->Branch("IsPileupFromSPDInMultBins",&fIsPileupFromSPDInMultBins,"IsPileupFromSPDInMultBins/O");
  fOutputTree->Branch("IsPhysSelected",&fIsPhysSelected,"IsPhysSelected/O");

  fOutputTree->ls(); 

 PostData(1,fOutputTree); 
 
} 

//_________________________________________________
void AliAnalysisTaskXeXeJPsiTree_Dimuon::UserExec(Option_t *)
{

  // mettere taglio in massa?
   
  fNMuons=0; 
  fNTracklets=-1;
  fNContributors=-1;
  fNDimu=0;
  fPercentV0M=-1.;	
  fPercentCL0=-1.;	
  fPercentCL1=-1.;	
  fPercentV0A-1.;	
  fPercentV0C=-1.;	
  fPercentZNA=-1.;	
  fPercentZNC=-1.;	
  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<300;i++){
    fPt[i]=999.;
    fE[i]=999.;
    fPx[i]=999; 
    fPy[i]=999; 
    fPz[i]=999; 
    fY[i]=999.; 
    fEta[i]=999.; 
    fMatchTrig[i]=999.; 
    fTrackChi2[i]=999.; 
    fMatchTrigChi2[i]=999.;
    fDCA[i]=999.;
    fCharge[i]=999;
    fRAtAbsEnd[i]=999;
  }
  for(Int_t i=0; i<3000;i++){  
    fDimuPt[i]=999.; 
    fDimuPx[i]=999.; 
    fDimuPy[i]=999.; 
    fDimuPz[i]=999.; 
    fDimuY[i]=999.; 
    fDimuMass[i]=999.;
    fDimuCharge[i]=999.;
    fDimuMatch[i]=0;
    for(Int_t k=0;k<2;k++) fDimuMu[i][k]=999;
  } 
 
//
// Execute analysis for current event
//
  fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  if ( ! fAODEvent ) {
    AliError ("AOD event not found. Nothing done!");
    return;
  }
  
  char hname[200];
  AliAODHeader *aodheader=dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
   TString firedtrigger = aodheader->GetFiredTriggerClasses();
   sprintf(fTrigClass,"%s",firedtrigger.Data());
   
  finpmask = aodheader->GetL0TriggerInputs(); 
   
 //   to apply physics selection
   UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
   fIsPhysSelected = fSelectMask & (AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonLikeLowPt7 | AliVEvent::kMuonSingleLowPt7 | AliVEvent::kMuonSingleHighPt7  | AliVEvent::kINT7inMUON  | AliVEvent::kINT7);               
   // warning: I added AliVEvent::kINT7 (MB in CENT cluster, while kINT7inMUON is MB in MUFAST)
   
   
   //printf("trigger = %s\n\n",fTrigClass);
  
  // centrality 
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PACentStudiesRun2
  Float_t PercV0M = 300; 
  Float_t PercCL0 = 300; 
  Float_t PercCL1 = 300; 
  Float_t PercV0A = 300; 
  Float_t PercV0C = 300; 
  Float_t PercZNA = 300;
  Float_t PercZNC = 300;
  
  AliMultSelection *MultSelection = 0x0; 
    MultSelection = (AliMultSelection * ) fAODEvent->FindListObject("MultSelection");
    if( !MultSelection) {
      //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
      AliWarning("AliMultSelection object not found!");
    }else{
      PercV0M = MultSelection->GetMultiplicityPercentile("V0M");
//       PercCL0 = MultSelection->GetMultiplicityPercentile("CL0");
//       PercCL1 = MultSelection->GetMultiplicityPercentile("CL1");
//       PercV0A = MultSelection->GetMultiplicityPercentile("V0A");
//       PercV0C = MultSelection->GetMultiplicityPercentile("V0C");
//       PercZNA = MultSelection->GetMultiplicityPercentile("ZNA");
//       PercZNC = MultSelection->GetMultiplicityPercentile("ZNC");
   }	 
  fPercentV0M = PercV0M;
  fPercentCL0 = PercCL0;
  fPercentCL1 = PercCL1;
  fPercentV0A = PercV0A;
  fPercentV0C = PercV0C;
  fPercentZNA = PercZNA;
  fPercentZNC = PercZNC;

//  fIsPileupFromSPDInMultBins = fAODEvent->IsPileupFromSPDInMultBins();
  
  AliAODVertex *PrimVertex =  fAODEvent->GetPrimaryVertex();
  fNContributors = PrimVertex->GetNContributors();
  //fNTracklets = PrimVertex->GetNumberOfTracklets();
  fVertex[0]=PrimVertex->GetX();
  fVertex[1]=PrimVertex->GetY();
  fVertex[2]=PrimVertex->GetZ();
   
    //
    // build dimuon object starting from single muons
    //
    //------------------------------------------------------------------------- 
    // uncomment if reading AliMuonAODs (in PbPb)
    //-------------------------------------------------------------------------
     Int_t numdimu = 0;
     Int_t nummu = 0;
//     TRefArray *mutracks = new TRefArray();
//     Int_t nmuons = fAODEvent->GetMuonTracks(mutracks);
//     for (Int_t i=0;i<nmuons;i++){
//       AliAODTrack *mu0=(AliAODTrack*)mutracks->At(i);
//       fCharge[i] = mu0->Charge();
//       fPt[i] = mu0->Pt();
//       fPx[i] = mu0->Px();
//       fPy[i] = mu0->Py();
//       fPz[i] = mu0->Pz();
//       fY[i]  = mu0->Y();
//       fEta[i]= mu0->Eta();
//       fE[i] = mu0->E();
// //      fDCA[i] = mu0->GetDCA();
// //      fTrackChi2[i] = mui->GetChi2()/(2.*mu0->GetNHit()-5.);
//       fMatchTrig[i]   = mu0->GetMatchTrigger();
//       fMatchTrigChi2[i]= mu0->GetChi2MatchTrigger();
//       fRAtAbsEnd[i]=mu0->GetRAtAbsorberEnd();
//       for(Int_t j=i+1;j<nmuons;j++){
//         AliAODTrack *mu1=(AliAODTrack*)mutracks->At(j);
//         AliAODDimuon *dimu=new AliAODDimuon(mu0,mu1);
   
    //----------------------------------------------------------
    // uncomment if reading AODs or old AliAOD.Muons (as old pA)
    //----------------------------------------------------------
    //	
        Int_t ntracks = fAODEvent->GetNumberOfTracks();   
        for (Int_t i=0;i<ntracks;i++){
	  
	  printf("i=%d \n",i);
          AliAODTrack *mu0=(AliAODTrack*)fAODEvent->GetTrack(i);
          fCharge[i] = mu0->Charge();
          fPt[i] = mu0->Pt();
          fPx[i] = mu0->Px();
          fPy[i] = mu0->Py();
          fPz[i] = mu0->Pz();
          fY[i]  = mu0->Y();
          fEta[i]= mu0->Eta();
          fE[i] = mu0->E();
//        fDCA[i] = mu0->GetDCA();
//        fTrackChi2[i] = mui->GetChi2()/(2.*mu0->GetNHit()-5.);
          fMatchTrig[i]   = mu0->GetMatchTrigger();
          fMatchTrigChi2[i]= mu0->GetChi2MatchTrigger();
          fRAtAbsEnd[i]=mu0->GetRAtAbsorberEnd();
	
          if (!mu0->IsMuonTrack()) continue;
          for(Int_t j=i+1;j<ntracks;j++){
	  printf("j=%d \n",j);
            AliAODTrack *mu1=(AliAODTrack*)fAODEvent->GetTrack(j);
            if (!mu1->IsMuonTrack()) continue;
           AliAODDimuon *dimu=new AliAODDimuon(mu0,mu1);
    //----------------------------------------------------------
	   
//	if(dimu->Mass()>1.8){
          fDimuMass[numdimu] = dimu->Mass();
          fDimuPt[numdimu] = dimu->Pt();
          fDimuPx[numdimu] = dimu->Px();
          fDimuPy[numdimu] = dimu->Py();
          fDimuPz[numdimu] = dimu->Pz();
          fDimuY[numdimu] = dimu->Y();
          fDimuCharge[numdimu]= dimu->Charge();
          fDimuMu[numdimu][0]=i;  fDimuMu[numdimu][1]=j;

//       printf("dimuon Pt %f\n",dimu->Pt());
//       printf("dimuon Px %f\n",dimu->Px());
//       printf("dimuon Py %f\n",dimu->Py());
//       printf("dimuon Pz %f\n",dimu->Pz());
//       printf("dimuon Y %f\n",dimu->Y());
//       printf("dimuon mass %f\n",dimu->Mass());
//       printf("dimuon charge %d\n",dimu->Charge());
       
        if(mu0->GetMatchTrigger()>1 || mu1->GetMatchTrigger()>1){
          fDimuMatch[numdimu]=1;
        }	
        if(mu0->GetMatchTrigger()>1 && mu1->GetMatchTrigger()>1){
          fDimuMatch[numdimu]=2;
        }	
        numdimu++;      
//       delete mu1;
//         delete dimu;
//	} // mass cut
         delete dimu;
//       delete mu1;	
      }
     nummu++;
//     delete mu0;
  }
  fNMuons =nummu;
  fNDimu=numdimu;     
  fOutputTree->Fill();
  PostData(1,fOutputTree);
  
}


//________________________________________________________________________
void AliAnalysisTaskXeXeJPsiTree_Dimuon::Terminate(Option_t *) 
{

 }



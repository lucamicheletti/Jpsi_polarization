/**************************************************************************
 * Copyright(c) 2013-2014, ALICE Experiment at CERN, All rights reserved. *
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

//#ifdef __ECLIPSE_IDE
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "TSystem.h"
#include "TChain.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"

extern AliAnalysisGrid* CreateAlienHandler_Grid(const char *, Bool_t);
AliAnalysisTaskPbPbJPsiTree_Dimuon* AddTaskPbPbJPsiTree_Dimuon_Grid(Int_t );
#endif

void RunPbPb_JPsi_Grid(const char *RunMode = "test",
                       Int_t RunNumber = 246980, //without 000
		                   Bool_t usePhysicsSelection = kTRUE,
                       //TString GridDir="/alice/data/2017/LHC17n",
                       //TString DataPattern="*muon_calo_pass2/AOD/*/AliAOD.root",
                       TString GridDir="/alice/data/2015/LHC15o",
                       TString DataPattern="*muon_calo_pass1/AOD197/*/AliAOD.root",
	                     TString AliPhysicsVersion="vAN-20180107-1",
	                     Bool_t gridMerge = kTRUE){
  // RunMode can be set to "test" to test running of the code
  // RunMode can be set to "full" to produce Trees
  // RunMode can be set to "terminate" to do the merging of the Trees

  //============================================================================
  // Load paths
  //============================================================================
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY");
  // lines to be included to set correctly the path to find my macros
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ROOTSYS/include");
  //============================================================================
  // Handler and analysis manager
  //============================================================================
  AliAnalysisManager *mgr = new AliAnalysisManager("Dimuon Tree creator");
  mgr -> SetDebugLevel(AliLog::kError);
  AliLog::SetGlobalLogLevel(AliLog::kDebug);

  AliAODInputHandler* aodH = new AliAODInputHandler();
  mgr -> SetInputEventHandler(aodH);

  AliAnalysisGrid *alienHandler = NULL;
  gROOT -> LoadMacro("CreateAlienHandler_Grid.C");
  alienHandler = CreateAlienHandler_Grid(RunMode, GridDir, DataPattern, RunNumber, AliPhysicsVersion, gridMerge);
  if (!alienHandler) return;
  mgr -> SetGridHandler(alienHandler);

  AliAnalysisDataContainer *cinput1 = NULL;
  cinput1 = mgr -> GetCommonInputContainer();
  //============================================================================
  // Physics selection
  //============================================================================
  gROOT -> ProcessLine(".L $ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask;
  if(usePhysicsSelection) physSelTask = AddTaskPhysicsSelection(kFALSE,kFALSE);    //second kTRUE to switch on the pile-up cut in PS
  //============================================================================
  // Read centrality information (https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PACentStudiesRun2)
  //============================================================================
  gROOT -> LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask *task = AddTaskMultSelection(kFALSE); // kFALSE for user mode: kTRUE for for calibration mode:
  task -> SetUseDefaultCalib(kTRUE); // to use the default calib for runs not yet calibrated
  //============================================================================
  // Add task to build Dimuon tree
  //============================================================================
  gROOT -> LoadMacro("AliAnalysisTaskPbPbJPsiTree_Dimuon.cxx++g");  //NEEDED to run LOCALLY!!!!
  gROOT -> LoadMacro("AddTaskPbPbJPsiTree_Dimuon_Grid.C");
  AliAnalysisTaskPbPbJPsiTree_Dimuon *taskTree = AddTaskPbPbJPsiTree_Dimuon_Grid(RunNumber);
  mgr -> AddTask(taskTree);
  //============================================================================
  // Init analysis on GRID
  //============================================================================
  if (!mgr -> InitAnalysis()) return;
  mgr -> PrintStatus();
  mgr -> StartAnalysis("grid");
}

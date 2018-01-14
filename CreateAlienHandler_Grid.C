#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliAnalysisAlien.h"
#endif

AliAnalysisGrid* CreateAlienHandler_Grid(const char *runMode, TString GridDir, TString DataPattern, Int_t RunNumber, TString AliPhysicsVersion, Bool_t gridMerge)
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.

  AliAnalysisAlien *plugin = new AliAnalysisAlien();

  //Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin -> SetRunMode(runMode);

  plugin -> SetNtestFiles(1); // num of test files in "test" mode

  // Set versions of used packages
  plugin -> SetAliPhysicsVersion(AliPhysicsVersion.Data());

  plugin -> AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  plugin -> SetDataPattern(DataPattern.Data());
  plugin -> SetGridDataDir(GridDir.Data()); // Data
  plugin -> SetRunPrefix("000");

  /* alternatively provide run number */
  plugin -> AddRunNumber(RunNumber);
  // Alternatively use run range
  //  plugin->SetRunRange(138653, 138666);

  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  char outdirname[250];
  sprintf(outdirname,"PbPb2015/LHC15o/muon_calo_pass1/NO_CUTS/%d",RunNumber); //without PS
  plugin -> SetGridWorkingDir(outdirname);  // NOTE: Change name here every new run!!!eclare alien output directory. Relative to working directory.
  plugin -> SetGridOutputDir("OutputTree"); // In this case will be $HOME/work/output
  plugin -> SetOutputToRunNo(kFALSE); // we want the run number as output subdirectory
  plugin -> SetDefaultOutputs(kTRUE);
 // plugin->SetMergeExcludes("AliAOD.Muons.root");

  plugin -> SetMergeViaJDL(gridMerge);

  // added by me
  plugin -> SetAnalysisSource("AliAnalysisTaskPbPbJPsiTree_Dimuon.cxx");
  plugin -> SetAdditionalLibs("AliAnalysisTaskPbPbJPsiTree_Dimuon.cxx AliAnalysisTaskPbPbJPsiTree_Dimuon.h AddTaskPbPbJPsiTree_Dimuon_Grid.C");

  // Declare the output file names separated by blanks.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  plugin -> SetOverwriteMode(kFALSE);

  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  Int_t nNoOfInputFiles = 30; //default is usually 0
  if(nNoOfInputFiles != 0)
    plugin -> SetSplitMaxInputFileNumber(nNoOfInputFiles);
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  //   plugin->SetMaxInitFailed(15);
  // Optionally resubmit threshold.
  plugin -> SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
  plugin -> SetTTL(60000); //default was 18000
  // Optionally set input format (default xml-single)
  plugin -> SetInputFormat("xml-single");
  return plugin;
}

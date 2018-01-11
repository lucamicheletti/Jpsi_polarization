AliAnalysisTaskPbPbJPsiTree_Dimuon *AddTaskPbPbJPsiTree_Dimuon_Grid(Int_t RunNumber){

//****************************************************************************************
// Add task class to fill a tree with dimuon infos
// Roberta
//****************************************************************************************

   printf("Creating Task for Muon/Dimuon Histos in PbPb\n");

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskPbPbJPsi_Dimuon", "No analysis manager to connect to.");
      return NULL;
   }
   TString fnameout_norun = "Tree_%d.root";
   TString fnameout;
   fnameout.Form(fnameout_norun.Data(), RunNumber);
   printf("Fnameout = %s\n",fnameout.Data());

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("ctree0",TTree::Class(),AliAnalysisManager::kOutputContainer,fnameout);

   AliAnalysisTaskPbPbJPsiTree_Dimuon *PbPbJPsiTask = new AliAnalysisTaskPbPbJPsiTree_Dimuon("AliAnalysisTaskPbPbJPsiTree_Dimuon");

   PbPbJPsiTask->SetBeamEnergy(5.44);  // define by hand the beam energy
   mgr -> AddTask(PbPbJPsiTask);

   mgr -> ConnectInput(PbPbJPsiTask,0,mgr->GetCommonInputContainer());
   mgr -> ConnectOutput(PbPbJPsiTask,1,coutput1);

   return PbPbJPsiTask;
}

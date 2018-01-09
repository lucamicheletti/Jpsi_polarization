AliAnalysisTaskXeXeJPsiTree_Dimuon *AddTaskXeXeJPsiTree_Dimuon_Grid(Int_t RunNumber){

//****************************************************************************************
// Add task class to fill a tree with dimuon infos
// Roberta
//****************************************************************************************

   printf("Creating Task for Muon/Dimuon Histos in XeXe\n");
    
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskXeXeJPsi_Dimuon", "No analysis manager to connect to.");
      return NULL;
   }   
   TString fnameout_norun = "Tree_%d.root";  
   TString fnameout;  
   fnameout.Form(fnameout_norun.Data(), RunNumber);
   printf("Fnameout = %s\n",fnameout.Data());
  
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("ctree0",TTree::Class(),AliAnalysisManager::kOutputContainer,fnameout);

   AliAnalysisTaskXeXeJPsiTree_Dimuon *XeXeJPsiTask = new AliAnalysisTaskXeXeJPsiTree_Dimuon("AliAnalysisTaskXeXeJPsiTree_Dimuon");
   
   XeXeJPsiTask->SetBeamEnergy(5.44);  // define by hand the beam energy
   mgr->AddTask(XeXeJPsiTask);
 
   mgr->ConnectInput(XeXeJPsiTask,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(XeXeJPsiTask,1,coutput1);
  
   return XeXeJPsiTask;
}

/*
Author:    Mohit Gola 10th July 2023
*/

#include "run.hh"

#include <fstream>
#include <sstream>

MyRunAction::MyRunAction()
{
  G4AnalysisManager::Instance();
}

MyRunAction::~MyRunAction()
{
  auto *man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
}

void MyRunAction::BeginOfRunAction(const G4Run *run)
{
  auto *man = G4AnalysisManager::Instance();
  auto *runManager = G4RunManager::GetRunManager();

  if (!fileOpened)
  {
    int run_num = 0;
    G4String filename;

    do
    {
      std::ostringstream oss;
      oss << name_template;
      if (run_num > 0)
        oss << run_num;
      oss << ".root";

      filename = oss.str();
      run_num++;
    }
    while (std::ifstream(filename).good());

    outputFilename = filename;
    man->OpenFile(outputFilename);

    // Master ntuple (created ONCE)
    master_ntupleId = man->CreateNtuple("Photons_Master", "Photons_Master");
    man->CreateNtupleIColumn(master_ntupleId, "Scanpoint_ID");
    man->CreateNtupleDColumn(master_ntupleId, "PosX_Initial");
    man->CreateNtupleDColumn(master_ntupleId, "PosY_Initial");
    man->CreateNtupleDColumn(master_ntupleId, "PosZ_Initial");
    man->FinishNtuple(master_ntupleId);

    fileOpened = true;
  }

  runId = run->GetRunID();

  std::ostringstream strRunId;
  strRunId << runId;

  scanpoint_ntupleId =
      man->CreateNtuple("Photons_" + strRunId.str(), "Photons");

  const G4UserEventAction *eventAction = runManager->GetUserEventAction();
  auto *myEventAction =
      dynamic_cast<MyEventAction *>(const_cast<G4UserEventAction *>(eventAction));

  myEventAction->SetScanpointNtupleID(scanpoint_ntupleId);
  myEventAction->ResetPhotNum();

  man->CreateNtupleIColumn(scanpoint_ntupleId, "Photon_ID");
  man->CreateNtupleIColumn(scanpoint_ntupleId, "Step_Number");
  man->CreateNtupleIColumn(scanpoint_ntupleId, "Step_Status");
  man->CreateNtupleDColumn(scanpoint_ntupleId, "PosX");
  man->CreateNtupleDColumn(scanpoint_ntupleId, "PosY");
  man->CreateNtupleDColumn(scanpoint_ntupleId, "PosZ");
  man->FinishNtuple(scanpoint_ntupleId);

  TotalNumAbsorbed = 0;
  TotalNumReflected = 0;
  TotalNumTransmitted = 0;
}

void MyRunAction::EndOfRunAction(const G4Run *)
{
  auto *man = G4AnalysisManager::Instance();
  auto *runManager = G4RunManager::GetRunManager();

  const G4UserEventAction *eventAction = runManager->GetUserEventAction();
  if (eventAction)
  {
    auto *myEventAction =
        dynamic_cast<MyEventAction *>(const_cast<G4UserEventAction *>(eventAction));

    if (myEventAction)
    {
      man->FillNtupleIColumn(master_ntupleId, 0, runId);
      man->FillNtupleDColumn(master_ntupleId, 1, myEventAction->GetPosX());
      man->FillNtupleDColumn(master_ntupleId, 2, myEventAction->GetPosY());
      man->FillNtupleDColumn(master_ntupleId, 3, myEventAction->GetPosZ());
      man->AddNtupleRow(master_ntupleId);

      myEventAction->ResetCounters();
    }
  }
}



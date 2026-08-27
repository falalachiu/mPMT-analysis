/*
Author:   Mohit Gola 10th July 2023
*/

#include <iostream>

#include "G4RunManager.hh"
#include "G4UIExecutive.hh"
#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4SteppingVerbose.hh"
#include "G4SteppingManager.hh"

#include "physics.hh"
#include "mysteppingaction.hh"
#include "WCSimTuningParameters.hh"
#include "generator.hh"
#include "MyPrimaryGeneratorMessenger.hh"

#include "G4VPhysicalVolume.hh"
#include "G4StepLimiterPhysics.hh"

#include "action.hh"

#include "constructInSitumPMT.hh"
#include "constructionPMT.hh"
#include "inSituPMT.hh"
using genmech = Laser;

int main(int argc, char **argv)
{
  G4RunManager *runManager = new G4RunManager();
  // G4SteppingVerbose::SetInstance(new G4SteppingVerbose);
  runManager->SetVerboseLevel(0);
  WCSimTuningParameters *tuningpars = new WCSimTuningParameters();
  MyPhysicsList *phys_list = new MyPhysicsList();
  runManager->SetUserInitialization(phys_list);
  //  G4VSteppingVerbose* steppingVerbose = new G4SteppingVerbose;
  //  G4VSteppingVerbose::SetInstance(steppingVerbose);

  enum DetConfiguration
  {
    wfm = 1,
    fwm = 2
  };
  G4int WCSimConfiguration = fwm;

  genmech *generator = new genmech();

  G4UImanager *UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/control/execute tuning_parameters.mac");
  UImanager->ApplyCommand("/gui/addMenu true");
  UImanager->ApplyCommand("/tracking/verbose 0");

  MyActionInitialization<genmech> *actionman = new MyActionInitialization<genmech>(
      generator,
      "../../projects/def-blairt2k/fchiu/outputs/output_mPMT");

  runManager->SetUserInitialization(actionman);
  // pmtConstruction *subtractionPMT = new pmtConstruction();
  // G4VSolid *PMTSolid = subtractionPMT->ConstructionPMT();

  ConstructInSitumPMT *InSitumPMT = new ConstructInSitumPMT();
  runManager->SetUserInitialization(InSitumPMT);

  runManager->Initialize();

  G4UIExecutive *ui = 0;

  if (argc == 1)
  {
    ui = new G4UIExecutive(argc, argv);
  }

  G4VisManager *visManager = new G4VisExecutive();
  visManager->Initialize();

  MyPrimaryGeneratorMessenger<genmech> *generatorMessenger = new MyPrimaryGeneratorMessenger<genmech>(generator);

  if (ui)
  {
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    //    delete ui;
  }
  else
  {
    G4String command = "/control/execute ";
    G4String filename = argv[1];
    UImanager->ApplyCommand(command + filename);
    G4cout << "FileName = " << filename << G4endl;
  }
  //  file_exists("tuning_parameters.mac");

  //  ui->SessionStart();

  delete ui;
  delete runManager;
  delete visManager;
  return 0;
}

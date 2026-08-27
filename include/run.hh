/*
Author:    Mohit Gola 10th July 2023
*/

#ifndef RUN_HH
#define RUN_HH

#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include <G4AnalysisManager.hh>
#include "G4RunManager.hh"

#include "event.hh"

class MySteppingAction;

class MyRunAction : public G4UserRunAction
{
public:
  MyRunAction();
  ~MyRunAction();

  void set_name_template(const G4String name) { name_template = name; }

  virtual void BeginOfRunAction(const G4Run *);
  virtual void EndOfRunAction(const G4Run *);

private:
  G4String name_template = "output";
  G4String outputFilename;
  G4bool fileOpened = false;
  G4int TotalNumAbsorbed;
  G4int TotalNumReflected;
  G4int TotalNumTransmitted;
  G4int master_ntupleId;
  G4int scanpoint_ntupleId;
  G4int runId;
};

#endif

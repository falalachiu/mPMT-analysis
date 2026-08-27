/*
Author:    Mohit Gola 10th July 2023
*/

#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"
#include "G4Event.hh"
#include "G4OpticalPhoton.hh"
#include <G4AnalysisManager.hh>
#include "G4DynamicParticle.hh"
#include "MyPrimaryGeneratorMessenger.hh"

const float A = 254.564663 * mm;
const float B = 254.110205 * mm;
const float E = 0.5 * (A + B);
const float C = 186.002389 * mm;

class Laser : public G4VUserPrimaryGeneratorAction
{
public:
  Laser();

  virtual void GeneratePrimaries(G4Event *);
  G4ThreeVector calculateNewPositionAndAngle(const G4ThreeVector &initialPosition, double initialAngle,
                                             G4ThreeVector &newPosition, double &newAngle);
  void SetAngle(G4double angle);
  void SetPAzimuthAngle(G4double aziangle);
  void SetPZenithAngle(G4double zenangle);
  void SetEnergy(G4double energy);
  void SetDiscRad(G4double discrad);
  void SetX(G4double xpos);
  void SetY(G4double ypos);
  void SetZ(G4double zpos);
  void SetSpread(G4double spread);

private:
  G4ParticleGun *fParticleGun;
  G4double angleDegrees;
  G4double discRadius = 0 * mm;
  G4double angleRadians;
  G4double pAzimuthAngle = 0.0;
  G4double pZenithAngle = 0.0;
  G4double particleEnergy;
  G4double xpos = 0.0;
  G4double ypos = 0.0;
  G4double zpos = 500.0 * mm;
  G4double spread = 0.0;
  MyPrimaryGeneratorMessenger<Laser> *fMessenger;
};

class PrecisionGen : public G4VUserPrimaryGeneratorAction
{
public:
  PrecisionGen();

  virtual void GeneratePrimaries(G4Event *);
  G4ThreeVector calculateNewPositionAndAngle(const G4ThreeVector &initialPosition, double initialAngle,
                                             G4ThreeVector &newPosition, double &newAngle);
  void SetAngle(G4double angle);
  void SetPAzimuthAngle(G4double aziangle);
  void SetPZenithAngle(G4double zenangle);
  void SetEnergy(G4double energy);
  void SetDiscRad(G4double discrad);
  void SetX(G4double xpos);
  void SetY(G4double ypos);
  void SetZ(G4double zpos);
  void SetSpread(G4double spread);

private:
  G4ParticleGun *fParticleGun;
  G4double angleDegrees;
  G4double discRadius = 1 * mm;
  G4double angleRadians;
  G4double pAzimuthAngle = 0.0;
  G4double pZenithAngle = 0.0;
  G4double particleEnergy;
  G4double xpos = 0.0;
  G4double ypos = 0.0;
  G4double zpos = 500.0 * mm;
  G4double spread = 0.0;
  MyPrimaryGeneratorMessenger<PrecisionGen> *fMessenger;
};

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
  MyPrimaryGenerator();
  ~MyPrimaryGenerator();

  virtual void GeneratePrimaries(G4Event *);
  G4ThreeVector calculateNewPositionAndAngle(const G4ThreeVector &initialPosition, double initialAngle,
                                             G4ThreeVector &newPosition, double &newAngle);
  void SetAngle(G4double angle);
  void SetPAzimuthAngle(G4double aziangle);
  void SetPZenithAngle(G4double zenangle);
  void SetEnergy(G4double energy);
  void SetDiscRad(G4double discrad);
  void SetX(G4double xpos) {};
  void SetY(G4double ypos) {};
  void SetZ(G4double zpos) {};
  void SetSpread(G4double spread) {};

private:
  G4ParticleGun *fParticleGun;
  G4double angleDegrees;
  G4double angleRadians;
  G4double pAzimuthAngle = 0.0;
  G4double pZenithAngle = 0.0;
  G4double discRadius = 250 * mm;
  G4double particleEnergy;
  MyPrimaryGeneratorMessenger<MyPrimaryGenerator> *fMessenger;
};

class NewGenerator : public G4VUserPrimaryGeneratorAction
{
public:
  NewGenerator();

  virtual void GeneratePrimaries(G4Event *);
  G4ThreeVector calculateNewPositionAndAngle(const G4ThreeVector &initialPosition, double initialAngle,
                                             G4ThreeVector &newPosition, double &newAngle);
  void SetAngle(G4double angle);
  void SetPAzimuthAngle(G4double aziangle);
  void SetPZenithAngle(G4double zenangle);
  void SetEnergy(G4double energy);
  void SetDiscRad(G4double discrad);
  void SetX(G4double xpos) {};
  void SetY(G4double ypos) {};
  void SetZ(G4double zpos) {};
  void SetSpread(G4double spread) {};

private:
  G4ParticleGun *fParticleGun;
  G4double angleDegrees;
  G4double angleRadians;
  G4double pAzimuthAngle = 0.0;
  G4double pZenithAngle = 0.0;
  G4double discRadius = 250 * mm;
  G4double particleEnergy;
  MyPrimaryGeneratorMessenger<NewGenerator> *fMessenger;
};

#endif

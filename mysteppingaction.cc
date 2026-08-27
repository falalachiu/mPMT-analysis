/*
Author:    Mohit Gola 10th July 2023
*/

#include "mysteppingaction.hh"

MySteppingAction::MySteppingAction()
{
  fEventAction = nullptr;
  photonPositions.clear();
}

MySteppingAction::~MySteppingAction()
{
}

void MySteppingAction::SetEventAction(MyEventAction *eventAction)
{
  fEventAction = eventAction;
}

double MySteppingAction::calculateIncidenceAngle(const G4ThreeVector &Momentum, const G4ThreeVector &Normal)
{
  G4ThreeVector FacetNormal;

  // Calculate alpha using the distribution p(alpha) = g(alpha; 0, sigma_alpha)*std::sin(alpha)
  // where g(alpha; 0, sigma_alpha) is a Gaussian distribution with mean 0 and standard deviation sigma_alpha.
  G4double alpha;
  G4double sigma_alpha = 0.0; // Set sigma_alpha value based on your requirements
  const G4double halfpi = 0.5 * CLHEP::pi;
  const G4double twopi = 2.0 * CLHEP::pi;

  if (sigma_alpha == 0.0)
  {
    FacetNormal = Normal;
  }
  else
  {
    double f_max = std::min(1.0, 4.0 * sigma_alpha);
    double phi, SinAlpha, CosAlpha, SinPhi, CosPhi, unit_x, unit_y, unit_z;
    G4ThreeVector tmpNormal;

    do
    {
      do
      {
        alpha = G4RandGauss::shoot(0.0, sigma_alpha);
      } while (G4UniformRand() * f_max > std::sin(alpha) || alpha >= halfpi);

      phi = G4UniformRand() * twopi;
      SinAlpha = std::sin(alpha);
      CosAlpha = std::cos(alpha);
      SinPhi = std::sin(phi);
      CosPhi = std::cos(phi);
      unit_x = SinAlpha * CosPhi;
      unit_y = SinAlpha * SinPhi;
      unit_z = CosAlpha;
      FacetNormal.setX(unit_x);
      FacetNormal.setY(unit_y);
      FacetNormal.setZ(unit_z);
      tmpNormal = Normal;
      FacetNormal.rotateUz(tmpNormal);
    } while (Momentum * FacetNormal >= 0.0);
  }
  // Calculate the incidence angle between Momentum and FacetNormal
  G4double PdotN = Momentum * FacetNormal;
  G4double magP = Momentum.mag();
  G4double magN = FacetNormal.mag();

  G4double incidenceAngle = CLHEP::pi - std::acos(PdotN / (magP * magN));
  return incidenceAngle;
}

void MySteppingAction::RecordAbsorption(MyEventAction *EventAction, G4Track *Track, G4String vol, bool opAbsorption, G4String postvol)
{
  // pmtPhysBulb and  pmtInnerPhysBulb
  if ((vol == "pmtPhysBulb" && postvol == "pmtInnerPhysBulb") || (vol == "pmtInnerPhysBulb" && postvol == "pmtPhysBulb") && !opAbsorption)
  {
    fEventAction->IncrementNumDetected();
    EventAction->RecordStep(2, Track->GetPosition(), 1);
  }
  // opAbsorption in phisCath (not at boundary) won't produce a pe, count it as absorbed in glass
  else if (vol == "totalPMT" || vol == "pmtPhysBulb" || vol == "pmtPhysReflector" || vol == "pmtPhysInnerTube" || vol == "pmtInnerPhysBulb" || vol == "innerReflector")
  {
    fEventAction->IncrementNumAbsorbed();
    EventAction->RecordStep(3, Track->GetPosition(), 1);
  }
  else if (vol == "pmtAbsorber")
  {
    fEventAction->IncrementNumAbsorbed();
    EventAction->RecordStep(4, Track->GetPosition(), 1);
  }
  else if (vol == "physWorld")
  {
    fEventAction->IncrementNumAbsorbed();
    EventAction->RecordStep(5, Track->GetPosition(), 1);
  }
  else if (vol == "Matrix")
  {
    fEventAction->IncrementNumAbsorbed();
    EventAction->RecordStep(6, Track->GetPosition(), 1);
  }
  else if (vol == "gelPhys")
  {
    fEventAction->IncrementNumAbsorbed();
    EventAction->RecordStep(7, Track->GetPosition(), 1);
  }
  else if (vol == "physDome")
  {
    fEventAction->IncrementNumAbsorbed();
    EventAction->RecordStep(8, Track->GetPosition(), 1);
  }
  else if (vol == "physCylinder")
  {
    fEventAction->IncrementNumAbsorbed();
    EventAction->RecordStep(9, Track->GetPosition(), 1);
  }
  else
  {
    bool found = false;
    G4String this_name;
    for (int pmtIndex = 0; pmtIndex < 20; pmtIndex++)
    {
      this_name = "PMTCopy" + G4String(std::to_string(pmtIndex));
      if (vol == this_name)
      {
        fEventAction->IncrementNumAbsorbed();
        EventAction->RecordStep(11 + pmtIndex, Track->GetPosition(), 1);
        found = true;
      }
    }
    if (not found)
    {
      std::ostringstream volerror;
      volerror << "This volume hasn't been accounted for: " << vol;
      throw std::runtime_error(volerror.str());
    }
  }
}

void MySteppingAction::UserSteppingAction(const G4Step *step)
{

  G4LogicalVolume *currentVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4String volume = currentVolume->GetName();

  G4Track *track = step->GetTrack();
  G4int trackID = track->GetTrackID();
  /////////////////////////////////////////////////////////////////
  if (track->GetCurrentStepNumber() == 2)
  {
    G4ThreeVector hitPosition = step->GetPreStepPoint()->GetPosition();
    //      G4cout<<"HIT POSITON == "<<hitPosition<<G4endl;
  }

  G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
  //  G4cout<<"PHOTON POSITON == "<<position<<G4endl;

  //////////////////////////incident angle//////////////////////////////////////

  G4ThreeVector momentum = step->GetTrack()->GetMomentumDirection();
  G4ThreeVector normal = step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetSolid()->SurfaceNormal(step->GetPreStepPoint()->GetPosition());

  // Calculate the incidence angle using the calculateIncidenceAngle function
  G4double incidenceAngle = calculateIncidenceAngle(momentum, normal);
  //  G4cout<<"THE INCIDENT ANGLE BEFORE LOOP:  "<<incidenceAngle<<G4endl;
  if (track->GetCurrentStepNumber() == 1)
  {
    G4double checkAngle = -(incidenceAngle - 0.5 * CLHEP::pi);
    //      fEventAction->SetPhotonAngle(checkAngle);
    fEventAction->SetPhotonAngle(incidenceAngle);

    fEventAction->SetPosX(position.getX());
    fEventAction->SetPosY(position.getY());
    fEventAction->SetPosZ(position.getZ());
  }
  ///////////////////////////////////////////////////////////////////////////////////
  //  if (volume == "pmtGlassLogic" || volume == "pmtInnerGlassLogic"){
  bool is_absorbed = false;
  bool is_other = false;

  if (track->GetDefinition()->GetParticleName() == "opticalphoton")
  {
    G4bool photonTransmitted = false;
    // Get the pre-step and post-step points of the step
    G4StepPoint *prePoint = step->GetPreStepPoint();
    G4StepPoint *postPoint = step->GetPostStepPoint();

    // Get the volume at the pre-step and post-step points
    G4VPhysicalVolume *preVolume = prePoint->GetPhysicalVolume();
    G4VPhysicalVolume *postVolume = postPoint->GetPhysicalVolume();

    // Check if the photon is interacting with an optical boundary
    if (preVolume && postVolume && preVolume != postVolume)
    {
      G4String preVolumeName = preVolume->GetName();
      G4String postVolumeName = postVolume->GetName();

      WCSimOpBoundaryProcessStatus boundaryStatus = Undefined;
      WCSimOpBoundaryProcess *opBoundary = NULL;

      // Get the particle definition for the current track
      G4ParticleDefinition *particleDefinition = track->GetDefinition();
      // Get the process manager for the particle definition
      G4ProcessManager *processManager = particleDefinition->GetProcessManager();
      // Get the boundary process for the current track
      G4int nProcesses = processManager->GetProcessListLength();
      //	  G4cout<<"NUMBER OF PROCESSES == "<<nProcesses<<G4endl;
      G4ProcessVector *processVector = processManager->GetProcessList();
      for (G4int i = 0; i < nProcesses; i++)
      {
        G4VProcess *process = (*processVector)[i];

        //	      G4cout<<"PROCESS NAMES == "<<process->GetProcessName()<<G4endl;
        if (process->GetProcessName() == "OpBoundary")
        {
          opBoundary = static_cast<WCSimOpBoundaryProcess *>(process);
          // G4cout<<"BOUNDARY PROCESS NAME == "<<opBoundary->GetProcessName()<<G4endl;
          if (opBoundary == nullptr)
          {
            std::cerr << "BOUNDARY IS NULL POINTER" << std::endl;
          }
          boundaryStatus = opBoundary->GetStatus();

          BoundaryMeta status = meta_status(boundaryStatus);

          if (preVolumeName == "mesh_grid_abyss" || postVolumeName == "mesh_grid_abyss")
          {
            G4double rand = G4UniformRand();
            if (rand > 0.50)
            {
              track->SetTrackStatus(fStopAndKill);
            }
          }
          else if (preVolumeName == "dynode_internal_logical" || postVolumeName == "dynode_internal_logical")
          {
            track->SetTrackStatus(fStopAndKill);
          }

          // NoRINDEX is recorded as absorption
          if (status == BNoRINDEX)
          {
            track->SetTrackStatus(fStopAndKill);
            // std::cout << "stepping from " << preVolume->GetName() << " to " << postVolume->GetName() << std::endl;
            // throw std::runtime_error("no r index");
          }
          if (status == BAbsorption) // && (volume == "pmtGlassLogic" || volume == "pmtInnerGlassLogic"))
          {
            // fEventAction->IncrementNumAbsorbed();
            if (status == BAbsorption)
            {
              RecordAbsorption(fEventAction, track, preVolumeName, false, postVolumeName);
            }
            else if (status == BNoRINDEX)
            {
              RecordAbsorption(fEventAction, track, postVolumeName, false);
            }
            photonPositions.push_back(position);
            is_absorbed = true;
          }

          else if (status == BTransmission) // && (volume == "pmtGlassLogic"|| volume == "pmtInnerGlassLogic"))
          {
            fEventAction->IncrementNumTransmitted();
            fEventAction->RecordStep(0, track->GetPosition(), 0);
          }

          else if (status == BReflection) // && (volume == "pmtGlassLogic" || volume == "pmtInnerGlassLogic"))
          {
            fEventAction->IncrementNumReflected();
            fEventAction->RecordStep(1, track->GetPosition(), 0);
          }

          else if (status == BOther)
          {
            if (track->GetTrackStatus() == fStopAndKill)
            {
              fEventAction->RecordStep(11, track->GetPosition(), 1);
            }
            else
            {
              fEventAction->RecordStep(11, track->GetPosition(), 0);
            }
            is_other = true;
          }
        }
      }
    }
  }
  else
  {
    // std::cout << track->GetDefinition()->GetParticleName();
  }
  //    }
  //	}

  ////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////Photon Absorption/////////////////////////////////
  const skDetCon *detectorConstruction = static_cast<const skDetCon *>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  const myDetectorConstruction *pmtDetCon = static_cast<const myDetectorConstruction *>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();

  if (track->GetTrackStatus() == fStopAndKill && !is_absorbed && !is_other)
  {
    G4String procname = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    if (!step->GetPostStepPoint()->GetPhysicalVolume())
    {
      // If photon leaves the world volume
      fEventAction->RecordStep(10, track->GetPosition(), 1);
    }
    else if (procname == "OpAbsorption")
    {
      // if the photon is not absorbed at a boundary but while it passes through material
      RecordAbsorption(fEventAction, track, step->GetPostStepPoint()->GetPhysicalVolume()->GetName(), true);
    }
    else if (procname == "Transportation")
    {
      RecordAbsorption(fEventAction, track, step->GetPostStepPoint()->GetPhysicalVolume()->GetName(), true);
    }
    else
    {
      std::ostringstream procerror;
      procerror << "This process hasn't been accounted for: " << procname;
      throw std::runtime_error(procerror.str());
    }
  }

  // I think this "if" block is unnecessary but haven't tested that
  if (currentVolume == fScoringVolume)
  {
    //      G4cout<< "ABSORBER VOLUME  == "<<fScoringVolume->GetName()<<G4endl;
    //      G4cout<<"STEP IS WITHIN THE TARGET LOGICAL VOLUME!!!"<<G4endl;
    //  G4Track* track = step->GetTrack();
    // is_absorbed

    if (false)
    {
      track->SetTrackStatus(fStopAndKill);
    }
  }
  else
  {
    G4double edep = step->GetTotalEnergyDeposit();
    fEventAction->AddEdep(edep);
    // step->GetTrack()->SetTrackStatus(fStopAndKill);
  }
  ///////////////////////////////////////////////////////////////////////////////////
}

G4bool MySteppingAction::IsTrackCounted(G4int trackID)
{
  return (trackStatus[trackID] == "Absorbed" || trackStatus[trackID] == "Reflected" || trackStatus[trackID] == "Transmitted");
}

const std::vector<G4ThreeVector> &MySteppingAction::GetPhotonPositions() const
{
  return photonPositions;
}

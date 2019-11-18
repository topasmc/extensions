//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id: TsTrackerHit.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file TsTrackerHit.hh
/// \brief Definition of the TsTrackerHit class

#ifndef TsTrackerHit_h
#define TsTrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"

/// Tracker hit class
///
/// It defines data members to store the trackID, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fEdep, fPos

class TsTrackerHit : public G4VHit
{
  public:
    TsTrackerHit();
    TsTrackerHit(const TsTrackerHit&);
    virtual ~TsTrackerHit();

    // operators
    const TsTrackerHit& operator=(const TsTrackerHit&);
    G4int operator==(const TsTrackerHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // Set methods
    void SetTrackID  (G4int track)      { fTrackID = track; };
    void SetEdep     (G4double de)      { fEdep = de; };
    void SetPos      (G4ThreeVector xyz){ fPos = xyz; };

    void SetIncidentEnergy (G4double e) {fIncidentEnergy = e;};

    void SetParticleName(G4String name) {fParticleName = name;}
    void SetParticleFlag(G4int flag) {fParticleFlag = flag;}

    // Get methods
    G4int GetTrackID() const     { return fTrackID; };
    G4double GetEdep() const     { return fEdep; };
    G4ThreeVector GetPos() const { return fPos; };

    G4double GetIncidentEnergy() const {return fIncidentEnergy;};
    G4String GetParticleName() const {return fParticleName;}
    G4int    GetParticleFlag() const {return fParticleFlag;}
  private:

      G4int         fTrackID;
      G4double      fEdep;
      G4ThreeVector fPos;
      G4String      fParticleName;
      G4int         fParticleFlag;

      G4double      fIncidentEnergy;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<TsTrackerHit> TsTrackerHitsCollection;

extern G4ThreadLocal G4Allocator<TsTrackerHit>* TsTrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* TsTrackerHit::operator new(size_t)
{
  if(!TsTrackerHitAllocator)
      TsTrackerHitAllocator = new G4Allocator<TsTrackerHit>;
  return (void *) TsTrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void TsTrackerHit::operator delete(void *hit)
{
  TsTrackerHitAllocator->FreeSingle((TsTrackerHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

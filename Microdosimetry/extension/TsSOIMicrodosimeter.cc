// Component for TsSOIMicrodosimeter
//
// ********************************************************************
// *                                                                  *
// * This  code implementation is an extension to developments by the *
// * TOPAS collaboration.                                             *
// * This extension  is an freely  available in  accordance  with the *
// * freeBSD license:                                                 *
// * Copyright (c) <2015>, <Harald Paganetti>                         *
// * All rights reserved.                                             *
// * Redistribution    and   use in   source and   binary    forms,   *
// * with or without modification, are permitted provided that the    *
// * following conditions are met:                                    *
// *                                                                  *
// *                                                                  *
// * 1. Redistributions of source code must retain the above          *
// * copyright notice, this                                           *
// * list of conditions and the following disclaimer.                 *
// * 2. Redistributions in binary form must reproduce the above       *
// * copyright notice, this list of conditions and the following      *
// * disclaimer in the documentation and/or other materials provided  *
// * with the distribution.                                           *
// *                                                                  *
// * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
// * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
// * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
// * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
// * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR             *
// * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
// * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT *
// * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF *
// * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED  *
// * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT      *
// * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING   *
// * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF   *
// * THE POSSIBILITY OF SUCH DAMAGE.                                  *
// *                                                                  *
// * The views and conclusions contained in the software and          *
// * documentation are those of the authors and should not be         *
// * interpreted as representing official policies, either expressed  *
// * or implied, of the FreeBSD Project.                              *
// *                                                                  *
// * Contacts: Jan Schuemann, jschuemann@mgh.harvard.edu              *
// *           Harald Paganetti, hpaganetti@mgh.harvard.edu           *
// *                                                                  *
// ********************************************************************
//

#include "TsSOIMicrodosimeter.hh"

#include "TsParameterManager.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcommand.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

//bool
#include "G4BooleanSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"


TsSOIMicrodosimeter::TsSOIMicrodosimeter(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			   TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


TsSOIMicrodosimeter::~TsSOIMicrodosimeter()
{
}


G4VPhysicalVolume* TsSOIMicrodosimeter::Construct()
{
	BeginConstruction();
	
	G4NistManager* man = G4NistManager::Instance();
    G4Material* Si = man->FindOrBuildMaterial("G4_Si");
    G4Material* SiO2 = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    G4Material* Al = man->FindOrBuildMaterial("G4_Al");
    
	//********************************************************************************************
	//-------------------------------------------------------------------------------------------
	//-------------------------- Detector Stats for Region Generation ---------------------------
	//-------------------------------------------------------------------------------------------
 
	G4double SVheight = 10.*micrometer;
	G4double insulationThickness = 1.*micrometer;
	G4double SiOoverlayerBottomThickness = 1.7*micrometer;
	G4double AlOverlayerThickness = 1.7*micrometer;
	G4double AlOverlayerwidth = 4.*micrometer;
	G4double SiOoverlayerTopThickness = 1.43*micrometer;
	G4double SiOoverlayerTopwidth = 10.5 * micrometer;
	G4double overLayerThickness = AlOverlayerThickness + SiOoverlayerTopThickness;
 
	G4double baseSiThickness = 300.*micrometer;
	G4double detectorHeight = (SVheight + baseSiThickness + insulationThickness + overLayerThickness);
	G4double maxElecRange = 2.*mm;
 
	G4double SVwidth = 30.*micrometer;
	G4double pitch = 20.*micrometer; //distance between the odd and even rows
 
	G4double bridgingWidth = 20.*micrometer;
	G4double bridgingLength = 15.*micrometer;
	G4double bridgingHeight = SVheight;
 
	G4double numberOfRows = 59;
	G4double numberOfColumns = (24*3);
	G4double SVareaWidth = (numberOfColumns * (SVwidth + bridgingWidth)) - 1.*bridgingWidth;
	G4double SVareaLength = numberOfRows*(SVwidth + pitch) - 1.*pitch;
	G4double bufferWidth = 100.*micrometer;
	G4double detectorWidth = SVareaWidth + bufferWidth;
	G4double detectorLength = SVareaLength + bufferWidth;
	
	//G4cout<<"SVareaLength="<<SVareaLength<<G4endl;
	//G4cout<<"SVareaWidth="<<SVareaWidth<<G4endl;

	//------------------------------------------------------------------------------------
	//---------------------------------- SOI Detector  -----------------------------------
	//------------------------------------------------------------------------------------	

    // Detector Envelope
    G4double HLX= SVareaLength/2.;
    G4double HLY= SVareaWidth /2.;
    G4double HLZ= detectorHeight/2.;

    G4Box* Detector = new G4Box(fName, HLX, HLY, HLZ);
    fEnvelopeLog = CreateLogicalVolume(Detector);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    // ----------------------------------- Si Base ---------------------------------------
    HLX= SVareaLength/2.;
    HLY= SVareaWidth /2.;
    HLZ= baseSiThickness/2.;
    G4double PosZ=7.065*um;
    G4String MaSi = "G4_Si";
    G4String subComponentName = "SiBase";

    G4Box* solidSiBase = new G4Box(subComponentName, HLX, HLY, HLZ);
	G4LogicalVolume* logicSiBase = CreateLogicalVolume(subComponentName, MaSi,solidSiBase);
	G4VPhysicalVolume* physicSiBase = CreatePhysicalVolume(subComponentName, logicSiBase, new G4RotationMatrix(), new G4ThreeVector(0,0,PosZ),fEnvelopePhys);

    // ------------------------------- SIO2 Insulation -----------------------------------
    HLX= SVareaLength/2.;
    HLY= SVareaWidth /2.;
    HLZ= insulationThickness/2.;
    PosZ=-143.435*um;
    G4String MaSiO2 = "G4_SILICON_DIOXIDE";
    subComponentName = "SiO2layer";

    G4Box* solidSiO2layer = new G4Box(subComponentName, HLX, HLY, HLZ);
	G4LogicalVolume* logicSiO2layer = CreateLogicalVolume(subComponentName, MaSiO2,solidSiO2layer);
	G4VPhysicalVolume* physicSiO2layer = CreatePhysicalVolume(subComponentName, logicSiO2layer, new G4RotationMatrix(), new G4ThreeVector(0,0,PosZ),fEnvelopePhys);
	
	// ------------------------------- Sensitive Volume ----------------------------------
	HLX= SVwidth/2.;
    HLY= SVwidth/2.;
    HLZ= SVheight/2.;
    PosZ=-148.935*um;
    subComponentName = "SensitiveVolume";
    G4double CenterXminSV = -SVareaLength/2.+SVwidth/2.;
    G4double CenterYminSV = -SVareaWidth/2. +SVwidth/2.;

    G4Box* solidSensitiveVolume = new G4Box(subComponentName, HLX, HLY, HLZ);
	G4LogicalVolume* logicSensitiveVolume = CreateLogicalVolume(subComponentName, MaSi,solidSensitiveVolume);
    for (int i = 0; i < 72; i++){
        for (int j = 0; j < 59; j++){

            G4double Px = CenterXminSV + j*50*um;
            G4double Py = CenterYminSV + i*50*um;
            G4int Noij =i*59+j;  
	        G4VPhysicalVolume* physicSensitiveVolume = CreatePhysicalVolume(subComponentName, Noij, true, logicSensitiveVolume, new G4RotationMatrix(), new G4ThreeVector(Px,Py,PosZ),fEnvelopePhys);

        }
    }
    
    // ---------------------------------- Bridge -----------------------------------------
    HLX= bridgingLength/2.;
    HLY= bridgingWidth/2.;
    HLZ= bridgingHeight/2.;
    PosZ=-148.935*um;
    subComponentName = "BridgeVolume";
    G4double CenterXminBridge = -SVareaLength/2.+bridgingLength/2.;
    G4double CenterYminBridge = -SVareaWidth/2. +bridgingWidth/2.+SVwidth ; 
    G4Box* solidBridge = new G4Box(subComponentName, HLX, HLY, HLZ);
	G4LogicalVolume* logicBridge = CreateLogicalVolume(subComponentName, MaSi,solidBridge);
	
	for (int i = 0; i < 72; i++){
        for (int j = 0; j < 59; j++){

            G4double Px = CenterXminBridge + j*50*um;
            G4double Py = CenterYminBridge + i*50*um;
            G4int Noij =i*59+j;  
	        G4VPhysicalVolume* physicBridge = CreatePhysicalVolume(subComponentName, Noij, true, logicBridge, new G4RotationMatrix(), new G4ThreeVector(Px,Py,PosZ),fEnvelopePhys);

        }
    }

    // ------------------------ SiO2 Over Layer on SV Top --------------------------------
    HLX= SVwidth/2.;
    HLY= SVwidth /2.;
    HLZ= SiOoverlayerBottomThickness/2.;
    PosZ=-154.785*um;
    subComponentName = "SiO2OverLayerTopSV";
    G4Box* solid1 = new G4Box("solid1", HLX, HLY, HLZ);
    //G4Box* solid2 = new G4Box("solid2", AlOverlayerwidth/2,  AlOverlayerwidth/2,  AlOverlayerThickness/2.);
    //G4SubtractionSolid* solidSiO2OverLayerTopSV  = new G4SubtractionSolid(subComponentName,solid1, solid2);
	G4LogicalVolume* logicSiO2OverLayerTopSV = CreateLogicalVolume(subComponentName, MaSiO2,solid1);
	
	for (int i = 0; i < 72; i++){
        for (int j = 0; j < 59; j++){

            G4double Px = CenterXminSV + j*50*um;
            G4double Py = CenterYminSV + i*50*um;
            G4int Noij =i*59+j;  
	        G4VPhysicalVolume* physicSiO2OverLayerTopSV = CreatePhysicalVolume(subComponentName, Noij, true, logicSiO2OverLayerTopSV, new G4RotationMatrix(), new G4ThreeVector(Px,Py,PosZ),fEnvelopePhys);

        }
    }

    	    
    //-------------------------------- Al Over Layer -------------------------------------
    HLX= AlOverlayerwidth/2.;
    HLY= AlOverlayerwidth/2.;
    HLZ= AlOverlayerThickness/2.;
    PosZ=-154.785*um;
    G4String MaAl = "G4_Al";
    subComponentName = "AlOverLayer";   
    G4Box* solidAlOverLayer = new G4Box(subComponentName, HLX, HLY, HLZ);
	G4LogicalVolume* logicAlOverLayer = CreateLogicalVolume(subComponentName, MaAl,solidAlOverLayer);
	
    for (int i = 0; i < 72; i++){
        for (int j = 0; j < 59; j++){

            G4double Px = CenterXminSV + j*50*um;
            G4double Py = CenterYminSV + i*50*um;
            G4int Noij =i*59+j;  
	        G4VPhysicalVolume* physicAlOverLayer = CreatePhysicalVolume(subComponentName, Noij, true, logicAlOverLayer, new G4RotationMatrix(), new G4ThreeVector(Px,Py,PosZ),fEnvelopePhys);

        }
    }
    
	    
	// ------------------------ SiO2 Over Layer on Bridge Top ----------------------------
    HLX= bridgingLength/2.;
    HLY= bridgingWidth/2.;
    HLZ= SiOoverlayerBottomThickness/2.;
    PosZ=-154.785*um;
    subComponentName = "SiO2OverLayerTopBridge";
    G4Box* solidSiO2OverLayerTopBridge = new G4Box(subComponentName, HLX, HLY, HLZ);
	G4LogicalVolume* logicSiO2OverLayerTopBridge = CreateLogicalVolume(subComponentName, MaSiO2,solidSiO2OverLayerTopBridge);
	
	for (int i = 0; i < 72; i++){
        for (int j = 0; j < 59; j++){

            G4double Px = CenterXminBridge + j*50*um;
            G4double Py = CenterYminBridge + i*50*um;
            G4int Noij =i*59+j;  
	        G4VPhysicalVolume* physicSiO2OverLayerTopBridge = CreatePhysicalVolume(subComponentName, Noij, true, logicSiO2OverLayerTopBridge, new G4RotationMatrix(), new G4ThreeVector(Px,Py,PosZ),fEnvelopePhys);

        }
    }
    
    // ------------------------ SiO2 Over Layer Top Top ----------------------------------
    HLX= SiOoverlayerTopwidth/2.;
    HLY= SiOoverlayerTopwidth/2.;
    HLZ= SiOoverlayerTopThickness/2.;
    PosZ=-156.35*um;
    subComponentName = "SiO2OverLayerTopTop";
    G4Box* solidSiO2OverLayerTopTop = new G4Box(subComponentName, HLX, HLY, HLZ);
	G4LogicalVolume* logicSiO2OverLayerTopTop = CreateLogicalVolume(subComponentName, MaSiO2,solidSiO2OverLayerTopTop);
		
	for (int i = 0; i < 72; i++){
        for (int j = 0; j < 59; j++){

            G4double Px = CenterXminSV + j*50*um;
            G4double Py = CenterYminSV + i*50*um;
            G4int Noij =i*59+j;  
	        G4VPhysicalVolume* physicSiO2OverLayerTopTop = CreatePhysicalVolume(subComponentName, Noij, true, logicSiO2OverLayerTopTop, new G4RotationMatrix(), new G4ThreeVector(Px,Py,PosZ),fEnvelopePhys);

        }
    }

	InstantiateChildren(fEnvelopePhys);

	return fEnvelopePhys;
}

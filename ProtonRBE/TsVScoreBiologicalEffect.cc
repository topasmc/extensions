// Extra Class for AllScorers
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

#include "TsVScoreBiologicalEffect.hh"

#include "TsParameterManager.hh"

TsVScoreBiologicalEffect::TsVScoreBiologicalEffect(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                         G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVBinnedScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{;}


TsVScoreBiologicalEffect::~TsVScoreBiologicalEffect()
{;}


void TsVScoreBiologicalEffect::PostConstructor()
{
    TsVBinnedScorer::PostConstructor();

    G4String* cellLines = fPm->GetStringVector(GetFullParmName("CellLines"));
    G4int cellLinesLength = fPm->GetVectorLength(GetFullParmName("CellLines"));

    G4String* structureNames = NULL;
    G4int structureNamesLength = 0;
    if (fPm->ParameterExists(GetFullParmName("RTStructures"))) {
        structureNames = fPm->GetStringVector(GetFullParmName("RTStructures"));
        structureNamesLength = fPm->GetVectorLength(GetFullParmName("RTStructures"));
    }

    if (cellLinesLength != structureNamesLength+1) {
        G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
        G4cerr << "Scorer \"" << GetName() << "\" has the wrong number of" << G4endl;
        G4cerr << "CellLines and/or RTStructures." << G4endl;
        exit(1);
    }

    // Construct a biological effect model for each cell line
    for (G4int i = 0; i < cellLinesLength; i++)
        if (fModelsByCellLine.find(cellLines[i]) == fModelsByCellLine.end())
            fModelsByCellLine[cellLines[i]] = ConstructModel(cellLines[i]);

    // Associate RTStructures with biological effect models
    for (G4int i = 0; i < structureNamesLength; i++) {
        G4int id = fComponent->GetStructureID(structureNames[i]);
        if (id == -1) {
            G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
            G4cerr << "Component: " << fComponent->GetNameWithCopyId() << " does not have stucture: " << structureNames[i] << G4endl;
            exit(1);
        }
        fRTStructureIDs.push_back(id);
        fRTStructureNames.push_back(structureNames[i]);
        fRTStructureCellLines.push_back(cellLines[i]);
        fModelsByStructureID[id] = fModelsByCellLine[cellLines[i]];
    }

    // Set default cell line
    fDefaultCellLine = cellLines[cellLinesLength-1];
    fModelsByStructureID[-1] = fModelsByCellLine[cellLines[cellLinesLength-1]];
}


TsVModelBiologicalEffect* TsVScoreBiologicalEffect::GetModelForVoxel(G4int index)
{
    G4int structureID = -1; // default model

    // reverse iterate through structures, to give priority to earlier structures
    for (std::vector<G4int>::reverse_iterator it_id = fRTStructureIDs.rbegin(); it_id != fRTStructureIDs.rend(); ++it_id)
    {
        if (fComponent->IsInNamedStructure(*it_id, index))
            structureID = *it_id;
    }

    return fModelsByStructureID[structureID];
}

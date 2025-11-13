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

//GEANT4 - Depth-of-Interaction enabled Positron emission tomography (PET) advanced example 

//Authors and contributors

// Author list to be updated, with names of co-authors and contributors from National Institute of Radiological Sciences (NIRS)

// Abdella M. Ahmed (1, 2), Andrew Chacon (1, 2), Harley Rutherford (1, 2),
// Hideaki Tashima (3), Go Akamatsu (3), Akram Mohammadi (3), Eiji Yoshida (3), Taiga Yamaya (3)
// Susanna Guatelli (2), and Mitra Safavi-Naeini (1, 2)

// (1) Australian Nuclear Science and Technology Organisation, Australia
// (2) University of Wollongong, Australia
// (3) National Institute of Radiological Sciences, Japan

//Implemetation of the doiPETAnalysis.cc class
//This implementation mimics readout (or digitizer) of the PET scanner. To mimic realistic PET detector, the signals are blurred. Blurring 
//parameters are given in inputParameters.txt file. Crystal based energy resolution and quantum efficiency has been applied. Deadtime (on 
//each detector block and axially multiplexed detector) is also applied before the event is rejected by the energy window. The units for 
//blurring parameters are in keV (for energy) and ns (nano sec) for dead time. If the units are different, exception will be thrown and the
//program quits. In this class, an ideal PMT is assumed to be placed at the corners of the crystal block. First, the ideal interaction position
//(obtained by G4) is used to determine the distance of the PMT from the interaction point. The signal on each PMT depends on the lateral (2D)
//distance from the PMTs. Light sharing method (reflector based) DOI identification method has been used. If the crystal ID is out of bound, 
//error message will be displayed and the event will be rejected. The output file is single based list-mode ASCII file and can be then be 
//processed into coinsidence list-mode data. As, an option, binary output method is also given.  
//Explanation is given for the methods provided. 

as
	//dist1y = dist3y, and dist2y = dist4y, so only take two of them or take the average
	//disty = ((dist1y + dist3y) + (dist2y+dist4y))/2;
	disty = dist1y + dist2y;

	//signalPMT1z = signalPMT2z, and signalPMT3z = signalPMT4z
	signalPMT1z = Edep * dist3z/(dist1z + dist3z);
	signalPMT3z = Edep * dist1z/(dist1z + dist3z);

	signalPMT2z = Edep * dist4z/(dist2z + dist4z);
	signalPMT4z = Edep * dist2z/(dist2z + dist4z);


	//signalPMT1y = signalPMT3y, and signalPMT2y = signalPMT4y
	signalPMT1y = Edep * dist2y/(dist1y + dist2y);
	signalPMT2y = Edep * dist1y/(dist1y + dist2y);

	signalPMT3y = Edep * dist4y/(dist3y + dist4y);
	signalPMT4y = Edep * dist3y/(dist3y + dist4y);

	//Calculate the signal on each PMT from the 'component' signal
	signalPMT1 = (signalPMT1z +  signalPMT1y)*0.5; 
	signalPMT2 = (signalPMT2z +  signalPMT2y)*0.5;
	signalPMT3 = (signalPMT3z +  signalPMT3y)*0.5;
	signalPMT4 = (signalPMT4z +  signalPMT4y)*0.5;


	signalZplus = (signalPMT3 + signalPMT4);
	signalZminus = (signalPMT1 + signalPMT2);
	signalYplus = (signalPMT2 + signalPMT4);
	signalYminus = (signalPMT1 + signalPMT3);


	//Position of interaction is calculated based on Anger logic method. 
	//To get the position by Anger calculation, the result should be multiplied by the dimenion of the total distance.
	PositionAngerZ = (signalZplus - signalZminus)/(signalZplus + signalZminus)*distz; 
	PositionAngerY = (signalYplus - signalYminus)/(signalYplus + signalYminus)*disty;


	//For detectors with reflector insertion (light sharing), the response is shifted depending on the reflector patter.
	//Here, it is assumed that the shift of the response is equal to half of the distance from the interaction position to the airgap in the lateral (transversal direction)

	//If reflector is only in the left side of the crystal, then response shift is to the right side (away from the reflector).
	//If reflector is only in the right side of the crystal, then response shift is to the left side (away from the reflector).

	//Response shift for 1st Layer
	if(i_doi == 0){
		//If reflector is only in one (left) side of the crystal, then response shifts to the right side (away from the reflector)
		if(ireflectorLayer1_Tangential[i_tan] == 1 && ireflectorLayer1_Tangential[i_tan + 1] == 0) PositionAngerY += (crystalPitch_tan/2)*shiftDis;

		//If reflector is only in one (right) side of the crystal, then response shifts to the left side (away from the reflector)
		if(ireflectorLayer1_Tangential[i_tan] == 0 && ireflectorLayer1_Tangential[i_tan + 1] == 1) PositionAngerY -= (crystalPitch_tan/2)*shiftDis;

		if(ireflectorLayer1_Axial[i_axial] == 1 && ireflectorLayer1_Axial [i_axial + 1] == 0) PositionAngerZ += (crystalPitch_axial/2)*shiftDis;
		if(ireflectorLayer1_Axial[i_axial] == 0 && ireflectorLayer1_Axial [i_axial + 1] == 1) PositionAngerZ -= (crystalPitch_axial/2)*shiftDis;
	}
	if(i_doi == 1){ //Response shift for 2nd Layer
		if(ireflectorLayer2_Tangential[i_tan] == 1 && ireflectorLayer2_Tangential[i_tan + 1] == 0) PositionAngerY += (crystalPitch_tan/2)*shiftDis;
		if(ireflectorLayer2_Tangential[i_tan] == 0 && ireflectorLayer2_Tangential[i_tan + 1] == 1) PositionAngerY -= (crystalPitch_tan/2)*shiftDis;

		if(ireflectorLayer2_Axial[i_axial] == 1 && ireflectorLayer2_Axial [i_axial + 1] == 0) PositionAngerZ += (crystalPitch_axial/2)*shiftDis;
		if(ireflectorLayer2_Axial[i_axial] == 0 && ireflectorLayer2_Axial [i_axial + 1] == 1) PositionAngerZ -= (crystalPitch_axial/2)*shiftDis;
	}
	if(i_doi == 2){ //Response shift for 3rd Layer
		if(ireflectorLayer3_Tangential[i_tan] == 1 && ireflectorLayer3_Tangential[i_tan + 1] == 0) PositionAngerY += (crystalPitch_tan/2)*shiftDis;
		if(ireflectorLayer3_Tangential[i_tan] == 0 && ireflectorLayer3_Tangential[i_tan + 1] == 1) PositionAngerY -= (crystalPitch_tan/2)*shiftDis;

		if(ireflectorLayer3_Axial[i_axial] == 1 && ireflectorLayer3_Axial [i_axial + 1] == 0) PositionAngerZ += (crystalPitch_axial/2)*shiftDis;
		if(ireflectorLayer3_Axial[i_axial] == 0 && ireflectorLayer3_Axial [i_axial + 1] == 1) PositionAngerZ -= (crystalPitch_axial/2)*shiftDis;
	}
	if(i_doi == 3){ //Response shift for 4th Layer
		if(ireflectorLayer4_Tangential[i_tan] == 1 && ireflectorLayer4_Tangential[i_tan + 1] == 0) PositionAngerY += (crystalPitch_tan/2)*shiftDis;
		if(ireflectorLayer4_Tangential[i_tan] == 0 && ireflectorLayer4_Tangential[i_tan + 1] == 1) PositionAngerY -= (crystalPitch_tan/2)*shiftDis;

		if(ireflectorLayer4_Axial[i_axial] == 1 && ireflectorLayer4_Axial [i_axial + 1] == 0) PositionAngerZ += (crystalPitch_axial/2)*shiftDis;
		if(ireflectorLayer4_Axial[i_axial] == 0 && ireflectorLayer4_Axial [i_axial + 1] == 1) PositionAngerZ -= (crystalPitch_axial/2)*shiftDis;
	}
   
   //Blur the 2D position (obtained by ANger Logic method) to include uncertainity of the PMT position response.
	if(isDOI_LUT){
		PositionAngerZ = G4RandGauss::shoot(PositionAngerZ,PMTblurring_axial/2.35);
		PositionAngerY = G4RandGauss::shoot(PositionAngerY,PMTblurring_tan/2.35);
	}
	//The main purpose of shifting the response is to be able to project the response of all the crytal elements into a 2D position histogram so that we can identify the DOI layer 
	//by comparing with a look-up-table which is prepared based on the reflector insertion. 

	//The crystal ID in 2D position histogram along the axial (z) direction. It can have values of: 0, 1, .. , 31, in 32x32 pixel position histogram
	crystalID_in2D_posHist_axial = (G4int)(PositionAngerZ/(crystalPitch_axial*0.5) + (G4double)numberOfPixel_axial*0.5);//Note! crystalPitch_axial*0.5 is the pitch for the 32x32 2D pixel space, and 0.5 is added for round off

	//The crystal ID in 2D position histogram along the tangential (y) direction. It can have values of: 0, 1, .. , 31, in 32x32 pixel position histogram
	crystalID_in2D_posHist_tan =   (G4int)(PositionAngerY/(crystalPitch_tan*0.5) + (G4double)numberOfPixel_tan * 0.5);
	
	//continuous crystal ID in the 2D position histogram. It will be from 0 to 1023 (in the case of 16x16x4 crystal array). 
	crystalID_in2D_posHist = crystalID_in2D_posHist_axial + crystalID_in2D_posHist_tan * numberOfPixel_tan;//32;


	//Now, lets find the crystal ID in 3D after applying Anger Logic calculation. NOTE that its value can be the same as the original crystal ID or not.

	//Crystal ID along the tangential diraction after Anger Logic calculation
	crystalIDNew_tan = (G4int)(crystalID_in2D_posHist_tan/2);

	//Crystal ID along the axial diraction after Anger Logic calculation
	crystalIDNew_axial = (G4int)(crystalID_in2D_posHist_axial/2);

	////Crystal ID along the DOI diraction after Anger Logic calculation
	if(crystalID_in2D_posHist>numberOfCrystal_tangential*numberOfCrystal_axial*numberOfCrystal_DOI) return;
	crystalIDNew_DOI = doi_table[crystalID_in2D_posHist];

	//If the crsytal ID is beyond the given the number of crystal in the detector, the following is excecuted and the event will be rejected
	if(crystalIDNew_tan < 0 || crystalIDNew_axial < 0 || crystalIDNew_DOI < 0 ||
		crystalIDNew_tan >= numberOfCrystal_tangential || crystalIDNew_axial >= numberOfCrystal_axial || crystalIDNew_DOI >= numberOfCrystal_DOI){
			return;
	}

	CrystalIDAfterAngerLogic(crystalIDNew_tan,crystalIDNew_axial,crystalIDNew_DOI);	
}

/////
void doiPETAnalysis::CrystalIDAfterAngerLogic(G4int i_tan, G4int i_axial, G4int i_doi){
	crystalID_tangential = i_tan;
	crystalID_axial = i_axial;
	DOI_ID = i_doi;
}

void doiPETAnalysis::book() 
{ 
  auto manager = G4AnalysisManager::Instance();
  
  //manager->SetVerboseLevel(2);
 
  G4bool fileOpen = manager->OpenFile(rootFileName);
  if (!fileOpen) {
    G4cout << "\n---> HistoManager::book(): cannot open " 
           << rootFileName << G4endl;
    return;
  }
  // Create directories  
  //manager->SetNtupleDirectoryName("ListModeData");

  manager->SetFirstNtupleId(1);

  if(getSinglesData){
	manager -> CreateNtuple("Singles", "Singles");
	fNtColId[0] = manager -> CreateNtupleIColumn("eventID");
	fNtColId[1] = manager -> CreateNtupleIColumn("blockID");
	//fNtColId[2] = manager -> CreateNtupleDColumn("crystalID_axial");
	//fNtColId[3] = manager -> CreateNtupleDColumn("crystalID_tangential");
	//fNtColId[4] = manager -> CreateNtupleDColumn("DOI_ID");
	fNtColId[2] = manager -> CreateNtupleDColumn("timeStamp");
	fNtColId[3] = manager -> CreateNtupleDColumn("totalEdep");

	//Interaction position of the photon with the detector
	fNtColId[4] = manager -> CreateNtupleDColumn("intPosX");
	fNtColId[5] = manager -> CreateNtupleDColumn("intPosY");
	fNtColId[6] = manager -> CreateNtupleDColumn("intPosZ");

	////source position (annihilation position)
	fNtColId[7] = manager -> CreateNtupleDColumn("spositionX");
	fNtColId[8] = manager -> CreateNtupleDColumn("spositionY");
	fNtColId[9] = manager -> CreateNtupleDColumn("spositionZ");

	manager -> FinishNtuple();	  
  }
  if(getCoincidenceData){
	  manager -> CreateNtuple("Coincidence", "Coincidence");
	fNtColId[0] = manager -> CreateNtupleIColumn("eventID0");
	fNtColId[1] = manager -> CreateNtupleIColumn("blockID0");
	fNtColId[2] = manager -> CreateNtupleIColumn("crystalID_axial0");
	fNtColId[3] = manager -> CreateNtupleIColumn("crystalID_tangential0");
	fNtColId[4] = manager -> CreateNtupleIColumn("DOI_ID0");
	fNtColId[5] = manager -> CreateNtupleDColumn("timeStamp0");
	fNtColId[6] = manager -> CreateNtupleDColumn("totalEdep0");

	fNtColId[7] = manager -> CreateNtupleIColumn("eventID1");
	fNtColId[8] = manager -> CreateNtupleIColumn("blockID1");
	fNtColId[9] = manager -> CreateNtupleIColumn("crystalID_axial1");
	fNtColId[10] = manager -> CreateNtupleIColumn("crystalID_tangential1");
	fNtColId[11] = manager -> CreateNtupleIColumn("DOI_ID1");
	fNtColId[12] = manager -> CreateNtupleDColumn("timeStamp1");
	fNtColId[13] = manager -> CreateNtupleDColumn("totalEdep1");

	//source position
	fNtColId[14] = manager -> CreateNtupleDColumn("spositionX");
	fNtColId[15] = manager -> CreateNtupleDColumn("spositionY");
	fNtColId[16] = manager -> CreateNtupleDColumn("spositionZ");

	manager -> FinishNtuple();

  }
  
  
  factoryOn = true;    
}
void doiPETAnalysis::FillListModeEvent()
{

  auto manager = G4AnalysisManager::Instance();
  if(getSinglesData){
	 	manager -> FillNtupleIColumn(1, fNtColId[0], G4int(eventID));
		manager -> FillNtupleIColumn(1, fNtColId[1], G4int(blockID));
		//manager -> FillNtupleDColumn(1, fNtColId[2], crystalID_axial);
		//manager -> FillNtupleDColumn(1, fNtColId[3], crystalID_tangential);
		//manager -> FillNtupleDColumn(1, fNtColId[4], DOI_ID);
		manager -> FillNtupleDColumn(1, fNtColId[2], timeStamp/s);// in second
		manager -> FillNtupleDColumn(1, fNtColId[3], totalEdep/keV); //in keV
		
		//Interaction position of the photon in the detector
		manager -> FillNtupleDColumn(1, fNtColId[4], intPosX); //mm
		manager -> FillNtupleDColumn(1, fNtColId[5], intPosY); //mm
		manager -> FillNtupleDColumn(1, fNtColId[6], intPosZ); //mm

		//
		//Add source position
		manager -> FillNtupleDColumn(1, fNtColId[7], spositionX);
		manager -> FillNtupleDColumn(1, fNtColId[8], spositionY);
		manager -> FillNtupleDColumn(1, fNtColId[9], spositionZ);

		manager -> AddNtupleRow(1);
  }

  if(getCoincidenceData){
	  //First Single
		manager -> FillNtupleIColumn(1, fNtColId[0], eventID0);
		manager -> FillNtupleIColumn(1, fNtColId[1], blockID0);
		manager -> FillNtupleIColumn(1, fNtColId[2], crystalID_axial0);
		manager -> FillNtupleIColumn(1, fNtColId[3], crystalID_tangential0);
		manager -> FillNtupleIColumn(1, fNtColId[4], DOI_ID0);
		manager -> FillNtupleDColumn(1, fNtColId[5], timeStamp0/s);
		manager -> FillNtupleDColumn(1, fNtColId[6], totalEdep0/keV);
	
	//Second Single
		manager -> FillNtupleIColumn(1, fNtColId[7], eventID1);
		manager -> FillNtupleIColumn(1, fNtColId[8], blockID1);
		manager -> FillNtupleIColumn(1, fNtColId[9], crystalID_axial1);
		manager -> FillNtupleIColumn(1, fNtColId[10], crystalID_tangential1);
		manager -> FillNtupleIColumn(1, fNtColId[11], DOI_ID1);
		manager -> FillNtupleDColumn(1, fNtColId[12], timeStamp1/s);
		manager -> FillNtupleDColumn(1, fNtColId[13], totalEdep1/keV);
	
	//Add source position
		manager -> FillNtupleDColumn(1, fNtColId[14], spositionX);
		manager -> FillNtupleDColumn(1, fNtColId[15], spositionY);
		manager -> FillNtupleDColumn(1, fNtColId[16], spositionZ);

		manager -> AddNtupleRow(1);
  }
    
}
void doiPETAnalysis::finish() 
{   
 if (factoryOn) 
   {
    auto manager = G4AnalysisManager::Instance();    
    manager -> Write();
    manager -> CloseFile();  
      
   // delete G4AnalysisManager::Instance();
    factoryOn = false;
   }
}

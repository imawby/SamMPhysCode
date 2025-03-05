#define anacode_cxx
#include "anacode.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void anacode::Loop()
{
//   In a ROOT session, you can do:
//      root> .L anacode.C
//      root> anacode t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch



  //////////////////////////////////////////////////////////////
  // Define signal/background trees for BDT
  //////////////////////////////////////////////////////////////
  // Create file to put our trees in
  TFile *file = new TFile("BDTTrees.root", "RECREATE");
  
  // Construct trees
  TTree *signalTree = new TTree("signalTree", "signalTree");
  TTree *backgroundTree = new TTree("backgroundTree", "backgroundTree");

  // Define variables to add to tree
  float exampleVariable = 999.f;
  float closestSepToParent = 999.f;
  int nHits = 0;
  float trackShowerScore = 999.f;
  int generation = 0;

  // Set the branches for the tree
  for (TTree *tree : {signalTree, backgroundTree})
  {
    tree->Branch("ExampleVariable", &exampleVariable, "ExampleVariable/F");
    tree->Branch("ClosestSepToParent", &closestSepToParent, "ClosestSepToParent/F");
    tree->Branch("nHits", &nHits, "nHits/I");
    tree->Branch("trackShowerScore", &trackShowerScore, "trackShowerScore/F");
    tree->Branch("generation", &generation, "generation/I");
  }
  
  if (fChain == 0) return;
   
  Long64_t nentries = fChain->GetEntriesFast();
   
  Long64_t nbytes = 0, nb = 0;
 
  //for (Long64_t jentry=0; jentry<nentries;jentry++)
   
  for (Long64_t jentry = 0; jentry < 1000; jentry++)   
    {
      Long64_t ientry = LoadTree(jentry);
      
      if (ientry < 0)
	break;
     
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
          
      if (Cut(ientry) < 0) continue;

      for (int recoIndex = 0; recoIndex < Event_NRecoPFPs; recoIndex++) //testing ground for functions
	{
	  if (recoIndex == 0)
	    {
	      continue;
	    }
	  std::cout << "recoIndex: " << recoIndex << std::endl;
	  std::vector<vector<double>> withinDist = findPointsWithinDist(recoIndex, 5);

	  std::cout << "PFP_SpacepointX->at(recoIndex).size(): " << PFP_SpacepointX->at(recoIndex).size() << ", withinDist[0].size(): " << withinDist[0].size() << std::endl;
	  
	  TVector3 PCA_full = findBestFitPCA(recoIndex, PFP_SpacepointX->at(recoIndex), PFP_SpacepointY->at(recoIndex), PFP_SpacepointZ->at(recoIndex));
	  TVector3 PCA_dist = findBestFitPCA(recoIndex, withinDist[0], withinDist[1], withinDist[2]);
	  TVector3 PtoC = GetPtoC(Track_EndX[recoIndex-1], Track_EndY[recoIndex-1], Track_EndZ[recoIndex-1], Track_StartX[recoIndex], Track_StartY[recoIndex], Track_StartZ[recoIndex]);
	  TVector3 TVectorChild = startEndChildVector(Track_StartX[recoIndex], Track_StartY[recoIndex], Track_StartZ[recoIndex], Track_EndX[recoIndex], Track_EndY[recoIndex], Track_EndZ[recoIndex]);
	  
	  double DCA_startEnd = DCA(TVectorChild, PtoC);
	  double DCA_PCAfull = DCA(PCA_full, PtoC);
	  double DCA_PCAdist = DCA(PCA_dist, PtoC);
	  double exDist_startEnd = extrapDist(TVectorChild, PtoC);
	  double exDist_PCAfull = extrapDist(PCA_full, PtoC);
	  double exDist_PCAdist = extrapDist(PCA_dist, PtoC);
	  
	  std::cout << "DCA startEnd: "  << DCA_startEnd << ", DCA PCA full: " << DCA_PCAfull << ", DCA PCA dist: " << DCA_PCAdist << std::endl;
	  std::cout << "extrapolated dist startEnd: " << exDist_startEnd << ", extrap dist PCA full: " << exDist_PCAfull << ", extrap dist PCA dist: " << exDist_PCAdist << std::endl;
	}
      
      /*
      std::cout << "-----------------------------" << std::endl;
      std::cout << "-----------------------------" << std::endl;
      std::cout << "jEntry: " << jentry << std::endl;
      std::cout << "-----------------------------" << std::endl;
      std::cout << "-----------------------------" << std::endl;     
      */
      
      /*
      // run the map filling function
      std::map<int, int> parentMap = mapParents(); //child->parent
      //PrintChildToParentRecoIDMap(parentMap);

      // For each entry in the tree (this should be a parent/child link) fill the tree
      bool isSignal = true;


      // Decide if background or signal
      for (int recoIndex_parent = 0; recoIndex_parent < Event_NRecoPFPs; recoIndex_parent++)
	{
	  std::cout << "recoID: " << PFP_RecoID[recoIndex_parent] << std::endl;
	  int PDG_parent = findPDG(PFP_RecoID[recoIndex_parent]);
	  std::cout << "parent PDG: " << PDG_parent << std::endl;
	  if (PDG_parent == 13 or PDG_parent == 211 or PDG_parent == -13 or PDG_parent == -211)
	    {
	      for (int recoIndex_child = 0; recoIndex_child < Event_NRecoPFPs; recoIndex_child++)
		{
		  if (parentMap[PFP_RecoID[recoIndex_child]] == PFP_RecoID[recoIndex_parent])
		    {
		      closestSepToParent = findClosestSep(PFP_RecoID[recoIndex_parent], PFP_RecoID[recoIndex_child]);
		      nHits = PFP_SpacepointX[recoIndex_child].size();
		      trackShowerScore = PFP_TrackShowerScore[recoIndex_child];
		      generation = findGeneration(PFP_RecoID[recoIndex_child], parentMap, 1);
		      
		      
		      signalTree->Fill();         //signal
		    }
		  else
		    {
		      if (PFP_RecoID[recoIndex_parent] != PFP_RecoID[recoIndex_child])
			{
			  closestSepToParent = findClosestSep(PFP_RecoID[recoIndex_parent], PFP_RecoID[recoIndex_child]);
			  nHits = PFP_SpacepointX[recoIndex_child].size();
			  trackShowerScore = PFP_TrackShowerScore[recoIndex_child];
			  
			  backgroundTree->Fill();     //background
			}
		    }
		}
	    }
	  else
	    {
	      continue;
	    }
	}
      
    }

  // Write trees to file and close file
  signalTree->Write();
  backgroundTree->Write();
  file->Close();
  */
      /////
      
    }
  ////
}

     

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::map<int, int> anacode::mapParents()                  //return recoID -> parentRecoID for the next reconstructed parent
{    
  std::map<int, std::vector<int>> simRecoMap;             //init. map simID -> recoID's
  std::vector<int> splitList;                             //init. list for split
  std::vector<int> deltaRayList;                          //init. delta ay list

  //make initial maps and lists from data  
  for (int i = 0; i < Event_NRecoPFPs; i++)                                                        
    {
      //use matchedMC_SimID to put in ALL the simID's of reconstructed particles
      if (IsDeltaRay(i))
	{
          deltaRayList.emplace_back(PFP_RecoID[i]);
	}
      else
	{
          simRecoMap[MatchedMC_SimID[i]].push_back(PFP_RecoID[i]);

	  if(simRecoMap[MatchedMC_SimID[i]].size() == 2)
	    splitList.emplace_back(MatchedMC_SimID[i]);      //if second entry under same simID, add entry to split list, if third etc then ignore	  
	}      
    }

  simRecoMap[0].emplace_back(-1);                          //add in neutrino simID->recoID

  std::map<int, int> childToParentRecoIDMap;
  std::map<int, int> childToParentSimIDMap;

  // Fill our childToParentRecoIDMap
  for (int i = 0; i < Event_NRecoPFPs; i++)                                                        //search through every recon particle
    {

      /////////////////////////////////////////////////////////////
      // First find the sim ID of the first reconstructed parent
      /////////////////////////////////////////////////////////////
      int recoParentSimID = 0; // default as a true primary
      int indexofCur = findSimArrayIndex(MatchedMC_SimID[i]);
      bool found = (indexofCur != -1);
      
      if (found)
	{	
	  if (IsDeltaRay(i))
	    {	      
	      if (simRecoMap.find(MC_SimID[indexofCur]) == simRecoMap.end())
		{
		  recoParentSimID = 0;
		}
	      else
		{
		  recoParentSimID = MC_SimID[indexofCur];
		}
	    }
	  else
	    {
              recoParentSimID = nextRecoParent(simRecoMap, MC_ParentSimID[indexofCur], MC_SimID[indexofCur]);
	    }
	}

      childToParentSimIDMap[PFP_RecoID[i]] = recoParentSimID; 

      /////////////////////////////////////////////////////////////
      // Now find the reco ID of the first reconstructed parent
      /////////////////////////////////////////////////////////////
      
      // first check if parent particle has been split in the recon
      bool isSplit = (std::find(splitList.begin(), splitList.end(), recoParentSimID) != splitList.end());

      // Now find the appropriate reco parent
      int recoParentRecoID = -1; // initialise variable first!
      
      if (!isSplit)
      {
          recoParentRecoID = simRecoMap.at(recoParentSimID).at(0);
      }
      else
      {
	  int closestParentRecoID = -1; // initialise variable first!
	  double closestSep = std::numeric_limits<double>::max();
          int nSpacePoints_child = PFP_SpacepointX->at(i).size();	  

          for (int coordIndex_child = 0; coordIndex_child < nSpacePoints_child; coordIndex_child++)                  //for each daughter spacepoint set, one loop
          {
	      double xPos_child = PFP_SpacepointX->at(i).at(coordIndex_child);                                       //coords of daughter particle
	      double yPos_child = PFP_SpacepointY->at(i).at(coordIndex_child);
	      double zPos_child = PFP_SpacepointZ->at(i).at(coordIndex_child);
		  
	      for (int parentRecoID : simRecoMap[recoParentSimID])
	      {
  		  int parentArrayIndex = findRecoArrayIndex(parentRecoID);

		  // if we couldn't find index (should be impossible)
		  if (parentArrayIndex == -1)
		    continue;
		  
		  int nSpacePoints_parent = PFP_SpacepointX->at(parentArrayIndex).size();
		
		  for (int coordIndex_parent = 0; coordIndex_parent < nSpacePoints_parent; coordIndex_parent++)                //for each parent spacepoint set, loop
		  {			  
		      double xPos_parent = PFP_SpacepointX->at(parentArrayIndex).at(coordIndex_parent);                           //coords of parent particle k (0-however many splits)
		      double yPos_parent = PFP_SpacepointY->at(parentArrayIndex).at(coordIndex_parent);
		      double zPos_parent = PFP_SpacepointZ->at(parentArrayIndex).at(coordIndex_parent);

		      double xDiff = xPos_child - xPos_parent;
		      double yDiff = yPos_child - yPos_parent;
		      double zDiff = zPos_child - zPos_parent;		    
			  
		      double sep = sqrt(pow(xDiff,2) + pow(yDiff,2) + pow(zDiff,2));

		      if (sep < closestSep)
		      {
			  closestSep = sep;
			  closestParentRecoID = parentRecoID;
		      }
		  }
	      }
	  }

	  // incase we can't find a closest parent
	  if (closestParentRecoID == -1)
	  {
  	      std::cout << "WARNING: couldn't find closest parent!" << std::endl;
	      recoParentRecoID = -1;
	  }
	  else
	  {
	      recoParentRecoID = closestParentRecoID;
	  }
      }
      childToParentRecoIDMap[PFP_RecoID[i]] = recoParentRecoID;                             //add to the map of all particles in the event's parents
  }
  return childToParentRecoIDMap;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int anacode::findSimArrayIndex(int simID)
{
  
    
  int index = -1; //initialise first!

  if (simID == -1)
    {
      return index; 
    }
  
  for (int i = 0; i < Event_NMCParticles; i++)
    {
        if (MC_SimID[i] == simID)
	  {
	    index = i;
	    break;
	  }
    }
  
    if (index == -1)
      {
	//std::cout << simID <<std::endl;
	//std::cout << Event_NMCParticles << std::endl; 
	std::cout << "WARNING: couldn't sim array index!" << std::endl;
      }
    
    return index;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int anacode::findRecoArrayIndex(int recoID)
{
    // Get array index of parent particle 
    int index = -1; //initialise first!

      if (recoID == -1)
    {
      return index; 
    }
		  
    for (int i = 0; i < Event_NRecoPFPs; i++)
    {
        if (PFP_RecoID[i] == recoID)
	{
	    index = i;
	    break;
	}
    }

    if (index == -1)
      std::cout << "WARNING: couldn't reco array index!" << std::endl;
    
    return index;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int anacode::nextRecoParent(std::map<int, std::vector<int>> simRecoMap, int trueParentSimID, int currentParticleSimID)
{
  if (simRecoMap.find(trueParentSimID) == simRecoMap.end()) //if empty at SimID in simRecoMap, go next
  {
      int index = findSimArrayIndex(trueParentSimID);
      bool found = (index != -1);

      // If we couldn't find the index then return the neutrino as the parent
      if (!found)
	return 0;

      //index = parentSimID's index in MC_SimID array 
      return nextRecoParent(simRecoMap, MC_ParentSimID[index], trueParentSimID);
  }
  else
  {
    return trueParentSimID;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool anacode::IsDeltaRay(int index)
{
  if (MatchedMC_IsSaved[index] == 1)
    return false;

  if (abs(MatchedMC_PDG[index]) != 13)
    return false;

  return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void anacode::PrintChildToParentSimIDMap(std::map<int,int> childToParentSimIDMap)
{
  std::cout << "///////////////////////////" << std::endl;
  std::cout << "// childToParentSimIDMap //" << std::endl;
  std::cout << "///////////////////////////" << std::endl;  
  
  for (auto entry : childToParentSimIDMap)
    std::cout << "Child RecoID: " << entry.first << ", Parent SimID: " << entry.second << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void anacode::PrintChildToParentRecoIDMap(std::map<int,int> childToParentRecoIDMap)
{
  std::cout << "////////////////////////////" << std::endl;
  std::cout << "// childToParentRecoIDMap //" << std::endl;
  std::cout << "////////////////////////////" << std::endl;
  
  for (auto entry : childToParentRecoIDMap)
    std::cout << "Child RecoID: " << entry.first << ", Parent RecoID: " << entry.second << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool anacode::isTrueParent(int parentSimID, int childSimID)
{
  if (MC_ParentSimID[findSimArrayIndex(childSimID)] == parentSimID)
    {
      return true;
    }
  else
    {
      return false;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double anacode::findClosestSep(int particleRecoIndex_1, int particleRecoIndex_2)
{
  double closestSep = std::numeric_limits<double>::max();
  int nSpacePoints_1 = PFP_SpacepointX->at(particleRecoIndex_1).size();
  int nSpacePoints_2 = PFP_SpacepointX->at(particleRecoIndex_2).size();
	      
  for (int coordIndex_1 = 0; coordIndex_1 < nSpacePoints_1; coordIndex_1++)
    {
      double xPos_1 = PFP_SpacepointX->at(particleRecoIndex_1).at(coordIndex_1);                                  
      double yPos_1 = PFP_SpacepointY->at(particleRecoIndex_1).at(coordIndex_1);
      double zPos_1 = PFP_SpacepointZ->at(particleRecoIndex_1).at(coordIndex_1);
	  
      int nSpacePoints_2 = PFP_SpacepointX->at(particleRecoIndex_2).size();
		
      for (int coordIndex_2 = 0; coordIndex_2 < nSpacePoints_2; coordIndex_2++)              
	{			  
	  double xPos_2 = PFP_SpacepointX->at(particleRecoIndex_2).at(coordIndex_2);                         
	  double yPos_2 = PFP_SpacepointY->at(particleRecoIndex_2).at(coordIndex_2);
	  double zPos_2 = PFP_SpacepointZ->at(particleRecoIndex_2).at(coordIndex_2);

	  double xDiff = xPos_1 - xPos_2;
	  double yDiff = yPos_1 - yPos_2;
	  double zDiff = zPos_1 - zPos_2;		    
			  
	  double sep = sqrt(pow(xDiff,2) + pow(yDiff,2) + pow(zDiff,2));

	  if (sep < closestSep)
	    {
	      closestSep = sep;
	    }
	}
    }

  return closestSep;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int anacode::findSimID(int recoID)
{
  int SimID = -1;
  if (recoID == -1)
    {
      return SimID;
    }
  else
    {
    return MatchedMC_SimID[findRecoArrayIndex(recoID)];
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int anacode::findPDG(int recoIndex)
{
  int PDG = -1;
  if (recoIndex == -1)
    {
      return PDG;
    }
  else
    {
      return MatchedMC_PDG[recoIndex];
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TVector3 anacode::findBestFitPCA(int recoParticleIndex, vector<double> SpacepointsX, vector<double> SpacepointsY, vector<double> SpacepointsZ)
{
  std::vector<double> eVals;
  std::vector<TVector3> eVecs;
  
  TPrincipal* principal = new TPrincipal(3, "D");
  
  const int nSpacepoints = SpacepointsX.size();

  for (int i=0; i < nSpacepoints; i++)
    {
      double coord[3] = {SpacepointsX.at(i), SpacepointsY.at(i), SpacepointsZ.at(i)};
      principal->AddRow(coord);
    }
	
    // PERFORM PCA
    principal->MakePrincipals();
    // GET EIGENVALUES AND EIGENVECTORS
    for (unsigned int i = 0; i < 3; ++i)
         eVals.push_back(principal->GetEigenValues()->GetMatrixArray()[i]);

    for (int i : {0, 3, 6})
    {
        const double eVec_x = principal->GetEigenVectors()->GetMatrixArray()[i];
        const double eVec_y = principal->GetEigenVectors()->GetMatrixArray()[i + 1];
        const double eVec_z = principal->GetEigenVectors()->GetMatrixArray()[i + 2];

        eVecs.push_back(TVector3(eVec_x, eVec_y, eVec_z));
    }
  
  // Make sure the PCA fit is sensible
  std::cout << "eVals.size : " << eVals.size() << std::endl;
  for (int i = 0; i < eVals.size(); i++)
    {
      if (std::isnan(eVals[i]))
	{
	  std::cout << "ERROR: eVal not found" << std::endl;
	  eVals[i] = -999;
	}
    }
  
  // Find largest eVal, and return its eVector (may be able to do it in the loop above but more clear to seperate)
  int eValMax_index = 0;
  double eValMax = 0;
  for (int i = 0; i < eVals.size(); i++)
    {
      if (eVals[i] < eValMax)
	{
	  continue;
	}
      eValMax = eVals[i];
      eValMax_index = i;
    }
  
  return eVecs[eValMax_index];
     
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TVector3 anacode::findMeanPoint(int recoIndex)
{
  const int nSpacepoints = PFP_SpacepointX->at(recoIndex).size();
  double xAdd = 0;
  double yAdd = 0;
  double zAdd = 0;
  TVector3 meanpoint;
  
  for (int i = 0; i < nSpacepoints; i++) //for set of points
    {
      xAdd += PFP_SpacepointX->at(recoIndex)[i];
      yAdd += PFP_SpacepointY->at(recoIndex)[i];
      zAdd += PFP_SpacepointZ->at(recoIndex)[i];
    }
  
  meanpoint.SetX(xAdd/nSpacepoints);
  meanpoint.SetY(yAdd/nSpacepoints);
  meanpoint.SetZ(zAdd/nSpacepoints);
  
  return meanpoint;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TVector3 anacode::GetPtoC(double ParentEndX, double ParentEndY, double ParentEndZ, double ChildStartX, double ChildStartY, double ChildStartZ)
{
  double PtoC_X = ChildStartX - ParentEndX;
  double PtoC_Y = ChildStartY - ParentEndY;
  double PtoC_Z = ChildStartZ - ParentEndZ;
  TVector3 PtoC;
  PtoC.SetX(PtoC_X);
  PtoC.SetY(PtoC_Y);
  PtoC.SetZ(PtoC_Z);
  return PtoC;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double anacode::DCA(TVector3 TVectorChild, TVector3 TVectorPtoC)
{
  TVector3 DCA_Vec;
  
  DCA_Vec = TVectorPtoC.Cross(TVectorChild);
  return DCA_Vec.Mag();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double anacode::extrapDist(TVector3 TVectorChild, TVector3 TVectorPtoC)
{
  double ExtrapDist = TVectorPtoC.Dot(TVectorChild);
  return abs(ExtrapDist);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double anacode::getAngle(TVector3 TVectorChild, TVector3 TVectorPtoC)
{
  double a = TVectorPtoC.Angle(TVectorChild);
  return a;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TVector3 anacode::startEndChildVector(double childStartX, double childStartY, double childStartZ, double childEndX, double childEndY, double childEndZ)
{
  double normaliser_X = abs(childStartX - childEndX);
  double normaliser_Y = abs(childStartY - childEndY);
  double normaliser_Z = abs(childStartZ - childEndZ);
  
  TVector3 ChildVector;
  ChildVector.SetX((childStartX - childEndX) / normaliser_X);
  ChildVector.SetY((childStartY - childEndY) / normaliser_Y);
  ChildVector.SetZ((childStartZ - childEndZ) / normaliser_Z);
  
  return ChildVector;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::vector<std::vector<double>> anacode::findPointsWithinDist(int recoIndex, double dist = 5)
{
  std::vector<double> withinDist_X;
  std::vector<double> withinDist_Y;
  std::vector<double> withinDist_Z;
  const int nSpacepoints = PFP_SpacepointX->at(recoIndex).size();
  
  for (int i = 0; i < nSpacepoints; i++)
    {
      double xDiff = PFP_SpacepointX->at(recoIndex)[i] - Track_StartX[recoIndex];
      double yDiff = PFP_SpacepointY->at(recoIndex)[i] - Track_StartY[recoIndex];
      double zDiff = PFP_SpacepointZ->at(recoIndex)[i] - Track_StartZ[recoIndex];
      double distFromTrackStart = sqrt(pow(xDiff,2) + pow(yDiff,2) + pow(zDiff,2));
      if (distFromTrackStart > dist)
	{
	  continue;
	}
      withinDist_X.emplace_back(PFP_SpacepointX->at(recoIndex)[i]);
      withinDist_Y.emplace_back(PFP_SpacepointY->at(recoIndex)[i]);
      withinDist_Z.emplace_back(PFP_SpacepointZ->at(recoIndex)[i]);
    }

  std::vector<vector<double>> withinSpacepoints;
  withinSpacepoints.emplace_back(withinDist_X);
  withinSpacepoints.emplace_back(withinDist_Y);
  withinSpacepoints.emplace_back(withinDist_Z);

  return withinSpacepoints;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

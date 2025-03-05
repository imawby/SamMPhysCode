//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 29 15:48:59 2024 by ROOT version 6.32.04
// from TTree ccnusel/CC nu selection
// found on file: ccnutree_SAM_0_v2.root
//////////////////////////////////////////////////////////

#ifndef anacode_h
#define anacode_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class anacode {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Event_Run;
   Int_t           Event_Subrun;
   Int_t           Event_Event;
   Int_t           Nu_True_PDG;
   Int_t           Nu_True_SimID;
   Int_t           Nu_True_IsNC;
   Int_t           Nu_True_Mode;
   Int_t           Nu_True_TargetZ;
   Double_t        Nu_True_Energy;
   Double_t        Nu_True_MomX;
   Double_t        Nu_True_MomY;
   Double_t        Nu_True_MomZ;
   Double_t        Nu_True_VertexX;
   Double_t        Nu_True_VertexY;
   Double_t        Nu_True_VertexZ;
   Int_t           Event_NMCParticles;
   Int_t           Event_NRecoPFPs;
   Double_t        Nu_Reco_VertexX;
   Double_t        Nu_Reco_VertexY;
   Double_t        Nu_Reco_VertexZ;
   Int_t           MC_SimID[1000];   //[Event_NMCParticles]
   Int_t           MC_ParentSimID[1000];   //[Event_NMCParticles]
   Int_t           PFP_RecoID[46];   //[Event_NRecoPFPs]
   Int_t           MatchedMC_PDG[46];   //[Event_NRecoPFPs]
   Int_t           MatchedMC_SimID[46];   //[Event_NRecoPFPs]
   Bool_t          MatchedMC_IsSaved[46];   //[Event_NRecoPFPs]
   Bool_t          MatchedMC_IsTruePrimary[46];   //[Event_NRecoPFPs]
   Int_t           MatchedMC_Generation[46];   //[Event_NRecoPFPs]
   Int_t           MatchedMC_ParentSimID[46];   //[Event_NRecoPFPs]
   Int_t           MatchedMC_ParentPDG[46];   //[Event_NRecoPFPs]
   Bool_t          PFP_CurrentReco_IsPrimary[46];   //[Event_NRecoPFPs]
   Int_t           PFP_CurrentReco_Generation[46];   //[Event_NRecoPFPs]
   Int_t           PFP_CurrentReco_ParentRecoID[46];   //[Event_NRecoPFPs]
   Double_t        PFP_TrackShowerScore[46];   //[Event_NRecoPFPs]
   Double_t        PFP_VertexX[46];   //[Event_NRecoPFPs]
   Double_t        PFP_VertexY[46];   //[Event_NRecoPFPs]
   Double_t        PFP_VertexZ[46];   //[Event_NRecoPFPs]
   Int_t           PFP_N2DHits[46];   //[Event_NRecoPFPs]
   Int_t           PFP_NSpacepoints[46];   //[Event_NRecoPFPs]
   vector<vector<double> > *PFP_SpacepointX;
   vector<vector<double> > *PFP_SpacepointY;
   vector<vector<double> > *PFP_SpacepointZ;
   Double_t        PFP_Completeness[46];   //[Event_NRecoPFPs]
   Double_t        PFP_Purity[46];   //[Event_NRecoPFPs]
   Bool_t          PFP_ShowerFitSuccess[46];   //[Event_NRecoPFPs]
   Bool_t          PFP_TrackFitSuccess[46];   //[Event_NRecoPFPs]
   Double_t        Track_StartX[46];   //[Event_NRecoPFPs]
   Double_t        Track_StartY[46];   //[Event_NRecoPFPs]
   Double_t        Track_StartZ[46];   //[Event_NRecoPFPs]
   Double_t        Track_EndX[46];   //[Event_NRecoPFPs]
   Double_t        Track_EndY[46];   //[Event_NRecoPFPs]
   Double_t        Track_EndZ[46];   //[Event_NRecoPFPs]
   Double_t        Track_Length[46];   //[Event_NRecoPFPs]
   Double_t        Track_DeflecAngleSD[46];   //[Event_NRecoPFPs]
   Double_t        Track_EvalRatio[46];   //[Event_NRecoPFPs]
   Double_t        Track_Concentration[46];   //[Event_NRecoPFPs]
   Double_t        Track_CoreHaloRatio[46];   //[Event_NRecoPFPs]
   Double_t        Track_Conicalness[46];   //[Event_NRecoPFPs]
   Double_t        Track_dEdxStart[46];   //[Event_NRecoPFPs]
   Double_t        Track_dEdxEnd[46];   //[Event_NRecoPFPs]
   Double_t        Track_dEdxEndRatio[46];   //[Event_NRecoPFPs]
   Int_t           Track_NTrajPoints;
   Int_t           Track_TrajPointStartIndex[46];   //[Event_NRecoPFPs]
   Int_t           Track_TrajPointEndIndex[46];   //[Event_NRecoPFPs]
   Double_t        Track_dEdX[1000];   //[Track_NTrajPoints]
   Double_t        Track_RR[1000];   //[Track_NTrajPoints]
   Double_t        Shower_StartX[46];   //[Event_NRecoPFPs]
   Double_t        Shower_StartY[46];   //[Event_NRecoPFPs]
   Double_t        Shower_StartZ[46];   //[Event_NRecoPFPs]
   Double_t        Shower_Length[46];   //[Event_NRecoPFPs]
   Double_t        Shower_OpeningAngle[46];   //[Event_NRecoPFPs]
   Double_t        Shower_InitialdEdx[46][3];   //[Event_NRecoPFPs]
   Int_t           Shower_BestPlane[46];   //[Event_NRecoPFPs]
   Double_t        Shower_Energy[46][3];   //[Event_NRecoPFPs]

   // List of branches
   TBranch        *b_Event_Run;   //!
   TBranch        *b_Event_Subrun;   //!
   TBranch        *b_Event_Event;   //!
   TBranch        *b_Nu_True_PDG;   //!
   TBranch        *b_Nu_True_SimID;   //!
   TBranch        *b_Nu_True_IsNC;   //!
   TBranch        *b_Nu_True_Mode;   //!
   TBranch        *b_Nu_True_TargetZ;   //!
   TBranch        *b_Nu_True_Energy;   //!
   TBranch        *b_Nu_True_MomX;   //!
   TBranch        *b_Nu_True_MomY;   //!
   TBranch        *b_Nu_True_MomZ;   //!
   TBranch        *b_Nu_True_VertexX;   //!
   TBranch        *b_Nu_True_VertexY;   //!
   TBranch        *b_Nu_True_VertexZ;   //!
   TBranch        *b_Event_NMCParticles;   //!
   TBranch        *b_Event_NRecoPFPs;   //!
   TBranch        *b_Nu_Reco_VertexX;   //!
   TBranch        *b_Nu_Reco_VertexY;   //!
   TBranch        *b_Nu_Reco_VertexZ;   //!
   TBranch        *b_MC_SimID;   //!
   TBranch        *b_MC_ParentSimID;   //!
   TBranch        *b_PFP_RecoID;   //!
   TBranch        *b_MatchedMC_PDG;   //!
   TBranch        *b_MatchedMC_SimID;   //!
   TBranch        *b_MatchedMC_IsSaved;   //!
   TBranch        *b_MatchedMC_IsTruePrimary;   //!
   TBranch        *b_MatchedMC_Generation;   //!
   TBranch        *b_MatchedMC_ParentSimID;   //!
   TBranch        *b_MatchedMC_ParentPDG;   //!
   TBranch        *b_PFP_CurrentReco_IsPrimary;   //!
   TBranch        *b_PFP_CurrentReco_Generation;   //!
   TBranch        *b_PFP_CurrentReco_ParentRecoID;   //!
   TBranch        *b_PFP_TrackShowerScore;   //!
   TBranch        *b_PFP_VertexX;   //!
   TBranch        *b_PFP_VertexY;   //!
   TBranch        *b_PFP_VertexZ;   //!
   TBranch        *b_PFP_N2DHits;   //!
   TBranch        *b_PFP_NSpacepoints;   //!
   TBranch        *b_PFP_SpacepointX;   //!
   TBranch        *b_PFP_SpacepointY;   //!
   TBranch        *b_PFP_SpacepointZ;   //!
   TBranch        *b_PFP_Completeness;   //!
   TBranch        *b_PFP_Purity;   //!
   TBranch        *b_PFP_ShowerFitSuccess;   //!
   TBranch        *b_PFP_TrackFitSuccess;   //!
   TBranch        *b_Track_StartX;   //!
   TBranch        *b_Track_StartY;   //!
   TBranch        *b_Track_StartZ;   //!
   TBranch        *b_Track_EndX;   //!
   TBranch        *b_Track_EndY;   //!
   TBranch        *b_Track_EndZ;   //!
   TBranch        *b_Track_Length;   //!
   TBranch        *b_Track_DeflecAngleSD;   //!
   TBranch        *b_Track_EvalRatio;   //!
   TBranch        *b_Track_Concentration;   //!
   TBranch        *b_Track_CoreHaloRatio;   //!
   TBranch        *b_Track_Conicalness;   //!
   TBranch        *b_Track_dEdxStart;   //!
   TBranch        *b_Track_dEdxEnd;   //!
   TBranch        *b_Track_dEdxEndRatio;   //!
   TBranch        *b_Track_NTrajPoints;   //!
   TBranch        *b_Track_TrajPointStartIndex;   //!
   TBranch        *b_Track_TrajPointEndIndex;   //!
   TBranch        *b_Track_dEdX;   //!
   TBranch        *b_Track_RR;   //!
   TBranch        *b_Shower_StartX;   //!
   TBranch        *b_Shower_StartY;   //!
   TBranch        *b_Shower_StartZ;   //!
   TBranch        *b_Shower_Length;   //!
   TBranch        *b_Shower_OpeningAngle;   //!
   TBranch        *b_Shower_InitialdEdx;   //!
   TBranch        *b_Shower_BestPlane;   //!
   TBranch        *b_Shower_Energy;   //!

  anacode(TTree *tree=0);
  virtual ~anacode();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual std::map<int, int>     mapParents();
  virtual int      findSimArrayIndex(int simID);
  virtual int      findRecoArrayIndex(int recoID);
  virtual int      nextRecoParent(std::map<int, std::vector<int>> simRecoMap, int trueParentSimID, int currentParticleSimID);
  virtual bool     IsDeltaRay(int index);
  virtual void     PrintChildToParentSimIDMap(std::map<int,int> childToParentSimIDMap);
  virtual void     PrintChildToParentRecoIDMap(std::map<int,int> childToParentRecoIDMap);
  virtual bool     isTrueParent(int parentSimID, int childSimID);
  virtual double   findClosestSep(int particleRecoIndex_1, int particleRecoIndex_2);
  virtual int      findSimID(int recoID);
  virtual int      findPDG(int recoIndex);
  virtual TVector3 findBestFitPCA(int recoParticleIndex, vector<double> SpacepointsX, vector<double> SpacepointsY, vector<double> SpacepointsZ);
  virtual TVector3 findMeanPoint(int recoIndex);
  virtual double   DCA(TVector3 TVectorChild, TVector3 TVectorPtoC);
  virtual double   extrapDist(TVector3 TVectorChild, TVector3 TVectorPtoC);
  virtual double   getAngle(TVector3 TVectorChild, TVector3 TVectorPtoC);
  virtual TVector3 GetPtoC(double ParentEndX, double ParentEndY, double ParentEndZ, double ChildStartX, double ChildStartY, double ChildStartZ);
  virtual TVector3 startEndChildVector(double childStartX, double childStartY, double childStartZ, double childEndX, double childEndY, double childEndZ);
  virtual std::vector<std::vector<double>> findPointsWithinDist(int recoIndex, double dist = 5);
  virtual bool     Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef anacode_cxx
anacode::anacode(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ccnutree_SAM_0_v2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ccnutree_SAM_0_v2.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ccnutree_SAM_0_v2.root:/ccnuselection");
      dir->GetObject("ccnusel",tree);

   }
   Init(tree);
}

anacode::~anacode()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t anacode::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t anacode::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void anacode::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PFP_SpacepointX = 0;
   PFP_SpacepointY = 0;
   PFP_SpacepointZ = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event_Run", &Event_Run, &b_Event_Run);
   fChain->SetBranchAddress("Event_Subrun", &Event_Subrun, &b_Event_Subrun);
   fChain->SetBranchAddress("Event_Event", &Event_Event, &b_Event_Event);
   fChain->SetBranchAddress("Nu_True_PDG", &Nu_True_PDG, &b_Nu_True_PDG);
   fChain->SetBranchAddress("Nu_True_SimID", &Nu_True_SimID, &b_Nu_True_SimID);
   fChain->SetBranchAddress("Nu_True_IsNC", &Nu_True_IsNC, &b_Nu_True_IsNC);
   fChain->SetBranchAddress("Nu_True_Mode", &Nu_True_Mode, &b_Nu_True_Mode);
   fChain->SetBranchAddress("Nu_True_TargetZ", &Nu_True_TargetZ, &b_Nu_True_TargetZ);
   fChain->SetBranchAddress("Nu_True_Energy", &Nu_True_Energy, &b_Nu_True_Energy);
   fChain->SetBranchAddress("Nu_True_MomX", &Nu_True_MomX, &b_Nu_True_MomX);
   fChain->SetBranchAddress("Nu_True_MomY", &Nu_True_MomY, &b_Nu_True_MomY);
   fChain->SetBranchAddress("Nu_True_MomZ", &Nu_True_MomZ, &b_Nu_True_MomZ);
   fChain->SetBranchAddress("Nu_True_VertexX", &Nu_True_VertexX, &b_Nu_True_VertexX);
   fChain->SetBranchAddress("Nu_True_VertexY", &Nu_True_VertexY, &b_Nu_True_VertexY);
   fChain->SetBranchAddress("Nu_True_VertexZ", &Nu_True_VertexZ, &b_Nu_True_VertexZ);
   fChain->SetBranchAddress("Event_NMCParticles", &Event_NMCParticles, &b_Event_NMCParticles);
   fChain->SetBranchAddress("Event_NRecoPFPs", &Event_NRecoPFPs, &b_Event_NRecoPFPs);
   fChain->SetBranchAddress("Nu_Reco_VertexX", &Nu_Reco_VertexX, &b_Nu_Reco_VertexX);
   fChain->SetBranchAddress("Nu_Reco_VertexY", &Nu_Reco_VertexY, &b_Nu_Reco_VertexY);
   fChain->SetBranchAddress("Nu_Reco_VertexZ", &Nu_Reco_VertexZ, &b_Nu_Reco_VertexZ);
   fChain->SetBranchAddress("MC_SimID", MC_SimID, &b_MC_SimID);
   fChain->SetBranchAddress("MC_ParentSimID", MC_ParentSimID, &b_MC_ParentSimID);
   fChain->SetBranchAddress("PFP_RecoID", PFP_RecoID, &b_PFP_RecoID);
   fChain->SetBranchAddress("MatchedMC_PDG", MatchedMC_PDG, &b_MatchedMC_PDG);
   fChain->SetBranchAddress("MatchedMC_SimID", MatchedMC_SimID, &b_MatchedMC_SimID);
   fChain->SetBranchAddress("MatchedMC_IsSaved", MatchedMC_IsSaved, &b_MatchedMC_IsSaved);
   fChain->SetBranchAddress("MatchedMC_IsTruePrimary", MatchedMC_IsTruePrimary, &b_MatchedMC_IsTruePrimary);
   fChain->SetBranchAddress("MatchedMC_Generation", MatchedMC_Generation, &b_MatchedMC_Generation);
   fChain->SetBranchAddress("MatchedMC_ParentSimID", MatchedMC_ParentSimID, &b_MatchedMC_ParentSimID);
   fChain->SetBranchAddress("MatchedMC_ParentPDG", MatchedMC_ParentPDG, &b_MatchedMC_ParentPDG);
   fChain->SetBranchAddress("PFP_CurrentReco_IsPrimary", PFP_CurrentReco_IsPrimary, &b_PFP_CurrentReco_IsPrimary);
   fChain->SetBranchAddress("PFP_CurrentReco_Generation", PFP_CurrentReco_Generation, &b_PFP_CurrentReco_Generation);
   fChain->SetBranchAddress("PFP_CurrentReco_ParentRecoID", PFP_CurrentReco_ParentRecoID, &b_PFP_CurrentReco_ParentRecoID);
   fChain->SetBranchAddress("PFP_TrackShowerScore", PFP_TrackShowerScore, &b_PFP_TrackShowerScore);
   fChain->SetBranchAddress("PFP_VertexX", PFP_VertexX, &b_PFP_VertexX);
   fChain->SetBranchAddress("PFP_VertexY", PFP_VertexY, &b_PFP_VertexY);
   fChain->SetBranchAddress("PFP_VertexZ", PFP_VertexZ, &b_PFP_VertexZ);
   fChain->SetBranchAddress("PFP_N2DHits", PFP_N2DHits, &b_PFP_N2DHits);
   fChain->SetBranchAddress("PFP_NSpacepoints", PFP_NSpacepoints, &b_PFP_NSpacepoints);
   fChain->SetBranchAddress("PFP_SpacepointX", &PFP_SpacepointX, &b_PFP_SpacepointX);
   fChain->SetBranchAddress("PFP_SpacepointY", &PFP_SpacepointY, &b_PFP_SpacepointY);
   fChain->SetBranchAddress("PFP_SpacepointZ", &PFP_SpacepointZ, &b_PFP_SpacepointZ);
   fChain->SetBranchAddress("PFP_Completeness", PFP_Completeness, &b_PFP_Completeness);
   fChain->SetBranchAddress("PFP_Purity", PFP_Purity, &b_PFP_Purity);
   fChain->SetBranchAddress("PFP_ShowerFitSuccess", PFP_ShowerFitSuccess, &b_PFP_ShowerFitSuccess);
   fChain->SetBranchAddress("PFP_TrackFitSuccess", PFP_TrackFitSuccess, &b_PFP_TrackFitSuccess);
   fChain->SetBranchAddress("Track_StartX", Track_StartX, &b_Track_StartX);
   fChain->SetBranchAddress("Track_StartY", Track_StartY, &b_Track_StartY);
   fChain->SetBranchAddress("Track_StartZ", Track_StartZ, &b_Track_StartZ);
   fChain->SetBranchAddress("Track_EndX", Track_EndX, &b_Track_EndX);
   fChain->SetBranchAddress("Track_EndY", Track_EndY, &b_Track_EndY);
   fChain->SetBranchAddress("Track_EndZ", Track_EndZ, &b_Track_EndZ);
   fChain->SetBranchAddress("Track_Length", Track_Length, &b_Track_Length);
   fChain->SetBranchAddress("Track_DeflecAngleSD", Track_DeflecAngleSD, &b_Track_DeflecAngleSD);
   fChain->SetBranchAddress("Track_EvalRatio", Track_EvalRatio, &b_Track_EvalRatio);
   fChain->SetBranchAddress("Track_Concentration", Track_Concentration, &b_Track_Concentration);
   fChain->SetBranchAddress("Track_CoreHaloRatio", Track_CoreHaloRatio, &b_Track_CoreHaloRatio);
   fChain->SetBranchAddress("Track_Conicalness", Track_Conicalness, &b_Track_Conicalness);
   fChain->SetBranchAddress("Track_dEdxStart", Track_dEdxStart, &b_Track_dEdxStart);
   fChain->SetBranchAddress("Track_dEdxEnd", Track_dEdxEnd, &b_Track_dEdxEnd);
   fChain->SetBranchAddress("Track_dEdxEndRatio", Track_dEdxEndRatio, &b_Track_dEdxEndRatio);
   fChain->SetBranchAddress("Track_NTrajPoints", &Track_NTrajPoints, &b_Track_NTrajPoints);
   fChain->SetBranchAddress("Track_TrajPointStartIndex", Track_TrajPointStartIndex, &b_Track_TrajPointStartIndex);
   fChain->SetBranchAddress("Track_TrajPointEndIndex", Track_TrajPointEndIndex, &b_Track_TrajPointEndIndex);
   fChain->SetBranchAddress("Track_dEdX", Track_dEdX, &b_Track_dEdX);
   fChain->SetBranchAddress("Track_RR", Track_RR, &b_Track_RR);
   fChain->SetBranchAddress("Shower_StartX", Shower_StartX, &b_Shower_StartX);
   fChain->SetBranchAddress("Shower_StartY", Shower_StartY, &b_Shower_StartY);
   fChain->SetBranchAddress("Shower_StartZ", Shower_StartZ, &b_Shower_StartZ);
   fChain->SetBranchAddress("Shower_Length", Shower_Length, &b_Shower_Length);
   fChain->SetBranchAddress("Shower_OpeningAngle", Shower_OpeningAngle, &b_Shower_OpeningAngle);
   fChain->SetBranchAddress("Shower_InitialdEdx", Shower_InitialdEdx, &b_Shower_InitialdEdx);
   fChain->SetBranchAddress("Shower_BestPlane", Shower_BestPlane, &b_Shower_BestPlane);
   fChain->SetBranchAddress("Shower_Energy", Shower_Energy, &b_Shower_Energy);
   Notify();
}

bool anacode::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void anacode::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t anacode::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef anacode_cxx

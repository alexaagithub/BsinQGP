#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TString.h>
#include <TNtuple.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <iostream>
#include <iomanip>
#include <cmath>

Double_t MUON_MASS = 0.10565837;
Double_t PION_MASS = 0.13957018;
Double_t KAON_MASS = 0.493677;
Double_t KSHORT_MASS = 0.497614;
Double_t KSTAR_MASS = 0.89594;
Double_t PHI_MASS = 1.019455;
Double_t JPSI_MASS = 3.096916;

Int_t MUON_PDGID = 13;
Int_t PION_PDGID = 211;
Int_t KAON_PDGID = 321;
Int_t KSTAR_PDGID = 313;
Int_t PHI_PDGID = 333;
Int_t JPSI_PDGID = 443;
Int_t BZERO_PDGID = 511;
Int_t BPLUS_PDGID = 521;
Int_t BSUBS_PDGID = 531;

#define MAX_XB 16384
#define MAX_MUON 4096
#define MAX_TRACK 8192
#define MAX_GEN 4096
#define MAX_BX 128
#define MAX_Vertices 4096

Double_t tk1mass[7] = {KAON_MASS, PION_MASS, PION_MASS,   KAON_MASS,  KAON_MASS,  KAON_MASS, PION_MASS};
Double_t tk2mass[7] = {0,         0,         PION_MASS,   PION_MASS,  PION_MASS,  KAON_MASS, PION_MASS};
Double_t midmass[7] = {0,         0,         KSHORT_MASS, KSTAR_MASS, KSTAR_MASS, PHI_MASS,  0};
//EvtInfo
Int_t      RunNo;
Int_t      EvtNo;
Int_t      LumiNo;
Int_t      Bsize;
Float_t   PVx;
Float_t   PVy;
Float_t   PVz;
Float_t   PVxE;
Float_t   PVyE;
Float_t   PVzE;
Float_t   PVnchi2;
Float_t   PVchi2;
Float_t   BSx;
Float_t   BSy;
Float_t   BSz;
Float_t   BSxErr;
Float_t   BSyErr;
Float_t   BSzErr;
Float_t   BSdxdz;
Float_t   BSdydz;
Float_t   BSdxdzErr;
Float_t   BSdydzErr;
Float_t   BSWidthX;
Float_t   BSWidthXErr;
Float_t   BSWidthY;
Float_t   BSWidthYErr;
Float_t pthatweight;
Int_t hiBin;
//DInfo
Int_t      Bindex[MAX_XB];
Int_t      Btype[MAX_XB];
Float_t   Bmass[MAX_XB];
Float_t   Bpt[MAX_XB];
Float_t   Beta[MAX_XB];
Float_t   Bphi[MAX_XB];
Float_t   By[MAX_XB];
Float_t   Balpha[MAX_XB];
Float_t   BvtxX[MAX_XB];
Float_t   BvtxY[MAX_XB];
Float_t   Bd0[MAX_XB];
Float_t   Bd0Err[MAX_XB];
Float_t   Bdxyz[MAX_XB];
Float_t   BdxyzErr[MAX_XB];
Float_t   Bchi2ndf[MAX_XB];
Float_t   Bchi2cl[MAX_XB];
Float_t   Bdtheta[MAX_XB];
Float_t   Blxy[MAX_XB];
Float_t   BlxyBS[MAX_XB];
Float_t   BlxyBSErr[MAX_XB];
Bool_t     Bmaxpt[MAX_XB];
Bool_t     Bmaxprob[MAX_XB];
Bool_t     Bbesttktkmass[MAX_XB];
Bool_t     BmaxptMatched[MAX_XB];
Bool_t     BmaxprobMatched[MAX_XB];
Bool_t     BbesttktkmassMatched[MAX_XB];

//BInfo.muonInfo
Float_t   Bmu1pt[MAX_XB];
Float_t   Bmu2pt[MAX_XB];
Float_t   Bmu1p[MAX_XB];
Float_t   Bmu2p[MAX_XB];
Float_t   Bmu1eta[MAX_XB];
Float_t   Bmu2eta[MAX_XB];
Float_t   Bmu1phi[MAX_XB];
Float_t   Bmu2phi[MAX_XB];
Float_t   Bmu1y[MAX_XB];
Float_t   Bmu2y[MAX_XB];
Float_t   Bmu1dzPV[MAX_XB];
Float_t   Bmu2dzPV[MAX_XB];
Float_t   Bmu1dxyPV[MAX_XB];
Float_t   Bmu2dxyPV[MAX_XB];
Float_t   Bmu1normchi2[MAX_XB];
Float_t   Bmu2normchi2[MAX_XB];
Float_t   Bmu1Chi2ndf[MAX_XB];
Float_t   Bmu2Chi2ndf[MAX_XB];
Int_t      Bmu1muqual[MAX_XB];
Int_t      Bmu2muqual[MAX_XB];
Bool_t     Bmu1TrackerMuArbitrated[MAX_XB];
Bool_t     Bmu2TrackerMuArbitrated[MAX_XB];
Bool_t     Bmu1isTrackerMuon[MAX_XB];
Bool_t     Bmu2isTrackerMuon[MAX_XB];
Bool_t     Bmu1isGlobalMuon[MAX_XB];
Bool_t     Bmu2isGlobalMuon[MAX_XB];
Bool_t     Bmu1TMOneStationTight[MAX_XB];
Bool_t     Bmu2TMOneStationTight[MAX_XB];
Int_t    pclusterCompatibilityFilter[MAX_XB];
Int_t    pprimaryVertexFilter[MAX_XB];
Int_t      Bmu1InPixelLayer[MAX_XB];
Int_t      Bmu2InPixelLayer[MAX_XB];
Int_t      Bmu1InStripLayer[MAX_XB];
Int_t      Bmu2InStripLayer[MAX_XB];
Int_t      Bmu1InTrackerLayer[MAX_XB];
Int_t      Bmu2InTrackerLayer[MAX_XB];
//BInfo.mumuInfo
Float_t   Bmumumass[MAX_XB];
Float_t   Bmumueta[MAX_XB];
Float_t   Bmumuphi[MAX_XB];
Float_t   Bmumuy[MAX_XB];
Float_t   Bmumupt[MAX_XB];
//BInfo.ujInfo
Float_t   Bujmass[MAX_XB];
Float_t   BujvProb[MAX_XB];
Float_t   Bujpt[MAX_XB];
Float_t   Bujeta[MAX_XB];
Float_t   Bujphi[MAX_XB];
Float_t   Bujy[MAX_XB];
Float_t   Bujlxy[MAX_XB];

//DInfo.trkInfo
Int_t   Btrk1Idx[MAX_XB];
Int_t   Btrk2Idx[MAX_XB];
Float_t   Btrk1Pt[MAX_XB];
Float_t   Btrk2Pt[MAX_XB];
Float_t   Btrk1D0[MAX_XB];
Float_t   Btrk2D0[MAX_XB];
Float_t   Btrk1Dz[MAX_XB];
Float_t   Btrk2Dz[MAX_XB];
Float_t   Btrk1Eta[MAX_XB];
Float_t   Btrk2Eta[MAX_XB];
Float_t   Btrk1Phi[MAX_XB];
Float_t   Btrk2Phi[MAX_XB];
Float_t   Btrk1PtErr[MAX_XB];
Float_t   Btrk2PtErr[MAX_XB];
Float_t   Btrk1EtaErr[MAX_XB];
Float_t   Btrk2EtaErr[MAX_XB];
Float_t   Btrk1PhiErr[MAX_XB];
Float_t   Btrk2PhiErr[MAX_XB];
Float_t   Btrk1Y[MAX_XB];
Float_t   Btrk2Y[MAX_XB];
Float_t   Btrk1Dxy[MAX_XB];
Float_t   Btrk2Dxy[MAX_XB];
Float_t   Btrk1D0Err[MAX_XB];
Float_t   Btrk2D0Err[MAX_XB];
Float_t   Btrk1PixelHit[MAX_XB];
Float_t   Btrk2PixelHit[MAX_XB];
Float_t   Btrk1StripHit[MAX_XB];
Float_t   Btrk2StripHit[MAX_XB];
Float_t   Btrk1nPixelLayer[MAX_XB];
Float_t   Btrk2nPixelLayer[MAX_XB];
Float_t   Btrk1nStripLayer[MAX_XB];
Float_t   Btrk2nStripLayer[MAX_XB];
Float_t   Btrk1Chi2ndf[MAX_XB];
Float_t   Btrk2Chi2ndf[MAX_XB];
Float_t   Btrk1MVAVal[MAX_XB];
Float_t   Btrk2MVAVal[MAX_XB];
Int_t      Btrk1Algo[MAX_XB];
Int_t      Btrk2Algo[MAX_XB];
Bool_t     Btrk1highPurity[MAX_XB];
Bool_t     Btrk2highPurity[MAX_XB];
Int_t      Btrk1Quality[MAX_XB];
Int_t      Btrk2Quality[MAX_XB];
//BInfo.tktkInfo
Float_t   Btktkmass[MAX_XB];
Float_t   BtktkmassKK[MAX_XB];
Float_t   BtktkvProb[MAX_XB];
Float_t   Btktkpt[MAX_XB];
Float_t   Btktketa[MAX_XB];
Float_t   Btktkphi[MAX_XB];
Float_t   Btktky[MAX_XB];
Float_t   Bdoubletmass[MAX_XB];
Float_t   Bdoubletpt[MAX_XB];
Float_t   Bdoubleteta[MAX_XB];
Float_t   Bdoubletphi[MAX_XB];
Float_t   Bdoublety[MAX_XB];

//BInfo.genInfo
Float_t   Bgen[MAX_XB];
Int_t      BgenIndex[MAX_XB];
Float_t   Bgenpt[MAX_XB];
Float_t   Bgeneta[MAX_XB];
Float_t   Bgenphi[MAX_XB];
Float_t   Bgeny[MAX_XB];

/*
Int_t   BMaxDoca[MAX_XB];
Double_t   Balpha[MAX_XB];*/
Float_t   BsvpvDistance[MAX_XB];
Float_t   BsvpvDisErr[MAX_XB];
Float_t   BsvpvDistance_2D[MAX_XB];
Float_t   BsvpvDisErr_2D[MAX_XB];
/*Int_t      kstar[MAX_XB]; 
Double_t   Btrk1MassHypo[MAX_XB];
Double_t   Btrk2MassHypo[MAX_XB];
Double_t   Btrkminpt[MAX_XB];
Double_t   Btrkmaxpt[MAX_XB];
Int_t      Btrkminptindex[MAX_XB];
Int_t      Btrkmaxptindex[MAX_XB];*/


void buildBranch(TTree* nt)
{
  //EvtInfo
  nt->Branch("RunNo",&RunNo);
  nt->Branch("EvtNo",&EvtNo);
  nt->Branch("LumiNo",&LumiNo);
  nt->Branch("Bsize",&Bsize);
  nt->Branch("PVx",&PVx);
  nt->Branch("PVy",&PVy);
  nt->Branch("PVz",&PVz);
  nt->Branch("PVxE",&PVxE);
  nt->Branch("PVyE",&PVyE);
  nt->Branch("PVzE",&PVzE);
  nt->Branch("PVnchi2",&PVnchi2);
  nt->Branch("PVchi2",&PVchi2);
  nt->Branch("BSx",&BSx);
  nt->Branch("BSy",&BSy);
  nt->Branch("BSz",&BSz);
  nt->Branch("BSxErr",&BSxErr);
  nt->Branch("BSyErr",&BSyErr);
  nt->Branch("BSzErr",&BSzErr);
  nt->Branch("BSdxdz",&BSdxdz);
  nt->Branch("BSdydz",&BSdydz);
  nt->Branch("BSdxdzErr",&BSdxdzErr);
  nt->Branch("BSdydzErr",&BSdydzErr);
  nt->Branch("BSWidthX",&BSWidthX);
  nt->Branch("BSWidthXErr",&BSWidthXErr);
  nt->Branch("BSWidthY",&BSWidthY);
  nt->Branch("BSWidthYErr",&BSWidthYErr);
  //BInfo
  nt->Branch("Bindex",Bindex,"Bindex[Bsize]/I");
  nt->Branch("Btype",Btype,"Btype[Bsize]/I");
  nt->Branch("Bmass",Bmass,"Bmass[Bsize]/D");
  nt->Branch("Bpt",Bpt,"Bpt[Bsize]/D");
  nt->Branch("Beta",Beta,"Beta[Bsize]/D");
  nt->Branch("Bphi",Bphi,"Bphi[Bsize]/D");
  nt->Branch("By",By,"By[Bsize]/D");
  nt->Branch("BvtxX",BvtxX,"BvtxX[Bsize]/D");
  nt->Branch("BvtxY",BvtxY,"BvtxY[Bsize]/D");
  nt->Branch("Bd0",Bd0,"Bd0[Bsize]/D");
  nt->Branch("Bd0Err",Bd0Err,"Bd0Err[Bsize]/D");
  nt->Branch("Bdxyz",Bdxyz,"Bdxyz[Bsize]/D");
  nt->Branch("BdxyzErr",BdxyzErr,"BdxyzErr[Bsize]/D");
  nt->Branch("Bchi2ndf",Bchi2ndf,"Bchi2ndf[Bsize]/D");
  nt->Branch("Bchi2cl",Bchi2cl,"Bchi2cl[Bsize]/D");
  nt->Branch("Bdtheta",Bdtheta,"Bdtheta[Bsize]/D");
  nt->Branch("Blxy",Blxy,"Blxy[Bsize]/D");
  nt->Branch("BlxyBS",BlxyBS,"BlxyBS[Bsize]/D");
  nt->Branch("BlxyBSErr",BlxyBSErr,"BlxyBSErr[Bsize]/D");
  nt->Branch("Bmaxpt",Bmaxpt,"Bmaxpt[Bsize]/O");
  nt->Branch("Bmaxprob",Bmaxprob,"Bmaxprob[Bsize]/O");
  nt->Branch("Bbesttktkmass",Bbesttktkmass,"Bbesttktkmass[Bsize]/O");
  nt->Branch("BmaxptMatched",BmaxptMatched,"BmaxptMatched[Bsize]/O");
  nt->Branch("BmaxprobMatched",BmaxprobMatched,"BmaxprobMatched[Bsize]/O");
  nt->Branch("BbesttktkmassMatched",BbesttktkmassMatched,"BbesttktkmassMatched[Bsize]/O");

  //BInfo.trkInfo
  nt->Branch("Btrk1Idx",Btrk1Idx,"Btrk1Idx[Bsize]/I");
  nt->Branch("Btrk2Idx",Btrk2Idx,"Btrk2Idx[Bsize]/I");
  nt->Branch("Btrk1Pt",Btrk1Pt,"Btrk1Pt[Bsize]/D");
  nt->Branch("Btrk2Pt",Btrk2Pt,"Btrk2Pt[Bsize]/D");
  nt->Branch("Btrk1Eta",Btrk1Eta,"Btrk1Eta[Bsize]/D");  
  nt->Branch("Btrk2Eta",Btrk2Eta,"Btrk2Eta[Bsize]/D");  
  nt->Branch("Btrk1Phi",Btrk1Phi,"Btrk1Phi[Bsize]/D");  
  nt->Branch("Btrk2Phi",Btrk2Phi,"Btrk2Phi[Bsize]/D");  
  nt->Branch("Btrk1PtErr",Btrk1PtErr,"Btrk1PtErr[Bsize]/D");  
  nt->Branch("Btrk2PtErr",Btrk2PtErr,"Btrk2PtErr[Bsize]/D");
  nt->Branch("Btrk1EtaErr",Btrk1EtaErr,"Btrk1EtaErr[Bsize]/D");
  nt->Branch("Btrk2EtaErr",Btrk2EtaErr,"Btrk2EtaErr[Bsize]/D");
  nt->Branch("Btrk1PhiErr",Btrk1PhiErr,"Btrk1PhiErr[Bsize]/D");
  nt->Branch("Btrk2PhiErr",Btrk2PhiErr,"Btrk2PhiErr[Bsize]/D");
  nt->Branch("Btrk1Y",Btrk1Y,"Btrk1Y[Bsize]/D");  
  nt->Branch("Btrk2Y",Btrk2Y,"Btrk2Y[Bsize]/D");  
  nt->Branch("Btrk1Dxy",Btrk1Dxy,"Btrk1Dxy[Bsize]/D");
  nt->Branch("Btrk2Dxy",Btrk2Dxy,"Btrk2Dxy[Bsize]/D");
  nt->Branch("Btrk1D0Err",Btrk1D0Err,"Btrk1D0Err[Bsize]/D");
  nt->Branch("Btrk2D0Err",Btrk2D0Err,"Btrk2D0Err[Bsize]/D");
  nt->Branch("Btrk1PixelHit",Btrk1PixelHit,"Btrk1PixelHit[Bsize]/D");
  nt->Branch("Btrk2PixelHit",Btrk2PixelHit,"Btrk2PixelHit[Bsize]/D");
  nt->Branch("Btrk1StripHit",Btrk1StripHit,"Btrk1StripHit[Bsize]/D");
  nt->Branch("Btrk2StripHit",Btrk2StripHit,"Btrk2StripHit[Bsize]/D");
  nt->Branch("Btrk1nPixelLayer",Btrk1nPixelLayer,"Btrk1nPixelLayer[Bsize]/D");
  nt->Branch("Btrk2nPixelLayer",Btrk2nPixelLayer,"Btrk2nPixelLayer[Bsize]/D");
  nt->Branch("Btrk1nStripLayer",Btrk1nStripLayer,"Btrk1nStripLayer[Bsize]/D");
  nt->Branch("Btrk2nStripLayer",Btrk2nStripLayer,"Btrk2nStripLayer[Bsize]/D");
  nt->Branch("Btrk1Chi2ndf",Btrk1Chi2ndf,"Btrk1Chi2ndf[Bsize]/D");
  nt->Branch("Btrk2Chi2ndf",Btrk2Chi2ndf,"Btrk2Chi2ndf[Bsize]/D");
  nt->Branch("Btrk1MVAVal",Btrk1MVAVal,"Btrk1MVAVal[Bsize]/D");
  nt->Branch("Btrk2MVAVal",Btrk2MVAVal,"Btrk2MVAVal[Bsize]/D");
  nt->Branch("Btrk1Algo",Btrk1Algo,"Btrk1Algo[Bsize]/I");
  nt->Branch("Btrk2Algo",Btrk2Algo,"Btrk2Algo[Bsize]/I");
  nt->Branch("Btrk1highPurity",Btrk1highPurity,"Btrk1highPurity[Bsize]/O");
  nt->Branch("Btrk2highPurity",Btrk2highPurity,"Btrk2highPurity[Bsize]/O");
  nt->Branch("Btrk1Quality",Btrk1Quality,"Btrk1Quality[Bsize]/I");
  nt->Branch("Btrk2Quality",Btrk2Quality,"Btrk2Quality[Bsize]/I");
  //BInfo.tktkInfo
  nt->Branch("Btktkmass",Btktkmass,"Btktkmass[Bsize]/D");
  nt->Branch("BtktkmassKK",BtktkmassKK,"BtktkmassKK[Bsize]/D");
  nt->Branch("BtktkvProb",BtktkvProb,"BtktkvProb[Bsize]/D");
  nt->Branch("Btktkpt",Btktkpt,"Btktkpt[Bsize]/D");
  nt->Branch("Btktketa",Btktketa,"Btktketa[Bsize]/D");
  nt->Branch("Btktkphi",Btktkphi,"Btktkphi[Bsize]/D");
  nt->Branch("Btktky",Btktky,"Btktky[Bsize]/D");
  nt->Branch("Bdoubletmass",Bdoubletmass,"Bdoubletmass[Bsize]/D");
  nt->Branch("Bdoubletpt",Bdoubletpt,"Bdoubletpt[Bsize]/D");
  nt->Branch("Bdoubleteta",Bdoubleteta,"Bdoubleteta[Bsize]/D");  
  nt->Branch("Bdoubletphi",Bdoubletphi,"Bdoubletphi[Bsize]/D");  
  nt->Branch("Bdoublety",Bdoublety,"Bdoublety[Bsize]/D");
  
  //BInfo.muonInfo
  nt->Branch("Bmu1pt",Bmu1pt,"Bmu1pt[Bsize]/D");
  nt->Branch("Bmu2pt",Bmu2pt,"Bmu2pt[Bsize]/D");
  nt->Branch("Bmu1p",Bmu1p,"Bmu1p[Bsize]/D");
  nt->Branch("Bmu2p",Bmu2p,"Bmu2p[Bsize]/D");
  nt->Branch("Bmu1eta",Bmu1eta,"Bmu1eta[Bsize]/D");
  nt->Branch("Bmu2eta",Bmu2eta,"Bmu2eta[Bsize]/D");
  nt->Branch("Bmu1phi",Bmu1phi,"Bmu1phi[Bsize]/D");
  nt->Branch("Bmu2phi",Bmu2phi,"Bmu2phi[Bsize]/D");
  nt->Branch("Bmu1y",Bmu1y,"Bmu1y[Bsize]/D");
  nt->Branch("Bmu2y",Bmu2y,"Bmu2y[Bsize]/D");
  nt->Branch("Bmu1dzPV",Bmu1dzPV,"Bmu1dzPV[Bsize]/D");
  nt->Branch("Bmu2dzPV",Bmu2dzPV,"Bmu2dzPV[Bsize]/D");
  nt->Branch("Bmu1dxyPV",Bmu1dxyPV,"Bmu1dxyPV[Bsize]/D");
  nt->Branch("Bmu2dxyPV",Bmu2dxyPV,"Bmu2dxyPV[Bsize]/D");
  nt->Branch("Bmu1normchi2",Bmu1normchi2,"Bmu1normchi2[Bsize]/D");
  nt->Branch("Bmu2normchi2",Bmu2normchi2,"Bmu2normchi2[Bsize]/D");
  nt->Branch("Bmu1Chi2ndf",Bmu1Chi2ndf,"Bmu1Chi2ndf[Bsize]/D");
  nt->Branch("Bmu2Chi2ndf",Bmu2Chi2ndf,"Bmu2Chi2ndf[Bsize]/D");
  nt->Branch("Bmu1muqual",Bmu1muqual,"Bmu1muqual[Bsize]/I");
  nt->Branch("Bmu1muqual",Bmu1muqual,"Bmu1muqual[Bsize]/I");
  nt->Branch("Bmu1TrackerMuArbitrated",Bmu1TrackerMuArbitrated,"Bmu1TrackerMuArbitrated[Bsize]/O");
  nt->Branch("Bmu2TrackerMuArbitrated",Bmu2TrackerMuArbitrated,"Bmu2TrackerMuArbitrated[Bsize]/O");
  nt->Branch("Bmu1isTrackerMuon",Bmu1isTrackerMuon,"Bmu1isTrackerMuon[Bsize]/O");
  nt->Branch("Bmu2isTrackerMuon",Bmu2isTrackerMuon,"Bmu2isTrackerMuon[Bsize]/O");
  nt->Branch("Bmu1isGlobalMuon",Bmu1isGlobalMuon,"Bmu1isGlobalMuon[Bsize]/O");
  nt->Branch("Bmu2isGlobalMuon",Bmu2isGlobalMuon,"Bmu2isGlobalMuon[Bsize]/O");
  nt->Branch("Bmu1TMOneStationTight",Bmu1TMOneStationTight,"Bmu1TMOneStationTight[Bsize]/O");
  nt->Branch("Bmu2TMOneStationTight",Bmu2TMOneStationTight,"Bmu2TMOneStationTight[Bsize]/O");
  nt->Branch("Bmu1InPixelLayer",Bmu1InPixelLayer,"Bmu1InPixelLayer[Bsize]/I");
  nt->Branch("Bmu2InPixelLayer",Bmu2InPixelLayer,"Bmu2InPixelLayer[Bsize]/I");
  nt->Branch("Bmu1InStripLayer",Bmu1InStripLayer,"Bmu1InStripLayer[Bsize]/I");
  nt->Branch("Bmu2InStripLayer",Bmu2InStripLayer,"Bmu2InStripLayer[Bsize]/I");
  nt->Branch("Bmu1InTrackerLayer",Bmu1InTrackerLayer,"Bmu1InTrackerLayer[Bsize]/I");
  nt->Branch("Bmu2InTrackerLayer",Bmu2InTrackerLayer,"Bmu2InTrackerLayer[Bsize]/I");
  nt->Branch("Bmumumass",Bmumumass,"Bmumumass[Bsize]/D");
  nt->Branch("Bmumueta",Bmumueta,"Bmumueta[Bsize]/D");
  nt->Branch("Bmumuphi",Bmumuphi,"Bmumuphi[Bsize]/D");
  nt->Branch("Bmumuy",Bmumuy,"Bmumuy[Bsize]/D");
  nt->Branch("Bmumupt",Bmumupt,"Bmumupt[Bsize]/D");
  nt->Branch("Bujmass",Bujmass,"Bujmass[Bsize]/D");
  nt->Branch("BujvProb",BujvProb,"BujvProb[Bsize]/D");
  nt->Branch("Bujpt",Bujpt,"Bujpt[Bsize]/D");
  nt->Branch("Bujeta",Bujeta,"Bujeta[Bsize]/D");
  nt->Branch("Bujphi",Bujphi,"Bujphi[Bsize]/D");
  nt->Branch("Bujy",Bujy,"Bujy[Bsize]/D");
  nt->Branch("Bujlxy",Bujlxy,"Bujlxy[Bsize]/D");

  //BInfo.genInfo
  nt->Branch("Bgen",Bgen,"Bgen[Bsize]/D");
  nt->Branch("BgenIndex",BgenIndex,"BgenIndex[Bsize]/I");
  nt->Branch("Bgenpt",Bgenpt,"Bgenpt[Bsize]/D");
  nt->Branch("Bgeny",Bgeny,"Bgeny[Bsize]/D");
  nt->Branch("Bgeneta",Bgeneta,"Bgeneta[Bsize]/D");
  nt->Branch("Bgenphi",Bgenphi,"Bgenphi[Bsize]/D");

  
  //nt->Branch("Balpha",Balpha,"Balpha[Bsize]/D");
  nt->Branch("BsvpvDistance",BsvpvDistance,"BsvpvDistance[Bsize]/F");
  nt->Branch("BsvpvDisErr",BsvpvDisErr,"BsvpvDisErr[Bsize]/F");
  /*nt->Branch("BsvpvDistance_2D",BsvpvDistance_2D,"BsvpvDistance_2D[Bsize]/D");
  nt->Branch("BsvpvDisErr_2D",BsvpvDisErr_2D,"BsvpvDisErr_2D[Bsize]/D");
  nt->Branch("BMaxDoca",BMaxDoca,"BMaxDoca[Bsize]/D");
  nt->Branch("Btrk1MassHypo",Btrk1MassHypo,"Btrk1MassHypo[Bsize]/D");
  nt->Branch("Btrk2MassHypo",Btrk2MassHypo,"Btrk2MassHypo[Bsize]/D");
  nt->Branch("Btrkminpt",Btrkminpt,"Btrkminpt[Bsize]/D");
  nt->Branch("Btrkmaxpt",Btrkmaxpt,"Btrkmaxpt[Bsize]/D");
  nt->Branch("Btrkminptindex",Btrkminptindex,"Btrkminptindex[Bsize]/I");
  nt->Branch("Btrkmaxptindex",Btrkmaxptindex,"Btrkmaxptindex[Bsize]/I");
  */

}

Int_t      Gsize;
Float_t  Gy[MAX_GEN];
Float_t   Geta[MAX_GEN];
Float_t   Gphi[MAX_GEN];
Float_t   Gpt[MAX_GEN];
Float_t   GpdgId[MAX_GEN];
Float_t   GisSignal[MAX_GEN];
Float_t   Gmu1pt[MAX_GEN];
Float_t   Gmu2pt[MAX_GEN];
Float_t   Gmu1p[MAX_GEN];
Float_t   Gmu2p[MAX_GEN];
Float_t   Gmu1eta[MAX_GEN];
Float_t   Gmu2eta[MAX_GEN];
Float_t   Gmu1phi[MAX_GEN];
Float_t   Gmu2phi[MAX_GEN];
Float_t   Gtk1pt[MAX_GEN];
Float_t   Gtk2pt[MAX_GEN];
Float_t   Gtk1eta[MAX_GEN];
Float_t   Gtk2eta[MAX_GEN];
Float_t   Gtk1phi[MAX_GEN];
Float_t   Gtk2phi[MAX_GEN];

void buildGenBranch(TTree* nt)
{
  nt->Branch("Gsize",&Gsize);
  nt->Branch("Gy",Gy,"Gy[Gsize]/D");
  nt->Branch("Geta",Geta,"Geta[Gsize]/D");
  nt->Branch("Gphi",Gphi,"Gphi[Gsize]/D");
  nt->Branch("Gpt",Gpt,"Gpt[Gsize]/D");
  nt->Branch("GpdgId",GpdgId,"GpdgId[Gsize]/D");
  nt->Branch("GisSignal",GisSignal,"GisSignal[Gsize]/D");
  nt->Branch("Gmu1eta",Gmu1eta,"Gmu1eta[Gsize]/D");
  nt->Branch("Gmu1phi",Gmu1phi,"Gmu1phi[Gsize]/D");
  nt->Branch("Gmu1pt",Gmu1pt,"Gmu1pt[Gsize]/D");
  nt->Branch("Gmu1p",Gmu1p,"Gmu1p[Gsize]/D");
  nt->Branch("Gmu2eta",Gmu2eta,"Gmu2eta[Gsize]/D");
  nt->Branch("Gmu2phi",Gmu2phi,"Gmu2phi[Gsize]/D");
  nt->Branch("Gmu2pt",Gmu2pt,"Gmu2pt[Gsize]/D");
  nt->Branch("Gmu2p",Gmu2p,"Gmu2p[Gsize]/D");
  nt->Branch("Gtk1pt",Gtk1pt,"Gtk1pt[Gsize]/D");
  nt->Branch("Gtk1eta",Gtk1eta,"Gtk1eta[Gsize]/D");
  nt->Branch("Gtk1phi",Gtk1phi,"Gtk1phi[Gsize]/D");
  nt->Branch("Gtk2pt",Gtk2pt,"Gtk2pt[Gsize]/D");
  nt->Branch("Gtk2eta",Gtk2eta,"Gtk2eta[Gsize]/D");
  nt->Branch("Gtk2phi",Gtk2phi,"Gtk2phi[Gsize]/D");
}

//EvtInfo
Int_t           EvtInfo_RunNo;
Int_t           EvtInfo_EvtNo;
Int_t           EvtInfo_BxNo;
Int_t           EvtInfo_LumiNo;
Int_t           EvtInfo_Orbit;
Bool_t          EvtInfo_McFlag;
Int_t           EvtInfo_nBX;
Int_t           EvtInfo_BXPU[MAX_BX];
Int_t           EvtInfo_nPU[MAX_BX];
Float_t         EvtInfo_trueIT[MAX_BX];
Double_t        EvtInfo_PVx;
Double_t        EvtInfo_PVy;
Double_t        EvtInfo_PVz;
Double_t        EvtInfo_PVxE;
Double_t        EvtInfo_PVyE;
Double_t        EvtInfo_PVzE;
Double_t        EvtInfo_PVnchi2;
Double_t        EvtInfo_PVchi2;
Double_t        EvtInfo_BSx;
Double_t        EvtInfo_BSy;
Double_t        EvtInfo_BSz;
Double_t        EvtInfo_BSxErr;
Double_t        EvtInfo_BSyErr;
Double_t        EvtInfo_BSzErr;
Double_t        EvtInfo_BSdxdz;
Double_t        EvtInfo_BSdydz;
Double_t        EvtInfo_BSdxdzErr;
Double_t        EvtInfo_BSdydzErr;
Double_t        EvtInfo_BSWidthX;
Double_t        EvtInfo_BSWidthXErr;
Double_t        EvtInfo_BSWidthY;
Double_t        EvtInfo_BSWidthYErr;
//MuonInfo
Int_t           MuonInfo_size;
Int_t           MuonInfo_index[MAX_MUON];
Int_t           MuonInfo_handle_index[MAX_MUON];
Int_t           MuonInfo_charge[MAX_MUON];
Double_t        MuonInfo_pt[MAX_MUON];
Double_t        MuonInfo_eta[MAX_MUON];
Double_t        MuonInfo_phi[MAX_MUON];
Double_t        MuonInfo_ptErr[MAX_MUON];
Double_t        MuonInfo_etaErr[MAX_MUON];
Double_t        MuonInfo_phiErr[MAX_MUON];
Bool_t          MuonInfo_isTrackerMuon[MAX_MUON];
Bool_t          MuonInfo_isGlobalMuon[MAX_MUON];
Int_t           MuonInfo_muqual[MAX_MUON];
Double_t        MuonInfo_iso_trk[MAX_MUON];
Double_t        MuonInfo_iso_ecal[MAX_MUON];
Double_t        MuonInfo_iso_hcal[MAX_MUON];
Int_t           MuonInfo_type[MAX_MUON];
Double_t        MuonInfo_n_matches[MAX_MUON];
Bool_t          MuonInfo_TMOneStationTight[MAX_MUON];
Bool_t          MuonInfo_TrackerMuonArbitrated[MAX_MUON];
Int_t           MuonInfo_geninfo_index[MAX_MUON];
Bool_t          MuonInfo_BfinderMuID[MAX_MUON];
Bool_t          MuonInfo_SoftMuID[MAX_MUON];
Bool_t          MuonInfo_isStandAloneMuon[MAX_MUON];
Int_t           MuonInfo_StandAloneMuon_charge[MAX_MUON];
Double_t        MuonInfo_StandAloneMuon_pt[MAX_MUON];
Double_t        MuonInfo_StandAloneMuon_eta[MAX_MUON];
Double_t        MuonInfo_StandAloneMuon_phi[MAX_MUON];
Double_t        MuonInfo_StandAloneMuon_d0[MAX_MUON];
Double_t        MuonInfo_StandAloneMuon_dz[MAX_MUON];
Double_t        MuonInfo_StandAloneMuon_dzPV[MAX_MUON];
Double_t        MuonInfo_StandAloneMuon_dxyPV[MAX_MUON];
Bool_t          MuonInfo_outerTrackisNonnull[MAX_MUON];
Bool_t          MuonInfo_innerTrackisNonnull[MAX_MUON];
Bool_t          MuonInfo_globalTrackisNonnull[MAX_MUON];
Int_t           MuonInfo_innerTrackQuality[MAX_MUON];
Double_t        MuonInfo_normchi2[MAX_MUON];
Int_t           MuonInfo_i_striphit[MAX_MUON];
Int_t           MuonInfo_i_pixelhit[MAX_MUON];
Int_t           MuonInfo_i_nStripLayer[MAX_MUON];
Int_t           MuonInfo_i_nPixelLayer[MAX_MUON];
Double_t        MuonInfo_i_chi2[MAX_MUON];
Double_t        MuonInfo_i_ndf[MAX_MUON];
Int_t           MuonInfo_fpbarrelhit[MAX_MUON];
Int_t           MuonInfo_fpendcaphit[MAX_MUON];
Double_t        MuonInfo_d0[MAX_MUON];
Double_t        MuonInfo_dz[MAX_MUON];
Double_t        MuonInfo_dzPV[MAX_MUON];
Double_t        MuonInfo_dxyPV[MAX_MUON];
Double_t        MuonInfo_g_chi2[MAX_MUON];
Double_t        MuonInfo_g_ndf[MAX_MUON];
Int_t           MuonInfo_g_striphit[MAX_MUON];
Int_t           MuonInfo_g_pixelhit[MAX_MUON];
Int_t           MuonInfo_nmuhit[MAX_MUON];
Bool_t          MuonInfo_isTriggered[MAX_MUON];
Int_t           MuonInfo_MuTrgMatchPathSize;

//TrackInfo
Int_t           TrackInfo_size;
Int_t           TrackInfo_index[MAX_TRACK];
Int_t           TrackInfo_handle_index[MAX_TRACK];
Int_t           TrackInfo_charge[MAX_TRACK];
Double_t        TrackInfo_pt[MAX_TRACK];
Double_t        TrackInfo_eta[MAX_TRACK];
Double_t        TrackInfo_phi[MAX_TRACK];
Double_t        TrackInfo_ptErr[MAX_TRACK];
Double_t        TrackInfo_etaErr[MAX_TRACK];
Double_t        TrackInfo_phiErr[MAX_TRACK];
Int_t           TrackInfo_striphit[MAX_TRACK];
Int_t           TrackInfo_pixelhit[MAX_TRACK];
Int_t           TrackInfo_nStripLayer[MAX_TRACK];
Int_t           TrackInfo_nPixelLayer[MAX_TRACK];
Int_t           TrackInfo_fpbarrelhit[MAX_TRACK];
Int_t           TrackInfo_fpendcaphit[MAX_TRACK];
Double_t        TrackInfo_chi2[MAX_TRACK];
Double_t        TrackInfo_ndf[MAX_TRACK];
Double_t        TrackInfo_d0[MAX_TRACK];
Double_t        TrackInfo_d0error[MAX_TRACK];
Double_t        TrackInfo_dzPV[MAX_TRACK];
Double_t        TrackInfo_dxyPV[MAX_TRACK];
Int_t           TrackInfo_geninfo_index[MAX_TRACK];
Int_t           TrackInfo_trackQuality[MAX_TRACK];
Bool_t          TrackInfo_highPurity[MAX_TRACK];
Double_t        TrackInfo_trkMVAVal[MAX_TRACK];
Int_t           TrackInfo_trkAlgo[MAX_TRACK];

//BInfo
Int_t           BInfo_uj_size;
Int_t           BInfo_uj_index[MAX_XB];
Double_t        BInfo_uj_mass[MAX_XB];
Double_t        BInfo_uj_pt[MAX_XB];
Double_t        BInfo_uj_eta[MAX_XB];
Double_t        BInfo_uj_phi[MAX_XB];
Double_t        BInfo_uj_vtxX[MAX_XB];
Double_t        BInfo_uj_vtxY[MAX_XB];
Double_t        BInfo_uj_vtxZ[MAX_XB];
Double_t        BInfo_uj_vtxXErr[MAX_XB];
Double_t        BInfo_uj_vtxYErr[MAX_XB];
Double_t        BInfo_uj_vtxZErr[MAX_XB];
Double_t        BInfo_uj_vtxYXErr[MAX_XB];
Double_t        BInfo_uj_vtxZXErr[MAX_XB];
Double_t        BInfo_uj_vtxZYErr[MAX_XB];
Double_t        BInfo_uj_vtxdof[MAX_XB];
Double_t        BInfo_uj_vtxchi2[MAX_XB];
Int_t           BInfo_uj_rfmu1_index[MAX_XB];
Int_t           BInfo_uj_rfmu2_index[MAX_XB];
Double_t        BInfo_uj_rfmu1_pt[MAX_XB];
Double_t        BInfo_uj_rfmu1_eta[MAX_XB];
Double_t        BInfo_uj_rfmu1_phi[MAX_XB];
Double_t        BInfo_uj_rfmu2_pt[MAX_XB];
Double_t        BInfo_uj_rfmu2_eta[MAX_XB];
Double_t        BInfo_uj_rfmu2_phi[MAX_XB];

Int_t           BInfo_size;
Int_t           BInfo_index[MAX_XB];
Double_t        BInfo_mass[MAX_XB];
Double_t        BInfo_pt[MAX_XB];
Double_t        BInfo_eta[MAX_XB];
Double_t        BInfo_phi[MAX_XB];
Double_t        BInfo_pxE[MAX_XB];
Double_t        BInfo_pyE[MAX_XB];
Double_t        BInfo_pzE[MAX_XB];
Double_t        BInfo_vtxX[MAX_XB];
Double_t        BInfo_vtxY[MAX_XB];
Double_t        BInfo_vtxZ[MAX_XB];
Double_t        BInfo_vtxXErr[MAX_XB];
Double_t        BInfo_vtxYErr[MAX_XB];
Double_t        BInfo_vtxZErr[MAX_XB];
Double_t        BInfo_vtxYXErr[MAX_XB];
Double_t        BInfo_vtxZXErr[MAX_XB];
Double_t        BInfo_vtxZYErr[MAX_XB];
Double_t        BInfo_vtxdof[MAX_XB];
Double_t        BInfo_vtxchi2[MAX_XB];
Int_t           BInfo_rfuj_index[MAX_XB];
Int_t           BInfo_rftk1_index[MAX_XB];
Int_t           BInfo_rftk2_index[MAX_XB];
Int_t           BInfo_type[MAX_XB];
Double_t        BInfo_rfmu1_pt[MAX_XB];
Double_t        BInfo_rfmu1_eta[MAX_XB];
Double_t        BInfo_rfmu1_phi[MAX_XB];
Double_t        BInfo_rfmu2_pt[MAX_XB];
Double_t        BInfo_rfmu2_eta[MAX_XB];
Double_t        BInfo_rfmu2_phi[MAX_XB];
Double_t        BInfo_rftk1_pt[MAX_XB];
Double_t        BInfo_rftk1_eta[MAX_XB];
Double_t        BInfo_rftk1_phi[MAX_XB];
Double_t        BInfo_rftk2_pt[MAX_XB];
Double_t        BInfo_rftk2_eta[MAX_XB];
Double_t        BInfo_rftk2_phi[MAX_XB];
Double_t        BInfo_tktk_mass[MAX_XB];
Double_t        BInfo_tktk_pt[MAX_XB];
Double_t        BInfo_tktk_eta[MAX_XB];
Double_t        BInfo_tktk_phi[MAX_XB];
Double_t        BInfo_tktk_vtxX[MAX_XB];
Double_t        BInfo_tktk_vtxY[MAX_XB];
Double_t        BInfo_tktk_vtxZ[MAX_XB];
Double_t        BInfo_tktk_vtxXErr[MAX_XB];
Double_t        BInfo_tktk_vtxYErr[MAX_XB];
Double_t        BInfo_tktk_vtxZErr[MAX_XB];
Double_t        BInfo_tktk_vtxYXErr[MAX_XB];
Double_t        BInfo_tktk_vtxZXErr[MAX_XB];
Double_t        BInfo_tktk_vtxZYErr[MAX_XB];
Double_t        BInfo_tktk_vtxdof[MAX_XB];
Double_t        BInfo_tktk_vtxchi2[MAX_XB];
Double_t        BInfo_tktk_rftk1_pt[MAX_XB];
Double_t        BInfo_tktk_rftk1_eta[MAX_XB];
Double_t        BInfo_tktk_rftk1_phi[MAX_XB];
Double_t        BInfo_tktk_rftk2_pt[MAX_XB];
Double_t        BInfo_tktk_rftk2_eta[MAX_XB];
Double_t        BInfo_tktk_rftk2_phi[MAX_XB];


Int_t           GenInfo_size;
Int_t           GenInfo_index[MAX_GEN];
Int_t           GenInfo_handle_index[MAX_GEN];
Double_t        GenInfo_pt[MAX_GEN];
Double_t        GenInfo_eta[MAX_GEN];
Double_t        GenInfo_phi[MAX_GEN];
Double_t        GenInfo_mass[MAX_GEN];
Int_t           GenInfo_pdgId[MAX_GEN];
Int_t           GenInfo_status[MAX_GEN];
Int_t           GenInfo_nMo[MAX_GEN];
Int_t           GenInfo_nDa[MAX_GEN];
Int_t           GenInfo_mo1[MAX_GEN];
Int_t           GenInfo_mo2[MAX_GEN];
Int_t           GenInfo_da1[MAX_GEN];
Int_t           GenInfo_da2[MAX_GEN];
Int_t           GenInfo_da3[MAX_GEN];
Int_t           GenInfo_da4[MAX_GEN];

void setBBranch(TTree *root)
{
   root->SetBranchAddress("EvtInfo.RunNo",&EvtInfo_RunNo);
   root->SetBranchAddress("EvtInfo.EvtNo",&EvtInfo_EvtNo);
   root->SetBranchAddress("EvtInfo.BxNo",&EvtInfo_BxNo);
   root->SetBranchAddress("EvtInfo.LumiNo",&EvtInfo_LumiNo);
   root->SetBranchAddress("EvtInfo.Orbit",&EvtInfo_Orbit);
   root->SetBranchAddress("EvtInfo.McFlag",&EvtInfo_McFlag);
   root->SetBranchAddress("EvtInfo.nBX",&EvtInfo_nBX);
   root->SetBranchAddress("EvtInfo.BXPU",&EvtInfo_BXPU);
   root->SetBranchAddress("EvtInfo.nPU",&EvtInfo_nPU);
   root->SetBranchAddress("EvtInfo.trueIT",&EvtInfo_trueIT);
   root->SetBranchAddress("EvtInfo.PVx",&EvtInfo_PVx);
   root->SetBranchAddress("EvtInfo.PVy",&EvtInfo_PVy);
   root->SetBranchAddress("EvtInfo.PVz",&EvtInfo_PVz);
   root->SetBranchAddress("EvtInfo.PVxE",&EvtInfo_PVxE);
   root->SetBranchAddress("EvtInfo.PVyE",&EvtInfo_PVyE);
   root->SetBranchAddress("EvtInfo.PVzE",&EvtInfo_PVzE);
   root->SetBranchAddress("EvtInfo.PVnchi2",&EvtInfo_PVnchi2);
   root->SetBranchAddress("EvtInfo.PVchi2",&EvtInfo_PVchi2);
   root->SetBranchAddress("EvtInfo.BSx",&EvtInfo_BSx);
   root->SetBranchAddress("EvtInfo.BSy",&EvtInfo_BSy);
   root->SetBranchAddress("EvtInfo.BSz",&EvtInfo_BSz);
   root->SetBranchAddress("EvtInfo.BSxErr",&EvtInfo_BSxErr);
   root->SetBranchAddress("EvtInfo.BSyErr",&EvtInfo_BSyErr);
   root->SetBranchAddress("EvtInfo.BSzErr",&EvtInfo_BSzErr);
   root->SetBranchAddress("EvtInfo.BSdxdz",&EvtInfo_BSdxdz);
   root->SetBranchAddress("EvtInfo.BSdydz",&EvtInfo_BSdydz);
   root->SetBranchAddress("EvtInfo.BSdxdzErr",&EvtInfo_BSdxdzErr);
   root->SetBranchAddress("EvtInfo.BSdydzErr",&EvtInfo_BSdydzErr);
   root->SetBranchAddress("EvtInfo.BSWidthX",&EvtInfo_BSWidthX);
   root->SetBranchAddress("EvtInfo.BSWidthXErr",&EvtInfo_BSWidthXErr);
   root->SetBranchAddress("EvtInfo.BSWidthY",&EvtInfo_BSWidthY);
   root->SetBranchAddress("EvtInfo.BSWidthYErr",&EvtInfo_BSWidthYErr);
   
   root->SetBranchAddress("MuonInfo.size",&MuonInfo_size);
   root->SetBranchAddress("MuonInfo.index",MuonInfo_index);
   root->SetBranchAddress("MuonInfo.handle_index",MuonInfo_handle_index);
   root->SetBranchAddress("MuonInfo.charge",MuonInfo_charge);
   root->SetBranchAddress("MuonInfo.pt",MuonInfo_pt);
   root->SetBranchAddress("MuonInfo.eta",MuonInfo_eta);
   root->SetBranchAddress("MuonInfo.phi",MuonInfo_phi);
   root->SetBranchAddress("MuonInfo.ptErr",MuonInfo_ptErr);
   root->SetBranchAddress("MuonInfo.etaErr",MuonInfo_etaErr);
   root->SetBranchAddress("MuonInfo.phiErr",MuonInfo_phiErr);
   root->SetBranchAddress("MuonInfo.isTrackerMuon",MuonInfo_isTrackerMuon);
   root->SetBranchAddress("MuonInfo.isGlobalMuon",MuonInfo_isGlobalMuon);
   root->SetBranchAddress("MuonInfo.muqual",MuonInfo_muqual);
   root->SetBranchAddress("MuonInfo.iso_trk",MuonInfo_iso_trk);
   root->SetBranchAddress("MuonInfo.iso_ecal",MuonInfo_iso_ecal);
   root->SetBranchAddress("MuonInfo.iso_hcal",MuonInfo_iso_hcal);
   root->SetBranchAddress("MuonInfo.type",MuonInfo_type);
   root->SetBranchAddress("MuonInfo.n_matches",MuonInfo_n_matches);
   root->SetBranchAddress("MuonInfo.TMOneStationTight",MuonInfo_TMOneStationTight);
   root->SetBranchAddress("MuonInfo.TrackerMuonArbitrated",MuonInfo_TrackerMuonArbitrated);
   root->SetBranchAddress("MuonInfo.geninfo_index",MuonInfo_geninfo_index);
   root->SetBranchAddress("MuonInfo.BfinderMuID",MuonInfo_BfinderMuID);
   root->SetBranchAddress("MuonInfo.SoftMuID",MuonInfo_SoftMuID);
   root->SetBranchAddress("MuonInfo.isStandAloneMuon",MuonInfo_isStandAloneMuon);
   root->SetBranchAddress("MuonInfo.StandAloneMuon_charge",MuonInfo_StandAloneMuon_charge);
   root->SetBranchAddress("MuonInfo.StandAloneMuon_pt",MuonInfo_StandAloneMuon_pt);
   root->SetBranchAddress("MuonInfo.StandAloneMuon_eta",MuonInfo_StandAloneMuon_eta);
   root->SetBranchAddress("MuonInfo.StandAloneMuon_phi",MuonInfo_StandAloneMuon_phi);
   root->SetBranchAddress("MuonInfo.StandAloneMuon_d0",MuonInfo_StandAloneMuon_d0);
   root->SetBranchAddress("MuonInfo.StandAloneMuon_dz",MuonInfo_StandAloneMuon_dz);
   root->SetBranchAddress("MuonInfo.StandAloneMuon_dzPV",MuonInfo_StandAloneMuon_dzPV);
   root->SetBranchAddress("MuonInfo.StandAloneMuon_dxyPV",MuonInfo_StandAloneMuon_dxyPV);
   root->SetBranchAddress("MuonInfo.outerTrackisNonnull",MuonInfo_outerTrackisNonnull);
   root->SetBranchAddress("MuonInfo.innerTrackisNonnull",MuonInfo_innerTrackisNonnull);
   root->SetBranchAddress("MuonInfo.globalTrackisNonnull",MuonInfo_globalTrackisNonnull);
   root->SetBranchAddress("MuonInfo.innerTrackQuality",MuonInfo_innerTrackQuality);
   root->SetBranchAddress("MuonInfo.normchi2",MuonInfo_normchi2);
   root->SetBranchAddress("MuonInfo.i_striphit",MuonInfo_i_striphit);
   root->SetBranchAddress("MuonInfo.i_pixelhit",MuonInfo_i_pixelhit);
   root->SetBranchAddress("MuonInfo.i_nStripLayer",MuonInfo_i_nStripLayer);
   root->SetBranchAddress("MuonInfo.i_nPixelLayer",MuonInfo_i_nPixelLayer);
   root->SetBranchAddress("MuonInfo.i_chi2",MuonInfo_i_chi2);
   root->SetBranchAddress("MuonInfo.i_ndf",MuonInfo_i_ndf);
   root->SetBranchAddress("MuonInfo.fpbarrelhit",MuonInfo_fpbarrelhit);
   root->SetBranchAddress("MuonInfo.fpendcaphit",MuonInfo_fpendcaphit);
   root->SetBranchAddress("MuonInfo.d0",MuonInfo_d0);
   root->SetBranchAddress("MuonInfo.dz",MuonInfo_dz);
   root->SetBranchAddress("MuonInfo.dzPV",MuonInfo_dzPV);
   root->SetBranchAddress("MuonInfo.dxyPV",MuonInfo_dxyPV);
   root->SetBranchAddress("MuonInfo.g_chi2",MuonInfo_g_chi2);
   root->SetBranchAddress("MuonInfo.g_ndf",MuonInfo_g_ndf);
   root->SetBranchAddress("MuonInfo.g_striphit",MuonInfo_g_striphit);
   root->SetBranchAddress("MuonInfo.g_pixelhit",MuonInfo_g_pixelhit);
   root->SetBranchAddress("MuonInfo.nmuhit",MuonInfo_nmuhit);
   root->SetBranchAddress("MuonInfo.isTriggered",MuonInfo_isTriggered);
   root->SetBranchAddress("MuonInfo.MuTrgMatchPathSize",&MuonInfo_MuTrgMatchPathSize);

   root->SetBranchAddress("TrackInfo.size",&TrackInfo_size);
   root->SetBranchAddress("TrackInfo.index",TrackInfo_index);
   root->SetBranchAddress("TrackInfo.handle_index",TrackInfo_handle_index);
   root->SetBranchAddress("TrackInfo.charge",TrackInfo_charge);
   root->SetBranchAddress("TrackInfo.pt",TrackInfo_pt);
   root->SetBranchAddress("TrackInfo.eta",TrackInfo_eta);
   root->SetBranchAddress("TrackInfo.phi",TrackInfo_phi);
   root->SetBranchAddress("TrackInfo.ptErr",TrackInfo_ptErr);
   root->SetBranchAddress("TrackInfo.etaErr",TrackInfo_etaErr);
   root->SetBranchAddress("TrackInfo.phiErr",TrackInfo_phiErr);
   root->SetBranchAddress("TrackInfo.striphit",TrackInfo_striphit);
   root->SetBranchAddress("TrackInfo.pixelhit",TrackInfo_pixelhit);
   root->SetBranchAddress("TrackInfo.nStripLayer",TrackInfo_nStripLayer);
   root->SetBranchAddress("TrackInfo.nPixelLayer",TrackInfo_nPixelLayer);
   root->SetBranchAddress("TrackInfo.fpbarrelhit",TrackInfo_fpbarrelhit);
   root->SetBranchAddress("TrackInfo.fpendcaphit",TrackInfo_fpendcaphit);
   root->SetBranchAddress("TrackInfo.chi2",TrackInfo_chi2);
   root->SetBranchAddress("TrackInfo.ndf",TrackInfo_ndf);
   root->SetBranchAddress("TrackInfo.d0",TrackInfo_d0);
   root->SetBranchAddress("TrackInfo.d0error",TrackInfo_d0error);
   root->SetBranchAddress("TrackInfo.dzPV",TrackInfo_dzPV);
   root->SetBranchAddress("TrackInfo.dxyPV",TrackInfo_dxyPV);
   root->SetBranchAddress("TrackInfo.geninfo_index",TrackInfo_geninfo_index);
   root->SetBranchAddress("TrackInfo.trackQuality",TrackInfo_trackQuality);
   root->SetBranchAddress("TrackInfo.highPurity",TrackInfo_highPurity);
   root->SetBranchAddress("TrackInfo.trkMVAVal",TrackInfo_trkMVAVal);
   root->SetBranchAddress("TrackInfo.trkAlgo",TrackInfo_trkAlgo);

   root->SetBranchAddress("BInfo.uj_size",&BInfo_uj_size);
   root->SetBranchAddress("BInfo.uj_index",BInfo_uj_index);
   root->SetBranchAddress("BInfo.uj_mass",BInfo_uj_mass);
   root->SetBranchAddress("BInfo.uj_pt",BInfo_uj_pt);
   root->SetBranchAddress("BInfo.uj_eta",BInfo_uj_eta);
   root->SetBranchAddress("BInfo.uj_phi",BInfo_uj_phi);
   root->SetBranchAddress("BInfo.uj_vtxX",BInfo_uj_vtxX);
   root->SetBranchAddress("BInfo.uj_vtxY",BInfo_uj_vtxY);
   root->SetBranchAddress("BInfo.uj_vtxZ",BInfo_uj_vtxZ);
   root->SetBranchAddress("BInfo.uj_vtxXErr",BInfo_uj_vtxXErr);
   root->SetBranchAddress("BInfo.uj_vtxYErr",BInfo_uj_vtxYErr);
   root->SetBranchAddress("BInfo.uj_vtxZErr",BInfo_uj_vtxZErr);
   root->SetBranchAddress("BInfo.uj_vtxYXErr",BInfo_uj_vtxYXErr);
   root->SetBranchAddress("BInfo.uj_vtxZXErr",BInfo_uj_vtxZXErr);
   root->SetBranchAddress("BInfo.uj_vtxZYErr",BInfo_uj_vtxZYErr);
   root->SetBranchAddress("BInfo.uj_vtxdof",BInfo_uj_vtxdof);
   root->SetBranchAddress("BInfo.uj_vtxchi2",BInfo_uj_vtxchi2);
   root->SetBranchAddress("BInfo.uj_rfmu1_index",BInfo_uj_rfmu1_index);
   root->SetBranchAddress("BInfo.uj_rfmu2_index",BInfo_uj_rfmu2_index);
   root->SetBranchAddress("BInfo.uj_rfmu1_pt",BInfo_uj_rfmu1_pt);
   root->SetBranchAddress("BInfo.uj_rfmu1_eta",BInfo_uj_rfmu1_eta);
   root->SetBranchAddress("BInfo.uj_rfmu1_phi",BInfo_uj_rfmu1_phi);
   root->SetBranchAddress("BInfo.uj_rfmu2_pt",BInfo_uj_rfmu2_pt);
   root->SetBranchAddress("BInfo.uj_rfmu2_eta",BInfo_uj_rfmu2_eta);
   root->SetBranchAddress("BInfo.uj_rfmu2_phi",BInfo_uj_rfmu2_phi);

   root->SetBranchAddress("BInfo.size",&BInfo_size);
   root->SetBranchAddress("BInfo.index",BInfo_index);
   root->SetBranchAddress("BInfo.mass",BInfo_mass);
   root->SetBranchAddress("BInfo.pt",BInfo_pt);
   root->SetBranchAddress("BInfo.eta",BInfo_eta);
   root->SetBranchAddress("BInfo.phi",BInfo_phi);
   root->SetBranchAddress("BInfo.pxE",BInfo_pxE);
   root->SetBranchAddress("BInfo.pyE",BInfo_pyE);
   root->SetBranchAddress("BInfo.pzE",BInfo_pzE);
   root->SetBranchAddress("BInfo.vtxX",BInfo_vtxX);
   root->SetBranchAddress("BInfo.vtxY",BInfo_vtxY);
   root->SetBranchAddress("BInfo.vtxZ",BInfo_vtxZ);
   root->SetBranchAddress("BInfo.vtxXErr",BInfo_vtxXErr);
   root->SetBranchAddress("BInfo.vtxYErr",BInfo_vtxYErr);
   root->SetBranchAddress("BInfo.vtxZErr",BInfo_vtxZErr);
   root->SetBranchAddress("BInfo.vtxYXErr",BInfo_vtxYXErr);
   root->SetBranchAddress("BInfo.vtxZXErr",BInfo_vtxZXErr);
   root->SetBranchAddress("BInfo.vtxZYErr",BInfo_vtxZYErr);
   root->SetBranchAddress("BInfo.vtxdof",BInfo_vtxdof);
   root->SetBranchAddress("BInfo.vtxchi2",BInfo_vtxchi2);
   root->SetBranchAddress("BInfo.rfuj_index",BInfo_rfuj_index);
   root->SetBranchAddress("BInfo.rftk1_index",BInfo_rftk1_index);
   root->SetBranchAddress("BInfo.rftk2_index",BInfo_rftk2_index);
   root->SetBranchAddress("BInfo.type",BInfo_type);
   root->SetBranchAddress("BInfo.rfmu1_pt",BInfo_rfmu1_pt);
   root->SetBranchAddress("BInfo.rfmu1_eta",BInfo_rfmu1_eta);
   root->SetBranchAddress("BInfo.rfmu1_phi",BInfo_rfmu1_phi);
   root->SetBranchAddress("BInfo.rfmu2_pt",BInfo_rfmu2_pt);
   root->SetBranchAddress("BInfo.rfmu2_eta",BInfo_rfmu2_eta);
   root->SetBranchAddress("BInfo.rfmu2_phi",BInfo_rfmu2_phi);
   root->SetBranchAddress("BInfo.rftk1_pt",BInfo_rftk1_pt);
   root->SetBranchAddress("BInfo.rftk1_eta",BInfo_rftk1_eta);
   root->SetBranchAddress("BInfo.rftk1_phi",BInfo_rftk1_phi);
   root->SetBranchAddress("BInfo.rftk2_pt",BInfo_rftk2_pt);
   root->SetBranchAddress("BInfo.rftk2_eta",BInfo_rftk2_eta);
   root->SetBranchAddress("BInfo.rftk2_phi",BInfo_rftk2_phi);
   root->SetBranchAddress("BInfo.tktk_mass",BInfo_tktk_mass);
   root->SetBranchAddress("BInfo.tktk_pt",BInfo_tktk_pt);
   root->SetBranchAddress("BInfo.tktk_eta",BInfo_tktk_eta);
   root->SetBranchAddress("BInfo.tktk_phi",BInfo_tktk_phi);
   root->SetBranchAddress("BInfo.tktk_vtxX",BInfo_tktk_vtxX);
   root->SetBranchAddress("BInfo.tktk_vtxY",BInfo_tktk_vtxY);
   root->SetBranchAddress("BInfo.tktk_vtxZ",BInfo_tktk_vtxZ);
   root->SetBranchAddress("BInfo.tktk_vtxXErr",BInfo_tktk_vtxXErr);
   root->SetBranchAddress("BInfo.tktk_vtxYErr",BInfo_tktk_vtxYErr);
   root->SetBranchAddress("BInfo.tktk_vtxZErr",BInfo_tktk_vtxZErr);
   root->SetBranchAddress("BInfo.tktk_vtxYXErr",BInfo_tktk_vtxYXErr);
   root->SetBranchAddress("BInfo.tktk_vtxZXErr",BInfo_tktk_vtxZXErr);
   root->SetBranchAddress("BInfo.tktk_vtxZYErr",BInfo_tktk_vtxZYErr);
   root->SetBranchAddress("BInfo.tktk_vtxdof",BInfo_tktk_vtxdof);
   root->SetBranchAddress("BInfo.tktk_vtxchi2",BInfo_tktk_vtxchi2);
   root->SetBranchAddress("BInfo.tktk_rftk1_pt",BInfo_tktk_rftk1_pt);
   root->SetBranchAddress("BInfo.tktk_rftk1_eta",BInfo_tktk_rftk1_eta);
   root->SetBranchAddress("BInfo.tktk_rftk1_phi",BInfo_tktk_rftk1_phi);
   root->SetBranchAddress("BInfo.tktk_rftk2_pt",BInfo_tktk_rftk2_pt);
   root->SetBranchAddress("BInfo.tktk_rftk2_eta",BInfo_tktk_rftk2_eta);
   root->SetBranchAddress("BInfo.tktk_rftk2_phi",BInfo_tktk_rftk2_phi);

   root->SetBranchAddress("GenInfo.size",&GenInfo_size);
   root->SetBranchAddress("GenInfo.index",&GenInfo_index);
   root->SetBranchAddress("GenInfo.handle_index",&GenInfo_handle_index);
   root->SetBranchAddress("GenInfo.pt",&GenInfo_pt);
   root->SetBranchAddress("GenInfo.eta",&GenInfo_eta);
   root->SetBranchAddress("GenInfo.phi",&GenInfo_phi);
   root->SetBranchAddress("GenInfo.mass",&GenInfo_mass);
   root->SetBranchAddress("GenInfo.pdgId",&GenInfo_pdgId);
   root->SetBranchAddress("GenInfo.status",&GenInfo_status);
   root->SetBranchAddress("GenInfo.nMo",&GenInfo_nMo);
   root->SetBranchAddress("GenInfo.nDa",&GenInfo_nDa);
   root->SetBranchAddress("GenInfo.mo1",&GenInfo_mo1);
   root->SetBranchAddress("GenInfo.mo2",&GenInfo_mo2);
   root->SetBranchAddress("GenInfo.da1",&GenInfo_da1);
   root->SetBranchAddress("GenInfo.da2",&GenInfo_da2);
   root->SetBranchAddress("GenInfo.da3",&GenInfo_da3);
   root->SetBranchAddress("GenInfo.da4",&GenInfo_da4);
}

//HltInfo
Int_t           Bf_HLT_Run;
ULong64_t       Bf_HLT_Event;
Int_t           Bf_HLT_LumiBlock;
void setHltBranch(TTree* hltroot)
{
  hltroot->SetBranchAddress("Run",&Bf_HLT_Run);
  hltroot->SetBranchAddress("Event",&Bf_HLT_Event);
  hltroot->SetBranchAddress("LumiBlock",&Bf_HLT_LumiBlock);
}

//hiEvtInfo
unsigned int       Bf_HiTree_Run;
unsigned long long Bf_HiTree_Evt;
unsigned int       Bf_HiTree_Lumi;
void setHiTreeBranch(TTree* hitreeroot)
{
  hitreeroot->SetBranchAddress("run",&Bf_HiTree_Run);
  hitreeroot->SetBranchAddress("evt",&Bf_HiTree_Evt);
  hitreeroot->SetBranchAddress("lumi",&Bf_HiTree_Lumi);
}

void readBranch(TTree* nt)
{
  //EvtInfo
  nt->SetBranchAddress("RunNo",&RunNo);
  nt->SetBranchAddress("EvtNo",&EvtNo);
  nt->SetBranchAddress("LumiNo",&LumiNo);
  nt->SetBranchAddress("Bsize",&Bsize);
  nt->SetBranchAddress("PVx",&PVx);
  nt->SetBranchAddress("PVy",&PVy);
  nt->SetBranchAddress("PVz",&PVz);
  nt->SetBranchAddress("PVxE",&PVxE);
  nt->SetBranchAddress("PVyE",&PVyE);
  nt->SetBranchAddress("PVzE",&PVzE);
  nt->SetBranchAddress("PVnchi2",&PVnchi2);
  nt->SetBranchAddress("PVchi2",&PVchi2);
  nt->SetBranchAddress("BSx",&BSx);
  nt->SetBranchAddress("BSy",&BSy);
  nt->SetBranchAddress("BSz",&BSz);
  nt->SetBranchAddress("BSxErr",&BSxErr);
  nt->SetBranchAddress("BSyErr",&BSyErr);
  nt->SetBranchAddress("BSzErr",&BSzErr);
  nt->SetBranchAddress("BSdxdz",&BSdxdz);
  nt->SetBranchAddress("BSdydz",&BSdydz);
  nt->SetBranchAddress("BSdxdzErr",&BSdxdzErr);
  nt->SetBranchAddress("BSdydzErr",&BSdydzErr);
  nt->SetBranchAddress("BSWidthX",&BSWidthX);
  nt->SetBranchAddress("BSWidthXErr",&BSWidthXErr);
  nt->SetBranchAddress("BSWidthY",&BSWidthY);
  nt->SetBranchAddress("BSWidthYErr",&BSWidthYErr);
  nt->SetBranchAddress("pthatweight",&pthatweight);
  nt->SetBranchAddress("hiBin",&hiBin);
//BInfo
  nt->SetBranchAddress("Bindex",Bindex);
  nt->SetBranchAddress("Btype",Btype);
  nt->SetBranchAddress("Bmass",Bmass);
  nt->SetBranchAddress("Bpt",Bpt);
  nt->SetBranchAddress("Beta",Beta);
  nt->SetBranchAddress("Bphi",Bphi);
  nt->SetBranchAddress("By",By);
  nt->SetBranchAddress("Balpha",Balpha);
  nt->SetBranchAddress("BvtxX",BvtxX);
  nt->SetBranchAddress("BvtxY",BvtxY);
  nt->SetBranchAddress("Bd0",Bd0);
  nt->SetBranchAddress("Bd0Err",Bd0Err);
  nt->SetBranchAddress("Bdxyz",Bdxyz);
  nt->SetBranchAddress("BdxyzErr",BdxyzErr);
  nt->SetBranchAddress("Bchi2ndf",Bchi2ndf);
  nt->SetBranchAddress("Bchi2cl",Bchi2cl);
  nt->SetBranchAddress("Bdtheta",Bdtheta);
  nt->SetBranchAddress("Blxy",Blxy);
  nt->SetBranchAddress("BlxyBS",BlxyBS);
  nt->SetBranchAddress("BlxyBSErr",BlxyBSErr);
  nt->SetBranchAddress("Bmaxpt",Bmaxpt);
  nt->SetBranchAddress("Bmaxprob",Bmaxprob);
  nt->SetBranchAddress("Bbesttktkmass",Bbesttktkmass);
  nt->SetBranchAddress("BmaxptMatched",BmaxptMatched);
  nt->SetBranchAddress("BmaxprobMatched",BmaxprobMatched);
  nt->SetBranchAddress("BbesttktkmassMatched",BbesttktkmassMatched);

  //BInfo.trkInfo
  nt->SetBranchAddress("Btrk1Idx",Btrk1Idx);
  nt->SetBranchAddress("Btrk2Idx",Btrk2Idx);
  nt->SetBranchAddress("Btrk1Pt",Btrk1Pt);
  nt->SetBranchAddress("Btrk2Pt",Btrk2Pt);
  nt->SetBranchAddress("Btrk1Eta",Btrk1Eta);  
  nt->SetBranchAddress("Btrk2Eta",Btrk2Eta);  
  nt->SetBranchAddress("Btrk1Phi",Btrk1Phi);  
  nt->SetBranchAddress("Btrk2Phi",Btrk2Phi);  
  nt->SetBranchAddress("Btrk1PtErr",Btrk1PtErr); 
  nt->SetBranchAddress("Btrk2PtErr",Btrk2PtErr);
  nt->SetBranchAddress("Btrk1EtaErr",Btrk1EtaErr);
  nt->SetBranchAddress("Btrk2EtaErr",Btrk2EtaErr);
  nt->SetBranchAddress("Btrk1PhiErr",Btrk1PhiErr);
  nt->SetBranchAddress("Btrk2PhiErr",Btrk2PhiErr);
  nt->SetBranchAddress("Btrk1Y",Btrk1Y);  
  nt->SetBranchAddress("Btrk2Y",Btrk2Y);  
  nt->SetBranchAddress("Btrk1D0", Btrk1D0);
  nt->SetBranchAddress("Btrk2D0", Btrk2D0);
  nt->SetBranchAddress("Btrk1Dz", Btrk1Dz);
  nt->SetBranchAddress("Btrk2Dz", Btrk2Dz);
  nt->SetBranchAddress("Btrk1Dxy",Btrk1Dxy);
  nt->SetBranchAddress("Btrk2Dxy",Btrk2Dxy);
  nt->SetBranchAddress("Btrk1D0Err",Btrk1D0Err);
  nt->SetBranchAddress("Btrk2D0Err",Btrk2D0Err);
  nt->SetBranchAddress("Btrk1PixelHit",Btrk1PixelHit);
  nt->SetBranchAddress("Btrk2PixelHit",Btrk2PixelHit);
  nt->SetBranchAddress("Btrk1StripHit",Btrk1StripHit);
  nt->SetBranchAddress("Btrk2StripHit",Btrk2StripHit);
  nt->SetBranchAddress("Btrk1nPixelLayer",Btrk1nPixelLayer);
  nt->SetBranchAddress("Btrk2nPixelLayer",Btrk2nPixelLayer);
  nt->SetBranchAddress("Btrk1nStripLayer",Btrk1nStripLayer);
  nt->SetBranchAddress("Btrk2nStripLayer",Btrk2nStripLayer);
  nt->SetBranchAddress("Btrk1Chi2ndf",Btrk1Chi2ndf);
  nt->SetBranchAddress("Btrk2Chi2ndf",Btrk2Chi2ndf);
  nt->SetBranchAddress("Btrk1MVAVal",Btrk1MVAVal);
  nt->SetBranchAddress("Btrk2MVAVal",Btrk2MVAVal);
  nt->SetBranchAddress("Btrk1Algo",Btrk1Algo);
  nt->SetBranchAddress("Btrk2Algo",Btrk2Algo);
  nt->SetBranchAddress("Btrk1highPurity",Btrk1highPurity);
  nt->SetBranchAddress("Btrk2highPurity",Btrk2highPurity);
  nt->SetBranchAddress("Btrk1Quality",Btrk1Quality);
  nt->SetBranchAddress("Btrk2Quality",Btrk2Quality);
  nt->SetBranchAddress("pclusterCompatibilityFilter", pclusterCompatibilityFilter); 

  nt->SetBranchAddress("pprimaryVertexFilter", pprimaryVertexFilter); 

  //BInfo.tktkInfo
  nt->SetBranchAddress("Btktkmass",Btktkmass);
  nt->SetBranchAddress("BtktkmassKK",BtktkmassKK);
  nt->SetBranchAddress("BtktkvProb",BtktkvProb);
  nt->SetBranchAddress("Btktkpt",Btktkpt);
  nt->SetBranchAddress("Btktketa",Btktketa);
  nt->SetBranchAddress("Btktkphi",Btktkphi);
  nt->SetBranchAddress("Btktky",Btktky);
  nt->SetBranchAddress("Bdoubletmass",Bdoubletmass);
  nt->SetBranchAddress("Bdoubletpt",Bdoubletpt);
  nt->SetBranchAddress("Bdoubleteta",Bdoubleteta);  
  nt->SetBranchAddress("Bdoubletphi",Bdoubletphi);  
  nt->SetBranchAddress("Bdoublety",Bdoublety);
  
  //BInfo.muonInfo
  nt->SetBranchAddress("Bmu1pt",Bmu1pt);
  nt->SetBranchAddress("Bmu2pt",Bmu2pt);
  nt->SetBranchAddress("Bmu1p",Bmu1p);
  nt->SetBranchAddress("Bmu2p",Bmu2p);
  nt->SetBranchAddress("Bmu1eta",Bmu1eta);
  nt->SetBranchAddress("Bmu2eta",Bmu2eta);
  nt->SetBranchAddress("Bmu1phi",Bmu1phi);
  nt->SetBranchAddress("Bmu2phi",Bmu2phi);
  nt->SetBranchAddress("Bmu1y",Bmu1y);
  nt->SetBranchAddress("Bmu2y",Bmu2y);
  nt->SetBranchAddress("Bmu1dzPV",Bmu1dzPV);
  nt->SetBranchAddress("Bmu2dzPV",Bmu2dzPV);
  nt->SetBranchAddress("Bmu1dxyPV",Bmu1dxyPV);
  nt->SetBranchAddress("Bmu2dxyPV",Bmu2dxyPV);
  nt->SetBranchAddress("Bmu1normchi2",Bmu1normchi2);
  nt->SetBranchAddress("Bmu2normchi2",Bmu2normchi2);
  nt->SetBranchAddress("Bmu1Chi2ndf",Bmu1Chi2ndf);
  nt->SetBranchAddress("Bmu2Chi2ndf",Bmu2Chi2ndf);
  nt->SetBranchAddress("Bmu1muqual",Bmu1muqual);
  nt->SetBranchAddress("Bmu1muqual",Bmu1muqual);
  nt->SetBranchAddress("Bmu1TrackerMuArbitrated",Bmu1TrackerMuArbitrated);
  nt->SetBranchAddress("Bmu2TrackerMuArbitrated",Bmu2TrackerMuArbitrated);
  nt->SetBranchAddress("Bmu1isTrackerMuon",Bmu1isTrackerMuon);
  nt->SetBranchAddress("Bmu2isTrackerMuon",Bmu2isTrackerMuon);
  nt->SetBranchAddress("Bmu1isGlobalMuon",Bmu1isGlobalMuon);
  nt->SetBranchAddress("Bmu2isGlobalMuon",Bmu2isGlobalMuon);
  nt->SetBranchAddress("Bmu1TMOneStationTight",Bmu1TMOneStationTight);
  nt->SetBranchAddress("Bmu2TMOneStationTight",Bmu2TMOneStationTight);
  nt->SetBranchAddress("Bmu1InPixelLayer",Bmu1InPixelLayer);
  nt->SetBranchAddress("Bmu2InPixelLayer",Bmu2InPixelLayer);
  nt->SetBranchAddress("Bmu1InStripLayer",Bmu1InStripLayer);
  nt->SetBranchAddress("Bmu2InStripLayer",Bmu2InStripLayer);
  nt->SetBranchAddress("Bmu1InTrackerLayer",Bmu1InTrackerLayer);
  nt->SetBranchAddress("Bmu2InTrackerLayer",Bmu2InTrackerLayer);

  nt->SetBranchAddress("Bmumumass",Bmumumass);
  nt->SetBranchAddress("Bmumueta",Bmumueta);
  nt->SetBranchAddress("Bmumuphi",Bmumuphi);
  nt->SetBranchAddress("Bmumuy",Bmumuy);
  nt->SetBranchAddress("Bmumupt",Bmumupt);
  nt->SetBranchAddress("Bujmass",Bujmass);
  nt->SetBranchAddress("BujvProb",BujvProb);
  nt->SetBranchAddress("Bujpt",Bujpt);
  nt->SetBranchAddress("Bujeta",Bujeta);
  nt->SetBranchAddress("Bujphi",Bujphi);
  nt->SetBranchAddress("Bujy",Bujy);
  nt->SetBranchAddress("Bujlxy",Bujlxy);

  //BInfo.genInfo
  nt->SetBranchAddress("Bgen",Bgen);
  nt->SetBranchAddress("BgenIndex",BgenIndex);
  nt->SetBranchAddress("Bgenpt",Bgenpt);
  nt->SetBranchAddress("Bgeny",Bgeny);
  nt->SetBranchAddress("Bgeneta",Bgeneta);
  nt->SetBranchAddress("Bgenphi",Bgenphi);

  
  //nt->SetBranchAddress("Balpha",Balpha,"Balpha[Bsize]/D");
  nt->SetBranchAddress("BsvpvDistance",BsvpvDistance);
  nt->SetBranchAddress("BsvpvDisErr",BsvpvDisErr);
  nt->SetBranchAddress("BsvpvDistance_2D",BsvpvDistance_2D);
  nt->SetBranchAddress("BsvpvDisErr_2D",BsvpvDisErr_2D);
  /*nt->SetBranchAddress("BMaxDoca",BMaxDoca,"BMaxDoca[Bsize]/D");
  nt->SetBranchAddress("Btrk1MassHypo",Btrk1MassHypo,"Btrk1MassHypo[Bsize]/D");
  nt->SetBranchAddress("Btrk2MassHypo",Btrk2MassHypo,"Btrk2MassHypo[Bsize]/D");
  nt->SetBranchAddress("Btrkminpt",Btrkminpt,"Btrkminpt[Bsize]/D");
  nt->SetBranchAddress("Btrkmaxpt",Btrkmaxpt,"Btrkmaxpt[Bsize]/D");
  nt->SetBranchAddress("Btrkminptindex",Btrkminptindex,"Btrkminptindex[Bsize]/I");
  nt->SetBranchAddress("Btrkmaxptindex",Btrkmaxptindex,"Btrkmaxptindex[Bsize]/I");
  */

}

void readGenBranch(TTree* nt)
{
  nt->SetBranchAddress("Gsize",&Gsize);
  nt->SetBranchAddress("Gy",Gy);
  nt->SetBranchAddress("Geta",Geta);
  nt->SetBranchAddress("Gphi",Gphi);
  nt->SetBranchAddress("Gpt",Gpt);
  nt->SetBranchAddress("GpdgId",GpdgId);
  nt->SetBranchAddress("GisSignal",GisSignal);
  nt->SetBranchAddress("Gmu1eta",Gmu1eta);
  nt->SetBranchAddress("Gmu1phi",Gmu1phi);
  nt->SetBranchAddress("Gmu1pt",Gmu1pt);
  nt->SetBranchAddress("Gmu1p",Gmu1p);
  nt->SetBranchAddress("Gmu2eta",Gmu2eta);
  nt->SetBranchAddress("Gmu2phi",Gmu2phi);
  nt->SetBranchAddress("Gmu2pt",Gmu2pt);
  nt->SetBranchAddress("Gmu2p",Gmu2p);
  nt->SetBranchAddress("Gtk1pt",Gtk1pt);
  nt->SetBranchAddress("Gtk1eta",Gtk1eta);
  nt->SetBranchAddress("Gtk1phi",Gtk1phi);
  nt->SetBranchAddress("Gtk2pt",Gtk2pt);
  nt->SetBranchAddress("Gtk2eta",Gtk2eta);
  nt->SetBranchAddress("Gtk2phi",Gtk2phi);
}


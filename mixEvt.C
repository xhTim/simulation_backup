#include "../common/headers.h"
#include "../common/funUtil.h"
#include "../common/VertexCompositeTree.h"

TH3D *hMvsPtvsRap;
TH3D *hMvsPtvsRap_NeuDir[nNeus][nNeus];
TH3D *hMvsPtvsRap_NeuDir_gg[nNeus][nNeus];

const int N_gen = 50000;
const double L_int = 7 * pow(10, 6); //mb^{-1}

const double sec_CohJpsi_0n0n = 1.796; // mb
const double sec_CohJpsi_0nXn = 0.460587;
const double sec_CohJpsi_XnXn = 0.180980;

const double sec_gg_0n0n = 12.243;
const double sec_gg_0nXn = 1.981;
const double sec_gg_XnXn = 0.562020;

Bool_t   goodMuPair(const TLorentzVector posFourMom, const Bool_t isPosTrig, const TLorentzVector negFourMom, const Bool_t isNegTrig);
Bool_t   goodMuPair(VertexCompositeTree& evtTree, const int icand);

void bookHistos();
void resetError();
void writeHistos(TString fileName = "test");

const std::vector<std::string>& generateFilenames(const std::string& prefix, int num=50) {
    static std::vector<std::string> filenames;

    filenames.clear();

    for (int i = 1; i <= num; ++i) {
        std::stringstream ss;
		ss << prefix << "dimuana_mc_" << i << ".root";
        filenames.push_back(ss.str());
    }

    return filenames;
}

void fillEvt(TH3D* hMvsPtvsRap_NeuCl, const std::string DirName, double sec, int N_generated=N_gen){

    const std::vector<std::string> &Filenames = generateFilenames(DirName);
    const auto& csTreeDir = "dimucontana_mc";

    // Extract the tree
	VertexCompositeTree csTree;
	if (!csTree.GetTree(Filenames, csTreeDir)) 
	{
		cout << "Invalid Correct-Sign tree!" << endl;
		return;
	}

	if(!init()) 
	{
		cout<<"Initialization failed !"<<endl;
		return;
	}

    for(Long64_t jentry = 1; jentry < csTree.GetEntries(); jentry++){

        if (jentry % (csTree.GetEntries() / 10) == 0)
			cout << "begin " << jentry << "th entry...." << endl;

		// Get the entry
		if (csTree.GetEntry(jentry) < 0) 
		{
			cout << "Invalid correct-sign entry!" << endl;
			return;
		}

        Int_t   nTrkHP       = csTree.NtrkHP();
        Bool_t passEvtSel = (nTrkHP == 2);
		if(!passEvtSel) continue;

        // Loop over the correct-sign candidate pairs
        for(UInt_t icand = 0; icand < csTree.candSize(); icand++){

            Float_t pt   = csTree.pT()[icand];
			Float_t mass = csTree.mass()[icand];
			Float_t y    = csTree.y()[icand];

			Double_t posPt     = csTree.chargeD1()[icand] > 0 ? csTree.pTD1()[icand]  : csTree.pTD2()[icand]; 
			Double_t posEta    = csTree.chargeD1()[icand] > 0 ? csTree.EtaD1()[icand] : csTree.EtaD2()[icand]; 
			Double_t posPhi    = csTree.chargeD1()[icand] > 0 ? csTree.PhiD1()[icand] : csTree.PhiD2()[icand]; 
			Bool_t   isPosTrig = csTree.chargeD1()[icand] > 0 ? csTree.trigMuon1()[trigIdx][icand] : csTree.trigMuon2()[trigIdx][icand];
			Double_t negPt     = csTree.chargeD1()[icand] < 0 ? csTree.pTD1()[icand]  : csTree.pTD2()[icand]; 
			Double_t negEta    = csTree.chargeD1()[icand] < 0 ? csTree.EtaD1()[icand] : csTree.EtaD2()[icand]; 
			Double_t negPhi    = csTree.chargeD1()[icand] < 0 ? csTree.PhiD1()[icand] : csTree.PhiD2()[icand]; 
			Bool_t   isNegTrig = csTree.chargeD1()[icand] < 0 ? csTree.trigMuon1()[trigIdx][icand] : csTree.trigMuon2()[trigIdx][icand];

            // soft muon selection
//			if( !csTree.softCand(icand) )   continue;

            // apply muon acceptance and trigger selections
			if( !goodMuPair(csTree, icand) ) continue;

            hMvsPtvsRap_NeuCl->Fill(y, pt, mass);
        }
    }

    double scalefactor = sec * L_int / N_generated;

    hMvsPtvsRap_NeuCl->Scale(scalefactor);

}

void mixEvt(){

    //-------------------------------
	bookHistos();
	//-------------------------------

    fillEvt(hMvsPtvsRap_NeuDir[0][0], "../rootfiles/CohJpsi_0n0n/", sec_CohJpsi_0n0n);
    fillEvt(hMvsPtvsRap_NeuDir[0][1], "../rootfiles/CohJpsi_0nXn/", sec_CohJpsi_0nXn);
    hMvsPtvsRap_NeuDir[0][1]->Scale(0.5);
    hMvsPtvsRap_NeuDir[1][0] = (TH3D *)hMvsPtvsRap_NeuDir[0][1]->Clone("hMvsPtvsRap_NeuDir1p0m");
    fillEvt(hMvsPtvsRap_NeuDir[1][1], "../rootfiles/CohJpsi_XnXn/", sec_CohJpsi_XnXn);

    fillEvt(hMvsPtvsRap_NeuDir_gg[0][0], "../rootfiles/Cohgg_0n0n/", sec_gg_0n0n);
    fillEvt(hMvsPtvsRap_NeuDir_gg[0][1], "../rootfiles/Cohgg_0nXn/", sec_gg_0nXn);
    hMvsPtvsRap_NeuDir_gg[0][1]->Scale(0.5);
    hMvsPtvsRap_NeuDir_gg[1][0] = (TH3D *)hMvsPtvsRap_NeuDir_gg[0][1]->Clone("hMvsPtvsRap_NeuDir1p0m_gg");
    fillEvt(hMvsPtvsRap_NeuDir_gg[1][1], "../rootfiles/Cohgg_XnXn/", sec_gg_XnXn);

    hMvsPtvsRap = (TH3D *)hMvsPtvsRap_NeuDir[0][0]->Clone("hMvsPtvsRap");
	hMvsPtvsRap->Add(hMvsPtvsRap_NeuDir_gg[0][0]);
	for (int ip = 0; ip < nNeus; ip++)
	{
		for (int im = 0; im < nNeus; im++)
		{
			hMvsPtvsRap_NeuDir[ip][im]->Add(hMvsPtvsRap_NeuDir_gg[ip][im]);
            if(!(ip==0&&im==0))
                hMvsPtvsRap->Add(hMvsPtvsRap_NeuDir[ip][im]);
		}
	}

	resetError();

	TString dirName = "jpsiHistos";
	system(Form("mkdir -p %s", dirName.Data()));

    TString fileName = "rawSig";

    cout<<"fileName: "<<fileName<<endl;

	writeHistos(Form("%s/%s", dirName.Data(), fileName.Data()));
}

void bookHistos(){

    const Int_t    mHistRapBins    = 60;
	const Double_t mHistRapLow     = -3;
	const Double_t mHistRapHi      = 3;

    const Int_t    mHistPtBins     = 600;
	const Double_t mHistPtLow      = 0;
	const Double_t mHistPtHi       = 6;

	const Int_t    mHistMassBins   = 300;
	const Double_t mHistMassLow    = 2;
	const Double_t mHistMassHi     = 5;

    hMvsPtvsRap = new TH3D("hMvsPtvsRap", "hMvsPtvsRap; Rapidity; p_{T} (GeV); M_{#mu#mu} (GeV)", mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);

    for(Int_t ip=0; ip<nNeus; ip++)
	{
		for(Int_t im=0; im<nNeus; im++)
		{
            hMvsPtvsRap_NeuDir[ip][im] = new TH3D(Form("hMvsPtvsRap_NeuDir%dp%dm", ip, im), "; Rapidity; p_{T} (GeV); M_{#mu#mu} (GeV)", mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
            hMvsPtvsRap_NeuDir_gg[ip][im] = new TH3D(Form("hMvsPtvsRap_NeuDir%dp%dm_gg", ip, im), "; Rapidity; p_{T} (GeV); M_{#mu#mu} (GeV)", mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
        }
    }
}

void resetError(){

	for (Int_t ibinx = 1; ibinx <= hMvsPtvsRap->GetNbinsX(); ibinx++)
	{
		for (Int_t ibiny = 1; ibiny <= hMvsPtvsRap->GetNbinsY(); ibiny++)
		{
			for (Int_t ibinz = 1; ibinz <= hMvsPtvsRap->GetNbinsZ(); ibinz++)
			{
				hMvsPtvsRap->SetBinError(ibinx, ibiny, ibinz, sqrt(hMvsPtvsRap->GetBinContent(ibinx, ibiny, ibinz)));
			}
		}
	}

	for (Int_t ip = 0; ip < nNeus; ip++)
	{
		for (Int_t im = 0; im < nNeus; im++)
		{
			for (Int_t ibinx = 1; ibinx <= hMvsPtvsRap_NeuDir[ip][im]->GetNbinsX(); ibinx++)
			{
				for (Int_t ibiny = 1; ibiny <= hMvsPtvsRap_NeuDir[ip][im]->GetNbinsY(); ibiny++)
				{
					for (Int_t ibinz = 1; ibinz <= hMvsPtvsRap_NeuDir[ip][im]->GetNbinsZ(); ibinz++)
					{
						hMvsPtvsRap_NeuDir[ip][im]->SetBinError(ibinx, ibiny, ibinz, sqrt(hMvsPtvsRap_NeuDir[ip][im]->GetBinContent(ibinx, ibiny, ibinz)));
					}
				}
			}
		}
	}
}

void writeHistos(TString fileName){

    TFile* fOut = new TFile(Form("%s.root", fileName.Data()), "recreate");

	fOut ->cd();

    hMvsPtvsRap->Write();

    for(Int_t ip=0; ip<nNeus; ip++)
	{
		for(Int_t im=0; im<nNeus; im++)
		{
			hMvsPtvsRap_NeuDir[ip][im]->Write();
			hMvsPtvsRap_NeuDir_gg[ip][im]->Write();
		}
	}

	fOut->Close();
}

Bool_t goodMuPair(VertexCompositeTree& evtTree, const int icand)
{
	Double_t mTrkPtTh1  = fTrkAcc  ->Eval(evtTree.EtaD1()[icand]);
	Double_t mTrkPtTh2  = fTrkAcc  ->Eval(evtTree.EtaD2()[icand]);
	Double_t mTrigPtTh1 = fTrigAcc ->Eval(evtTree.EtaD1()[icand]);
	Double_t mTrigPtTh2 = fTrigAcc ->Eval(evtTree.EtaD2()[icand]);

	if(evtTree.pTD1()[icand] < mTrkPtTh1 || evtTree.pTD2()[icand] < mTrkPtTh2) return kFALSE;

	Bool_t isTrigAcc1 = kFALSE, isTrigAcc2 = kFALSE;
	if(evtTree.trigMuon1()[trigIdx][icand] && evtTree.pTD1()[icand] >= mTrigPtTh1) isTrigAcc1 = kTRUE;
	if(evtTree.trigMuon2()[trigIdx][icand] && evtTree.pTD2()[icand] >= mTrigPtTh2) isTrigAcc2 = kTRUE;

//	if(!isTrigAcc1 && !isTrigAcc2) return kFALSE;

	return kTRUE;
}

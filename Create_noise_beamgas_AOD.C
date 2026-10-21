#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Background analysis (noise+beam-gas) using AO2D's
//-------------------------------------------------------

void Create_scan_noise_beamgas(Int_t Fill, Int_t scan_type, Int_t scan, const char *rate_name, const char *sep_type, const char *intensity_type)
{
    gStyle->SetOptFit(1);
    // first get the files and trees
    // Rate of ALL bc (selections applied later!)
	// --> create rate file name
    char *rate_file_name = new char[kg_string_size];
	if (scan_type == 1) sprintf(rate_file_name,
			"../Fill-%d/AllBcRawRate_%s_x_Scan_%d.root", g_vdm_Fill,rate_name,scan);
	if (scan_type == 2) sprintf(rate_file_name,
			"../Fill-%d/AllBcRawRate_%s_y_Scan_%d.root", g_vdm_Fill,rate_name,scan);

    // --> create separation file name
	char *sep_file_name = new char[kg_string_size];
	if (scan_type == 1) sprintf(sep_file_name,"../Fill-%d/%sSep_x_Scan_%d.root",g_vdm_Fill,sep_type,scan);
	if (scan_type == 2) sprintf(sep_file_name,"../Fill-%d/%sSep_y_Scan_%d.root",g_vdm_Fill,sep_type,scan);

    // --> create intensity file name
	char *intensity_file_name = new char[kg_string_size];
	sprintf(intensity_file_name,"../Fill-%d/%s_Scan_%d.root",g_vdm_Fill,intensity_type,scan);

	// --> open files
	TFile *rate_file = new TFile(rate_file_name);
	TFile *sep_file = new TFile(sep_file_name);
    TFile *intensity_file = new TFile(intensity_file_name);

    // --> get the trees
	TTree *rate_tree = (TTree *) rate_file->Get("Rate");
	TTree *sep_tree = (TTree *) sep_file->Get("Separations");
    TTree *intensity_tree = (TTree *) intensity_file->Get("Bunch_Intensity_typeA");

    // Next step: prepare variables to be fitted
	// --> get number of separations
	Int_t n_sep = FindNumberSeparations(scan_type, scan);

    // --> reserve space for rates and separations
	Double_t *rate = new Double_t[n_sep];
	Double_t *rate_error = new Double_t[n_sep];  
	Double_t *sep = new Double_t[n_sep];

    // --> set branches
	rate_tree->ResetBranchAddresses();
	rate_tree->SetBranchAddress("rate",rate);
	rate_tree->SetBranchAddress("rate_error",rate_error);
	sep_tree->ResetBranchAddresses();
	sep_tree->SetBranchAddress("separation",sep);

    // get bunch crossing information & bucket info
	// -- number of bc for each type
    Int_t nIBC = GetNumberInteractingBunchCrossings(); //number of interacting bc
	Int_t nABC = GetNumberBunchCrossingsTypeA(); //number of non-interacting type-A bc
	Int_t nCBC = GetNumberBunchCrossingsTypeC(); //number of non-interacting type-C bc
	Int_t nEBC; //number of empty BC
	if(strcmp(rate_name, "aodzn") == 0) nEBC = 3564 - 2*nIBC - nABC - nCBC;
	if(strcmp(rate_name, "aodft0vtx") == 0) nEBC = 3564 - nIBC - nABC - nCBC;
	//
	Int_t *bunches_I = new Int_t [nIBC];
	Int_t *bunches_A = new Int_t [nABC];
	Int_t *bunches_C = new Int_t [nCBC];
	Int_t *bunches_E = new Int_t [nEBC];

    //get bunch indices for all types
	GetBunchIndices(bunches_I);
	GetBunchIndices_typeA(bunches_A);
	GetBunchIndices_typeC(bunches_C);

    // Next step: get intensity related information for type_A
	// --> get intensity information (only 1 bunch is full in this case)
    Double_t *bunch_intensity_1 = new Double_t[nABC];
    Double_t cf_dcct_1;
    intensity_tree->ResetBranchAddresses();
    intensity_tree->SetBranchAddress("cf_dcct_1",&cf_dcct_1);
	intensity_tree->SetBranchAddress("bunch_intensity_1",bunch_intensity_1);

    // Now reserve space for output
    //Double_t noise;
    //Double_t noise_err;
    Double_t *noise = new Double_t[n_sep];
    Double_t *noise_err = new Double_t[n_sep];
	Double_t beam_gas_A;
	Double_t beam_gas_A_err;

	// create file and tree for output
    // --> file
	char *h_file_name = new char[kg_string_size];
	if (scan_type == 1) sprintf(h_file_name,"../Fill-%d/NoiseBeamgas_%s_%sSep_x_Scan_%d.root",
			g_vdm_Fill,rate_name,sep_type,scan);
	if (scan_type == 2) sprintf(h_file_name,"../Fill-%d/NoiseBeamgas_%s_%sSep_y_Scan_%d.root",
			g_vdm_Fill,rate_name,sep_type,scan);  
    cout<<"FILE NAME"<<h_file_name;
	TFile *hFile = new TFile(h_file_name,"recreate");
    // --> tree
	TTree *h_tree = new TTree("NoiseBeamgasPar","NoiseBeamgasPar");
    //h_tree->Branch("noise",&noise,"noise/D");
    //h_tree->Branch("noise_err",&noise_err,"noise_err/D");
    //Consider noise vs sep as an array; same values for all BC indices
    char txt_tmp[kg_string_size];
    sprintf(txt_tmp,"noise[%d]/D",n_sep);
    h_tree->Branch("noise",noise,txt_tmp);
    sprintf(txt_tmp,"noise_err[%d]/D",n_sep);
    h_tree->Branch("noise_err",noise_err,txt_tmp);
    h_tree->Branch("beam_gas_A",&beam_gas_A,"beam_gas_A/D");
    h_tree->Branch("beam_gas_A_err",&beam_gas_A_err,"beam_gas_A_err/D");
    
    //Graph RA vs NA
	TGraphErrors* grRvsN = new TGraphErrors();

    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// NOISE ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    Int_t index_E = 0; //index of empty bc
    //Define the total rate for each scan; it is an array, dimension = # of separations
	Double_t *RateAvg        = new Double_t[n_sep];
	Double_t *RateAvg_error  = new Double_t[n_sep];
    Double_t noisePol0;
	//Initialize
	for(Int_t i=0; i<n_sep; i++){
		RateAvg[i]=0;
		RateAvg_error[i]=0;
	}
	Double_t *SepTot = new Double_t[n_sep];
    for (Int_t bc=0; bc<3564; bc++) // loop over bunches
    {
        //---------------------
        //select only empty BC:
        bool skipBc_notE = false;
		for(Int_t k=0; k<nIBC; k++) if(bc == bunches_I[k]){skipBc_notE = true; break;}  //skip bc if interacting
		for(Int_t k=0; k<nABC; k++) if(bc == bunches_A[k]){skipBc_notE = true; break;}  //skip bc if type A
		for(Int_t k=0; k<nCBC; k++) if(strcmp(rate_name, "aodft0vtx") == 0 && bc == bunches_C[k]){skipBc_notE = true; break;}  //skip bc if type C
		for(Int_t k=0; k<nCBC; k++) if(strcmp(rate_name, "aodzn") == 0 && bc == bunches_C[k]-31){skipBc_notE = true; break;}  //skip bc if it's 31 bc before type C
		for(Int_t k=0; k<nIBC; k++) if(strcmp(rate_name, "aodzn") == 0 && bc == bunches_I[k]-31){skipBc_notE = true; break;}  //skip bc if it's 31 bc before an interacting
		if(skipBc_notE) continue;
        bunches_E[index_E] = bc;

        rate_tree->GetEntry(bc);
		sep_tree->GetEntry(bc);

        for (Int_t i=0; i<n_sep; i++) // loop over steps
	    {
            if(index_E == 0) SepTot[i] = sep[i];
            RateAvg[i]+=rate[i];
			RateAvg_error[i] += rate_error[i]*rate_error[i];
        }

        index_E++;
    }
    if(index_E != nEBC && scan_type == 1) cout<<"ERROR: Number of empty does not match in x scan, scan "<<scan<<" : nEBC = "<<nEBC<<" ; index_E = "<<index_E<<endl;
	if(index_E != nEBC && scan_type == 2) cout<<"ERROR: Number of empty does not match in y scan, scan "<<scan<<" : nEBC = "<<nEBC<<" ; index_E = "<<index_E<<endl;

    for(Int_t i=0; i<n_sep; i++){
		RateAvg[i]= RateAvg[i]/nEBC;
		if(RateAvg[i] == 0) RateAvg_error[i] = 0;
		else RateAvg_error[i] = TMath::Sqrt(RateAvg_error[i])/nEBC;
        noise[i] = RateAvg[i];
        noise_err[i] = RateAvg_error[i];
        //cout<<noise[i]<<" "<<noise_err[i]<<endl;
	}

    Double_t sep_min = 0;
    Double_t sep_max = 0;
    if (SepTot[n_sep-1]>SepTot[0]) { 
        sep_min = SepTot[0] - 0.5*TMath::Abs(SepTot[0]-SepTot[1]);
        sep_max = SepTot[n_sep-1] + 0.5*TMath::Abs(SepTot[n_sep-1]-SepTot[n_sep-2]);
    } else {
        sep_max = SepTot[0] + 0.5*TMath::Abs(SepTot[0]-SepTot[1]);
        sep_min = SepTot[n_sep-1] - 0.5*TMath::Abs(SepTot[n_sep-1]-SepTot[n_sep-2]);
    }

    TGraphErrors *gr = new TGraphErrors(n_sep, SepTot, RateAvg, NULL, RateAvg_error);
	TF1 *fitFuncPol0 = new TF1("fitFuncPol0", "pol0", sep_min, sep_max);
	if(strcmp(rate_name, "aodzn") == 0){
		TFitResultPtr r = gr->Fit("fitFuncPol0", "Q0RS");
		//gStyle->SetOptFit(1);
        noisePol0 = fitFuncPol0->GetParameter(0);
        //noisePol0_err = fitFuncPol0->GetParError(0);
	}
	if(strcmp(rate_name, "aodft0vtx") == 0){
		gr -> Fit("gaus");
		//gStyle->SetOptFit(1);
        noisePol0 = 0;
        //noisePol0_err = 0;
	}

    const char* cName1 = Form("Fill%i_%s_%s_RawRateEmpty_vs_%sSep", g_vdm_Fill, sep_type, rate_name, sep_type);
	const char* cName2;
	if(scan_type == 1) cName2 = Form("%s_scan%i_x", cName1, scan);
	if(scan_type == 2) cName2 = Form("%s_scan%i_y", cName1, scan);
	TCanvas* c1 = new TCanvas(cName2, cName2, 800, 600);
	c1->cd();
    gPad->SetLeftMargin(0.15);
	gr->SetMarkerStyle(8);
	gr->GetXaxis()->SetTitle("Separation [mm]");
	if(strcmp(rate_name, "aodzn") == 0) gr->GetYaxis()->SetTitle("#LT R_{empty}^{ZNC} #GT [Hz]");
	if(strcmp(rate_name, "ft0vtx") == 0) gr->GetYaxis()->SetTitle("#LT R_{empty}^{MTVX} #GT [Hz]");
	gr->SetTitle(cName2);
	gr->Draw("ape");
	if(strcmp(rate_name, "aodzn") == 0){
		fitFuncPol0->SetLineStyle(2);
		fitFuncPol0->Draw("same");
	}
	if(strcmp(rate_name, "aodzn") == 0) c1->Print(Form("../Fill-%d/BkgPlotsZNC/%s.pdf", g_vdm_Fill, cName2));
	if(strcmp(rate_name, "ft0vtx") == 0) c1->Print(Form("../Fill-%d/BkgPlotsMTVX/%s.pdf", g_vdm_Fill, cName2));



    /////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// BEAM-GAS A ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    TF1 *fitFuncPol0A[nABC];
    for (Int_t k=0; k<3564; k++) // loop over bunches
    {
        //------------------------------------------------------------------------------
        //select only typeA BC:
        Int_t bc_index_A;
        Bool_t skipBc_notA = kTRUE;
		for(Int_t bcA=0; bcA<nABC; bcA++){if(k == bunches_A[bcA]){skipBc_notA = false; bc_index_A = bcA; break;}}  //skip bc if not type A
		if(skipBc_notA) continue;

        // get info
		rate_tree->GetEntry(k);
		sep_tree->GetEntry(k);
        intensity_tree->GetEntry(bc_index_A);

        //NEW: subtract noise as a distribution from rate
        Double_t rate_A_noNoise[n_sep];
        Double_t rate_error_A_noNoise[n_sep];
        for(Int_t i=0; i<n_sep; i++){
            if(rate[i] == 0)
            {
                rate_A_noNoise[i] = 0;
                rate_error_A_noNoise[i] = 0;
            }
            else
            {
                rate_A_noNoise[i] = rate[i] - noise[i];
                rate_error_A_noNoise[i] = TMath::Sqrt(rate_error[i]*rate_error[i] + noise_err[i]*noise_err[i]);
            }
            if(scan==0 && scan_type==1) cout<<rate[i]<<" "<<noise[i]<<" "<<rate_A_noNoise[i]<<endl;
	    }

        // do fit
        fitFuncPol0A[bc_index_A] = new TF1(Form("fitFuncPol0A_%d", bc_index_A), "pol0", sep_min, sep_max);

        // define a TGraph to perform the fit
	    TGraphErrors *grA = new TGraphErrors(n_sep, sep, rate_A_noNoise, NULL, rate_error_A_noNoise);
        // fit and check
        TVirtualFitter::SetDefaultFitter("Minuit");
        //fitFuncPol0A[bc_index_A]->SetParameter(0,0.2);
        TFitResultPtr r;
        if(scan==0 && scan_type==1)
            r = grA->Fit(Form("fitFuncPol0A_%d", bc_index_A), "RS");
        else
            r = grA->Fit(Form("fitFuncPol0A_%d", bc_index_A), "Q0RS");

        //gStyle->SetOptFit(1);
        //chi2_dof = fitFuncPol0A[bc_index_A]->GetChisquare()/((Double_t)fitFuncPol0A[bc_index_A]->GetNDF());

        //TF1* f = grA->GetFunction("fitFuncPol0A");
        Double_t N_A = cf_dcct_1 * bunch_intensity_1[bc_index_A];
        Double_t R_A = fitFuncPol0A[bc_index_A]->GetParameter(0);
        //Double_t R_A = f->GetParameter(0);
        Double_t  R_A_err;
        //if(fitFuncPol0A[bc_index_A]->GetParameter(0) == 0) R_A_err = 0;
        if(R_A == 0) R_A_err = 0;
        else R_A_err = fitFuncPol0A[bc_index_A]->GetParError(0);

        grRvsN->SetPoint(grRvsN->GetN(), N_A, R_A);
        grRvsN->SetPointError(grRvsN->GetN()-1, 0, R_A_err);
        if(scan==0 && scan_type==1) cout<<"NA="<<N_A<<"; RA="<<R_A<<endl;

        const char* cNameBC0 = Form("Fill%i_%s_Raw_%s", Fill, sep_type, rate_name);
        const char* cNameBC;
        if(scan_type == 1) cNameBC = Form("TypeARawRate_%s_scan%i_x_iBC%d", cNameBC0, scan, k);
        if(scan_type == 2) cNameBC = Form("TypeARawRate_%s_scan%i_x_iBC%d", cNameBC0, scan, k);
        TCanvas* c0 = new TCanvas(cNameBC, cNameBC, 800, 600);
        grA->SetMarkerStyle(20);
        grA->SetTitle(cNameBC);
        grA->GetXaxis()->SetTitle("Separation [mm]");
        grA->GetYaxis()->SetTitle("(Rate_{A} - #LT Rate_{empty} #GT)) [Hz]");
        grA->Draw("AP");
        fitFuncPol0A[bc_index_A]->SetParameters(r->GetParams());
        fitFuncPol0A[bc_index_A]->SetParErrors(r->GetErrors());
        grA->GetListOfFunctions()->Remove(fitFuncPol0A[bc_index_A]); // optional cleanup
        grA->Fit(fitFuncPol0A[bc_index_A], "R");
        fitFuncPol0A[bc_index_A]->SetLineColor(kGreen+1);
        fitFuncPol0A[bc_index_A]->Draw("SAME");
        if(strcmp(rate_name, "aodzn") == 0) c0->Print(Form("../Fill-%d/BkgPlotsZNC/TypeA_fits/%s_AOD.pdf", g_vdm_Fill, cNameBC));
		
    }// end loop over bunches

    // fit rate vs intensity
    if(strcmp(rate_name, "aodzn") == 0){
        TF1 *fitFuncPol1 = new TF1("fitFuncPol1", "pol1", 2e6, 4e6);
        //cout<<noisePol0<<endl;
        //fitFuncPol1->FixParameter(0, noisePol0);
        fitFuncPol1->FixParameter(0, 0); //noise already subtracted in by bin
        fitFuncPol1->SetParameter(1, 7.62087e-12);
        fitFuncPol1->SetLineColor(kGreen+1);
        grRvsN->Fit("fitFuncPol1");
        //gStyle->SetOptFit(1112);
        //gStyle->SetOptFit(1);

        //save result for each scan number and type
        const char* cName0 = Form("Fill%i_%s_Raw_%s", Fill, sep_type, rate_name);
        const char* cName3;
        if(scan_type == 1) cName3 = Form("RaVsNa_%s_scan%i_x", cName0, scan);
        if(scan_type == 2) cName3 = Form("RaVsNa_%s_scan%i_y", cName0, scan);
        TCanvas* c2 = new TCanvas(cName3, cName3, 800, 600);
        c2->cd();
        grRvsN->SetMarkerStyle(8);
        grRvsN->GetXaxis()->SetTitle("N_{A}");
        grRvsN->GetYaxis()->SetTitle("(R_{A}^{ZNC} - #LT R_{empty}^{ZNC} #GT) [Hz]");
        const char* cTitle2; 
        if(scan_type == 1) cTitle2 = Form("ZNC rate vs intensity, fill %i, scan %i, x", Fill, scan);
        if(scan_type == 2) cTitle2 = Form("ZNC rate vs intensity, fill %i, scan %i, y", Fill, scan);
        grRvsN->SetTitle(cTitle2);
        grRvsN->GetYaxis()->SetRangeUser(0.1, 0.35);
        grRvsN->Draw("ape");
        c2->Print(Form("../Fill-%d/BkgPlotsZNC/%s_AOD.pdf", g_vdm_Fill, cName3));

        beam_gas_A = fitFuncPol1->GetParameter(1);
        beam_gas_A_err = fitFuncPol1->GetParError(1);
        h_tree->Fill();
	
    }
    if(strcmp(rate_name, "aodft0vtx") == 0)
    {
        //save result for each scan number and type
        const char* cName0 = Form("Fill%i_%s_Raw_%s", Fill, sep_type, rate_name);
        const char* cName3;
        if(scan_type == 1) cName3 = Form("RaVsNa_%s_scan%i_x", cName0, scan);
        if(scan_type == 2) cName3 = Form("RaVsNa_%s_scan%i_y", cName0, scan);
        TCanvas* c2 = new TCanvas(cName3, cName3, 800, 600);
        c2->cd();
        grRvsN->SetMarkerStyle(8);
        grRvsN->GetXaxis()->SetTitle("N_{A}");
        grRvsN->GetYaxis()->SetTitle("(R_{A}^{MTVX} - #LT R_{empty}^{MTVX} #GT) [Hz]");
        const char* cTitle2; 
        if(scan_type == 1) cTitle2 = Form("MTVX rate vs intensity, fill %i, scan %i, x", Fill, scan);
        if(scan_type == 2) cTitle2 = Form("MTVX rate vs intensity, fill %i, scan %i, y", Fill, scan);
        grRvsN->SetTitle(cTitle2);
        grRvsN->Draw("ape");
        c2->Print(Form("../Fill-%d/BkgPlotsMTVX/%s_AOD.pdf", g_vdm_Fill, cName3));

        beam_gas_A = 0;
        beam_gas_A_err = 0;
        h_tree->Fill();
    }

    // save
	hFile->cd();
	h_tree->SetDirectory(hFile);
	h_tree->Write();
	hFile->Close();
    
    // clean up
	delete [] h_file_name;    
	delete [] rate_file_name;
	delete [] sep_file_name;
	delete [] rate;
	delete [] rate_error;  
	delete [] sep;
    delete [] bunch_intensity_1;
}

void Create_noise_beamgas_AOD(
    Int_t Fill, const char *rate_name, const char *sep_type, const char *intensity_type
)
{
    // initialize
	Set_input_file_names(Fill);
	Set_pointers_to_input_files_and_trees();
	gROOT->SetBatch(kTRUE); // not to display canvases

    // create files for all scans
	for (Int_t i=0;i<g_n_Scans_in_Fill;i++)
	{
		Create_scan_noise_beamgas(Fill, 1, i, rate_name, sep_type, intensity_type); //x-scans
		Create_scan_noise_beamgas(Fill, 2, i, rate_name, sep_type, intensity_type); //y-scans
	} 

	return;
}

#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Create a background corrected rate file for a given
// trigger, background satCorr and scan
// It includes satellites, beam-gas and noise
//-------------------------------------------------------

void Create_one_bkgd_corrected_rate_file_FULL(
		Int_t scan_type, Int_t scan, const char *rate_name, const char *intensity_type, const char *sat_corr_name
		)
// scan_type: 1 => x-scan; 2 => y-scan
// sat_corr_name -> satellite fraction, multiplicative factor BC by BC and sep by sep
// noise -> common baseline which is subtracted
// beam_gas_A -> Beam-gas from A side, factor thet has to be multiplied by the charge of the A-side bunch
// beam_gas_C -> Same for C-side, mostly 0 due to filling scheme
// ?afterpulses? -> N, mu, sigma -> gaussian shape of afterpulses ?? 
{
	// create names for files
	char *file_name_rate = new char[kg_string_size];  
	char *file_name_sat = new char[kg_string_size];
	char *file_name_noise_beamgas = new char[kg_string_size];
	char *file_name_new_rate = new char[kg_string_size];

	if (scan_type == 1)
	{
		sprintf(file_name_rate, "../Fill-%d/RawRate_%s_x_Scan_%d.root",g_vdm_Fill,rate_name,scan);
		sprintf(file_name_sat, "../Fill-%d/BkgdCorr_%s_x_Scan_%d.root",g_vdm_Fill,sat_corr_name,scan);
		sprintf(file_name_noise_beamgas, "../Fill-%d/NoiseBeamgas_%s_NomSep_x_Scan_%d.root",g_vdm_Fill,rate_name,scan);
		sprintf(file_name_new_rate, "../Fill-%d/BkgdCorrRate_%s_x_Scan_%d.root",g_vdm_Fill,rate_name,scan);
	}
	if (scan_type == 2)
	{
		sprintf(file_name_rate, "../Fill-%d/RawRate_%s_y_Scan_%d.root",g_vdm_Fill,rate_name,scan);
		sprintf(file_name_sat, "../Fill-%d/BkgdCorr_%s_y_Scan_%d.root",g_vdm_Fill,sat_corr_name,scan);
		sprintf(file_name_noise_beamgas, "../Fill-%d/NoiseBeamgas_%s_NomSep_y_Scan_%d.root",g_vdm_Fill,rate_name,scan);
		sprintf(file_name_new_rate, "../Fill-%d/BkgdCorrRate_%s_y_Scan_%d.root",g_vdm_Fill,rate_name,scan);
	}

	// open files and get trees
	TFile *RateFile = new TFile(file_name_rate);
	TFile *SatFile = new TFile(file_name_sat);
	TFile *NoiseBeamgasFile = new TFile(file_name_noise_beamgas);
	TTree *rate_tree = (TTree *) RateFile->Get("Rate");
	TTree *sat_tree = (TTree *) SatFile->Get("BkgdCorr");  
	TTree *noise_beamgas_tree = (TTree *) NoiseBeamgasFile->Get("NoiseBeamgasPar");  

	// set branch adresses
	Int_t n_separations = FindNumberSeparations(scan_type, scan);
	Double_t *rate = new Double_t[n_separations];
	Double_t *rate_error = new Double_t[n_separations];  
	rate_tree->SetBranchAddress("rate",rate);
	rate_tree->SetBranchAddress("rate_error",rate_error);  
	cout<<"0TH CHECK"<<endl; 
	Double_t *satCorr = new Double_t[n_separations];
	Double_t *satCorr_error = new Double_t[n_separations];  
	sat_tree->SetBranchAddress("correction",satCorr);
	sat_tree->SetBranchAddress("correction_error",satCorr_error);
	cout<<"FIRST CHECK"<<endl;
	cout<<file_name_noise_beamgas<<endl;
	Double_t *noise = new Double_t[n_separations];
	Double_t *noise_err = new Double_t[n_separations];
	Double_t beam_gas_A;
	Double_t beam_gas_A_err;
	noise_beamgas_tree->SetBranchAddress("noise",noise);
	noise_beamgas_tree->SetBranchAddress("noise_err",noise_err);
	noise_beamgas_tree->SetBranchAddress("beam_gas_A",&beam_gas_A);
	noise_beamgas_tree->SetBranchAddress("beam_gas_A_err",&beam_gas_A_err);
	cout<<"SECOND CHECK"<<endl;
	Double_t beam_gas_C = 0;
	Double_t beam_gas_C_err = 0;

	// --------------------
	// get intensity file and tree
	//---------------------
	// -- number of bc
	Int_t nIBC = GetNumberInteractingBunchCrossings();

	// open the correct file
	char IntFileName[120];
	sprintf(IntFileName,"../Fill-%d/%s_Scan_%d.root",g_vdm_Fill,intensity_type,scan);
	TFile *IntensityFile = new TFile(IntFileName);

	// get the info from the tree
	Double_t *bunch_intensity_1 = new Double_t[nIBC];
	Double_t *bunch_intensity_2 = new Double_t[nIBC];
	Double_t cf_dcct_1;
	Double_t cf_dcct_2;

	TTree *intensity_tree = (TTree *) IntensityFile->Get("Bunch_Intensity");
	intensity_tree->ResetBranchAddresses();
	intensity_tree->SetBranchAddress("cf_dcct_1",&cf_dcct_1);
	intensity_tree->SetBranchAddress("cf_dcct_2",&cf_dcct_2);    
	intensity_tree->SetBranchAddress("bunch_intensity_1",bunch_intensity_1);
	intensity_tree->SetBranchAddress("bunch_intensity_2",bunch_intensity_2);
	intensity_tree->GetEntry(0);

	cout<<"THIRD CHECK"<<endl;

	//---------------------
	// prepare new tree
	//---------------------
	TFile *NewRateFile = new TFile(file_name_new_rate,"recreate");
	TTree *new_rate_tree = new TTree("Rate","Rate");
	Double_t *new_rate = new Double_t[n_separations];
	Double_t *new_rate_error = new Double_t[n_separations];
	char *txt_tmp = new char[kg_string_size];
	sprintf(txt_tmp,"rate[%d]/D",n_separations);
	new_rate_tree->Branch("rate",new_rate,txt_tmp);
	sprintf(txt_tmp,"rate_error[%d]/D",n_separations);
	new_rate_tree->Branch("rate_error",new_rate_error,txt_tmp);

	cout<<"FOURTH CHECK"<<endl;

	// loop over input data to compute new rate
	// -- loop over bc
	for(Int_t k=0;k<nIBC;k++)
	{
		// get entries
		rate_tree->GetEntry(k);
		sat_tree->GetEntry(k);
		intensity_tree->GetEntry(k);
		noise_beamgas_tree->GetEntry(0); //Only one entry for all BCs

		Double_t NA = cf_dcct_1 * bunch_intensity_1[k];
		Double_t NC = cf_dcct_2 * bunch_intensity_2[k];

		// compute new rate
		for(Int_t i=0; i<n_separations; i++) // loop over steps
		{
			cout<<"RawRate="<<rate[i]<<" Noise="<<noise[i]<<" B-Gas="<<beam_gas_A<<" Sat="<<satCorr[i]<<endl;
			Double_t rate_noise_beam_gas = rate[i] - noise[i] - beam_gas_A*NA -  beam_gas_C*NC; //subtract beam-gas int on both sides
			Double_t rate_noise_beam_gas_err = TMath::Sqrt(rate_error[i]*rate_error[i] + noise_err[i]*noise_err[i] + NA*NA*beam_gas_A_err*beam_gas_A_err + NC*NC*beam_gas_C_err*beam_gas_C_err);
			new_rate[i] = rate_noise_beam_gas*satCorr[i]; //including satellites
			new_rate_error[i] = BkgdCorrectedRateError(rate_noise_beam_gas,rate_noise_beam_gas_err,satCorr[i],satCorr_error[i]);
		}

		// store the info
		new_rate_tree->Fill();
	}

	cout<<"FIFTH CHECK"<<endl;


	// save  tree
	NewRateFile->cd();
	new_rate_tree->SetDirectory(NewRateFile);
	new_rate_tree->Write();
	NewRateFile->Close();

	cout<<"SIXTH CHECK"<<endl;


	// free memory
	delete [] file_name_rate;
	delete [] file_name_sat;  
	delete [] file_name_new_rate;
	delete [] txt_tmp;  
	delete [] rate;
	delete [] rate_error;
	delete [] new_rate;
	delete [] new_rate_error;
	delete [] satCorr;
	delete [] satCorr_error;
	return;
}

//-------------------------------------------------------
// Create background corrected  rate files
//-------------------------------------------------------

void Create_bkgd_corrected_rate_file_FULL(Int_t Fill, const char *rate_name, const char *intensity_type, const char *sat_corr_name
	)
{
	// get name of files and set pointers to trees
	Set_input_file_names(Fill);
	Set_pointers_to_input_files_and_trees();

	// create files for all scans
	for(Int_t i=0;i<g_n_Scans_in_Fill;i++)
	{
		Create_one_bkgd_corrected_rate_file_FULL(1, i, rate_name, intensity_type, sat_corr_name); //x-scans
		Create_one_bkgd_corrected_rate_file_FULL(2, i, rate_name, intensity_type, sat_corr_name); //y-scans
	}

	return;
}

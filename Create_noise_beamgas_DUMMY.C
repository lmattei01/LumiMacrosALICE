#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Background analysis (noise+beam-gas) using AO2D's
//-------------------------------------------------------

void Create_scan_noise_beamgas(Int_t Fill, Int_t scan_type, Int_t scan, const char *rate_name, const char *sep_type, const char *intensity_type)
{

    // Now reserve space for output
    Double_t noise;
    Double_t noise_err;
	Double_t beam_gas_A;
	Double_t beam_gas_A_err;

	// create file and tree for output
    // --> file
	char *h_file_name = new char[kg_string_size];
	if (scan_type == 1) sprintf(h_file_name,"../Fill-%d/NoiseBeamgas_%s_%sSep_x_Scan_%d.root",
			g_vdm_Fill,rate_name,sep_type,scan);
	if (scan_type == 2) sprintf(h_file_name,"../Fill-%d/NoiseBeamgas_%s_%sSep_y_Scan_%d.root",
			g_vdm_Fill,rate_name,sep_type,scan);  
	TFile *hFile = new TFile(h_file_name,"recreate");
    // --> tree
	TTree *h_tree = new TTree("NoiseBeamgasPar","NoiseBeamgasPar");
    h_tree->Branch("noise",&noise,"noise/D");
    h_tree->Branch("noise_err",&noise_err,"noise_err/D");
    h_tree->Branch("beam_gas_A",&beam_gas_A,"beam_gas_A/D");
    h_tree->Branch("beam_gas_A_err",&beam_gas_A_err,"beam_gas_A_err/D");

    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// NOISE ///////////////////////////////////////
    ///////////////////////////////////////// & /////////////////////////////////////////
    ///////////////////////////////////// BEAM-GAS A ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
	if(strcmp(rate_name, "aodft0vtx") == 0){
        noise = 0;
        noise_err = 0;
        beam_gas_A = 0;
        beam_gas_A_err = 0;
        h_tree->Fill();
	}
    else{
        noise = 0;
        noise_err = 0;
        beam_gas_A = 0;
        beam_gas_A_err = 0;
        h_tree->Fill();
        //cout<<"ERROR: "<<rate_name<<" Not dummy!!"<<endl;
        //cout<<"run .L Create_noise_beamgas_AOD.C istead!"<<endl;
    }
    
    // save
	hFile->cd();
	h_tree->SetDirectory(hFile);
	h_tree->Write();
	hFile->Close();

    // clean up
	delete [] h_file_name;    
}

void Create_noise_beamgas_DUMMY(
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

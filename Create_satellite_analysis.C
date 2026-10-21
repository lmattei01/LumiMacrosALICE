#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"
#include <iostream>
#include <fstream>
#include <vector>

//-------------------------------------------------------
// Macro to create satellite analysis
// It requires:
// 1) Running QA_print_timestamps.C in order to get the 
//    timestamps of the VdM scan 
// 2) Producing an AO2D from O2Physics task, which 
//    contains for each bc the timestamp, timeZNA, timeZNC
//-------------------------------------------------------


void Create_satellite_analysis(Int_t Fill, Int_t scan, const char*TrigType){

    //Options for TrigType: MTVX, 1ZNC, 1ZNCTypeA, 1ZNCTypeC31

    // initialize
	Set_input_file_names(Fill);
	Set_pointers_to_input_files_and_trees();

    //Get interacting BCs and type-A BCs
    Int_t nIBC = GetNumberInteractingBunchCrossings(); //number of interacting bc
	Int_t nABC = GetNumberBunchCrossingsTypeA(); //number of non-interacting type-A bc
	Int_t nCBC = GetNumberBunchCrossingsTypeC(); //number of non-interacting type-C bc

    Int_t *bunches_I = new Int_t [nIBC];
	Int_t *bunches_A = new Int_t [nABC];
	Int_t *bunches_C = new Int_t [nCBC];
    //get bunch indices for all types
	GetBunchIndices(bunches_I);
	GetBunchIndices_typeA(bunches_A);
	GetBunchIndices_typeC(bunches_C);

    //Get vdM timestamps
    char *file_sep_vdm_name = new char[kg_string_size];
    sprintf(file_sep_vdm_name,"../Fill-%d/Timestamps_Scan_%d.txt",g_vdm_Fill, scan);
    ifstream file_sep_vdm(file_sep_vdm_name);

    vector<uint64_t> tstamps_start, tstamps_end, tstamps_mean; 
    uint64_t tstamp_start, tstamp_end, tstamp_mean;

    while (file_sep_vdm >> tstamp_start >> tstamp_end >> tstamp_mean) {   // read 3 columns per line
        tstamps_start.push_back(tstamp_start);
        tstamps_end.push_back(tstamp_end);
        tstamps_mean.push_back(tstamp_mean);
    }
    file_sep_vdm.close();

    Int_t nTimestamps = tstamps_start.size();
    //nTimestamps = 10;
    Int_t n_separations = FindNumberSeparations(1, scan); //n_separations along x
    //Int_t n_separations = 21;
    cout<<"Number of sepration in each scan: "<<n_separations<<endl;

    char *ao2d_name = new char[kg_string_size];
    TFile *fAo2d;
    TTree *zdc_tree;

    if(strcmp(TrigType, "MTVX") == 0){
        sprintf(ao2d_name,"../Fill-%d/AO2D_HY_NEW_MTVX.root",g_vdm_Fill);
        fAo2d = new TFile(ao2d_name, "READ");
        if(Fill == 10802) zdc_tree = (TTree *) fAo2d->Get("DF_2423894932451520/O2mybcaod");
        else if(Fill == 10782) zdc_tree = (TTree *) fAo2d->Get("DF_2423538867473184/O2mybcaod");
    }
    else if(strcmp(TrigType, "1ZNC") == 0 || strcmp(TrigType, "1ZNCTypeA") == 0 || strcmp(TrigType, "1ZNCTypeC31") == 0 ){
        //sprintf(ao2d_name,"../Fill-%d/AO2D_HY_NEW_1ZNC.root",g_vdm_Fill); //previously using 1ZNC trigger, not enough info
        sprintf(ao2d_name,"../Fill-%d/AO2D_HY_NEW_1ZED.root",g_vdm_Fill); //using 1ZED trigger to have all info
        fAo2d = new TFile(ao2d_name, "READ");
        if(Fill == 10802) zdc_tree = (TTree *) fAo2d->Get("DF_2423894932451520/O2mybcaod");
        else if(Fill == 10782) zdc_tree = (TTree *) fAo2d->Get("DF_2423538867473184/O2mybcaod");
    }
    else{
        cout<<"Trigger name not valid!"<<endl;
        exit(1);
    }

    Int_t nEntries = zdc_tree->GetEntries();
    Int_t bcId = 0;
    uint64_t timestamp = 0;
    Float_t timeZNA = 0;
    Float_t timeZNC = 0;
    Float_t amplZNA = 0;
    Float_t amplZNC = 0;

    zdc_tree->ResetBranchAddresses();
    zdc_tree->SetBranchAddress("fBCid", &bcId);
    zdc_tree->SetBranchAddress("fTimestamp", &timestamp);
    zdc_tree->SetBranchAddress("fTimeZNA", &timeZNA);
    zdc_tree->SetBranchAddress("fTimeZNC", &timeZNC);
    zdc_tree->SetBranchAddress("fAmplitudeZNA", &amplZNA);
    zdc_tree->SetBranchAddress("fAmplitudeZNC", &amplZNC);

    //uint64_t tstamp_start_0_x = 1751688982; //start of scan 0 x: 1751688982000
    //uint64_t tstamp_end_1_y = 1751692031;//end   of scan 1 y: 1751692031000

    //The region limit for satellite rejection is defined here
    Double_t region_lim_pos = 1.5;
    Double_t region_lim_neg = -4;
    TLine *l1 = new TLine(region_lim_neg,region_lim_neg,region_lim_pos,region_lim_neg);
    TLine *l2 = new TLine(region_lim_pos,region_lim_neg,region_lim_pos,region_lim_pos);
    TLine *l3 = new TLine(region_lim_pos,region_lim_pos,region_lim_neg,region_lim_pos);
    TLine *l4 = new TLine(region_lim_neg,region_lim_pos,region_lim_neg,region_lim_neg);
    l1->SetLineColor(kRed); l1->SetLineWidth(2); 
    l2->SetLineColor(kRed); l2->SetLineWidth(2); 
    l3->SetLineColor(kRed); l3->SetLineWidth(2);
    l4->SetLineColor(kRed); l4->SetLineWidth(2);

    //TH for output
    TH2D *hZDCTimes[nTimestamps][nIBC];
    TH1D *hZNCTimes[nTimestamps][nIBC];
    TH1D *hZNCnoZNATimes[nTimestamps][nIBC];
    TH1D *hZNAnoZNCTimes[nTimestamps][nIBC];



    //Double_t counts_ratio[nTimestamps];
    TH1D *hCountsRatio_vs_x[nIBC];
    TH1D *hCountsRatio_vs_y[nIBC];


    TH1D *hQAZDC = new TH1D("hQAZDC", "ZDC QA", 4, -0.5, 3.5);
    hQAZDC->GetXaxis()->SetBinLabel(1, "ZNC && ZNA");
    hQAZDC->GetXaxis()->SetBinLabel(2, "ZNC && !ZNA");
    hQAZDC->GetXaxis()->SetBinLabel(3, "!ZNC && ZNA");
    hQAZDC->GetXaxis()->SetBinLabel(4, "!ZNC && !ZNA");

    // prepare new trees (for scan x and scan y)
    char *file_name_corr_x = new char[kg_string_size];
    sprintf(file_name_corr_x, "../Fill-%d/BkgdCorr_%s_x_Scan_%d.root",g_vdm_Fill,TrigType,scan);
    TFile *CorrFile_x = new TFile(file_name_corr_x, "RECREATE");
    TTree *corr_tree_x = new TTree("BkgdCorr","BkgdCorr in x");
    Double_t *correction_x = new Double_t[n_separations];
	Double_t *correction_x_error = new Double_t[n_separations];  
	corr_tree_x->Branch("correction", correction_x, Form("correction[%d]/D", n_separations));
    corr_tree_x->Branch("correction_error", correction_x_error, Form("correction_error[%d]/D", n_separations));

    char *file_name_corr_y = new char[kg_string_size];
    sprintf(file_name_corr_y, "../Fill-%d/BkgdCorr_%s_y_Scan_%d.root",g_vdm_Fill,TrigType,scan);
    TFile *CorrFile_y = new TFile(file_name_corr_y, "RECREATE");
    TTree *corr_tree_y = new TTree("BkgdCorr","BkgdCorr in x");
    Double_t *correction_y = new Double_t[n_separations];
	Double_t *correction_y_error = new Double_t[n_separations];  
	corr_tree_y->Branch("correction", correction_y, Form("correction[%d]/D", n_separations));
    corr_tree_y->Branch("correction_error", correction_y_error, Form("correction_error[%d]/D", n_separations));
    
    cout<<"Reading input ao2d ..."<<endl;

    for(Int_t iTstamp=0; iTstamp<nTimestamps; iTstamp++){
        cout<<"Processing separation index "<<iTstamp<<endl;
        uint64_t tstamp_min = tstamps_start[iTstamp]; 
        uint64_t tstamp_max = tstamps_end[iTstamp];

        Double_t counts_insideTH2 = 0;
        Double_t counts_outsideTH2 = 0;

        Double_t counts_insideTH1 = 0;
        Double_t counts_outsideTH1 = 0;

        for(Int_t k=0; k<nIBC; k++){
            hZDCTimes[iTstamp][k] = new TH2D(Form("hZDCTimes_sep%d_bc%d", iTstamp, k), Form("t_{ZNA} vs t_{ZNC}, iBC%d (%d), sep idx #%d", k, bunches_I[k], iTstamp), 800, -50, 50, 800, -50, 50);
            hZNCTimes[iTstamp][k] = new TH1D(Form("hZNCTimes_sep%d_bc%d", iTstamp, k), Form("t_{ZNC}, iBC%d (%d), sep idx #%d", k, bunches_I[k], iTstamp), 800, -50, 50);
            hZNCnoZNATimes[iTstamp][k] = new TH1D(Form("hZNCnoZNATimes_sep%d_bc%d", iTstamp, k), Form("t_{ZNC} (!t_{ZNA}), iBC%d (%d), sep idx #%d", k, bunches_I[k], iTstamp), 400, -50, 50);
            hZNAnoZNCTimes[iTstamp][k] = new TH1D(Form("hZNAnoZNCTimes_sep%d_bc%d", iTstamp, k), Form("t_{ZNA} (!t_{ZNC}), iBC%d (%d), sep idx #%d", k, bunches_I[k], iTstamp), 400, -50, 50);
        }

        /////////////////////////////////////////////////////////////////////////
        //                         FILLING HISTOGRAMS
        /////////////////////////////////////////////////////////////////////////
        for(int bc=0; bc<zdc_tree->GetEntries(); bc++){
            zdc_tree->GetEntry(bc);

            //select interacting Bc only, save index
            Int_t bc_index_int;
            bool skipBc = true;
            for(Int_t k=0; k<nIBC; k++){
                if(bcId == bunches_I[k]){
                    bc_index_int = k;
                    skipBc = false; 
                    break;
                }  //skip bc if non-interacting
            }
            
            if(skipBc) continue;
            if(timestamp <= tstamp_max*1e3 && timestamp >= tstamp_min*1e3){ //also adjusting measure units
                
                hZDCTimes[iTstamp][bc_index_int]->Fill(timeZNC, timeZNA);
                
                hZNCTimes[iTstamp][bc_index_int]->Fill(timeZNC);
                if(timeZNA<45 && timeZNA>40) hZNCnoZNATimes[iTstamp][bc_index_int]->Fill(timeZNC); //if t(ZNA) at infinity, it is set at a dummy value of 42
                if(timeZNC<45 && timeZNC>40) hZNAnoZNCTimes[iTstamp][bc_index_int]->Fill(timeZNA); //if t(ZNC) at infinity, it is set at a dummy value of 42
            }
        }
        
    }

    cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    cout<<"Now computing correction factor plots ..."<<endl;
    
    for(Int_t k=0; k<nIBC; k++){
        cout<<" Processing BC "<<k<<" , in orbit "<<bunches_I[k]<<endl;
        hCountsRatio_vs_x[k] = new TH1D(Form("hCountsRatio_vs_x_iBC%d", k), Form("CountsRatio vs x, iBC%d (%d)", k, bunches_I[k]), nTimestamps, 0, nTimestamps);
        hCountsRatio_vs_x[k]->SetLineColor(kBlue);
        hCountsRatio_vs_y[k] = new TH1D(Form("hCountsRatio_vs_y_iBC%d", k), Form("CountsRatio vs y, iBC%d (%d)", k, bunches_I[k]), nTimestamps, 0, nTimestamps);
        hCountsRatio_vs_y[k]->SetLineColor(kRed);

        for(Int_t iTstamp=0; iTstamp<nTimestamps; iTstamp++){

            Int_t bin_low_x_TOT = hZDCTimes[iTstamp][k]->GetXaxis()->FindBin(-20);
            Int_t bin_high_x_TOT = hZDCTimes[iTstamp][k]->GetXaxis()->FindBin(20);
            Int_t bin_low_y_TOT = hZDCTimes[iTstamp][k]->GetYaxis()->FindBin(-20);
            Int_t bin_high_y_TOT = hZDCTimes[iTstamp][k]->GetYaxis()->FindBin(20);

            Int_t bin_low_x_REG = hZDCTimes[iTstamp][k]->GetXaxis()->FindBin(region_lim_neg);
            Int_t bin_high_x_REG = hZDCTimes[iTstamp][k]->GetXaxis()->FindBin(region_lim_pos);
            Int_t bin_low_y_REG = hZDCTimes[iTstamp][k]->GetYaxis()->FindBin(region_lim_neg);
            Int_t bin_high_y_REG = hZDCTimes[iTstamp][k]->GetYaxis()->FindBin(region_lim_pos);
            
            double counts_TH2_TOT = hZDCTimes[iTstamp][k]->Integral(bin_low_x_TOT, bin_high_x_TOT, bin_low_y_TOT, bin_high_y_TOT);
            double counts_TH2_inside = hZDCTimes[iTstamp][k]->Integral(bin_low_x_REG, bin_high_x_REG, bin_low_y_REG, bin_high_y_REG);

            Int_t bin_low_TOT = hZNCnoZNATimes[iTstamp][k]->GetXaxis()->FindBin(-20);
            Int_t bin_high_TOT = hZNCnoZNATimes[iTstamp][k]->GetXaxis()->FindBin(20);
            
            Int_t bin_low_REG = hZNCnoZNATimes[iTstamp][k]->GetXaxis()->FindBin(region_lim_neg);
            Int_t bin_high_REG = hZNCnoZNATimes[iTstamp][k]->GetXaxis()->FindBin(region_lim_pos);

            double counts_TH1C_TOT = hZNCnoZNATimes[iTstamp][k]->Integral(bin_low_TOT, bin_high_TOT);
            double counts_TH1C_inside = hZNCnoZNATimes[iTstamp][k]->Integral(bin_low_REG, bin_high_REG);
            double counts_TH1C_outside = counts_TH1C_TOT - counts_TH1C_inside;

            Int_t bin_low_TOT_A = hZNAnoZNCTimes[iTstamp][k]->GetXaxis()->FindBin(-20);
            Int_t bin_high_TOT_A = hZNAnoZNCTimes[iTstamp][k]->GetXaxis()->FindBin(20);
            
            Int_t bin_low_REG_A = hZNAnoZNCTimes[iTstamp][k]->GetXaxis()->FindBin(region_lim_neg);
            Int_t bin_high_REG_A = hZNAnoZNCTimes[iTstamp][k]->GetXaxis()->FindBin(region_lim_pos);

            double counts_TH1A_TOT = hZNAnoZNCTimes[iTstamp][k]->Integral(bin_low_TOT_A, bin_high_TOT_A);
            double counts_TH1A_inside = hZNAnoZNCTimes[iTstamp][k]->Integral(bin_low_REG_A, bin_high_REG_A);
            double counts_TH1A_outside = counts_TH1A_TOT - counts_TH1A_inside;
    
            Double_t pass;
            if(strcmp(TrigType, "1ZNC") == 0 && g_vdm_Fill==10802) pass = counts_TH2_inside + counts_TH1C_inside - counts_TH1A_outside; //OO subtract satellites
            if(strcmp(TrigType, "1ZNC") == 0 && g_vdm_Fill==10782) pass = counts_TH2_inside + counts_TH1A_inside; //pO subtract satellites
            if(strcmp(TrigType, "MTVX") == 0) pass = counts_TH2_inside;
            Double_t tot;
            if(strcmp(TrigType, "1ZNC") == 0 && g_vdm_Fill==10802) tot = counts_TH2_TOT + counts_TH1C_TOT; //all that have t(ZNC) defined
            if(strcmp(TrigType, "1ZNC") == 0 && g_vdm_Fill==10782) tot = counts_TH2_TOT + counts_TH1A_TOT; //all that have t(ZNC) defined
            if(strcmp(TrigType, "MTVX") == 0) tot = counts_TH2_TOT; //I require both times for TVX=FT0A&&FT0C
            Double_t ratio = pass/tot;
            Double_t ratio_err = TMath::Sqrt(ratio*(1-ratio) / tot); //binomial error
            
            //Fill debug histograms
            if(iTstamp < n_separations){ //the n_separations first separation indices are for the x scan
                hCountsRatio_vs_x[k]->SetBinContent(iTstamp+1, ratio);
                hCountsRatio_vs_x[k]->SetBinError(iTstamp+1, ratio_err);
            }
            else{
                hCountsRatio_vs_y[k]->SetBinContent(iTstamp+1, ratio);
                hCountsRatio_vs_y[k]->SetBinError(iTstamp+1, ratio_err);
            }

            //Fill TTree variables
            if(iTstamp < n_separations){
                int idx_sep = iTstamp;
                correction_x[idx_sep] = ratio;
                correction_x_error[idx_sep] = ratio_err;
            }
            else{
                int idx_sep = iTstamp - n_separations;
                correction_y[idx_sep] = ratio;
                correction_y_error[idx_sep] = ratio_err;
            }

            //Print debug histograms for 1! BC; default is BC = 0, but it can be changed in this code
            Int_t index_to_print = 0;
            uint64_t tstamp_min = tstamps_start[iTstamp]; 
            uint64_t tstamp_max = tstamps_end[iTstamp];
            if(k == index_to_print)
            {
                const char* cName5 = Form("Satellites_ZNCnotZNA_timestamps_%llu-%llu", tstamp_min, tstamp_max);
                TCanvas* c5 = new TCanvas(cName5, cName5, 800, 700);
                c5->cd();
                c5->SetLogy();
                gPad->SetLeftMargin(0.15);
                //hZNCnoZNATimes[iTstamp]->SetStats(0);
                hZNCnoZNATimes[iTstamp][k]->GetXaxis()->SetRangeUser(-20,20);
                hZNCnoZNATimes[iTstamp][k]->GetXaxis()->SetTitle("ZNC arrival time (ns)");
                hZNCnoZNATimes[iTstamp][k]->SetMarkerStyle(8);
                hZNCnoZNATimes[iTstamp][k]->Draw("hist");
                int bin_max = hZNCnoZNATimes[iTstamp][k]->GetMaximumBin();
                double peak_val = hZNCnoZNATimes[iTstamp][k]->GetBinContent(bin_max);
                double ymin = gPad->GetUymin();
                double ymax = gPad->GetUymax();
                TLine *lL = new TLine(region_lim_neg,ymin,region_lim_neg,peak_val); lL->SetLineColor(kRed);  lL->Draw("same");
                TLine *lR = new TLine(region_lim_pos,ymin,region_lim_pos,peak_val); lR->SetLineColor(kRed);  lR->Draw("same");
                c5->Update();
                c5->Print(Form("../Fill-%d/Satellites%sforBC%d/SatellitesZNCnotZNA/%s.pdf", g_vdm_Fill, TrigType, k, cName5));
                
                const char* cName1 = Form("Satellites_ZED_sep%d_ibc%d", iTstamp, k);
                TCanvas* c1 = new TCanvas(cName1, cName1, 800, 700);
                c1->cd();
                gPad->SetRightMargin(0.15);
                gPad->SetLeftMargin(0.15);
                gStyle->SetPalette(kRainBow);
                hZDCTimes[iTstamp][k]->SetStats(0);
                hZDCTimes[iTstamp][k]->GetXaxis()->SetRangeUser(-20,20);
                hZDCTimes[iTstamp][k]->GetYaxis()->SetRangeUser(-20,20);
                hZDCTimes[iTstamp][k]->GetXaxis()->SetTitle("ZNC arrival time (ns)");
                hZDCTimes[iTstamp][k]->GetYaxis()->SetTitle("ZNA arrival time (ns)");
                hZDCTimes[iTstamp][k]->Draw("colz");
                l1->Draw("same");
                l2->Draw("same");
                l3->Draw("same");
                l4->Draw("same");
                c1->Update();
                c1->Print(Form("../Fill-%d/Satellites%sforBC%d/%s.pdf", g_vdm_Fill, TrigType, k, cName1));

                
                const char* cName3 = Form("Satellites_ZNAnotZNC_timestamps_%llu-%llu", tstamp_min, tstamp_max);
                TCanvas* c3 = new TCanvas(cName3, cName3, 800, 700);
                c3->cd();
                c3->SetLogy();
                gPad->SetLeftMargin(0.15);
                hZNAnoZNCTimes[iTstamp][k]->GetXaxis()->SetRangeUser(-20,40);
                hZNAnoZNCTimes[iTstamp][k]->GetXaxis()->SetTitle("ZNC arrival time (ns)");
                hZNAnoZNCTimes[iTstamp][k]->Draw("hist");
                int bin_max_C = hZNAnoZNCTimes[iTstamp][k]->GetMaximumBin();
                double ymin_C = gPad->GetUymin();
                double peak_val_C = hZNAnoZNCTimes[iTstamp][k]->GetBinContent(bin_max_C);
                TLine *lLC = new TLine(region_lim_neg,ymin_C,region_lim_neg,peak_val_C); lLC->SetLineColor(kRed);  lLC->Draw("same");
                TLine *lRC = new TLine(region_lim_pos,ymin_C,region_lim_pos,peak_val_C); lRC->SetLineColor(kRed);  lRC->Draw("same");
                c3->Update();
                c3->Print(Form("../Fill-%d/Satellites%sforBC%d//SatellitesZNAnotZNC/%s.pdf", g_vdm_Fill, TrigType, k, cName1));

            }

        }

        //Fill tree for each entry, which is the bc index
        corr_tree_x->Fill();
        corr_tree_y->Fill();

        //Output histograms for every bc
        TCanvas* c2 = new TCanvas("c2", Form("Scan %d, iBC%d, Counts in region / total", scan, k), 800, 700);
        c2->cd();
        gPad->SetLeftMargin(0.15);
        hCountsRatio_vs_x[k]->SetStats(0);
        hCountsRatio_vs_x[k]->SetMarkerStyle(8);
        hCountsRatio_vs_x[k]->SetMarkerColor(kBlue);
        hCountsRatio_vs_x[k]->SetLineColor(kBlue);
        hCountsRatio_vs_x[k]->GetYaxis()->SetRangeUser(0.0, 1.02);
        hCountsRatio_vs_x[k]->SetTitle(Form("%s trig., scan%d, iBC%d (%d), CountsRatio vs sep", TrigType, scan, k, bunches_I[k]));
        hCountsRatio_vs_x[k]->GetXaxis()->SetTitle("idx");
        hCountsRatio_vs_x[k]->GetYaxis()->SetTitle("Counts in / (Counts in + Counts out)");
        hCountsRatio_vs_x[k]->Draw("EP");
        hCountsRatio_vs_y[k]->SetMarkerStyle(8);
        hCountsRatio_vs_y[k]->SetMarkerColor(kRed);
        hCountsRatio_vs_y[k]->SetLineColor(kRed);
        hCountsRatio_vs_y[k]->Draw("EPSAME");
        TLegend *Lx= new TLegend(0.25, 0.82, 0.4, 0.88);
        Lx->SetLineWidth(0);
        Lx->SetFillColor(0);
        Lx->AddEntry(hCountsRatio_vs_x[k], "Scan x");
        Lx->Draw("SAME");
        TLegend *Ly= new TLegend(0.65, 0.82, 0.8, 0.88);
        Ly->SetLineWidth(0);
        Ly->SetFillColor(0);
        Ly->AddEntry(hCountsRatio_vs_y[k], "Scan y");
        Ly->Draw("SAME");
        c2->Update();
        c2->SaveAs(Form("../Fill-%d/Satellites%svsBC/HistCounts_scan%d_%s_iBC%d.pdf", g_vdm_Fill, TrigType, scan, TrigType, k));
        if(k<=1) c2->SaveAs(Form("../Fill-%d/Satellites%svsBC/HistCounts_scan%d_%s_iBC%d.root", g_vdm_Fill, TrigType, scan, TrigType, k));
    }

    //Save trees to files (depending on direction)
    CorrFile_x->cd();
    corr_tree_x->Write();
    CorrFile_x->Close();

    CorrFile_y->cd();
    corr_tree_y->Write();
    CorrFile_y->Close();

}
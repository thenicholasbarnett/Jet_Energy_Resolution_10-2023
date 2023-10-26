// this code takes in HiForestAOD_11.root, HiForestAOD_186.root, and HiForestAOD_187.root
// and prints out the number of total events in them as well as the number of passed events in them
// where passed events are events with |vz|<15 and have passed the following triggers:
// pPAprimaryVertexFilter, HBHENoiseFilterResultRun2Loose, pBeamScrapingFilter, HLT_HIAK4CaloJet80_v1
// this scripts also prints out how long it takes to run

// imports
#include "TROOT.h"
#include "TClass.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"

// PLOTTING FUNCTIONS

void ploth1d_1(TH1D *h, TString xtitle, TString ytitle, TString htitle){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle(htitle);
    TH1D *h_c = (TH1D*)h->Clone(htitle);
    h_c->Draw("e1p");
    h_c->SetMarkerStyle(20);
    h_c->SetMarkerColor(kBlack);
    h_c->SetLineColor(kBlack);
    h_c->GetXaxis()->CenterTitle(true);
    h_c->GetYaxis()->SetTitle(ytitle);
    h_c->GetYaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        c->SetLogy();
        h_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    }
    if(xtitle != "pt"){h_c->GetXaxis()->SetTitle(xtitle);}
}

void ploth2d_1(TH2D *h, TString xtitle, TString ytitle, TString htitle, TString option){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle(htitle);
    TH2D *h_c = (TH2D*)h->Clone(htitle);
    if(option == ""){h_c->Draw("e1p");}
    if(option != ""){h_c->Draw(option);}
    h_c->SetMarkerStyle(20);
    h_c->SetMarkerColor(kBlack);
    h_c->SetLineColor(kBlack);
    h_c->GetXaxis()->CenterTitle(true);
    h_c->GetYaxis()->SetTitle(ytitle);
    h_c->GetYaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        c->SetLogy();
        h_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    }
    if(xtitle != "pt"){h_c->GetXaxis()->SetTitle(xtitle);}

}

void ploth1d_2(TH1D *h0, TString h0name, TH1D *h1, TString h1name, TString xtitle, TString ytitle, TString htitle){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle(htitle);
    // h0
    TH1D *h0_c = (TH1D*)h0->Clone(h1name);
    h0_c->Draw("e1p");
    h0_c->SetMarkerStyle(20);
    h0_c->SetMarkerColor(kBlack);
    h0_c->SetLineColor(kBlack);
    h0_c->SetTitle(htitle);
    h0_c->GetXaxis()->CenterTitle(true);
    h0_c->GetYaxis()->SetTitle(ytitle);
    h0_c->GetYaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        h0_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");    
        c->SetLogy();}
    if(xtitle != "pt"){h0_c->GetXaxis()->SetTitle(xtitle);}
    // h1
    TH1D *h1_c = (TH1D*)h1->Clone(h1name);
    h1_c->Draw("e1psame");
    h1_c->SetMarkerStyle(21);
    h1_c->SetMarkerColor(kRed);
    h1_c->SetLineColor(kRed);
    // legend
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(h0_c,h0name,"pl");
    l->AddEntry(h1_c,h1name,"pl");
    l->Draw("same");
}

void ploth1d_3(TH1D *h0, TString h0name, TH1D *h1, TString h1name, TH1D *h2, TString h2name, TString xtitle, TString ytitle, TString htitle){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle(htitle);
    // h0
    TH1D *h0_c = (TH1D*)h0->Clone(h1name);
    h0_c->Draw("e1p");
    h0_c->SetMarkerStyle(20);
    h0_c->SetMarkerColor(kBlack);
    h0_c->SetLineColor(kBlack);
    h0_c->SetTitle(htitle);
    h0_c->GetXaxis()->CenterTitle(true);
    h0_c->GetYaxis()->SetTitle(ytitle);
    h0_c->GetYaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        h0_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");    
        c->SetLogy();}
    if(xtitle != "pt"){h0_c->GetXaxis()->SetTitle(xtitle);}
    // h1
    TH1D *h1_c = (TH1D*)h1->Clone(h1name);
    h1_c->Draw("e1psame");
    h1_c->SetMarkerStyle(21);
    h1_c->SetMarkerColor(kRed);
    h1_c->SetLineColor(kRed);
    // h1
    TH1D *h2_c = (TH1D*)h2->Clone(h2name);
    h2_c->Draw("e1psame");
    h2_c->SetMarkerStyle(21);
    h2_c->SetMarkerColor(kBlue);
    h2_c->SetLineColor(kBlue);
    // legend
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(h0_c,h0name,"pl");
    l->AddEntry(h1_c,h1name,"pl");
    l->AddEntry(h2_c,h2name,"pl");
    l->Draw("same");
}

void plotg1d_1(TGraph *h, TString xtitle, TString ytitle, TString htitle){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle(htitle);
    TGraph *h_c = (TGraph*)h->Clone(htitle);
    h_c->Draw("ALP");
    h_c->SetTitle("");
    h_c->SetMarkerStyle(20);
    h_c->SetMarkerColor(kBlack);
    h_c->SetLineColor(kBlack);
    h_c->GetXaxis()->CenterTitle(true);
    h_c->GetYaxis()->SetTitle(ytitle);
    h_c->GetYaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        c->SetLogy();
        h_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    }
    if(xtitle != "pt"){h_c->GetXaxis()->SetTitle(xtitle);}
}

void plotg1d_2(TGraph *h0, TString h0name, TGraph *h1, TString h1name, TString xtitle, TString ytitle, TString htitle){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle(htitle);
    TGraph *h0_c = (TGraph*)h0->Clone(htitle);
    TGraph *h1_c = (TGraph*)h1->Clone(htitle);
    // h0
    h0_c->Draw("apl");
    h0_c->SetTitle("");
    h0_c->SetMarkerStyle(20);
    h0_c->SetMarkerColor(kBlack);
    h0_c->SetLineColor(kBlack);
    h0_c->GetXaxis()->CenterTitle(true);
    h0_c->GetYaxis()->SetTitle(ytitle);
    h0_c->GetYaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        c->SetLogy();
        h0_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    }
    if(xtitle != "pt"){h0_c->GetXaxis()->SetTitle(xtitle);}
    // h1
    h1_c->Draw("pl same");
    h1_c->SetMarkerStyle(21);
    h1_c->SetMarkerColor(kRed);
    h1_c->SetLineColor(kRed);
    // legend
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(h0_c,h0name,"pl");
    l->AddEntry(h1_c,h1name,"pl");
    l->Draw("same");

}

void normalizeh(TH1D *h){
    double a = h->Integral();
    h->Scale(1/a);
}

void getthesliceinfo(TH1D *h, int in, Double_t sl[4][11]){
    sl[0][in] = h->GetStdDev();
    sl[1][in] = h->GetMean();
    sl[2][in] = h->GetStdDevError();
    sl[3][in] = h->GetMeanError();
}

// the script all runs in this function
void JER_A()
{
    // determining the inital time of the script
    clock_t ti = clock();

    // creating the histograms of interest
    TH1::SetDefaultSumw2();

    // creating some binning parameters
    // _h1dN is the Nth set of binning parameters for 1d hists of _
    double vzh1d0[3] = {40,-20,20};
    double pth1d0[3] = {100,20,1000};
    double etah1d0[3] = {50,-5,5};
    double etah1d1[3] = {25,-1.7,1.7};
    double ratio1d0[3] = {50,0,1};
    // making bins for some pT histograms
    const Int_t bin1 = 5;
    double pth1d1[bin1+1] = {60,80,100,120,200,1000};
    // making a binning for the jer vs refpt plot
    const Int_t bin2 = 12;
    double pth1d2[bin2+1] = {80,90,100,110,120,130,140,150,180,220,300,500};
    // making binning parameters for the normalized jer values
    double jerh1d0[3] = {10000,-1.2,1.2};

    // unweighted parameters
    // TH1D *hvzuw = new TH1D("vz_unweighted","",vzh1d0[0],vzh1d0[1],vzh1d0[2]);
    // TH1D *hpthatuw = new TH1D("hpthat_unweighted","",pth1d0[0],pth1d0[1],pth1d0[2]);
    // TH1D *hjtptuw = new TH1D("hjtpt_uw","",pth1d0[0],pth1d0[1],pth1d0[2]);
    // TH1D *hgenptuw = new TH1D("hgenpt_uw","",pth1d0[0],pth1d0[1],pth1d0[2]);
    // TH1D *hjeruw = new TH1D("hjeruw","",jerh1d0[0],jerh1d0[1],jerh1d0[2]);

    // established parameters
    TH1D *hvz = new TH1D("vz","",vzh1d0[0],vzh1d0[1],vzh1d0[2]);
    TH1D *hpthat = new TH1D("hpthat","",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hgenpt = new TH1D("hgenpt","",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hjtpt = new TH1D("hjtpt","",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hrefpt = new TH1D("hrefpt","",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hjer = new TH1D("hjer","",bin2,pth1d2);
    // TH1D *hjer = new TH1D("hjer","",jerh1d0[0],jerh1d0[1],jerh1d0[2]);
    TH2D *hjerref = new TH2D("hjerref","",pth1d0[0],pth1d0[1],pth1d0[2],jerh1d0[0],jerh1d0[1],jerh1d0[2]);
    TH2D *hjerjt = new TH2D("hjerjt","",pth1d0[0],pth1d0[1],pth1d0[2],jerh1d0[0],jerh1d0[1],jerh1d0[2]);

    // newish parameters
    TH1D *hrawpt = new TH1D("hrawpt","",bin2,pth1d2);
    TH1D *htrackMax = new TH1D("htrackMax","",bin2,pth1d2);
    TH1D *hratio = new TH1D("hratio","",ratio1d0[0],ratio1d0[1],ratio1d0[2]);
    TH1D *hjteta = new TH1D("hjteta","",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *hjteta_uc = new TH1D("hjteta_uc","",etah1d0[0],etah1d0[1],etah1d0[2]);
    

    // making slice plots
    TH1D *hslice0 = new TH1D("slice 0","",jerh1d0[0],jerh1d0[1],jerh1d0[2]);
    TH1D *hslice1 = new TH1D("slice 1","",jerh1d0[0],jerh1d0[1],jerh1d0[2]);
    TH1D *hslice2 = new TH1D("slice 2","",jerh1d0[0],jerh1d0[1],jerh1d0[2]);
    TH1D *hslice3 = new TH1D("slice 3","",jerh1d0[0],jerh1d0[1],jerh1d0[2]);
    TH1D *hslice4 = new TH1D("slice 4","",jerh1d0[0],jerh1d0[1],jerh1d0[2]);
    TH1D *hslice5 = new TH1D("slice 5","",jerh1d0[0],jerh1d0[1],jerh1d0[2]);
    TH1D *hslice6 = new TH1D("slice 6","",jerh1d0[0],jerh1d0[1],jerh1d0[2]);
    TH1D *hslice7 = new TH1D("slice 7","",jerh1d0[0],jerh1d0[1],jerh1d0[2]);
    TH1D *hslice8 = new TH1D("slice 8","",jerh1d0[0],jerh1d0[1],jerh1d0[2]);
    TH1D *hslice9 = new TH1D("slice 9","",jerh1d0[0],jerh1d0[1],jerh1d0[2]);
    TH1D *hslice10 = new TH1D("slice 10","",jerh1d0[0],jerh1d0[1],jerh1d0[2]);

    // a_ and b_ are the number of total events and accepted events respectively
    // n is the number of files that were looked at
    // tw is the total weight to normalize histograms with
    // tu is the total unweighted to normalize hists with
    Float_t a_=0, b_=0, n=0, tw=0, tu=0, nc=0, oc=0;

    // open files.txt to see the names of the files
    // ifstream myfile("public/filenamesAll.txt");
    // ifstream myfile("public/filenames50.txt");
    // ifstream myfile("public/filenames100.txt");
    ifstream myfile("files1.txt");
    // ifstream myfile("files.txt");
    string filename;

    // loop over the files by file names
    while (getline(myfile, filename))
    {
        // adding one to the number of files looked at
        n+=1;

        // turning the string filename into a TString so it passes through TFile
        TString q = filename;
        // TString q = "/eos/cms/store/group/phys_heavyions/jviinika/ppMC2017_5p02TeV_ptHat15_Dijet_Pythia8CP5_RunIIpp5Spring18DR_AOD_94X_allFiles_2023-01-05/0000/"+filename;

        cout<<q<<endl;
        
        // pointing fi0 and fi1 to the files holding the data of interest
        TFile *fi = TFile::Open(q,"read");
        fi->cd();

        // declaring variables
        int pPApVF, HBHENFRR2L, pBSF, HLT_HIAKCJ80v1;
        float vz, w;
        Float_t pthat;
        const Int_t nm = 200000;
        Float_t genpt[nm];
        Float_t jtpt[nm];
        Float_t refpt[nm];
        Float_t rawpt[nm];
        Float_t trackMax[nm];
        Float_t jteta[nm];

        // getting HltTree from fi, and the appropriate branches
        TTree *t0 = (TTree*)fi->Get("skimanalysis/HltTree");
        t0->SetBranchAddress("pPAprimaryVertexFilter",&pPApVF);
        t0->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHENFRR2L);
        t0->SetBranchAddress("pBeamScrapingFilter",&pBSF);

        // getting HiTree from fi, and the appropriate branches
        TTree *t1 = (TTree*)fi->Get("hiEvtAnalyzer/HiTree");
        t1->SetBranchAddress("vz",&vz);
        t1->SetBranchAddress("weight",&w);
        t1->SetBranchAddress("pthat",&pthat);

        // getting HltTree from fi, and the appropriate branches
        TTree *t2 = (TTree*)fi->Get("hltanalysis/HltTree");
        t2->SetBranchAddress("HLT_HIAK4CaloJet80_v1",&HLT_HIAKCJ80v1);

        // getting ak4PFJetAnalyzer from fi, and the appropriate branches
        TTree *t3 = (TTree*)fi->Get("ak4PFJetAnalyzer/t");
        t3->SetBranchAddress("genpt",genpt);
        t3->SetBranchAddress("jtpt",jtpt);
        t3->SetBranchAddress("refpt",refpt);
        t3->SetBranchAddress("rawpt",rawpt);
        t3->SetBranchAddress("trackMax",trackMax);
        t3->SetBranchAddress("jteta",jteta);

        // for loop for events in trees t0, t1, and t2
        for(unsigned int i=0; i<t0->GetEntries(); i++)
        // for(unsigned int i=0; i<20; i++)
        {
            // getting the entries
            t0->GetEntry(i);
            t1->GetEntry(i);
            t2->GetEntry(i);
            t3->GetEntry(i);

            // adding one to total events for every event
            a_+=1;

            // making trackmax/rawpt
            Double_t ratio = trackMax[0]/rawpt[0];

            // some print statements to see the trackMax and rawpt values
            // if(i<100){
            //     cout << "trackMax[0] is " << trackMax[0] << endl;
            //     cout << "rawpt[0] is " << trackMax[0] << endl;
            //     cout << "trackMax[0]/rawpt[0] is " << ratio << endl;
            // }

            // only events with |vz|<15 and all the triggers of interest are passed
            // if((TMath::Abs(vz)<15)&&(pPApVF==1)&&(HBHENFRR2L==1)&&(pBSF==1)&&(HLT_HIAKCJ80v1==1)&&(pthat>15)&&(refpt[0]>80)&&(TMath::Abs(jteta[0])<1.6)){
            if((TMath::Abs(vz)<15)&&(pPApVF==1)&&(HBHENFRR2L==1)&&(pBSF==1)&&(HLT_HIAKCJ80v1==1)&&(pthat>15)&&(refpt[0]>80)){
            // if((TMath::Abs(vz)<15)&&(pPApVF==1)&&(HBHENFRR2L==1)&&(pBSF==1)&&(HLT_HIAKCJ80v1==1)&&(pthat>15)&&(0.01<ratio)&&(ratio<0.98)){
                
                // event count without new cut
                // oc +=1;

                hjteta_uc->Fill(jteta[0],w);

                //event count with new cut
                // nc +=1;

                if(TMath::Abs(jteta[0])<1.6){
                    
                    // adding one to passed events iff all the conditionals are true
                    b_+=1;

                    // calculating the total weight by adding each weight to it
                    tw+=w;
                    tu+=1;

                    // filling the vz, w, and pthat histograms
                    // hvzuw->Fill(vz);
                    hvz->Fill(vz,w);
                    // hpthatuw->Fill(pthat);
                    hpthat->Fill(pthat,w);

                    // filling genpt, jtpt, and refpt
                    hgenpt->Fill(genpt[0],w);
                    // hgenptuw->Fill(genpt[0]);
                    hjtpt->Fill(jtpt[0],w);
                    // hjtptuw->Fill(jtpt[0]);
                    hrefpt->Fill(refpt[0],w);

                    // filling jet energy resolution
                    double jer = (jtpt[0] - refpt[0])/refpt[0];
                    hjer->Fill(jer,w);
                    // hjeruw->Fill(jer);
                    hjerref->Fill(refpt[0],jer,w);
                    hjerjt->Fill(jtpt[0],jer,w);
                    hjteta->Fill(jteta[0],w);

                    // filling hists
                    if((rawpt[0]>80)&&(rawpt[0]<500)&&(trackMax[0]>80)&&(trackMax[0]<500)){
                        hrawpt->Fill(rawpt[0],w);
                        htrackMax->Fill(trackMax[0],w);
                        hratio->Fill(ratio,w);
                    }

                    // filling the slices
                    // slice 0
                    if((refpt[0]>80) && (refpt[0]<90)){hslice0->Fill(jer,w);}
                    // slice 1
                    if((refpt[0]>90) && (refpt[0]<100)){hslice1->Fill(jer,w);}
                    // slice 2
                    if((refpt[0]>100) && (refpt[0]<110)){hslice2->Fill(jer,w);}
                    // slice 3
                    if((refpt[0]>110) && (refpt[0]<120)){hslice3->Fill(jer,w);}
                    // slice 4
                    if((refpt[0]>120) && (refpt[0]<130)){hslice4->Fill(jer,w);}
                    // slice 5
                    if((refpt[0]>130) && (refpt[0]<140)){hslice5->Fill(jer,w);}
                    // slice 6
                    if((refpt[0]>140) && (refpt[0]<150)){hslice6->Fill(jer,w);}
                    // slice 7
                    if((refpt[0]>150) && (refpt[0]<180)){hslice7->Fill(jer,w);}
                    // slice 8
                    if((refpt[0]>180) && (refpt[0]<220)){hslice8->Fill(jer,w);}
                    // slice 9
                    if((refpt[0]>220) && (refpt[0]<300)){hslice9->Fill(jer,w);}
                    // slice 10
                    if((refpt[0]>300) && (refpt[0]<500)){hslice10->Fill(jer,w);}

                    // PRINT STATEMENTS
                    // int lol= a_;
                    // if(lol%50000==0){
                    // cout << jer*w << " = jer * w, pthat is " << pthat << " and jer*w/tw = " << jer*w/tw << endl;
                    // cout << jer << " is jer, w is " << w << endl;
                    //     // cout << "genpt size is " << genpt.GetSize() << endl;
                    //     // cout << "jtpt size is " << jtpt.GetSize() << endl;
                    //     // cout << "jtpt[0] is " << jtpt[0] << endl;
                    //     // cout << "the number of t3 events is " << t3->GetEntries() << endl;
                    // }
                }   
            }
        }
    }

    // some print statements
    // float ps = 100*nc/oc;
    // cout << ps << " percent of the event counts survided the 0.01<trackMax/rawpt<0.98 cut" << endl;
    // cout << ps << " percent of the event counts survided the new cuts" << endl;

    TH1D *hjteta_c = (TH1D*)hjteta->Clone("hjteta_c");
    double hjteta_uc_w = hjteta_uc->Integral();
    hjteta_c->Scale(1/hjteta_uc_w);

    // normalizing the weighted histograms with the total weights
    // hvz->Scale(1/tw);
    // hpthat->Scale(1/tw);
    // hgenpt->Scale(1/tw);
    // hjtpt->Scale(1/tw);
    // hjer->Scale(1/tw);
    // hrefpt->Scale(1/tw);
    hjerref->Scale(1/tw);
    hjerjt->Scale(1/tw);
    // htrackMax->Scale(1/tw);

    // normalizing the unweighted histograms
    // hvzuw->Scale(1/tu);
    // hpthatuw->Scale(1/tu);
    // hgenptuw->Scale(1/tu);
    // hjtptuw->Scale(1/tu);
    // hjeruw->Scale(1/tu);
    // hratio->Scale(1/tw);
    // hjteta_uc->Scale(1/tw);

    // normalizing with the integral function
    normalizeh(hvz);
    normalizeh(hpthat);
    normalizeh(hgenpt);
    normalizeh(hjtpt);
    normalizeh(hjer);
    // normalizeh(hrefpt);
    // normalizeh(hjerref);
    // normalizeh(hjerjt);
    normalizeh(htrackMax);
    normalizeh(hrawpt);
    normalizeh(htrackMax);
    normalizeh(hratio);
    normalizeh(hjteta);
    normalizeh(hjteta_uc);

    // normalizing the slice histograms
    normalizeh(hslice0);
    normalizeh(hslice1);
    normalizeh(hslice2);
    normalizeh(hslice3);
    normalizeh(hslice4);
    normalizeh(hslice5);
    normalizeh(hslice6);
    normalizeh(hslice7);
    normalizeh(hslice8);
    normalizeh(hslice9);
    normalizeh(hslice10);

    // // double checking the hists are normalized with print statements 
    // cout << "hslice0 integral is " << hslice0->Integral() << endl;
    // cout << "hslice3 integral is " << hslice3->Integral() << endl;
    // cout << "jteta integral is " <<hjteta->Integral() << endl;
    // cout << "uncut jteta integral is " << hjteta_uc->Integral() << endl;
    // cout << "trackMax/rawpt integral is " << hratio->Integral() << endl;
    // cout << htrackMax->Integral() << endl;

    // making the 2d plot of standard deviation and mean of jer per bin as a function of refpt
    Double_t SlicesStdDevs[bin2-1], SlicesMeans[bin2-1], SlicesStdDevErrs[bin2-1], SlicesMeanErrs[bin2-1], sliceinfo[4][bin2-1];

    // using a function to get all the important info
    getthesliceinfo(hslice0, 0, sliceinfo);
    getthesliceinfo(hslice1, 1, sliceinfo);
    getthesliceinfo(hslice2, 2, sliceinfo);
    getthesliceinfo(hslice3, 3, sliceinfo);
    getthesliceinfo(hslice4, 4, sliceinfo);
    getthesliceinfo(hslice5, 5, sliceinfo);
    getthesliceinfo(hslice6, 6, sliceinfo);
    getthesliceinfo(hslice7, 7, sliceinfo);
    getthesliceinfo(hslice8, 8, sliceinfo);
    getthesliceinfo(hslice9, 9, sliceinfo);
    getthesliceinfo(hslice10, 10, sliceinfo);

    for(unsigned int i=0; i<11; i++){
        SlicesStdDevs[i] = sliceinfo[0][i];
        SlicesMeans[i] = sliceinfo[1][i];
        SlicesStdDevErrs[i] = sliceinfo[2][i];
        SlicesMeanErrs[i] = sliceinfo[3][i];
    }

    // plotting the stddev and mean as a function of refpt
    // getting x axis for plots
    double xax[bin2-1], xaxerr[bin2-1];
    for(unsigned int i=0; i<bin2-1; i++){
        double xmid = (pth1d2[i]+pth1d2[i+1])/2;
        xax[i] = xmid;
        xaxerr[i] = xmid-pth1d2[i];
    }

    // making graphs out of the slice info
    TGraph *gslicestddev0 = new TGraphErrors(bin2-1,xax,SlicesStdDevs,xaxerr,SlicesStdDevErrs);
    TGraph *gslicemean0 = new TGraphErrors(bin2-1,xax,SlicesMeans,xaxerr,SlicesMeanErrs);
    gslicemean0->SetMinimum(0);
    gslicemean0->SetMaximum(0.11);
    gslicemean0->GetXaxis()->SetLimits(80,500);
    TGraph *gslicestddev1 = new TGraph(bin2-1,xax,SlicesStdDevs);
    TGraph *gslicemean1 = new TGraph(bin2-1,xax,SlicesMeans);

    // PLOTTING
    //
    // getting rid of boxes of stats in plots
    gStyle->SetOptStat(0);

    // important parameter plots
    // ploth1d_1(hrawpt, "rawpt", "Probability", "rawpt");
    // ploth1d_1(htrackMax, "trackMax", "Probability", "trackMax");
    // ploth1d_1(hratio, "trackMax/rawpt", "Probability", "trackMax/rawpt");
    // ploth1d_1(hjteta, "jteta", "Probability", "jteta");
    // ploth1d_1(hjteta_uc, "jteta uncut", "Probability", "jteta uncut");
    // ploth1d_2(hjteta_uc, "uncut jteta", hjteta, "cut jteta", "#eta", "Probability", "#eta^{jet} before and after #eta cut, individually normalized");
    // ploth1d_2(hjteta_uc, "uncut jteta", hjteta_c, "cut jteta", "#eta", "Probability", "#eta^{jet} before and after #eta cut, normalized with precut scaling");

    // important graph plots
    // plotg1d_1(gslicestddev0, "P_{T}^{ref}", "Std. Dev. of JER", "Standard Deviation of JER vs P_{T}");
    // plotg1d_1(gslicemean0, "P_{T}^{ref}", "Mean of JER", "Mean of JER vs pT");
    plotg1d_2(gslicemean0, "Mean", gslicestddev0, "Std. Dev.", "P_{T}^{ref} [GeV/c]", "JER Probability", "Mean & Std. Dev. of JER vs pT");

    // TH1D *hslicestddev = new TH1D("hslicestddev","",bin2,pth1d2);
    // TH1D *hslicemean = new TH1D("hslicemean","",bin2,pth1d2);
    // for(unsigned int i=0; i<(bin2-1); i++){
    //     double s = (pth1d2[i]+pth1d2[i+1])/2;
    //     hslicestddev->Fill(s);
    //     hslicestddev->Fill(s);
    // }

    // // determining hjer min and max
    // int hjermaxbin = hjer->GetMaximumBin();
    // double hjermaxc = hjer->GetYaxis()->GetBinCenter(hjermaxbin);
    // double hjermaxu = hjer->GetYaxis()->GetBinUpEdge(hjermaxbin);
    // double hjermaxl = hjer->GetYaxis()->GetBinLowEdge(hjermaxbin);
    // int hjerminbin = hjer->GetMinimumBin();
    // double hjerminc = hjer->GetYaxis()->GetBinCenter(hjermaxbin);
    // double hjerminl = hjer->GetYaxis()->GetBinLowEdge(hjermaxbin);
    // double hjerminu = hjer->GetYaxis()->GetBinUpEdge(hjermaxbin);

    // PRINT STATEMENTS
    // // expressing passed/total as a percent
    // float p = b_/a_*100;
    //
    // // printing out the info of interest
    // cout<<"the following values were found by looking through "<< n << " root files" << endl;
    // cout<<"number of total events: "<< a_ << endl;
    // cout<<"number of passed events: "<< b_ << endl;
    // cout<<"percent of passed events: "<< p <<"%"<< endl;

    // common plots
    //
    ploth2d_1(hjerref, "P^{ref}_{T} [GeV/c]", "JER", "Jet Energy Resolution vs refpt", "COLZ");
    // ploth1d_1(hjer, "Jet Energy Resolution ", "Probability", "Jet Energy Resolution  ");
    
    // other plots
    //
    // ploth2d_1(hjerjt, "Probability", "P_{T} [GeV/c]", "Jet Energy Resolution","COLZ");
    // jer
    // TH1D *hjer_c1 = (TH1D*)hjer->Clone("Jet Energy Resolution");
    // hjer_c1->SetAxisRange(hjerminl*1.2, hjermaxu*1.2,"Y");
    // cout << hjerminu << " is the upper bin edge of jerminbin, and jermaxbin bin upper edge is " << hjermaxu << endl;
    // cout << hjerminl << " is the lower bin edge of jerminbin, and jermaxbin bin lower edge is " << hjermaxl << endl;
    // cout << hjerminc << " is the bin center of jerminbin, and jermaxbin bin center is " << hjermaxc << endl;
    // ploth1d_1(hjeruw, "#frac{P_{T}^{ref} - P_{T}^{gen}}{P_{T}^{ref}}", "Probability", "Jet Energy Resolution unweighted");
    // ploth1d_1(hjer_c1, "#frac{P_{T}^{ref} - P_{T}^{gen}}{P_{T}^{ref}}", "Probability", "Jet Energy Resolution");
    // vz
    // ploth1d_2(hvzuw, "Unweighted Vz", hvz, "Weighted Vz", "Position [cm]", "Probability", "Normalized Vertex Z Position");
    // pthat
    // ploth1d_2(hpthatuw, "Unweighted pthat", hpthat, "Weighted pthat", "pt", "Probability", "Normalized pthat");
    // genpt and jtpt
    // ploth1d_2(hjtptuw, "Unweighted jtpt", hjtpt, "Weighted jtpt", "pt", "Probability", "Normalized jtpt");
    // ploth1d_2(hgenptuw, "Unweighted genpt", hgenpt, "Weighted genpt", "pt", "Probability", "Normalized genpt");
    // ploth1d_2(hjtpt, "jtpt", hgenpt, "genpt", "pt", "Probability", "genpt & jtpt Overlay");
    // ploth1d_3(hgenpt, "genpt", hjtpt, "jtpt", hrefpt, "refpt", "pt", "Probability", "genpt, jtpt, and refpt Overlay");
    // ploth1d_2(hjtpt, "jtpt", hrefpt, "refpt", "pt", "Probability", "refpt & jtpt Overlay");
    
    // graph plots
    // plotg1d_1(gslicestddev1, "P_{T}^{ref}", "Std. Dev. of JER", "Standard Deviation of JER vs pT");
    // plotg1d_1(gslicemean1, "P_{T}^{ref}", "Mean of JER", "Mean of JER vs pT");
    
    // plotting the slices
    //
    ploth1d_1(hslice0, "JER", "Probability", "80 < P_{T} < 90");
    ploth1d_1(hslice1, "JER", "Probability", "90 < P_{T} < 100");
    ploth1d_1(hslice2, "JER", "Probability", "100 < P_{T} < 110");
    ploth1d_1(hslice3, "JER", "Probability", "110 < P_{T} < 120");
    ploth1d_1(hslice4, "JER", "Probability", "120 < P_{T} < 130");
    ploth1d_1(hslice5, "JER", "Probability", "130 < P_{T} < 140");
    ploth1d_1(hslice6, "JER", "Probability", "140 < P_{T} < 150");
    ploth1d_1(hslice7, "JER", "Probability", "150 < P_{T} < 180");
    ploth1d_1(hslice8, "JER", "Probability", "180 < P_{T} < 220");
    ploth1d_1(hslice9, "JER", "Probability", "220 < P_{T} < 300");
    ploth1d_1(hslice10, "JER", "Probability", "300 < P_{T} < 500");

    // determining the final time of the script
    clock_t tf = clock();
    // printing out the time between ti and tf in seconds
    cout<<"this script took "<<double(tf-ti)/(double)CLOCKS_PER_SEC<<" seconds to run"<<endl;
}

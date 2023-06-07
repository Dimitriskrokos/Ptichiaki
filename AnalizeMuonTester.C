#define AnalizeMuonTester_cxx
#include "AnalizeMuonTester.h"
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>


void AnalizeMuonTester::Loop() 
{
    fChain->SetBranchStatus("*",0);  // disable all branches
    fChain->SetBranchStatus("mmOnTrack_MuonLink",1);  // activate branchname
    fChain->SetBranchStatus("muons_eta",1);  // activate branchname
    fChain->SetBranchStatus("mmOnTrackResidualTrack",1);  // activate branchnames
    fChain->SetBranchStatus("mmOnTrackResidualTrackMS", 1);
    fChain->SetBranchStatus("mmOnTrackNStrips", 1);    
    fChain->SetBranchStatus("mmOnTrack_multiplet", 1);
    fChain->SetBranchStatus("mmOnTrack_gas_gap", 1);
    fChain->SetBranchStatus("mmOnTrackLocalPos_x", 1);
    
    TFile *rfile  = new TFile( "OutputHistograms.root","RECREATE");
    
    
    
    TH1F *h_muons_eta = new TH1F("muons_eta","muons eta ", 60,-3,3);
    TH1F *h_mmOnTrackNStrips = new TH1F("mmOnTrackNStrips","mmOnTrackNStrips ", 100,0,100);
    TH1F *h_muons_theta = new TH1F("h_muons_theta","muons theta ", 100,0,200);
    TH1F *hist1 = new TH1F("hist1","Strips between 8 and 11 degrees",100,0,60);
    TH1F *hist2 = new TH1F("hist2","Strips between 11 and 14 degrees",100,0,60);
    TH1F *hist3 = new TH1F("hist3","Strips between 14 and 17 degrees",100,0,60);
    TH1F *hist4 = new TH1F("hist4","Strips between 17 and 20 degrees",100,0,60);
    TH1F *hist5 = new TH1F("hist5","Strips between 20 and 23 degrees",100,0,60);
    TH1F *hist6 = new TH1F("hist6","Strips between 23 and 26 degrees",100,0,60);
    TH1F *hist7 = new TH1F("hist7","Strips between 26 and 29 degrees",100,0,60);
    TH1F *hist8 = new TH1F("hist8","Strips between 29 and 32 degrees",100,0,60);
    TH1F *h_x1_multiplet_1_station_1 = new TH1F("h_x1_multiplet_1_station_1","Hits at layers 1 at multiplet 1 station 1",100,0,5);
    TH1F *h_x1_multiplet_1_station_2 = new TH1F("h_x1_multiplet_1_station_2","Hits at layers 1 at multiplet 1 station 2",100,0,5);
    TH1F *h_x1_multiplet_2_station_1 = new TH1F("h_x1_multiplet_2_station_1","Hits at layers 1 at multiplet 2 station 1",100,0,5);
    TH1F *h_x1_multiplet_2_station_2 = new TH1F("h_x1_multiplet_2_station_2","Hits at layers 1 at multiplet 2 station 2",100,0,5);
    TH1F *h_x2_multiplet_1_station_1 = new TH1F("h_x2_multiplet_1_station_1","Hits at layers 2 at multiplet 1 station 1",100,0,5);    
    TH1F *h_x2_multiplet_1_station_2 = new TH1F("h_x2_multiplet_1_station_2","Hits at layers 2 at multiplet 1 station 2",100,0,5);
    TH1F *h_x2_multiplet_2_station_1 = new TH1F("h_x2_multiplet_2_station_1","Hits at layers 2 at multiplet 2 station 1",100,0,5);
    TH1F *h_x2_multiplet_2_station_2 = new TH1F("h_x2_multiplet_2_station_2","Hits at layers 2 at multiplet 2 station 2",100,0,5);    
    TH1F *h_x3_multiplet_1_station_1 = new TH1F("h_x3_multiplet_1_station_1","Hits at layers 3 at multiplet 1 station 1",100,0,5);       
    TH1F *h_x3_multiplet_1_station_2 = new TH1F("h_x3_multiplet_1_station_2","Hits at layers 3 at multiplet 1 station 2",100,0,5);       
    TH1F *h_x3_multiplet_2_station_1 = new TH1F("h_x3_multiplet_2_station_1","Hits at layers 3 at multiplet 2 station 1",100,0,5);       
    TH1F *h_x3_multiplet_2_station_2 = new TH1F("h_x3_multiplet_2_station_2","Hits at layers 3 at multiplet 2 station 2",100,0,5);       
    TH1F *h_x4_multiplet_1_station_1 = new TH1F("h_x4_multiplet_2_station_2","Hits at layers 4 at multiplet 1 station 1",100,0,5);       
    TH1F *h_x4_multiplet_1_station_2 = new TH1F("h_x4_multiplet_2_station_2","Hits at layers 4 at multiplet 1 station 2",100,0,5);       
    TH1F *h_x4_multiplet_2_station_1 = new TH1F("h_x4_multiplet_2_station_2","Hits at layers 4 at multiplet 2 station 1",100,0,5);       
    TH1F *h_x4_multiplet_2_station_2 = new TH1F("h_x4_multiplet_2_station_2","Hits at layers 4 at multiplet 2 station 2",100,0,5);       
    
    
    /*//TCanvas* canvas1 = new TCanvas("canvas1", "Histogram Canvas1", 1400, 1000);
    TCanvas* canvas2 = new TCanvas("canvas2", "1 station", 1400, 1000);
    TCanvas* canvas3 = new TCanvas("canvas3", "2 station", 1400, 1000);    
    //canvas1->Divide(3, 3);
    canvas2->Divide(2, 4);
    canvas3->Divide(2, 4);
    */
double Ntotal = 0.0;
double Npositive = 0.0;
double efficiency = 0.0;
    
    
    
///////////////////////     Pseudorabitity to theta angle   ///////////////////////////////////////////    
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast(); 
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) //analysis here
   {
   	Long64_t ientry = LoadTree(jentry);
   	
   	if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        //Analysis progress
        //if(jentry%100==0)std::cout << "Analysed a total of: " << jentry << " events" << std::endl;
        
        
         Float_t weight =1.;
         
         //ANALYSIS AND HISTROGRAM FILLING
         for(int i = 0; i < mmOnTrackNStrips->size(); i++)
         {
         	//cout <<mmOnTrackNStrips->size()<<"   "<< i <<endl;
         	int index_mu=mmOnTrack_MuonLink->at(i);
         	double eta_mu=muons_eta->at(index_mu);
		//cout << eta_mu << endl;
		double theta_mu = (2*TMath::ATan(TMath::Exp(-eta_mu)))*180/TMath::Pi(); 
		if(theta_mu >32.0)break;
		if(theta_mu <8.0)break;
		h_mmOnTrackNStrips->Fill(mmOnTrackNStrips->at(i));
		if(theta_mu >8.0 and theta_mu <11.0)   hist1 ->Fill(mmOnTrackNStrips->at(i)); 
		if(theta_mu >=11.0 and theta_mu <14.0) hist2 -> Fill(mmOnTrackNStrips->at(i));
		if(theta_mu >=14.0 and theta_mu <17.0) hist3 -> Fill(mmOnTrackNStrips->at(i));
		if(theta_mu >=17.0 and theta_mu <20.0) hist4 -> Fill(mmOnTrackNStrips->at(i));
		if(theta_mu >=20.0 and theta_mu <23.0) hist5 -> Fill(mmOnTrackNStrips->at(i));
		if(theta_mu >=23.0 and theta_mu <26.0) hist6 -> Fill(mmOnTrackNStrips->at(i));
		if(theta_mu >=26.0 and theta_mu <29.0) hist7 -> Fill(mmOnTrackNStrips->at(i));
		if(theta_mu >=29.0 and theta_mu <=32.0)hist8 -> Fill(mmOnTrackNStrips->at(i));
		 
		
		     	
		h_muons_theta -> Fill(theta_mu);         
         
         
         }
         
         
///////////////////////////////////     efficiency   //////////////////////////////////////////////////////////     
                 
         int muon_index[10];
         double Ntotal = 0.0;
         double Npositive = 0.0;
         double efficiency = 0.0;
         for(int m = 0; m < muons_eta->size(); m++)
         {
		double mu_eta=muons_eta->at(m);
	        double mu_pt=muons_pt->at(m);
	        double mu_phi=muons_phi->at(m);
		if(fabs(mu_eta)<1.3 || fabs(mu_eta)>2.7)continue;
		double number_of_mm_hits=0;
		int i = 0;

		
		 
		while(i-1 < mmOnTrackResidualTrack->size())
		{
			if(i >= mmOnTrackResidualTrack->size())
			{
				number_of_mm_hits += 1.0;//
				break;
			}
			Ntotal += 1.0;//The total hits before cuts
			int index_mu=mmOnTrack_MuonLink->at(i);
			i++;
			if(index_mu !=m)continue;
			double residual = mmOnTrackResidualTrack->at(i);
			if(residual <5) number_of_mm_hits += 1.0;
			i++;
		}
                Npositive += number_of_mm_hits;
		
		
         
         }
         
         efficiency = Npositive/Ntotal;
         cout<< "The efficiency value is " << efficiency <<endl; 
         
         
         
       
         
         
         
////////////////////////////////The maximun and minimum value of h_muons_eta
     
     /*double xMax = -std::numeric_limits<double>::infinity(); // Initialize with negative infinity
     int numBins = h_muons_theta->GetNbinsX();
     for (int i = 1; i <= numBins; ++i) {
        double binContent = h_muons_theta->GetBinContent(i);
        double binCenter = h_muons_theta->GetBinCenter(i);
        if (binContent > 0 && binCenter > xMax) {
           xMax = binCenter;
        }
     }

     cout << "The max value is: " << xMax << "degree" << endl;

   
    double xMin = std::numeric_limits<double>::infinity(); // Initialize with positive infinity

  
    for (int i = 1; i <= numBins; ++i) {
       double binContent = h_muons_theta->GetBinContent(i);
       double binCenter = h_muons_theta->GetBinCenter(i);
       if (binContent > 0 && binCenter < xMin) {
          xMin = binCenter;
       }
    }

cout << "The minimum value is: " << xMin << " degree" << endl;
*/

   
//////////////////////////////////////////////  Histograms draw    /////////////////////////////////////////////////////
 

/*////h_muons_theta ->Draw(); 
canvas1->cd(1);
h_mmOnTrackNStrips ->Draw();
   
canvas1->cd(2);
hist1->Draw();
    
canvas1->cd(3);
hist2->Draw();
     
canvas1->cd(4);
hist3->Draw();
   
canvas1->cd(5);
hist4->Draw();  
   
canvas1->cd(6);
hist5->Draw();      
   
canvas1->cd(7);
hist6->Draw(); 
   
canvas1->cd(8);
hist7->Draw();
   
canvas1->cd(9);
hist8->Draw(); 

//////////////  station 1

canvas2->cd(1);
h_x1_multiplet_1_station_1->Draw();


canvas2->cd(2);
h_x2_multiplet_1_station_1->Draw();

canvas2->cd(3);
h_x3_multiplet_1_station_1->Draw();

canvas2->cd(4);
h_x4_multiplet_1_station_1->Draw();


canvas2->cd(5);
h_x1_multiplet_2_station_1->Draw();

canvas2->cd(6);
h_x2_multiplet_2_station_1->Draw();

canvas2->cd(7);
h_x3_multiplet_2_station_1->Draw();

canvas2->cd(8);
h_x4_multiplet_2_station_1->Draw();
   
/////////// station 2
canvas3->cd(1);
h_x1_multiplet_1_station_2->Draw();


canvas3->cd(2);
h_x2_multiplet_1_station_2->Draw();

canvas3->cd(3);
h_x3_multiplet_1_station_2->Draw();

canvas3->cd(4);
h_x4_multiplet_1_station_2->Draw();


canvas3->cd(5);
h_x1_multiplet_2_station_2->Draw();

canvas3->cd(6);
h_x2_multiplet_2_station_2->Draw();

canvas3->cd(7);
h_x3_multiplet_2_station_2->Draw();

canvas3->cd(8);
h_x4_multiplet_2_station_2->Draw();





// At this point we write our histograms to the file and we close the file to finish
hist1->Write();
hist2->Write();
hist3->Write();
hist4->Write();
hist5->Write();
hist6->Write();
hist7->Write();
hist8->Write();
h_mmOnTrackNStrips->Write();
//////////////////////////// station 1
h_x1_multiplet_1_station_1->Write();
h_x2_multiplet_1_station_1->Write();
h_x3_multiplet_1_station_1->Write();
h_x4_multiplet_1_station_1->Write();
h_x1_multiplet_2_station_1->Write();
h_x2_multiplet_2_station_1->Write();
h_x3_multiplet_2_station_1->Write();
h_x4_multiplet_2_station_1->Write();
///////////////////////////// station 2
h_x1_multiplet_1_station_2->Write();
h_x2_multiplet_1_station_2->Write();
h_x3_multiplet_1_station_2->Write();
h_x4_multiplet_1_station_2->Write();
h_x1_multiplet_2_station_2->Write();
h_x2_multiplet_2_station_2->Write();
h_x3_multiplet_2_station_2->Write();
h_x4_multiplet_2_station_2->Write();




*/


//rfile->Close()    
}
}

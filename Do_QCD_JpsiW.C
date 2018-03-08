
#include "TH1.h"
#include "TFile.h"
#include <iostream>
#include <sstream>
#include "THStack.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "Riostream.h"
#include "TROOT.h"
#include "TClass.h"
#include "TTree.h"

//Define supporting functions
Double_t Difference_err(Double_t a, Double_t b);
Double_t Quotient_err(Double_t a, Double_t a_err, Double_t b, Double_t b_err);
void PrintDiff(TH1D **d, TH1D **mc, Double_t p, Double_t p_mc, int Hist_p, Double_t m, Double_t m_mc, int Hist_m, int ErrBin, Double_t *ErrArray);

//Main function
void Do_QCD_JpsiW(){
   // All backgrounds are calculated using MC, except QCD, which uses the ABCD method

   const char * status;
   try{

      status = "  //Declare Files and Histograms";
      cout << "" << endl;

      //input data file
      TFile df("/afs/cern.ch/user/d/dbertsch/work/Jpsi/Code/JpsiW/Output/out_data1.root");

      //input mc file
      TFile mcf("/afs/cern.ch/user/d/dbertsch/work/Jpsi/Code/JpsiW/Output/mc_all.root");

      //create output file
      TFile *ntup = new TFile("Output/ABCD_JpsiW.root","RECREATE");

      //create log file
      //  ofstream myfile;
      //  myfile.open ("QCD_out.txt");

      //h31=mu_plus, h32=mu_minus
      TH1D *data    = (TH1D*)df.Get("h01")->Clone("Data_Wmt");
      TH1D *data_p  = (TH1D*)df.Get("h31")->Clone("Data_Wmt_p");
      TH1D *data_m  = (TH1D*)df.Get("h32")->Clone("Data_Wmt_p");
      
      Float_t bins[] = { 25, 35, 45, 55, 70, 100, 200 };
      // Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1; // or just = 6
      Int_t  binnum = 6;

      //QCD background, which comes from data driven ABCD method, must be calculated
      TH1D *ABp  = new TH1D("ABp","Binned fake-factors A/B #mu+", binnum, bins);
      TH1D *ABm  = new TH1D("ABm","Binned fake-factors A/B #mu-", binnum, bins);
      TH1D *Dp   = new TH1D("Dp","Binned D (data-MC) #mu+",80,0,200);
      TH1D *Dm   = new TH1D("Dm","Binned D (data-MC) #mu-",80,0,200);
      TH1D *QCDp = new TH1D("QCDp","Binned QCD #mu+",80,0,200);
      TH1D *QCDm = new TH1D("QCDm","Binned QCD #mu-",80,0,200);

      //muon pt independent plots
      TH1D *Dp_I_data = new TH1D("Dp_I_data","All D (data) #mu+",80,0,200);
      TH1D *Dm_I_data = new TH1D("Dm_I_data","All D (data) #mu-",80,0,200);
      TH1D *Dp_I_MC   = new TH1D("Dp_I_MC","All D (MC) #mu+",80,0,200);
      TH1D *Dm_I_MC   = new TH1D("Dm_I_MC","All D (MC) #mu-",80,0,200);
      TH1D *Dp_I      = new TH1D("Dp_I","All D (data-MC) #mu+",80,0,200);
      TH1D *Dm_I      = new TH1D("Dm_I","All D (data-MC) #mu-",80,0,200);
      TH1D *QCDp_I    = new TH1D("QCDp_I","Independent QCD #mu+",80,0,200);
      TH1D *QCDm_I    = new TH1D("QCDm_I","Independent QCD #mu-",80,0,200);

      ABp->Sumw2();
      ABm->Sumw2();

      // create an array of histograms to be able to refer to them
      // by number later inside a loop
      TH1D *d[1000];
      TH1D *mc[1000];

      //Error Calculation and Yeild Tally
      Double_t p = 0;
      Double_t m = 0;
      Double_t p_mc = 0;
      Double_t m_mc = 0;
      Double_t ApTotal = 0;
      Double_t AmTotal = 0;
      Double_t BpTotal = 0;
      Double_t BmTotal = 0;
      Double_t DpTotal = 0;
      Double_t DmTotal = 0;
      
      // To store the Diff errors
      // plus:     1 = all events, 25   = bin 25 ...   35,   45,   55,   70,   100
      // minus: 1001 = all events, 1025 = bin 25 ... 1035, 1045, 1055, 1070, 10100
      Double_t A_err[2000] = {0};
      Double_t B_err[2000] = {0};
      Double_t D_err[2000] = {0};
      
      //first get all the histograms
      //Region A_data
      int BinStart = 1;
      int BinEnd = 6;
      int BinIncrement = 1;

      int Hist_p = 33;
      int Hist_m = 34;

      TString tempString_p;
      TString tempString_m;

      status = "  //Get A data";
      cout << "  //Get A data" << endl;

      // get all
      d[177] = (TH1D*)df.Get("h177");
      d[178] = (TH1D*)df.Get("h178");
      mc[177] = (TH1D*)mcf.Get("h177w");
      mc[178] = (TH1D*)mcf.Get("h178w");

      // get mu_pt binned version
      
      
      for(int b=BinStart;b<=BinEnd;b+=BinIncrement)
      {
         // First get data
         // Get hist number in char form
         tempString_p = "h";
         tempString_p += (long)Hist_p;
         tempString_m = "h";
         tempString_m += (long)Hist_m;

         d[Hist_p] = (TH1D*)df.Get(tempString_p);
         d[Hist_m] = (TH1D*)df.Get(tempString_m);

         // Then get mc
         tempString_p = "h";
         tempString_p += (long)Hist_p;
         tempString_p += "w";
         tempString_m = "h";
         tempString_m += (long)Hist_m;
         tempString_m += "w";

         mc[Hist_p] = (TH1D*)mcf.Get(tempString_p);
         mc[Hist_m] = (TH1D*)mcf.Get(tempString_m);

         Hist_p+=4;
         Hist_m+=4;
      }// End A Loop

      //Region B_data
      status = "  //Get B data";
      cout << "  //Get B data" << endl;
      Hist_p = 35;
      Hist_m = 36;

      // get all
      d[179] = (TH1D*)df.Get("h179");
      d[180] = (TH1D*)df.Get("h180");
      mc[179] = (TH1D*)mcf.Get("h179w");
      mc[180] = (TH1D*)mcf.Get("h180w");

      // get mu_pt binned version
      for(int b=BinStart;b<=BinEnd;b+=BinIncrement)
      {
         // First get data
         // Get hist number in char form
         tempString_p = "h";
         tempString_p += (long)Hist_p;
         tempString_m = "h";
         tempString_m += (long)Hist_m;

         d[Hist_p] = (TH1D*)df.Get(tempString_p);
         d[Hist_m] = (TH1D*)df.Get(tempString_m);
         
         // Then get mc
         tempString_p = "h";
         tempString_p += (long)Hist_p;
         tempString_p += "w";
         tempString_m = "h";
         tempString_m += (long)Hist_m;
         tempString_m += "w";

         mc[Hist_p] = (TH1D*)mcf.Get(tempString_p);
         mc[Hist_m] = (TH1D*)mcf.Get(tempString_m);
         
         Hist_p+=4;
         Hist_m+=4;
      }// End B Loop

      //Region D_data
      status = "  //Get D data";
      cout << "  //Get D data" << endl;
      Hist_p = 400;
      Hist_m = 401;

      // get all
      d[175] = (TH1D*)df.Get("h175");
      d[176] = (TH1D*)df.Get("h176");
      mc[175] = (TH1D*)mcf.Get("h175w");
      mc[176] = (TH1D*)mcf.Get("h176w");

      // get mu_pt binned version
      for(int b=BinStart;b<=BinEnd;b+=BinIncrement)
      {
         // First get data
         // Get hist number in char form
         tempString_p = "h";
         tempString_p += (long)Hist_p;
         tempString_m = "h";
         tempString_m += (long)Hist_m;
         
         d[Hist_p] = (TH1D*)df.Get(tempString_p);
         d[Hist_m] = (TH1D*)df.Get(tempString_m);
         
         // Then get mc
         tempString_p = "h";
         tempString_p += (long)Hist_p;
         tempString_p += "w";
         tempString_m = "h";
         tempString_m += (long)Hist_m;
         tempString_m += "w";
         
         mc[Hist_p] = (TH1D*)mcf.Get(tempString_p);
         mc[Hist_m] = (TH1D*)mcf.Get(tempString_m);
         
         Hist_p+=2;
         Hist_m+=2;
      }// End D_data Loop

      //Subtract off the weighted MC histograms from the data histograms
      //Region A
      status = "  //Subtract A";
      cout << endl << "  //Subtract A" << endl;

      //mu_pt independent version
      Hist_p = 177;
      Hist_m = 178;
      
      p = d[Hist_p]->Integral() ;
      p_mc = mc[Hist_p]->Integral();
      m = d[Hist_m]->Integral();
      m_mc = mc[Hist_m]->Integral();

      cout << "mu_pt independent version." << endl;
      PrintDiff(d, mc, p, p_mc, Hist_p, m, m_mc, Hist_m, 1, A_err);

      //mu_pt binned version
      Hist_p = 33;
      Hist_m = 34;
      for(int b=BinStart;b<=BinEnd;b+=BinIncrement)
      {
         p = d[Hist_p]->Integral() ;
         p_mc = mc[Hist_p]->Integral();
         m = d[Hist_m]->Integral();
         m_mc = mc[Hist_m]->Integral();
      
         cout << "  //  bin: " << b << endl;
         PrintDiff(d, mc, p, p_mc, Hist_p, m, m_mc, Hist_m, b, A_err);
      
         Hist_p+=4;
         Hist_m+=4;
      }// End A Loop

      //Region B
      status = "  //Subtract B";
      cout << "  //Subtract B" << endl;

      //mu_pt independent version
      Hist_p = 179;
      Hist_m = 180;
      
      p = d[Hist_p]->Integral() ;
      p_mc = mc[Hist_p]->Integral();
      m = d[Hist_m]->Integral();
      m_mc = mc[Hist_m]->Integral();
      
      cout << "mu_pt independent version." << endl;
      PrintDiff(d, mc, p, p_mc, Hist_p, m, m_mc, Hist_m, 1, B_err);

      //mu_pt binned version
      Hist_p = 35;
      Hist_m = 36;
      for(int b=BinStart;b<=BinEnd;b+=BinIncrement)
      {
         p = d[Hist_p]->Integral() ;
         p_mc = mc[Hist_p]->Integral();
         m = d[Hist_m]->Integral();
         m_mc = mc[Hist_m]->Integral();
         
         cout << "  //  bin: " << b << endl;
         PrintDiff(d, mc, p, p_mc, Hist_p, m, m_mc, Hist_m, b, B_err);

         Hist_p+=4;
         Hist_m+=4;
      }// End B Loop

      //Region D
      status = "  //Subtract D";
      cout << "  //Subtract D" << endl;

      //mu_pt independent version
      Hist_p = 175;
      Hist_m = 176;
      
      p = d[Hist_p]->Integral() ;
      p_mc = mc[Hist_p]->Integral();
      m = d[Hist_m]->Integral();
      m_mc = mc[Hist_m]->Integral();
      
      cout << "mu_pt independent version." << endl;
      PrintDiff(d, mc, p, p_mc, Hist_p, m, m_mc, Hist_m, 1, D_err);
      
      //save un-scaled region D
      Dp_I->Add(d[Hist_p]);
      Dm_I->Add(d[Hist_m]);

      //mu_pt binned version
      Hist_p = 400;
      Hist_m = 401;
      for(int b=BinStart;b<=BinEnd;b+=BinIncrement)
      {
      
         p = d[Hist_p]->Integral() ;
         p_mc = mc[Hist_p]->Integral();
         m = d[Hist_m]->Integral();
         m_mc = mc[Hist_m]->Integral();
         
         cout << "  //  bin: " << b << endl;
         PrintDiff(d, mc, p, p_mc, Hist_p, m, m_mc, Hist_m, b, D_err);

         Hist_p+=2;
         Hist_m+=2;
      }// End D Loop

      //Calculate the A/B fake-factors, and the (A/B)*D QCD fraction
      status = "  //calculate fake-factors";
      cout << "  //calculate fake-factors" << endl;

      //mu_pt independent version
      cout << "mu_pt independent version." << endl;

      Double_t A_p = d[177]->Integral();
      Double_t B_p = d[179]->Integral();
      Double_t A_m = d[178]->Integral();
      Double_t B_m = d[180]->Integral();
      Double_t ABratio_p = A_p/B_p;
      Double_t ABratio_m = A_m/B_m;
      cout << "A_p: " << A_p << " +/- " << A_err[1] << endl;
      cout << "A_m: " << A_m << " +/- " << A_err[1001] << endl;
      cout << "B_p: " << B_p << " +/- " << B_err[1] << endl;
      cout << "B_m: " << B_m << " +/- " << B_err[1001] << endl;
      cout << "A/B_p: " << ABratio_p << " +/- " << Quotient_err(A_p, A_err[1], B_p, B_err[1]) << endl;
      cout << "A/B_m: " << ABratio_m << " +/- " << Quotient_err(A_m, A_err[1001], B_m, B_err[1001]) << endl << endl;
      
      //Scale D histogram by the appropriate fake factor
      d[175]->Scale(ABratio_p);
      d[176]->Scale(ABratio_m);
      
      //save scaled data-MC for each region D in a histogram
      QCDp_I->Add(d[175]);
      QCDm_I->Add(d[176]);

      //mu_pt binned version
      int Hist_A_p = 33;
      int Hist_A_m = 34;
      int Hist_B_p = 35;
      int Hist_B_m = 36;
      int Hist_D_p = 400;
      int Hist_D_m = 401;
      int binIndex = 1; //first bin to fill

      for (int b=BinStart;b<=BinEnd;b+=BinIncrement){

         A_p = d[Hist_A_p]->Integral();
         B_p = d[Hist_B_p]->Integral();
         A_m = d[Hist_A_m]->Integral();
         B_m = d[Hist_B_m]->Integral();
         ABratio_p = A_p/B_p;
         ABratio_m = A_m/B_m;
         
         if (fabs(ABratio_p) > 10000) ABratio_p = 1;
         if (fabs(ABratio_m) > 10000) ABratio_m = 1;
         
         cout << "bin: " << b << endl;
         cout << "A_p: " << A_p << " +/- " << A_err[b] << endl;
         cout << "A_m: " << A_m << " +/- " << A_err[b+1000] <<  endl;
         cout << "B_p: " << B_p << " +/- " << B_err[b] <<  endl;
         cout << "B_m: " << B_m << " +/- " << B_err[b+1000] <<  endl;
         cout << "A/B_p: " << ABratio_p << " +/- " << Quotient_err(A_p, A_err[b], B_p, B_err[b]) << endl;
         cout << "A/B_m: " << ABratio_m << " +/- " << Quotient_err(A_m, A_err[b+1000], B_m, B_err[b+1000]) << endl << endl;
         
         ApTotal+=A_p;
         AmTotal+=A_m;
         BpTotal+=B_p;
         BmTotal+=B_m;
         DpTotal+=d[Hist_D_p]->Integral();
         DmTotal+=d[Hist_D_m]->Integral();
         
         //save fake factors in a histogram
         ABp->AddBinContent(binIndex,ABratio_p);
         ABm->AddBinContent(binIndex,ABratio_m);
         
         //save un-scaled region D
         Dp->Add(d[Hist_D_p]);
         Dm->Add(d[Hist_D_m]);
         
         //Scale each D histogram by the appropriate fake factor
         d[Hist_D_p]->Scale(ABratio_p);
         d[Hist_D_m]->Scale(ABratio_m);
         
         //save scaled data-MC for each region D in a histogram
         QCDp->Add(d[Hist_D_p]);
         QCDm->Add(d[Hist_D_m]);
         
         Hist_A_p += 4;
         Hist_A_m += 4;
         Hist_B_p += 4;
         Hist_B_m += 4;
         Hist_D_p += 2;
         Hist_D_m += 2;
         
         binIndex++;
      }
      
      cout << "Ap Bins Sum: " << ApTotal << endl;
      cout << "Am Bins Sum: " << AmTotal << endl;
      cout << "Bp Bins Sum: " << BpTotal << endl;
      cout << "Bm Bins Sum: " << BmTotal << endl;
      cout << "Dp Bins Sum: " << DpTotal << endl;
      cout << "Dm Bins Sum: " << DmTotal << endl << endl;

      cout << "Writing files" << endl;

      ABp->Write();
      ABm->Write();
      Dp->Write();
      Dm->Write();
      QCDp->Write();
      QCDm->Write();

      //for debugging
//      Dp_I_data->Write();
//      Dm_I_data->Write();
//      Dp_I_MC->Write();
//      Dm_I_MC->Write();
      
      Dp_I->Write();
      Dm_I->Write();
      QCDp_I->Write();
      QCDm_I->Write();
      
//      TCanvas c1("c1", "First canvas", 400, 400);
//      c1.Divide(2,2);
//      c1.cd(1);
//      Dp_I->DrawCopy();
//      c1.cd(2);
//      Dm_I->DrawCopy();
//      c1.cd(3);
//      QCDp_I->DrawCopy();
//      c1.cd(4);
//      QCDm_I->DrawCopy();
//      //c1.Save("plots.pdf");
      
      //Close output file
      ntup->Close();
      //   myfile.close();

   }// end try
   
//error handling
catch( char * e){
   cout << "Loop: CatchError: " << status << endl;
   throw e;
}
catch(...){
   cout << "Loop: CatchError: " << status << endl;
   throw;
}
   
}//end function Do_QCD()


Double_t Difference_err(Double_t a, Double_t b){
   Double_t a_err = TMath::Sqrt(a);
   Double_t b_err = TMath::Sqrt(b);
   return TMath::Sqrt(((a_err)*(a_err))+((b_err)*(b_err)));
}

Double_t Quotient_err(Double_t a, Double_t a_err, Double_t b, Double_t b_err){
   return fabs(a/b)*TMath::Sqrt(((a_err/a)*(a_err/a))+((b_err/b)*(b_err/b)));
}

void PrintDiff(TH1D **d, TH1D **mc, Double_t p, Double_t p_mc, int Hist_p, Double_t m, Double_t m_mc, int Hist_m, int ErrBin, Double_t *ErrArray){

   cout << "Using hist: " << Hist_p << "," << " data_p: " << p << ", mc_p: " << p_mc << " -> ";
   d[Hist_p]->Add(d[Hist_p],mc[Hist_p],1,-1);
   cout << "Diff: " << d[Hist_p]->Integral() << " +/- " << Difference_err(p,p_mc) << endl;
   ErrArray[ErrBin] = Difference_err(p,p_mc);
   
   cout << "Using hist: " << Hist_m << "," << " data_m: " << m << ", mc_m: " << m_mc << " -> ";
   d[Hist_m]->Add(d[Hist_m],mc[Hist_m],1,-1);
   cout << "Diff: " << d[Hist_m]->Integral() << " +/- " << Difference_err(m,m_mc) << endl << endl;
   ErrBin+=1000;
   ErrArray[ErrBin] = Difference_err(m,m_mc);
   
}






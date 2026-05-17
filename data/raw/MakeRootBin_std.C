#include <string>
#define File "Spectra_A.root"
#define N 10
#define prefix "a_series/A"
void MakeASCII(TH1D* h){
  std::string name = h->GetName();
  name = "c_series_ascii/"+name+".txt";
  std::ofstream out(name);
  for(int i = 0; i< h->GetXaxis()->GetNbins();i++){
     out<<i<<std::setw(15)<<h->GetBinContent(i)<<std::endl;
  }

}
TH1D* ReadFile(std::string file,std::string name)
{
  ifstream is;
  is.open (file, ios::binary );
  is.seekg (0, ios::end);
  int length = is.tellg();
  is.seekg (0, ios::beg);
   int * buffer = new  int [length/4+1000];
  std::cout<<file<<"  "<<length<<std::endl;
  is.read ((char*)buffer,length);
  is.close();
  std::string hname = name+"_hist";//0.19632445+N*0.10451050

 //  TH1D* h = new TH1D(hname.c_str(),"",4096,0.19632445,0.10451050*4096+0.19632445);
   TH1D* h = new TH1D(hname.c_str(),"",4096,0.19637099,0.19637099+0.10501098*4096);
   for(int i =0; i < 4096;i++) {h->SetBinContent(i,buffer[i])/*std::cout<<i<<"  "<<buffer[i]<<std::endl*/;}
   return h;
}

void MakeRootBin_std(void)
{
  int length;
  int * buffer;
 
  TFile* f = new TFile(File,"recreate");


 TH1D* h_cnt;
 TH1D* h_cnt2;
 h_cnt = new TH1D("Temporal","",1000,0,1000);
 h_cnt2 = new TH1D("Temporal_low","",1000,0,1000);
 TH1D* h_tot;
 TH1D* h_tot1;
 TH1D* h_tot2;

 for(int i = 1;i<N;i++)
 {
  stringstream ss;
  ss << i;
  string str = ss.str();

  stringstream ss2;
  ss2 << i+1;
  string str2 = ss2.str();
//calib_pb1.txt  MakeRootBin_std.C  тигель1.txt  тигель2.txt
  std::string n1 = prefix+str2;
  std::string n2 = prefix+str;
  std::cout<<n1<<"  "<<n2<<std::endl;
  std::string n_h = "spectrum_"+str2;
  std::string n_h2 = "spectrum_"+str;
  TH1D* h1 = ReadFile(n1,n_h);
  TH1D* h2 = ReadFile(n2,n_h2);
 
  if(i==1) {
      h_tot = (TH1D*) h1->Clone("total_spectrum");
      h_tot1 = (TH1D*) h1->Clone("total_spectrum1");
      h_tot2 = (TH1D*) h1->Clone("total_spectrum2");
      h_tot->Reset();
      h_tot1->Reset();
      h_tot2->Reset();
  }
   h1->Add(h2,-1);
//  h->Add(ReadFile(n2.c_str()),-1);
//  h1->SetName("calib");
//  h1->Write();
  h1->Write();
  h_cnt->SetBinContent(i,h1->Integral(1427,2760));
  h_cnt2->SetBinContent(i,h1->Integral(0,100));
  h_tot->Add(h1);
  if(i<96)  h_tot1->Add(h1);
  if(i>95 && i < 191)  h_tot2->Add(h1);
 // h3->Write();
 // MakeASCII(h1);
  MakeASCII(h1);
//  MakeASCII(h3);
 }
  h_tot->Write();
  h_tot1->Write();
  h_tot2->Write();

  h_cnt->Write();
//  TH1D* h_cal = ReadFile("a_series/CAM2.DAT","calib");
//  h_cal->Write();
//  f->Close();
}  

 

double ConvolutionIntegral(double,TH1D*,TH1D*);


const int numtrials = 1e6;

void LandauGaussian()
{
  TF1* funG = new TF1("funG", "gaus", -10, 20);
  TF1* funL = new TF1("funL", "landau", -10, 20);

  TH1D* histG = new TH1D("histG", "", 300, -10, 20);
  TH1D* histL = new TH1D("histL", "", 300, -10, 20);
  TH1D* histCombine = new TH1D("histCombine", "", 300, -10, 20);

  funG->SetParameter(0,1.0);
  funG->SetParameter(1,0.0);
  funG->SetParameter(2,1.0);
  funL->SetParameter(0,1.0);
  funL->SetParameter(1,0.0);
  funL->SetParameter(2,1.0);

  for(int i = 0; i<numtrials; ++i)
    {
      double G = funG -> GetRandom();
      double L = funL -> GetRandom();
      double C = G+L;
      histG -> Fill(G);
      histL -> Fill(L);
      histCombine -> Fill(C);
    }

  TH1D* histConvolve = new TH1D("histConvolve", "", 300, -10, 20);

  for ( int i = 0; i < histConvolve->GetNbinsX(); ++i)
    {
      double t = histConvolve->GetBinCenter(i+1);
      double convolution = ConvolutionIntegral(t,histG,histL);
      histConvolve->SetBinContent(i+1,convolution);
    }


  TCanvas* c1 = new TCanvas("c1", "");


  histG -> Draw();
  c1 -> Print("figures/lc_histG.png");
  histL -> Draw();
  c1 -> Print("figures/lc_histL.png");
  histCombine -> Draw();
  c1 -> Print("figures/lc_histCombine.png");
  histConvolve -> Draw();
  c1 -> Print("figures/lc_histConvole.png");



  TF1 *fgumbel = new TF1("fgumbel","([0]/sqrt(6.28))*TMath::Exp(-0.5*((x-[1])/[2] + TMath::Exp(-(x-[1]/[2]))))",-10,20);
  fgumbel->SetParameter(0,numtrials/100.0);
  fgumbel->SetParameter(1,0.0);
  fgumbel->SetParameter(2,1.0);

  TF1 *funcG = new TF1("funcG","gaus",-10, 10);
  TF1 *funcL = new TF1("funcL","landau",-10, 20);

  histG -> Fit("funcG","R");
  c1 -> Print("figures/lc_histG_fit.png");
  histL -> Fit("funcL","R");
  c1 -> Print("figures/lc_histL_fit.png");
  histCombine -> Fit("fgumbel", "R");
  c1 -> Print("figures/lc_histCombine_fit_gumbel.png");
  histConvolve -> Fit("fgumbel", "R");
  c1 -> Print("figures/lc_histConvolve_fit_gumbel.png");

  // ---

  TF1 *gengaus = new TF1("gengaus","[0]*exp(-pow(fabs(x-[1]),[3])/[2])",-10,20);
  gengaus->SetParameter(0,fgumbel->GetParameter(0));
  gengaus->SetParameter(1,fgumbel->GetParameter(1));
  gengaus->SetParameter(2,fgumbel->GetParameter(2));
  gengaus->SetParameter(3,2.2);

  TF1 *skewgaus = new TF1("skewgaus","[0]*exp(-pow(fabs(x-[1]),2)/[2])*(0.5+0.5*(TMath::Erf([3]*(x-[1])/[2])))",-10,20);
  skewgaus->SetParameter(0,fgumbel->GetParameter(0));
  skewgaus->SetParameter(1,fgumbel->GetParameter(1));
  skewgaus->SetParameter(2,fgumbel->GetParameter(2));
  skewgaus->SetParameter(3,0.1);

  TF1 *skewgengaus = new TF1("skewgengaus","[0]*exp(-pow(fabs(x-[1]),[3])/[2])*(0.5+0.5*(TMath::Erf([4]*(x-[1])/[2])))",-10,20);
  skewgengaus->SetParameter(0,fgumbel->GetParameter(0));
  skewgengaus->SetParameter(1,fgumbel->GetParameter(1));
  skewgengaus->SetParameter(2,fgumbel->GetParameter(2));
  skewgengaus->SetParameter(3,2.2);
  skewgengaus->SetParameter(4,0.1);

  //histCombine -> Fit("skewgaus", "R");
  // c1 -> Print("figures/lc_histCombine_skewgaus.png");

  histCombine -> Fit("skewgengaus", "R");
  c1 -> Print("figures/lc_histCombine_skewgengaus.png");

  histConvolve -> Fit("skewgengaus", "R");
  c1 -> Print("figures/lc_histConvolve_skewgengaus.png");

  TF1* landaugaus = new TF1("landaugaus","landau(0)*gaus(3)",-10,20);
  landaugaus->SetParameter(0,funL->GetParameter(0));
  landaugaus->SetParameter(1,funL->GetParameter(1));
  landaugaus->SetParameter(2,funL->GetParameter(2));
  landaugaus->SetParameter(3,funG->GetParameter(0));
  landaugaus->SetParameter(4,funG->GetParameter(1));
  landaugaus->SetParameter(5,funG->GetParameter(2));

  histCombine -> Fit("landaugaus", "R");
  c1 -> Print("figures/lc_histCombine_landaugaus.png");

  histConvolve -> Fit("landaugaus", "R");
  c1 -> Print("figures/lc_histConvolve_landaugaus.png");

}

double ConvolutionIntegral(double t, TH1D* hf, TH1D* hg)
{
  // limitation is that f and g need to have same range and number of bins...
  int nbins = hf->GetNbinsX();
  double lowerlimit = hg->GetBinLowEdge(1);
  double upperlimit = hg->GetBinLowEdge(nbins+2);
  double integral = 0;
  for ( int i = 0; i < nbins; ++i )
    {
      double tau = hf->GetBinCenter(i+1);
      double dtau = hf->GetBinWidth(i+1);
      double t_tau = t - tau;
      double f_tau = hf->GetBinContent(hf->FindBin(tau));
      double g_t_tau = hg->GetBinContent(hg->FindBin(t_tau));
      double integrand = f_tau * g_t_tau;
      bool inbounds = ( t_tau > lowerlimit && t_tau < upperlimit );
      if ( inbounds ) integral += integrand*dtau;
    }
  return integral;
}

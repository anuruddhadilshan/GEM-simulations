//First simulation of the stage 2 simulations.
//Creates GEM signal events using the APV function. Allow the parameters such as MaxADC, t0, and tau to be distributed randomly (use a Landau dis for Max ADC and gaus distribution for t0 and tau?)
//No pedestal noise and time jitter added. Also NO random background events included at this point.
//Aplies Time Deconvolution to those events. 
//Display a few events with the results and make the basic set of histograms to study the performance of the Time Deconvolution. Ideally we should see very good results here as this is the most simplest case for time deconvolution (I think).

// These are the Fit Values obtained by fitting the six time sample data from GMn SBS-4
const Double_t t0_Mean {3.173};
const Double_t t0_Sigma {8.363};
const Double_t maxADC_MPV {180.3};
const Double_t maxADC_Sigma {66.4};
const Double_t tau_Mean {56};
const Double_t tau_Sigma {24.42};

const Double_t timeSampleValues[6] {12,36,60,84,108,132};
const Int_t nTimeSamples {6};

//Function to generate a simulated APV signal with No pedestal noise added.
void sim_APVSignal_noPed(Double_t signal_noPed[6],Double_t& l_t0,Double_t& l_tau,Double_t& l_maxADC)
{
	TRandom2* rand_t0 = new TRandom2(0);
	Double_t t0 {rand_t0->Gaus(t0_Mean,t0_Sigma)};
	TRandom2* rand_tau = new TRandom2(0);
	Double_t tau {rand_tau->Gaus(tau_Mean,tau_Sigma)};
	TRandom2* rand_maxADC = new TRandom2(0);
	Double_t maxADC {rand_maxADC->Landau(maxADC_MPV,maxADC_Sigma)};

	l_t0 = t0;
	l_tau = tau;
	l_maxADC = maxADC;

	for (Int_t i_timesamp=0; i_timesamp<nTimeSamples; ++i_timesamp)
	{
		Double_t arg {(timeSampleValues[i_timesamp]-t0)/tau};
		signal_noPed[i_timesamp] = maxADC*arg*TMath::Exp(-arg+1);
	}
}
 
//Function to generate a simulated APV signal WITH pedestal noise added.
void sim_APVSignal_withPed(Double_t signal_withPed[6])
{
	TRandom2* rand_t0 = new TRandom2(0);
	Double_t t0 {rand_t0->Gaus(t0_Mean,t0_Sigma)};
	TRandom2* rand_tau = new TRandom2(0);
	Double_t tau {rand_tau->Gaus(tau_Mean,tau_Sigma)};
	TRandom2* rand_maxADC = new TRandom2(0);
	Double_t maxADC {rand_maxADC->Landau(maxADC_MPV,maxADC_Sigma)};
	TRandom2* rand_ped_noise = new TRandom2(0);

	for (Int_t i_timesamp=0; i_timesamp<nTimeSamples; ++i_timesamp)
	{
		Double_t ped_noise = rand_ped_noise->Gaus(0,20); //Pedestal noise with a gaussian distribution with mean 0 and std.dev 20.
		Double_t arg {(timeSampleValues[i_timesamp]-t0)/tau};
		signal_withPed[i_timesamp] = maxADC*arg*TMath::Exp(-arg+1)+ped_noise;
	}
}

void make_pdf_for_displayed_events(const Int_t nevents_display,TGraph* gr[nevents_display],const TString output_file_name)
{
	TCanvas* c[nevents_display];

	TString pdffilename = output_file_name;
	TString openfilename = pdffilename+"(";
	TString closefilename = pdffilename+")";

	Double_t lmargin=0.15;
  	Double_t rmargin=0.15;
    Double_t bmargin=0.15;
    Double_t tmargin=0.09;

	for (Int_t ievent_display = 0; ievent_display < nevents_display; ++ievent_display)
	{
		TString event = "Event Number: "+std::to_string(ievent_display);
		c[ievent_display] = new TCanvas(event, event);
		c[ievent_display]->SetGrid();
		gr[ievent_display]->SetMarkerStyle(21);
		gr[ievent_display]->SetMarkerColor(1);
		gr[ievent_display]->SetMarkerSize(1);
		gr[ievent_display]->SetTitle(event);
		gr[ievent_display]->GetXaxis()->SetTitle("Time (ns)");
		gr[ievent_display]->GetYaxis()->SetTitle("ADC");
		gr[ievent_display]->Draw("AP");
		//c[ievent_display]->BuildLegend();
		
		if(ievent_display == 0) c[ievent_display]->Print(openfilename,"Title:"+event);
		else	if (ievent_display == nevents_display-1) c[ievent_display]->Print(closefilename,"Title:"+event);
		else c[ievent_display]->Print(pdffilename,"Title:"+event);
	}
}

void timedeconv_stage2_sim1(const Int_t nevents_simulate=10000,const Int_t nevents_display=5)
{
	
	Double_t apv_signal_noPed[6] {0.};
	Double_t l_t0 {0.};
	Double_t l_tau {0.};
	Double_t l_maxADC {0.};

	Double_t apv_signal_withPed[6] {0.};
	Double_t l_ped_noise {0.};

	TH1D* h1_t0 = new TH1D("h1_t0","t0;t0(ns)",300,-150,150);
	TH1D* h1_tau = new TH1D("h1_tau","tau;tau(ns)",600,-100,500);
	TH1D* h1_maxADC = new TH1D("h1_maxADC","maxADC;ADC",350,-500,3000);
	
	TGraph* gr_SimAPVSigNoPed[nevents_display];
	TGraph* gr_SimAPVSigWithPed[nevents_display];

	Int_t ievent_display {0};


	for (Int_t ievent_simulate=0; ievent_simulate<nevents_simulate; ++ievent_simulate) 
	{
		sim_APVSignal_noPed(apv_signal_noPed,l_t0,l_tau,l_maxADC);
		sim_APVSignal_withPed(apv_signal_withPed);

		if (ievent_display<nevents_display && ievent_simulate%9==0) // Just a way to randomly select some simulated events to be displayed.
		{
			gr_SimAPVSigNoPed[ievent_display] = new TGraph(nTimeSamples,timeSampleValues,apv_signal_noPed);
			gr_SimAPVSigWithPed[ievent_display] = new TGraph(nTimeSamples,timeSampleValues,apv_signal_withPed);
			++ievent_display;
		}
			
		h1_t0->Fill(l_t0);
		h1_tau->Fill(l_tau);
		h1_maxADC->Fill(l_maxADC);
		
		if (ievent_simulate%1000==0) std::cout << "Simulating Event Number: " << ievent_simulate <<'\n';
	}

	const TString displayFileName_noPedSimEvnts {"APV_Signals_with_No_Pedestal_Noise.pdf"};
	make_pdf_for_displayed_events(nevents_display,gr_SimAPVSigNoPed,displayFileName_noPedSimEvnts);
	const TString displayFileName_withPedSimEvnts {"APV_Signals_with_Pedestal_Noise.pdf"};
	make_pdf_for_displayed_events(nevents_display,gr_SimAPVSigWithPed,displayFileName_withPedSimEvnts);

	gStyle->SetOptFit(1);

	TCanvas* a1 = new TCanvas();
	h1_t0->Fit("gaus");
	h1_t0->Draw();

	TCanvas* a2 = new TCanvas();
	h1_tau->Fit("gaus");
	h1_tau->Draw();

	TCanvas* a3 = new TCanvas();
	h1_maxADC->Fit("landau");
	h1_maxADC->Draw();

}

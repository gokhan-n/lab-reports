{
  
  std::vector<std::string> filenames = {
        "ba133 1s.txt",
        "ba133 10s.txt",
        "cs137 1s.txt",
        "cs137 10s.txt",
       
     };

  std::vector<std::string> titles = {
        "Barium-133, 1 second",
        "Barium-133, 10 seconds",
        "Caesium-137, 1 second",
        "Caesium-137, 10 seconds"
       
     };

  std::vector<std::vector<double>> histoRanges = {
    {0,10},{0,40},{0,20},{0,200}
  };
  std::vector<std::vector<double>> fitRanges = {
    {0.0,8},{0.0,30},{0.0,15},{0.0,120}
  };
  float histogramBins [4] = {10,40,20,100};

  std::vector<float> means;
  std::vector<float> sigmas;
  float kai2Gauss[4];
  float kai2Poisson[4];



  for (int i = 0; i < filenames.size(); ++i) {
    
    std::string gaussianTitle = "Gaussian(red line) Fit and Poisson(blue line) fit for " + titles[i];
    std::string poissonTitle = "Poisson Fit for " + titles[i];


    TH1F *histoGauss = new TH1F(filenames[i].c_str(),gaussianTitle.c_str(),histogramBins[i], histoRanges[i][0], histoRanges[i][1]);
    histoGauss->GetXaxis()->SetTitle("Number of Counts");
    histoGauss->GetYaxis()->SetTitle("Weight Of Counts");
    TH1F *histoPoisson = new TH1F(filenames[i].c_str(),poissonTitle.c_str(),histogramBins[i], histoRanges[i][0], histoRanges[i][1]);  
    histoPoisson->GetXaxis()->SetTitle("Number of Counts");
    histoPoisson->GetYaxis()->SetTitle("Weight Of Counts");
  
    std::ifstream infile(filenames[i]);
    float value;

    while (infile >> value) {
          histoGauss->Fill(value);
          histoPoisson->Fill(value);
    }
    means.push_back(histoGauss->GetMean());
    sigmas.push_back(histoGauss->GetStdDev());
    

    //gStyle->SetOptFit(1111);

    TF1 *gauss = new TF1("gauss","gaus");
    TF1 *poisson = new TF1("poisson","[0]*TMath::Poisson(x,[1])",0,histoRanges[i][1]);
    //gauss->SetParNames("Gauss_Constant","Gauss_Mean","Gauss_Sigma");
    poisson->SetParNames("Poisson_Scale","Poisson_Mean");

    gauss->SetParameters(histoGauss->GetMaximum(),histoGauss->GetMean(), histoGauss->GetRMS() );
    histoGauss->Fit("gauss","S",0,fitRanges[i][1]);
    poisson->SetParameters(histoPoisson->GetMaximum(),histoPoisson->GetMean(),0,fitRanges[i][1]);
    poisson->SetLineColor(kBlue);
    histoGauss->Fit("poisson","S+");

    float kaigauss = gauss->GetChisquare();
    float ndfgauss = gauss->GetNDF();
    float kaipoisson = poisson->GetChisquare();
    float ndfpoisson = poisson->GetNDF();
    kai2Poisson[i] = kaipoisson/ndfpoisson;
    kai2Gauss[i] = kaigauss/ndfgauss;

    cout<< "kai2/ndf value for gaussian fit for"<< titles[i].c_str()<<"is "<<kai2Gauss[i]<<endl;
    cout<< "kai2/ndf value for poisson fit for"<< titles[i].c_str()<<"is "<<kai2Poisson[i]<<endl;


    /*poisson->SetParameters(histoPoisson->GetMaximum(),histoPoisson->GetMean());
    histoPoisson->Fit("poisson", "R");
*/
    TCanvas *c1 = new TCanvas();
    histoGauss->Draw();
    /*TCanvas *c2 = new TCanvas();
    histoPoisson->Draw();
  */
    }

    /*for(int j = 0; j < 4; j++) {
        cout<< "means "<<means[j]<<"sigmas "<<sigmas[j]<< endl;
    } 
*/
    float x[4], y[4];
    for (int j = 0; j < 4; j++) {
        x[j] = means[j];
        float yValue = sqrt(means[j]) / sigmas[j];
        y[j]= yValue;
    }
    TCanvas *c3 = new TCanvas();
    TGraphErrors secondGraph(4, x, y);
    secondGraph.SetTitle("Sqrt(mean) / Standart Deviation as a function of mean");
    secondGraph.GetXaxis()->SetTitle("mean");
    secondGraph.GetYaxis()->SetTitle("sqrt(mean) / Standart Deviation");
    secondGraph.SetMinimum(0);   
    secondGraph.SetMaximum(1.2);

    TF1 *fitFunction = new TF1("fitFunction", "[0] + [1]*x",0,120); 
    secondGraph.Fit("fitFunction", "R");
    secondGraph.Draw("A*");
      
    fitFunction->Draw("same");  

    // PART 2

    TF1 *n0 = new TF1("n0", "[1]*[0]*TMath::Exp(-[0]*x)");
    TF1 *n1 = new TF1("n1", "[1]*[0]*[0]*x*TMath::Exp(-[0]*x)",2,60);
    n0->SetParName(0, "alpha");
    n1->SetParName(0, "alpha");

    TCanvas *c4 = new TCanvas();

    

    TH1F *lastHisto = new TH1F("lastHisto","n=0",20, 0.0, 43.0);
    lastHisto->GetXaxis()->SetTitle("Time Difference Between Successive Pulses(s)");
    lastHisto->GetYaxis()->SetTitle("number of occurances of Time Differences");
    lastHisto->SetTitle("Histogram for n=0");

    

    std::vector<std::string> filenames2 = {
        "occurance.txt",
     
     };
    std::ifstream infile2(filenames2[0]);
    float value2;
    int number_of_events;
    float total_time;
    while (infile2 >> value2) {
          lastHisto->Fill(value2);
          number_of_events += 1;
          total_time += value2;
    }

    float theoreticalAlpha = number_of_events/total_time;

    n0->SetParameters(lastHisto->GetMean(),lastHisto->GetMaximum());
    n0->SetRange(0.1 ,40);
    lastHisto->Fit("n0");
    lastHisto->Draw();

    TCanvas *c5 = new TCanvas();
    TH1F *lastHisto2 = new TH1F("lastHisto","n=0",20, 0.0, 60.0);
    lastHisto2->GetXaxis()->SetTitle("Time Difference Between Successive Pulses(s)");
    lastHisto2->GetYaxis()->SetTitle("number of occurances of Time Differences");
    lastHisto2->SetTitle("Histogram for n=1");


    std::vector<std::string> filenames3 = {
        "occurances2.txt",
     
     };

    std::ifstream infile3(filenames3[0]);
    float value3;
    while (infile3 >> value3) {
          lastHisto2->Fill(value3);
    }

    n1->SetParameters(0.1,120);
    lastHisto2->Fit("n1","R");
    lastHisto2->Draw();
    gStyle->SetOptFit(1111);

    cout << "theoretical value of alpha calculated is "<<theoreticalAlpha<<endl;




 
  
}


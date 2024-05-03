{
   
    std::vector<std::string> filenames = {
        "fixed.txt",
        "u1down.txt",
        "u1up.txt",
        "u3down.txt",
        "u3up.txt",
        "tdown.txt",
        "tup.txt",
        "random1.txt",
        "random2.txt",
        "random3.txt",


     };

     std::vector<std::string> titles = {
        "Fixed",
        "U1 Down",
        "U1 Up",
        "U3 Down",
        "U3 Up",
        "T Down",
        "T Up",
        "Random 1",
        "Random 2",
        "Random 3"
    };

    std::vector<std::vector<double>> peakRanges = {
        {13.9,15.1}, {19.2, 20.3},{24.1, 25.1},{29.1, 30.4},
        {14.1,15.3}, {18.8, 20.0},{24.0, 25.0},{29.2, 30.5},
        {14.3,15.2}, {18.8, 20.18},{24.03, 25.1},{29.1, 30.4},
        {14.2,15.13}, {19.1, 20.28},{24.1, 25.1},{29.13, 30.43},
        {14.05,15.11}, {19.12, 20.3},{24.1, 25.16},{29.13, 30.4},
        {13.98,15.16}, {18.9, 20.21},{24.14, 25.1},{29.12, 30.43},
        {13.93,15.12}, {18.82, 20.24},{24.02, 25.08},{28.95, 30.32},
        {13.9,15.02}, {18.9, 20.2},{24.08, 25.08},{29.04, 30.38},
        {13.88,15.13}, {19.16, 20.28},{24.1, 25.1},{29.1, 30.4},
        {13.93,15.16}, {19.09, 20.19},{24.04, 25.09},{29.14, 30.28},

   
          
          
    };

    std::vector<double> peakValues; 
    std::vector<double> peakDifferences; 
    std::vector<double> sigmaValues;
    std::vector<double> peakSigmas;
    std::vector<float> weightedPeaks;
    std::vector<float> weightedSigmas;
    



     for (int i = 0; i < filenames.size()+1; ++i) {
         TGraphErrors *gr = new TGraphErrors(filenames[i].c_str());
         TCanvas *c = new TCanvas();
         gr->Draw("AP");
         gr->SetTitle((titles[i] + "; U2 ; I").c_str());

         for (int j = 0; j < 4; ++j) {
            TF1 *gauss = new TF1(Form("gauss_%d_%d", i, j), "gaus", peakRanges[j][0], peakRanges[j][1]);
            gr->Fit(gauss, "R");
            gauss->Draw("same");
            double peak = gauss->GetParameter(1);
            double sigma = gauss->GetParameter(2);
            cout << "graph "<< i <<  "'s  " << j << " peak is ="<< peak <<"and its sigma is="<< sigma << "\n";
            peakValues.push_back(peak);            
            sigmaValues.push_back(sigma);

        }

         c->Draw();


}
    for (int n = 0; n < 40; ++n) { // calculating the difference between each peak value for each dataset 
        
        float difference = peakValues[n+1] - peakValues[n];
        if (difference>0){
            peakDifferences.push_back(difference);
        }
    }



    for (int k = 0; k < 40; ++k) { // getting sigma values for each peak
        if (k != 39 || k % 4 != 0){
            float first = sigmaValues[k]*sigmaValues[k];
            float second = sigmaValues[k+1]*sigmaValues[k+1];
            float final = sqrt(first + second);
            peakSigmas.push_back(final);

        }

}
    for (int q=0;q<10;q++){//calculating weighted averages
        int index = 3*q;

        float weightedupper = 0; 
        float weightedlower = 0; 

        for(int p=index;p<(index+3);p++){
            weightedupper += peakDifferences[p]/(peakSigmas[p]*peakSigmas[p]);
            weightedlower += 1/(peakSigmas[p]*peakSigmas[p]);

            //weighted average for peak differences
            float weightedPeak = weightedupper / weightedlower;
            weightedPeaks.push_back(weightedPeak);

             
            //weighted average for sigma values
            float weightedSigma = sqrt(1 / weightedlower);
            weightedSigmas.push_back(weightedSigma);

        }
    }
  
    
    for (int w=0;w<10;w++){
        cout<< "weighted peaks are" << weightedPeaks[w]<<"\n";
        cout<< "weighted sigma are" << weightedSigmas[w]<<"\n";

    }

    
    TH1F *histo = new TH1F("histo","Histogram",250, 4, 6);
    for (int e=0; e<30; e++) histo->Fill(peakDifferences[e]);
    histo->Draw();

    TF1 *gaussian = new TF1("gaussian", "gaus", 4, 6);

    histo->Fit("gaussian", "R");
    gaussian->SetLineColor(kRed);
    gaussian->Draw("same");


    cout << "Fit parameters:" << endl;
    cout << "Final result: " << gaussian->GetParameter(1) << endl;
    cout << "Standard deviation: " << gaussian->GetParameter(2) << endl;

}

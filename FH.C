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
        {13.2,15.8}, {17.8, 20.8},{23, 26},{28.7, 31.3},
        {12.8,15.8}, {18, 21},{22.3, 25.8},{28.3, 31.7},
        {12.8,15.8}, {17.9, 20.9},{22.5, 25.5},{27.3, 30.6},   
        {13.8,16.8}, {18.2, 21.2},{22.3, 25.3},{28, 31.8},   
        {12.9,15.9}, {18, 21},{23, 26},{28.8, 31.4},   
        {13.3,16.3}, {18.2, 21.2},{22.7, 25.7},{28.5, 31.4},   
        {13.4,16.4}, {17.8, 20.8},{22.9, 25.9},{28.2, 31.5},   
        {13.2,16.2}, {18.1, 21.1},{23.3, 26.3},{228.7, 31.6},   
        {13.1,16.1}, {18.2, 21.2},{22.7, 25.7},{28.4, 31.5},   
        {13.1,16.2}, {18.1, 21.2},{23, 26},{28.5, 31.5},   

   
          
          
    };

    std::vector<double> peakValues; 
    std::vector<double> sigmaValues;
    std::vector<double> peakSigmas;

    



     for (int i = 0; i < filenames.size(); ++i) {
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

    float peakDifferences[30] = {4.9776, 5.0419, 5.1249, 5.0140, 5.0480, 5.1246, 4.9793, 5.0770, 5.1597,
    4.9942, 5.0686, 5.1626, 4.9995, 5.0672, 5.1429, 5.0213, 5.0847, 5.1874,
    4.9057, 5.0261, 5.0073, 4.9834, 5.0416, 5.0812, 4.9873, 5.0358, 5.1427,
    4.9623, 5.0404, 5.1169}; 

}
    for (int m = 0; m < 40; ++m) {
        cout << "number" << m << "peak value is"<< peakValues[m] << "\n";

}
    for (int k = 0; k < 40; ++k) {
        if (k != 39 || k % 4 != 0){
            float first = sigmaValues[k]*sigmaValues[k];
            float second = sigmaValues[k+1]*sigmaValues[k+1];
            float result = sqrt(first + second);
            peakSigmas.push_back(result);

        }
        cout << "number" << k << "peak sigma value is"<< peakSigmas[k] << "\n";


}

    
}

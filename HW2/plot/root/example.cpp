/// @file example.cpp
/// @brief example root script to read ascii file and 
/// make a plot
/// @author pistack (Junho Lee)
/// @date 2021. 11. 17.

{
TFile *file = new TFile("example.root", "recreate");
TTree *tree = new TTree("example", "tree title");
// canvas for zeta
TCanvas *plot_zeta = new TCanvas("#zeta", "Zeta", 
800, 600);

Double_t t; // time
Double_t zeta; // zeta

tree -> Branch("t", &t);
tree -> Branch("#zeta", &zeta);

std::ifstream infile("../txtfile/zeta100.txt");
std::string line;
std::getline(infile, line); // discard first line
while(std::getline(infile, line))
{
    std::istringstream iss(line);
    iss >> t >> zeta;
    tree -> Fill();
}

tree -> Write(); // write tree to file
plot_zeta -> cd();
std::size_t n = tree -> Draw("t:#zeta", "", "gOff"); // zeta
TGraph *graph_zeta = new TGraph(n,tree->GetV1(),tree->GetV2());
graph_zeta -> SetMarkerStyle(21);
graph_zeta -> SetMarkerSize(0.5);
graph_zeta -> SetTitle("#zeta; t; #zeta");
(graph_zeta -> GetXaxis()) -> Set(n, 0, 3.25);
graph_zeta -> Draw("alp");
file -> Close(); // close file 
}


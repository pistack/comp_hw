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
// canvas for theta
TCanvas *plot_theta = new TCanvas("#theta", "Theta", 
800, 600);
// canvas for trajectory
TCanvas *plot_traj = new TCanvas("Trajectory", "Trajectory",
800, 600);

Double_t t; // time
Double_t zeta; // zeta
Double_t theta; // theta
Double_t x; // x coordinate
Double_t y; // y coordinate

tree -> Branch("t", &t);
tree -> Branch("#zeta", &zeta);
tree -> Branch("#theta", &theta);
tree -> Branch("x", &x);
tree -> Branch("y", &y);

std::ifstream infile("../txtfile/zeta100.txt");
std::string line;
std::getline(infile, line); // discard first line
while(std::getline(infile, line))
{
    std::istringstream iss(line);
    iss >> t >> zeta >> theta;
    x = zeta*std::cos(theta);
    y = zeta*std::sin(theta);
    tree -> Fill();
}

tree -> Write(); // write tree to file
plot_zeta -> cd();
std::size_t n = tree -> Draw("t:#zeta", "", "gOff"); // zeta
TGraph *graph_zeta = new TGraph(n,tree->GetV1(),tree->GetV2());
graph_zeta -> SetMarkerStyle(21);
graph_zeta -> SetMarkerSize(0.5);
graph_zeta -> SetTitle("#zeta; t; #zeta");
(graph_zeta -> GetXaxis()) -> Set(n, 0, 10);
graph_zeta -> Draw("alp");
plot_theta -> cd();
tree -> Draw("#theta:t", "", "gOff"); // theta
TGraph *graph_theta = new TGraph(n,tree->GetV1(),tree->GetV2());
graph_theta -> SetMarkerStyle(21);
graph_theta -> SetMarkerSize(0.5);
graph_theta -> SetTitle("#theta; t; #theta");
(graph_theta -> GetXaxis()) -> Set(n, 0, 10);
graph_theta -> Draw("alp");
plot_traj -> cd();
tree -> Draw("y:x", "", "gOff"); // trajectory
TGraph *graph_traj = new TGraph(n,tree->GetV1(),tree->GetV2());
graph_traj -> SetMarkerStyle(21);
graph_traj -> SetMarkerSize(0.5);
graph_traj -> SetTitle("Trajectory; x; y");
(graph_traj -> GetXaxis()) -> Set(n, -1.125, 1.125);
graph_traj -> Draw("alp");
file -> Close(); // close file 
}


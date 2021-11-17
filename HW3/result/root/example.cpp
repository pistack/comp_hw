/// @file example.cpp
/// @brief example root scripts to read ascii file and 
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

std::ifstream infile("../txtfile/zeta4_2.txt");
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
tree -> Draw("#zeta:t", ""); // zeta
plot_theta -> cd();
tree -> Draw("#theta:t", ""); // theta
plot_traj -> cd();
tree -> Draw("y:x", ""); // trajectory
file -> Close(); // close file 
}


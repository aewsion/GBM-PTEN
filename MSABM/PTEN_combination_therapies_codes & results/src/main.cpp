#include <iostream>
#include <fstream>
#include <random>

#include "environment.h"

extern std::random_device rd;

int main(int argc, char **argv){
    std::string folder = argv[1]; // save folder
    std::string set = argv[2]; // set
    int N = std::stoi(argv[3]); // number of simulations
    double simTime = std::stod(argv[4]); // simulation time
    int if_responder = std::stoi(argv[5]); // rand = 1 for responders 
                                           // rand = 0 for non-responders
    int pTEN = std::stoi(argv[6]);
    double dose_C = std::stod(argv[7]);
    double dose_I = std::stod(argv[8]);
    double dose_A = std::stod(argv[9]);
    
    double M0RecProb = 0.0000035; // basal recruitment probability of M0 macrophages
    double TCProl = 14.0; // tumor cell mature and division time
    double K1 = 0.02; // Michaelis constant for the Hill function of CSF1
    double K_CSF1RI = 5.0; // Michaelis constant for the Hill function of CSF1R_I
    double K_A = 5.7; // Michaelis constant for the Hill function of CSF1R_I accumulation
    double p_TC = 0.0016; // tumor cell basal proliferation probability
    if(pTEN == 0){
        p_TC = 0.0016;
    }
    double p_D = 0.0001; // tumor cell basal death probability
    int pha = 2; // Maximum phagocytosis number of each M1 macrophage before having a rest

    // save output information to a txt file
    std::string str = "mkdir -p "+folder+"/set_" + set;
    const char *command = str.c_str();
    std::system(command);   
    std::streambuf* coutBuf = std::cout.rdbuf();
    std::ofstream of(folder+"/set_"+set+"/out"+set+".txt");
    std::streambuf* fileBuf = of.rdbuf();
    std::cout.rdbuf(fileBuf);

    // set random key parameters
    std::mt19937 g(rd());
    std::normal_distribution<double> Dis_n(0,1); // normal distribution
    double p1 = Dis_n(g);
    double p2 = Dis_n(g);
    double p3 = Dis_n(g); 
    double b_C = 0.675; // adjustment coefficient of Hill function H_C
    double b_CI = 0.45; // adjustment coefficient of Hill function H_I
    double a_ERK = 7.15; // coefficient associated with ERK to promote tumor growth
    double b_G = 1.75/12.69 * 0.5;

    if(if_responder == 0){
        K_CSF1RI = K_CSF1RI*(0.9);
        b_CI = b_CI*(0.6);
        a_ERK = a_ERK*(0.9);
    }

    if(K_CSF1RI < 0){
        K_CSF1RI = 0;
    }
    if(b_CI < 0){
        b_CI = 0;
    }
    if(a_ERK < 0){
        a_ERK = 0;
    }

    // record running time
    double start = clock();

    for (int i = 0; i < N; i++) {
        std::cout << "************************************" << std::endl;
        std::cout << "simulation: "<< i << " start" << std::endl;
        std::cout << "b_C: " << b_C << std::endl;
        std::cout << "b_CI: " << b_CI << std::endl;
        std::cout << "a_ERK: " << a_ERK << std::endl;
        std::cout << "pTEN: " << pTEN << std::endl;
        std::cout << "dose_CSF1R_I: " << dose_C << std::endl;
        std::cout << "dose_IGF1R_I: " << dose_I << std::endl;
        std::cout << "dose_AKT_I: " << dose_A << std::endl;
        // create virtual TME, step size = 0.5 hour
        Environment model(0.5, folder, std::stoi(set), M0RecProb, TCProl,
                     K1, K_CSF1RI, K_A, p_TC, p_D, pha, rand,
                     b_C, b_CI, a_ERK, b_G, pTEN, dose_C, dose_I, dose_A);

        // start simulation             
        model.simulate(simTime);

        std::cout << "simulation: "<< i << " end" <<std::endl;
        std::cout << "************************************" << std::endl;
    }

    double stop = clock();

    // running time
    std::cout << "Duration: " << (stop-start)/CLOCKS_PER_SEC << std::endl;    
    of.flush();
    of.close();
    std::cout.rdbuf(coutBuf);
    
    std::cout << "Duration: "+folder+"/set_"+set <<" "<< (stop-start)/CLOCKS_PER_SEC << std::endl;

    return 0;
}

#ifndef MACROPHAGE_H
#define MACROPHAGE_H

#include <array>
#include <vector>

#include "cellgrids.h"
#include "diffusibles.h"

class Macrophage{
public:
    Macrophage(std::array<int, 2> loc, int initial, std::string State, double k1, 
                double k_CSF1RI, double k_A, int pha, double b_ci, double b_c, double b_g); // create TAM agent
    void migration(CellGrids &cg, Diffusibles &diff); // simualte migration of TAMs
    void differentiaion(CellGrids &cg); // simualte differentiaion of M0
    int transformation(CellGrids &cg, Diffusibles &diff); // simualte transformation of M1 & M2
    int phagocytosis(CellGrids &cg, Diffusibles &diff); // simualte phagocytosis of M1
    void simulaion(double tstep, CellGrids &cg, Diffusibles &diff, double depletion); // simulate TAM agent behaviors

    std::string state; //M0, M1, M2, dead
    std::array<int, 2> location; // loaction in the cell grid

    double life_span; // lifespan of the TAMs
    double age; // age

    int CSF1R = 1; // 1 means exsiting active CSF1R on TAMs
    int max_phagocytose; // maximum phagocytosis number of M1 macrophages before having a rest
    int now_phagocytose; // current phagocytosis number of the M1 macrophage

    double b_C; // adjustment coefficient of Hill function H_C
    double b_G;
    double b_CI; // adjustment coefficient of Hill function H_I

    int trans; // 1 means TAMs have transformation potential at the current location
               // 0 means TAMs cannot transform at the current location
    int trans_num; // control transformation & differetiation interval
    double rest_time; // rest time after reaching maximum phagocytosis number

private:
    double D_M; // Diffusion coefficient of M0, M1 and M2 macrophages
    double a_MC; // Chemotaxis coefficient of CSF1 for TAMs
    double a_MG; // Chemotaxis coefficient of Gal9_en for TAMs

    double p_M01; // Probability of M0 macrophage differentiate into M1
    double p_M02; // Probability of M0 macrophage differentiate into M2

    double K1; // Michaelis constant for the Hill function of CSF1
    double K2; // Michaelis constant for the Hill function of Gal9_en
    double K_CSF1RI; // Michaelis constant for the Hill function of CSF1R_I
    double K_Gal9I; // Michaelis constant for the Hill function of Gal9_I
    double K_A; // Michaelis constant for the Hill function of CSF1R_I accumulation

    int Pha; // phagocytosis rate
    double p_M21; // Probability of M2 macrophage polarize into M1
    double depletion; // Dead probability of TAMs while using drugs
};

#endif // MACROPHAGE_H

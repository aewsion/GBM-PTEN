#ifndef TUMOR_H
#define TUMOR_H

#include "cellgrids.h"
#include "macrophage.h"
#include "diffusibles.h"

class Tumor{
public:
    Tumor(std::array<int, 2> loc, double prolTime, int index, int initial, double p_tc, double p_d,
                double EGFRI, double IGF1RI, double AKTI, double ERKI, double GSK3bI, 
                double IRF1I, double Gal9I, double A_ERK, int TC_max, int pTEN); // Tumor cell agent
    void proliferation(CellGrids &cg, std::vector<Tumor> &tc_list, Diffusibles &diff, double prolTime); // simulate prolifetation of tumor cell
    void pro_prob(Diffusibles &diff, std::vector<Tumor> &tc_list); // calculate proliferation probability of tumor cell
    void migration(Diffusibles &diff, CellGrids &cg); // simulate migration of tumor cell
    void ODE_solution(Diffusibles &diff, double tstep); // odes representing signal pathways in tumor cell
    void simulation(double tstep, CellGrids &cg, std::vector<Tumor> &tc_list, Diffusibles &diff, double prolTime); // simulate tumor cells' hebaviors

    std::string state; // tumor cell state: active or quiescent
    int idx; // index
    std::array<int,2> location; // (i,j) location

    double age; // age
    double mature; // mature time
    double div_time; // time after the latest division
    int maxDiv; // maximum number of division
    int nDiv; // current number of division

    int pten; // pten = 0 for PTEN-null
              // pten = 1 for PTEN_wild

    double cell_cycle; // cell cycle
    double life_span; // lifespan
    
    int if_migra; // used to control tumor cell migration interval
    double prob; // tumor cell proliferate probability restricted by maximum capacity
    double prob_0; // proliferation probability without maximum capacity restriction

    // protein concentrations in tumor cell signaling pathways 
    double EGFR; 
    double IGF1R; 
    double AKT;
    double ERK; 
    double GSK3b;
    double IRF1;
    double Gal9;

    double a_ERK; // coefficient associated with ERK to promote tumor growth rate
    double a_AKT; // coefficient associated with AKT to promote tumor growth rate 

private:
    double tc_max; // tumor cell maximum capacity number in the TME

    double p_TC; // tumor cell basal proliferation probability
    double p_D; // tumor cell basal death probability

    double a_TI; // chemotaxis coefficient of IGF1 for tumor cells
    double a_TE; // chemotaxis coefficient of EGF for tumor cells
    double D_TC; // diffusion coefficient of tumor cells

    // maximum protein concentrations in tumor cell signaling pathways 
    double EGFR_max;
    double ERK_max;
    double AKT_max;
    double IGF1R_max;
    double GSK3b_max;
    double IRF1_max;
    double Gal9_max; 

    double dt; // dt for solving odes
    double time; // record simulation time

    // Michaelis constant for Hill functions in tumor cell
    // proliferation probability calculating
    double K4, K5;

    // parameters in tumor cell signal pathways
    double V1, K11, K12, K13, d1;
    double V21, V22, V2, K21, K22, K23, d2; 
    double V3, K31, K32, K33, d3;
    double V4, K41, K42, K43, K44, d4;
    double V5, V51, K51, K52, d5;
    double V6, K61, d6;
    double V7, K71, d7;
    double n1, n2, n3, n4, n5, n6;
    double K_ERKI, K_GSK3bI, K_IRF1I; 
};

#endif // TUMOR_H

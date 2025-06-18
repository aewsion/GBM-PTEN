#ifndef DIFFUSIBLES_H
#define DIFFUSIBLES_H

#include "cellgrids.h"

class Diffusibles {
public:
    Diffusibles(double DX);
    void diffusion(CellGrids &cg, double tstep, double dose_C);

    // cytokines
    double CSF1[100][100];
    double EGF[100][100];
    double IGF1[100][100];
    double Gal9_en[100][100]; // the concentration of Gal-9 in the extracellular microenvironment
    double I[100][100]; // save integral parts in IGF1 concentration caculation

    // drugs on TAMs
    double CSF1R_I[100][100];
    double Gal9_I[100][100];
    // drugs on TCs
    double IGF1R_I[100][100];
    double EGFR_I[100][100];
    double AKT_I[100][100];
    double ERK_I[100][100];
    double GSK3b_I[100][100];
    double IRF1_I[100][100];

    // signal pathway 
    double AKT[100][100];
    double ERK[100][100];
    double IGF1R[100][100];
    double EGFR[100][100];
    double GSK3b[100][100];
    double IRF1[100][100];
    double Gal9[100][100]; // intracellular Gal-9 activity

    // tumor cell proliferation probability
    // only used for recording, no calculation here
    double prob[100][100];
    double prob_0[100][100];

    double dx; // dx = grid size = 15um
    double dt; // dt for calculating cytokines
    double dt_drug; // dt for calculating IGF1 and drugs
    double time; // record simulating time (hour)

    // only used for recording continuous drug treatment
    double drug_time; // accumulative CSF1R_I using time (second)
    double drug_feed_time; // drug feeding time (second)
    double drug_suspend_time; // drug suspending time (second)

    double S_A; // basal rate for calculating I
    double A_0; // initial rate for calculating IGF1

    // maximum cytokines carrying capacity in the TME
    double max_CSF1;
    double max_IGF1;
    double max_EGF; 
    double max_Gal9_en;

    int if_CSF1R_I, if_IGF1R_I, if_EGFR_I; // if feed drug

private:
    // diffusion coefficient of cytokines
    double D_IGF1;
    double D_CSF1;
    double D_EGF;
    double D_Gal9_en;

    double D_CSF1R_I; // diffusion coefficient of CSF1R_I
    double D_IGF1R_I; // diffusion coefficient of IGF1R_I
    double D_EGFR_I; // diffusion coefficient of EGFR_I 
    double D_AKT_I;
    double D_ERK_I;
    double D_GSK3b_I;
    double D_IRF1_I;
    double D_Gal9_I;

    double q_CSF1R_I; // vascular permeability of CSF1R_I
    double q_IGF1R_I; // vascular permeability of IGF1R_I
    double q_EGFR_I; // vascular permeability of IGF1R_I 
    double q_AKT_I;
    double q_ERK_I;
    double q_GSK3b_I;
    double q_IRF1_I;
    double q_Gal9_I;

    double S_CSF1; // CSF1 secretion rate of tumor cell
    double S_IGF1; // IGF1 secretion rate of M2
    double S_EGF; // EGF secretion rate of M2
    double S_Gal9_en;

    // degradation rate of cytokines
    double d_CSF1;
    double d_IGF1;
    double d_EGF;
    double d_Gal9_en;

    // natural decay rate of drugs
    double d_CSF1R_I;
    double d_IGF1R_I;
    double d_EGFR_I; // if necessary
    double d_AKT_I;
    double d_ERK_I;
    double d_GSK3b_I;
    double d_IRF1_I;
    double d_Gal9_I;

    double u_CSF1R_I; // TAMs uptake rate of CSF1R_I 
    double u_IGF1R_I; // tumor cells uptake rate of IGF1R_I 
    double u_EGFR_I; // tumor cells uptake rate of EGFR_I
    double u_AKT_I;
    double u_ERK_I;
    double u_GSK3b_I;
    double u_IRF1_I;
    double u_Gal9_I;

    double K3; // Michaelis constant for the Hill function of Gal-9 in the signaling pathways
};

#endif // DIFFUSIBLES_H

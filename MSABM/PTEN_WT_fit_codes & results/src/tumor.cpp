#include <iostream>
#include <random>

#include "tumor.h"

extern std::random_device rd;

Tumor::Tumor(std::array<int, 2> loc, double prolTime, int index, int initial, double p_tc, double p_d,
                double EGFRI, double IGF1RI, double AKTI, double ERKI, double GSK3bI, 
                double IRF1I, double Gal9I, double A_ERK, int TC_max, int pTEN){
    // representation, unit, references of parameters
    // are listed in Supplementary Table 'Tumor cells' part 
    // and 'Signaling pathways' part 
    state = "alive";
    idx = index;
    location = loc;

    age = 0;
    mature = prolTime;
    div_time = 0;
    maxDiv = 8;
    
    // initial = 1 represents intializeCells before 
    // tumor growth simulation, not the initialization
    // of virtual TME
    std::mt19937 g(rd()); 
    std::uniform_real_distribution<>  ageStart(0.0,24*0.1);
    if(initial == 1){ 
        age = ageStart(g);
        mature = prolTime * 3;  
        div_time = age;
        maxDiv = 12;
    }

    nDiv = 0;

    cell_cycle = mature;
    life_span = 24.0 * 8;

    if_migra = 0;
    prob = 0.0;
    prob_0 = 0.0;

    tc_max = TC_max;
    
    p_D = p_d;
    p_TC = p_tc;

    a_TI = 2600; 
    a_TE = 2600;
    D_TC = 1e-9; 

    dt = 0.00069; // 1 min = 1/1440 day
    time = 0;

    K4 = 0.9502;
    K5 = 0.95; // 0.015

    a_AKT = 44;
    a_ERK = A_ERK;

    EGFR_max = 1;
    ERK_max = 1;
    AKT_max = 1;
    IGF1R_max = 1;
    GSK3b_max = 1;
    IRF1_max = 1;
    Gal9_max = 1; 
    
    EGFR = EGFRI;
    IGF1R = IGF1RI;
    AKT = AKTI;
    ERK = ERKI;
    GSK3b = GSK3bI;
    IRF1 = IRF1I;
    Gal9 = Gal9I;
    
    pten = pTEN;
    
    V1 = 26.0784;
    V2 = 3.9216;
    V21 = 45.0980;
    V22 = 4.1177;
    V3 = 10.7843;
    V4 = 47.8431;
    V5 = 6.4706;
    V51 = 7.6471;
    V6 = 23.7255;
    V7 = 28.0392;

    d1 = 30.5882;
    d2 = 24.3137;
    d3 = 9.6079;
    d4 = 5.6863;
    d5 = 2.1569;
    d6 = 40.7843;
    d7 = 2.7451;

    K11 = 0.0588;
    K12 = 0.4353;
    K13 = 0.2353;
    K21 = 0.5412;
    K22 = 2.0588;
    K23 = 1.2941;
    K31 = 1.4471;
    K32 = 2.5059;
    K33 = 0.1961;
    K41 = 0.0706;
    K42 = 1e-5;
    K43 = 0.1177;
    K44 = 0.4314;
    K51 = 1.2941;
    K52 = 0.3039;
    K61 = 0.4824;
    K71 = 0.7882;

    n1 = 6; 
    n2 = 2;
    n3 = 1;
    n4 = 2;
    n5 = 4;
    n6 = 5;

    K_ERKI = 0.05; 
    K_GSK3bI = 0.05;
    K_IRF1I = 0.5;
}

// calculate proliferation probability of tumor cell
void Tumor::pro_prob(Diffusibles &diff, std::vector<Tumor> &tc_list){

    // Hill functions
    double H1 = ERK/(K4 + ERK);
    double H2 = AKT/(K5 + AKT);

    // current number of tumor cells
    double tc_now = tc_list.size();

    // calculate probability of proliferation
    if(state == "quiescent"){
        prob = 1.0;
    }
    else{
        prob = p_TC * (1 + a_ERK*H1 + a_AKT*H2)*(1 - tc_now*1.0/tc_max); 
    }
    
    // probability without restriction of maximum capacity number
    // used for testing
    prob_0 = p_TC * (1 + a_ERK*H1 + a_AKT*H2);
    
    // location
    int i = location[0];
    int j = location[1];
    
    // probability should <= 1.0
    if(prob >= 1.0){
        prob = 1.0;
    }

    // record probabilities
    diff.prob[i][j] = prob;
    diff.prob_0[i][j] = prob_0;   
}

// simulate prolifetation of tumor cell
void Tumor::proliferation(CellGrids &cg, std::vector<Tumor> &tc_list, Diffusibles &diff, double prolTime){
    
    // (i, j) loaction
    int i = location[0];
    int j = location[1];

    // eight latent proliferation directions
    // Moore neighborhood
    int z = 0;
    std::vector<int> ix;
    std::vector<int> jx;
    std::mt19937 g(rd());
    std::uniform_real_distribution<> Dis(0.0,1.0);
    for(int il=-1; il<2; il++){
        for(int jl=-1; jl<2; jl++){
            ix.push_back(il);
            jx.push_back(jl);
            z++;
        }
    }

    double probs[z];
    double sum=0;
    for(int q=0; q<z; q++){
        probs[q] = (1 - cg.allcells[i+ix[q]][j+jx[q]]);    
        sum += probs[q];
    }

    // simulate proliferation
    double gg = Dis(g);
    if(gg > prob){
        return;
    }
    else if((sum == 0 || tc_list.size() >= tc_max) && gg <= prob){
        state = "quiescent";
        return;
    }
    else if(sum > 0 && tc_list.size() < tc_max && gg <= prob){
        state = "alive";

        // normalization
        double norm_probs[z];
        for(int q=0;q<z;q++){norm_probs[q] = probs[q]/(sum);} 
        for(int q=1; q<z; q++){norm_probs[q] = norm_probs[q] + norm_probs[q-1];} 

        // choose proliferation direction
        double p = Dis(rd);
        int choice = 0;
        for(double norm_prob : norm_probs){
            if(p > norm_prob){choice++;}
        }

        // proliferation, update cell grid and cell list
        int ni = i + ix[choice];
        int nj = j + jx[choice];
    
        cg.allcells[ni][nj] = 1 ;
        cg.tc[ni][nj] = 1;
        if(cell_cycle < prolTime){
            tc_list.push_back(Tumor({ni,nj}, cell_cycle*2, tc_list.size(), 0, p_TC, p_D, 
                                        EGFR/2, IGF1R/2, AKT/2, ERK/2, GSK3b/2, IRF1/2, Gal9/2,
                                        a_ERK, tc_max, pten));
        }
        else{
            tc_list.push_back(Tumor({ni,nj}, prolTime, tc_list.size(), 0, p_TC, p_D, 
                                        EGFR/2, IGF1R/2, AKT/2, ERK/2, GSK3b/2, IRF1/2, Gal9/2,
                                        a_ERK, tc_max, pten));
        }

        // division of proteins
        diff.AKT[ni][nj] = AKT/2;
        diff.ERK[ni][nj] = ERK/2;
        diff.IGF1R[ni][nj] = IGF1R/2;
        diff.EGFR[ni][nj] = EGFR/2;
        diff.GSK3b[ni][nj] = GSK3b/2;
        diff.IRF1[ni][nj] = IRF1/2;
        diff.Gal9[ni][nj] = Gal9/2;
        EGFR = EGFR/2;
        IGF1R = IGF1R/2;
        AKT = AKT/2;
        ERK = ERK/2;
        GSK3b = GSK3b/2;
        IRF1 = IRF1/2;
        Gal9 = Gal9/2;
        nDiv++; // current number of divisions 

        // reset
        prob = 0;
        div_time = 0;
    }   
}

// simulate migration of tumor cell
void Tumor::migration(Diffusibles &diff, CellGrids &cg){

    // (i, j) location
    int i = location[0];
    int j = location[1];

    // five latent choices of migration
    std::array<int, 5> ix={0,1,0,-1,0};
    std::array<int, 5> jx={1,0,-1,0,0};

    // EGF and IGF1 concentration from diffusion equation
    double I[5]={0,0,0,0,0};
    double E[5]={0,0,0,0,0};
    for(int q=0; q<5; q++){
        E[q] = diff.EGF[i+ix[q]][j+jx[q]];
        I[q] = diff.IGF1[i+ix[q]][j+jx[q]];
    }
    std::mt19937 g(rd());
    std::uniform_real_distribution<> Dis(0.0,1.0);
    
    double h = diff.dx; // 15um = grid size
    double k = 2250*0.05; // choose a suitable k

    // calculate probability of latent migration choices
    double probs[5]={0,0,0,0,0};
    probs[4] = (1 - 4*k*D_TC/pow(h,2) - k*a_TI*(I[0]+I[1]+I[2]+I[3]-4*I[4])/pow(h,2) - k*a_TE*(E[0]+E[1]+E[2]+E[3]-4*E[4])/pow(h,2)); // remain stationary
    probs[0] = (1 - cg.allcells[i+ix[0]][j+jx[0]]) * (k*D_TC/pow(h,2) - k*a_TI*(I[2]-I[0])/(4*pow(h,2)) - k*a_TE*(E[2]-E[0])/(4*pow(h,2))); // move up
    probs[1] = (1 - cg.allcells[i+ix[1]][j+jx[1]]) * (k*D_TC/pow(h,2) - k*a_TI*(I[3]-I[1])/(4*pow(h,2)) - k*a_TE*(E[3]-E[1])/(4*pow(h,2))); // move right
    probs[2] = (1 - cg.allcells[i+ix[2]][j+jx[2]]) * (k*D_TC/pow(h,2) - k*a_TI*(I[0]-I[2])/(4*pow(h,2)) - k*a_TE*(E[0]-E[2])/(4*pow(h,2))); // move down
    probs[3] = (1 - cg.allcells[i+ix[3]][j+jx[3]]) * (k*D_TC/pow(h,2) - k*a_TI*(I[1]-I[3])/(4*pow(h,2)) - k*a_TE*(E[1]-E[3])/(4*pow(h,2))); // move left
    
    // identify migration direction
    double sum=0;
    for(int q=0; q<5; q++){
        sum += probs[q];
    }
    if(sum == 0){
        return;
    }
    double norm_probs[5];
    for(int q=0;q<5;q++){norm_probs[q] = probs[q]/(sum);} // normalization
    for(int q=1; q<5; q++){norm_probs[q] = norm_probs[q] + norm_probs[q-1];} 

    double p = Dis(rd); 
    int choice = 0;
    // choose migration direction
    for(double norm_prob : norm_probs){
        if(p > norm_prob){choice++;}
    }
    int ni = i + ix[choice];
    int nj = j + jx[choice];

    // migration and update cell grid
    location[0] = ni;
    location[1] = nj;
    cg.allcells[i][j]=0;
    cg.m0[i][j]=0;
    cg.m1[i][j]=0;
    cg.m2[i][j]=0;
    cg.tc[i][j]=0;
    cg.allcells[ni][nj]=1;
    cg.m1[ni][nj]=0;
    cg.m2[ni][nj]=0;
    cg.m0[ni][nj]=0;
    cg.tc[ni][nj]=1; 

    // update signal pathways protein
    diff.AKT[i][j] = 0;
    diff.ERK[i][j] = 0;
    diff.IGF1R[i][j] = 0;
    diff.EGFR[i][j] = 0; 
    diff.GSK3b[i][j] = 0;
    diff.IRF1[i][j] = 0;
    diff.Gal9[i][j] = 0;
    diff.AKT[ni][nj] = AKT;
    diff.ERK[ni][nj] = ERK;
    diff.IGF1R[ni][nj] = IGF1R;
    diff.EGFR[ni][nj] = EGFR; 
    diff.GSK3b[ni][nj] = GSK3b;
    diff.IRF1[ni][nj] = IRF1;
    diff.Gal9[ni][nj] = Gal9;
}

// odes representing signal pathways in tumor cell
void Tumor::ODE_solution(Diffusibles &diff, double tstep){
    time = time + tstep;
    
    double t = 30; // tstep = 0.5h = 30 * 1 min
    double maxDif = 0;
    double c = 0;

    // loaction
    int i = location[0];
    int j = location[1];

    // cytokines EGF and IGF1
    double E = diff.EGF[i][j] / diff.max_EGF; // / 5;
    double I = diff.IGF1[i][j] / diff.max_IGF1; // / 3000;

    // drugs from diffusion equations
    double EGFR_I = diff.EGFR_I[i][j];
    double IGF1R_I = diff.IGF1R_I[i][j];
    double AKT_I = diff.AKT_I[i][j];
    double ERK_I = diff.ERK_I[i][j];
    double GSK3b_I = diff.GSK3b_I[i][j];
    double IRF1_I = diff.IRF1_I[i][j];

    double GSK_I = 0; // only used in signal pathway parameters estimation
    
    // forward Eular method
    for(int q=0; q<t; q++){
        maxDif = 0;
        double EGF_R_0 = EGFR;
        double ERK_0 = ERK;
        double IGF1_R_0 = IGF1R;
        double AKT_0 = AKT;
        double GSK3b_0 = GSK3b;
        double IRF1_0 = IRF1;
        double Gal9_0 = Gal9;
        EGFR = EGFR + dt * ( V1 * E/(K11+E) * 1/(1 + ERK/K12) * 1/(1 + EGFR_I/K13) * (EGFR_max - EGFR) - d1 * EGFR);

        ERK = ERK + dt * ( V2 * (1 + V21*pow(EGFR, n1))/(pow(K21, n1)+pow(EGFR, n1)) 
                            * (1 + V22*IGF1R/(K22+IGF1R)) * 1 / (1 + AKT/K23) * (1/(1 + ERK_I/K_ERKI)) * (ERK_max -ERK) - d2 * ERK );

        IGF1R = IGF1R + dt * ( V3 * I/(K31 + I) * 1/(1 + ERK/K32) * 1/(1 + IGF1R_I/K33) * (IGF1R_max - IGF1R) - d3 * IGF1R );

        AKT = AKT + dt * ( V4/(1 + pten/K41) * EGFR/(K42 + EGFR) * IGF1R/(K43 + IGF1R)
                             * (1/(1 + pow(AKT_I, n2)/pow(K44, n2))) * (AKT_max - AKT) - d4 * AKT );

        GSK3b = GSK3b + dt * ( V5 * pow(AKT, n3)/(pow(K51, n3) + pow(AKT, n3)) * (1 + V51*pow(GSK_I, n4)/(pow(K52, n4) + pow(GSK_I, n4)))
                                * 1/(1 + GSK3b_I/K_GSK3bI) * (GSK3b_max - GSK3b) - d5 * GSK3b );

        IRF1 = IRF1 + dt * ( V6/(1 + IRF1_I/K_IRF1I) * (IRF1_max - IRF1) - d6 * IRF1 * (1/(1 + pow(GSK3b, n5)/pow(K61, n5))) );

        Gal9 = Gal9 + dt * ( V7 * pow(IRF1, n6)/(pow(K71, n6) + pow(IRF1, n6)) * (Gal9_max - Gal9) - d7 * Gal9 );

        if(EGF_R_0 > 0){
            c = (EGFR - EGF_R_0)/EGF_R_0;
            if(c > maxDif){maxDif = c;}
        }
        if(ERK_0 > 0){
            c = (ERK - ERK_0)/ERK_0;
            if(c > maxDif){maxDif = c;}
        }
        if(IGF1_R_0 > 0){
            c = (IGF1R - IGF1_R_0)/IGF1_R_0;
            if(c > maxDif){maxDif = c;}
        }
        if(AKT > 0){
            c = (AKT - AKT_0)/AKT_0;
            if(c > maxDif){maxDif = c;}
        }
        if(GSK3b > 0){
            c = (GSK3b - GSK3b_0)/GSK3b_0;
            if(c > maxDif){maxDif = c;}
        }
        if(IRF1 > 0){
            c = (IRF1 - IRF1_0)/IRF1_0;
            if(c > maxDif){maxDif = c;}
        }
        if(Gal9 > 0){
            c = (Gal9 - Gal9_0)/Gal9_0;
            if(c > maxDif){maxDif = c;}
        }
        // break condition
        if(maxDif < 0.00001 && q > 5){break;}
    }

    // record protein concentration
    diff.AKT[i][j] = AKT;
    diff.ERK[i][j] = ERK;
    diff.IGF1R[i][j] = IGF1R;
    diff.EGFR[i][j] = EGFR; 
    diff.GSK3b[i][j] = GSK3b;
    diff.IRF1[i][j] = IRF1;
    diff.Gal9[i][j] = Gal9;
}

// simulate tumor cells' hebaviors
void Tumor::simulation( double tstep, CellGrids &cg, std::vector<Tumor> &tc_list, Diffusibles &diff, double prolTime){ 
    if(state == "dead"){
        return;
    }

    if(state == "alive"){
        age = age + tstep;
        div_time = div_time + tstep;
        if_migra++;
    }
    
    std::mt19937 g(rd());
    std::uniform_real_distribution<> Dis(0.0, 1.0);
    // unnatural death according to a fixed probability or natural death   
    if((Dis(g) < p_D && state == "alive" && age >= mature) || age >= life_span){
        state = "dead";
        return;
    }
    
    // migration
    if(if_migra == 48){
        migration(diff, cg);
        if_migra = 0;
    }
    
    // signal pathway   
    ODE_solution(diff, tstep);

    // see if the tumor cell can proliferate
    if(age >= mature && div_time >= cell_cycle && nDiv < maxDiv){
        pro_prob(diff, tc_list);
        proliferation(cg, tc_list, diff, prolTime);
        if(nDiv > 7){
            std::cout<<"nDiv = "<<nDiv<<std::endl;
        }       
    }   
}


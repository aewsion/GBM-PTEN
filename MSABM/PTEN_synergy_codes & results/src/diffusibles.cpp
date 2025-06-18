#include <iostream>

#include "diffusibles.h"

Diffusibles::Diffusibles(double DX) {
    // representation, unit, references of parameters
    // are listed in Supplementary Table Cytokines part & Drugs part
    dx = DX;
    dt = 2;
    dt_drug = 2;
    time = 0;
    
    D_CSF1 = 3e-7; 
    D_IGF1 = 3e-7;
    D_EGF = 3e-7;
    D_Gal9_en = 3e-7;
    
    D_CSF1R_I = 2e-10;
    D_IGF1R_I = 2e-10; 
    D_EGFR_I = 2e-10;
    D_AKT_I = 2e-10;
    D_ERK_I = 2e-10;
    D_GSK3b_I = 2e-10;
    D_IRF1_I = 2e-10;
    D_Gal9_I = 2e-10;
    
    q_CSF1R_I = 0.8; 
    q_IGF1R_I = 0.8;
    q_EGFR_I = 0.8;
    q_AKT_I = 0.8;
    q_ERK_I = 0.8;
    q_GSK3b_I = 0.8;
    q_IRF1_I = 0.8;
    q_Gal9_I = 0.8;
    
    S_CSF1 = 2.5e-14; // 2.5e-26 mol per second per cell = 2.5e-14 mol/L/s
    S_Gal9_en = 2.5e-14;
    S_IGF1 = 5e-15; // 2.5e-26 mol per second per cell; M2 cell diameter = 20um
    S_EGF = 5e-15; // same as S_IGF1
    
    d_CSF1 = 5e-4; 
    d_Gal9_en = 5e-4;
    d_IGF1 = 2.5e-4; 
    d_EGF = 2.5e-4; 
    
    u_CSF1R_I = 2e-7;
    u_IGF1R_I = 2e-7;
    u_EGFR_I = 2e-7;
    u_AKT_I = 2e-7;
    u_ERK_I = 2e-7;
    u_GSK3b_I = 2e-7;
    u_IRF1_I = 2e-7;
    u_Gal9_I = 2e-7;
    
    d_CSF1R_I = 2e-7; 
    d_IGF1R_I = 2e-7;
    d_EGFR_I = 2e-7;
    d_AKT_I = 2e-7;
    d_ERK_I = 2e-7;
    d_GSK3b_I = 2e-7;
    d_IRF1_I = 2e-7;
    d_Gal9_I = 2e-7;
    
    drug_time = 0;
    drug_feed_time = (0)*24*60*60; // continuous_CSF1R_I_treatment_case here
    // set drug_feed_time = 0-200 days for other continuous treatment
    drug_suspend_time = (205)*24*60*60;
    
    S_A = 1.2905e-7*3;
    A_0 = 0.0;
    
    max_IGF1 = 1.2e-13;
    max_EGF = 6e-13;
    max_CSF1 = 4e-12;
    max_Gal9_en = 4e-12; // flexible
    
    K3 = 0.7;
    
    if_CSF1R_I = 0;
    if_IGF1R_I = 0;
    if_EGFR_I = 0;   
    
    for(int i=0; i<100; ++i){
        for(int j=0; j<100; ++j){
            CSF1[i][j] = 0;
            EGF[i][j] = 0;
            IGF1[i][j] = 0;
            Gal9_en[i][j] = 0;
            I[i][j] = 0;
            
            CSF1R_I[i][j] = 0;
            Gal9_I[i][j] = 0;
            IGF1R_I[i][j] = 0;
            EGFR_I[i][j] = 0;
            AKT_I[i][j] = 0;
            ERK_I[i][j] = 0;
            GSK3b_I[i][j] = 0;
            IRF1_I[i][j] = 0;
            
            AKT[i][j] = 0;
            ERK[i][j] = 0;
            IGF1R[i][j] = 0;
            EGFR[i][j] = 0;
            GSK3b[i][j] = 0;
            IRF1[i][j] = 0;
            Gal9[i][j] = 0;
            
            prob[i][j] = 0;
            prob_0[i][j] = 0;
        }
    }
}

// diffusion
void Diffusibles::diffusion(CellGrids &cg, double tstep, double dose_c,
                            double dose_i, double dose_a) {
    
    // if feeding CSF1R_I
    double S_a = 0;
    if(time * 60 * 60 >= drug_feed_time){
        drug_time = drug_time + tstep * 60 * 60;
        S_a = S_A;
    }
    //std::cout<<"drug_time = "<<drug_time<<std::endl;
    
    // break condition
    double maxDif;
    double c;
    
    // steps of forward Eular method
    double t = 60 * 60 * tstep / dt; 
    double t_drug = 60 * 60 * tstep / dt_drug;
    
    // finite difference & forward Eular method
    for (int q = 0; q < t; q++) {
        for(int i=1; i<99; ++i){
            for(int j=1; j<99; ++j){
                if(cg.vas[i][j] == 1){
                    IGF1R_I[i][j] = IGF1R_I[i][j] + q_IGF1R_I * (dose_i - IGF1R_I[i][j]);
                    AKT_I[i][j] = AKT_I[i][j] + q_AKT_I * (dose_a - AKT_I[i][j]);
                }      
            }
        }
        maxDif = 0;
        for (int i = 1; i < 99; i++) {
            for (int j = 1; j < 99; j++) {
                double CSF1_0 = CSF1[i][j];
                double Gal9_en_0 = Gal9_en[i][j];       
                double EGF_0 = EGF[i][j];  
                double IGF1_0 = IGF1[i][j]; 
                double IGF1R_I_0 = IGF1R_I[i][j];  
                double AKT_I_0 = AKT_I[i][j]; 
                
                CSF1[i][j] = CSF1[i][j] + dt * ((S_CSF1*cg.tc[i][j]) * (1 - CSF1[i][j]/max_CSF1)) - dt * d_CSF1 * CSF1[i][j]
                + (dt * D_CSF1 / (dx * dx)) * (CSF1[i + 1][j] + CSF1[i - 1][j]+ CSF1[i][j + 1] + CSF1[i][j - 1]- 4 * CSF1[i][j]);
                
                Gal9_en[i][j] = Gal9_en[i][j] + dt * ((S_Gal9_en*cg.tc[i][j]) * Gal9[i][j]/(K3 + Gal9[i][j])) - dt * d_Gal9_en * Gal9_en[i][j]
                + (dt * D_Gal9_en / (dx * dx)) * (Gal9_en[i + 1][j] + Gal9_en[i - 1][j]+ Gal9_en[i][j + 1] + Gal9_en[i][j - 1]- 4 * Gal9_en[i][j]);
                
                EGF[i][j] = EGF[i][j] + dt * ((S_EGF*cg.m2[i][j]) * (1 - EGF[i][j]/max_EGF)) - dt * d_EGF * EGF[i][j]
                + (dt * D_EGF / (dx * dx)) * (EGF[i + 1][j] + EGF[i - 1][j] + EGF[i][j + 1] + EGF[i][j - 1] - 4 * EGF[i][j]);
                
                IGF1[i][j] = IGF1[i][j] + dt * ((S_IGF1*cg.m2[i][j]) * (1 - IGF1[i][j]/max_IGF1)) - dt * d_IGF1 * IGF1[i][j]
                + (dt * D_IGF1 / (dx * dx)) * (IGF1[i + 1][j] + IGF1[i - 1][j]+ IGF1[i][j + 1] + IGF1[i][j - 1]- 4 * IGF1[i][j]); 
                
                IGF1R_I[i][j] = IGF1R_I[i][j] + ((dt_drug) * D_IGF1R_I / (dx * dx))
                    * (IGF1R_I[i + 1][j] + IGF1R_I[i - 1][j]+ IGF1R_I[i][j + 1] + IGF1R_I[i][j - 1]- 4 * IGF1R_I[i][j])
                    - (dt_drug) * (d_IGF1R_I + cg.tc[i][j] * u_IGF1R_I) * IGF1R_I[i][j];
                    
                AKT_I[i][j] = AKT_I[i][j] + ((dt_drug) * D_AKT_I / (dx * dx))
                    * (AKT_I[i + 1][j] + AKT_I[i - 1][j]+ AKT_I[i][j + 1] + AKT_I[i][j - 1]- 4 * AKT_I[i][j])
                    - (dt_drug) * (d_AKT_I + cg.tc[i][j] * u_AKT_I) * AKT_I[i][j];   
                    
                // check positivity
                if(EGF[i][j] < 0 || CSF1[i][j] < 0 || Gal9_en[i][j] < 0 ||
                    IGF1[i][j] < 0 || IGF1R_I[i][j] < 0 || AKT_I[i][j] < 0){
                    std::cout << "EGF " << EGF[i][j] << std::endl;
                    std::cout << "CSF1 " << CSF1[i][j] << std::endl;
                    std::cout << "Gal9_en " << Gal9_en[i][j] << std::endl;
                    std::cout << "IGF1 " << IGF1[i][j] << std::endl;
                    std::cout << "IGF1R_I " << IGF1R_I[i][j] << std::endl;
                    std::cout << "AKT_I " << AKT_I[i][j] << std::endl;
                    throw std::runtime_error("Diffusibles::diffusion error 1"); 
                }
                

                /*
                Gal9_I[i][j] = Gal9_I[i][j] + ((dt_drug) * D_Gal9_I / (dx * dx)) 
                    * (Gal9_I[i + 1][j] + Gal9_I[i - 1][j]+ Gal9_I[i][j + 1] + Gal9_I[i][j - 1]- 4 * Gal9_I[i][j])
                    - (dt_drug) * (d_Gal9_I + (cg.m0[i][j]+cg.m1[i][j]+cg.m2[i][j]) * u_Gal9_I) * Gal9_I[i][j];
                    

                    
                EGFR_I[i][j] = EGFR_I[i][j] + ((dt_drug) * D_EGFR_I / (dx * dx))
                    * (EGFR_I[i + 1][j] + EGFR_I[i - 1][j]+ EGFR_I[i][j + 1] + EGFR_I[i][j - 1]- 4 * EGFR_I[i][j])
                    - (dt_drug) * (d_EGFR_I + cg.tc[i][j] * u_EGFR_I) * EGFR_I[i][j];
                    

                        
                ERK_I[i][j] = ERK_I[i][j] + ((dt_drug) * D_ERK_I / (dx * dx))
                    * (ERK_I[i + 1][j] + ERK_I[i - 1][j]+ ERK_I[i][j + 1] + ERK_I[i][j - 1]- 4 * ERK_I[i][j])
                    - (dt_drug) * (d_ERK_I + cg.tc[i][j] * u_ERK_I) * ERK_I[i][j];
        
                GSK3b_I[i][j] = GSK3b_I[i][j] + ((dt_drug) * D_GSK3b_I / (dx * dx))
                    * (GSK3b_I[i + 1][j] + GSK3b_I[i - 1][j]+ GSK3b_I[i][j + 1] + GSK3b_I[i][j - 1]- 4 * GSK3b_I[i][j])
                    - (dt_drug) * (d_GSK3b_I + cg.tc[i][j] * u_GSK3b_I) * GSK3b_I[i][j];
                    
                IRF1_I[i][j] = IRF1_I[i][j] + ((dt_drug) * D_IRF1_I / (dx * dx))
                    * (IRF1_I[i + 1][j] + IRF1_I[i - 1][j]+ IRF1_I[i][j + 1] + IRF1_I[i][j - 1]- 4 * IRF1_I[i][j])
                    - (dt_drug) * (d_IRF1_I + cg.tc[i][j] * u_IRF1_I) * IRF1_I[i][j];
                                    
                 */
                if(CSF1_0 > 0){
                    c = (CSF1[i][j] - CSF1_0)/CSF1_0;
                    if(c > maxDif){maxDif=c;}
                }
                if(Gal9_en_0 > 0){
                    c = (Gal9_en[i][j] - Gal9_en_0)/Gal9_en_0;
                    if(c > maxDif){maxDif=c;}
                }
                if(EGF_0 > 0){
                    c = (EGF[i][j] - EGF_0)/EGF_0;
                    if(c > maxDif){maxDif=c;}
                }
                if(IGF1_0 > 0){
                    c = (IGF1[i][j] - IGF1_0)/IGF1_0;
                    if(c > maxDif){maxDif=c;}
                }
                if(IGF1R_I_0 > 0){
                    c = (IGF1R_I[i][j] - IGF1R_I_0)/IGF1R_I_0;
                    if(c > maxDif){maxDif=c;}
                }
                if(AKT_I_0 > 0){
                    c = (AKT_I[i][j] - AKT_I_0)/AKT_I_0;
                    if(c > maxDif){maxDif=c;}
                }
            }            
        }
        
        // zero flux boundary condition
        for(int i = 1; i < 99; i++){
            CSF1[i][0] = CSF1[i][2]; 
            CSF1[i][99] = CSF1[i][97];
            Gal9_en[i][0] = Gal9_en[i][2]; 
            Gal9_en[i][99] = Gal9_en[i][97];
            EGF[i][0] = EGF[i][2];
            EGF[i][99] = EGF[i][97];
            IGF1[i][0] = IGF1[i][2];
            IGF1[i][99] = IGF1[i][97];
            Gal9_I[i][0] = Gal9_I[i][2];
            Gal9_I[i][99] = Gal9_I[i][97];
            IGF1R_I[i][0] = IGF1R_I[i][2];
            IGF1R_I[i][99] = IGF1R_I[i][97];
            EGFR_I[i][0] = EGFR_I[i][2];
            EGFR_I[i][99] = EGFR_I[i][97];
            AKT_I[i][0] = AKT_I[i][2];
            AKT_I[i][99] = AKT_I[i][97];
            ERK_I[i][0] = ERK_I[i][2];
            ERK_I[i][99] = ERK_I[i][97];
            GSK3b_I[i][0] = GSK3b_I[i][2];
            GSK3b_I[i][99] = GSK3b_I[i][97];
            IRF1_I[i][0] = IRF1_I[i][2];
            IRF1_I[i][99] = IRF1_I[i][97];
        }
        for(int j = 1; j < 99; j++){
            CSF1[0][j] = CSF1[2][j];
            CSF1[99][j] = CSF1[97][j];
            Gal9_en[0][j] = Gal9_en[2][j];
            Gal9_en[99][j] = Gal9_en[97][j];
            EGF[0][j] = EGF[2][j];
            EGF[99][j] = EGF[97][j];
            IGF1[0][j] = IGF1[2][j];
            IGF1[99][j] = IGF1[97][j];
            Gal9_I[0][j] = Gal9_I[2][j];
            Gal9_I[99][j] = Gal9_I[97][j];
            IGF1R_I[0][j] = IGF1R_I[2][j];
            IGF1R_I[99][j] = IGF1R_I[97][j];
            EGFR_I[0][j] = EGFR_I[2][j];
            EGFR_I[99][j] = EGFR_I[97][j];
            AKT_I[0][j] = AKT_I[2][j];
            AKT_I[99][j] = AKT_I[97][j];
            ERK_I[0][j] = ERK_I[2][j];
            ERK_I[99][j] = ERK_I[97][j];
            GSK3b_I[0][j] = GSK3b_I[2][j];
            GSK3b_I[99][j] = GSK3b_I[97][j];
            IRF1_I[0][j] = IRF1_I[2][j];
            IRF1_I[99][j] = IRF1_I[97][j];
        }
        
        // break condition
        if(maxDif < 0.00001 && q > 5){
            std::cout<<"diff break q1 = "<<q<<std::endl;
            break;}      
    }
    
    // finite difference & forward Eular method
    for (int q = 0; q < t_drug; q++) {
        // feed drug 
        for(int i=1; i<99; ++i){
            for(int j=1; j<99; ++j){
                if(cg.vas[i][j] == 1){
                    CSF1R_I[i][j] = CSF1R_I[i][j] + q_CSF1R_I * (dose_c - CSF1R_I[i][j]);
                }      
            }
        }
        
        for (int i = 1; i < 99; i++) {
            for (int j = 1; j < 99; j++) {
                // integral parts in IGF1 concentration caculation
                I[i][j] = I[i][j] + dt_drug * S_a * CSF1R_I[i][j];
                
                CSF1R_I[i][j] = CSF1R_I[i][j] + ((dt_drug) * D_CSF1R_I / (dx * dx)) 
                    * (CSF1R_I[i + 1][j] + CSF1R_I[i - 1][j]+ CSF1R_I[i][j + 1] + CSF1R_I[i][j - 1]- 4 * CSF1R_I[i][j])
                    - (dt_drug) * (d_CSF1R_I + (cg.m0[i][j]+cg.m1[i][j]+cg.m2[i][j]) * u_CSF1R_I) * CSF1R_I[i][j];
                    
            }
        }
        
        // zero flux boundary condition
        for(int i = 1; i < 99; i++){
            CSF1R_I[i][0] = CSF1R_I[i][2];
            CSF1R_I[i][99] = CSF1R_I[i][97];
        }
        for(int j = 1; j < 99; j++){
            CSF1R_I[0][j] = CSF1R_I[2][j];
            CSF1R_I[99][j] = CSF1R_I[97][j];
        }    
    }
    
    // set zero, only for recording
    // no calculation here
    for(int i=0; i<100; ++i){
        for(int j=0; j<100; ++j){
            AKT[i][j] = 0;
            ERK[i][j] = 0;
            IGF1R[i][j] = 0;
            EGFR[i][j] = 0;
            GSK3b[i][j] = 0;
            IRF1[i][j] = 0;
            Gal9[i][j] = 0;
            prob[i][j] = 0;
            prob_0[i][j] = 0;
        }
    }
    
    // record simulation time
    time = time + tstep;          
}








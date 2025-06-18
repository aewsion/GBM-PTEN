       for j=1:N
         a=1:N;
         a(a==j)=[];
        [rho_RG,pval_RG] = partialcorr(r_RG_O,P_O(:,j),P_O(:,a));
        [rho_K,pval_K] = partialcorr(r_K_O(~isnan(r_K_O)),P_O(~isnan(r_K_O),j),P_O(~isnan(r_K_O),a));
        PCC_RG(j)=rho_RG;
        PCC_K(j)=rho_K;
        Pvalue_RG(j)=pval_RG;
        Pvalue_K(j)=pval_K;
       end
       for i=1:1:200
           if r_K(i)>=0
               a=0;
           else
               r_K(i)=0;
           end
       end
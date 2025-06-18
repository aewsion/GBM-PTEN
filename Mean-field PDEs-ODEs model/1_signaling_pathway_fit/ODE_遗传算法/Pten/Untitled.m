global Ptest1
load('pa.mat')
Ptest1=Ptest1;
 Conditions_ODES_2;
    tick1={'EGFR','ERK','IGF1R','AKT','GSK','IRF1','Gal9'};
    b1=[147/255,224/255,255/255];
    b2=[153/255,77/255,82/255];
        figure,
        plot(1:7,Y1,'-o','color',b1,'MarkerFaceColor',b1);hold on;
        xticks(1:1:7)
        set(gca,'Xticklabel',tick1,'FontWeight','Bold','FontSize',10);
        set(gca,'XTickLabelRotation',45);
        axis([1 7 0 1 ]);

  
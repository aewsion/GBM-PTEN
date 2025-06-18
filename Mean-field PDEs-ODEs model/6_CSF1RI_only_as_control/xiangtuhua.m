   load('D:\Matlab代码存放\2018胶质瘤\总程序\5_15画图\肿瘤模拟图1\CSF1RI_1_a_ERK_1_a_A_0.05_grow.mat','u1','u2','u3','u4','u5','u6','u7','u8','u9','u10','C_T_mean','M1_mean','M2_mean');
   C_T_mean1=C_T_mean;M1_mean1=M1_mean;M2_mean1=M2_mean;
   [aaa1,aaa2]=min(C_T_mean1(25*24+1:200*24+1,1));
   aaa3=25*24+aaa2;
   aaa4=find(C_T_mean1(25*24+1:200*24+1,1)<=0.01);
   aaa5=25*24+aaa4(1,1);
   load('D:\Matlab代码存放\2018胶质瘤\总程序\5_15画图\肿瘤模拟图1\CSF1RI_1_a_ERK_1_a_A_2.5_grow.mat','u1','u2','u3','u4','u5','u6','u7','u8','u9','u10','C_T_mean','M1_mean','M2_mean');
   C_T_mean2=C_T_mean;M1_mean2=M1_mean;M2_mean2=M2_mean;
   [bbb1,bbb2]=min(C_T_mean2(25*24+1:200*24+1,1));
   bbb3=25*24+bbb2;
    


%%%画图一
   figure,
   for i=1:601
   plot(M1_mean1(i,1),C_T_mean1(i,1),'.-','color',[250/255-i/601*(250-137)/255,227/255-i/601*(227-157)/255,113/255-i/601*(113-192)/255],'LineWidth',2.5);
   hold on
   end
   plot(M1_mean1(601,1),C_T_mean1(601,1),'*-','color',[137/255,157/255,192/255],'LineWidth',10);hold on,
   plot(M1_mean1(601:4801,1),C_T_mean1(601:4801,1),'--','color',[137/255,157/255,192/255],'LineWidth',2.5);hold on,
%    for i=601:aaa5
%    plot(M1_mean1(i,1),C_T_mean1(i,1),'.-','color',[250/255-(i-601)/(aaa5-601)*(250-147)/255,227/255-(i-601)/(aaa5-601)*(227-224)/255,113/255-(i-601)/(aaa5-601)*(113-255)/255],'LineWidth',2.5);
%    hold on
%    end
%    for i=aaa5:4801
%    plot(M1_mean1(i,1),C_T_mean1(i,1),'b .-','LineWidth',2.5);
%    hold on
%    end
   xlim([0 1.1]),ylim([0 1]);
   xlabel('Average M1 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
   title('Positive patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);
%    colorbar
   
   figure,
   for i=1:601
   plot(M1_mean2(i,1),C_T_mean2(i,1),'.-','color',[137/255-i/601*(137-217)/255,157/255-i/601*(157-104)/255,192/255-i/601*(192-49)/255],'LineWidth',2.5);
   hold on
   end
   plot(M1_mean2(601,1),C_T_mean2(601,1),'*-','color',[217/255,104/255,49/255],'LineWidth',10);hold on,
   plot(M1_mean2(601:bbb3,1),C_T_mean2(601:bbb3,1),'--','color',[217/255,104/255,49/255],'LineWidth',2.5);hold on,
   for i=bbb3:4801
   plot(M1_mean2(i,1),C_T_mean2(i,1),'.-','color',[217/255-(i-bbb3)/(4801-bbb3)*(217-6)/255,104/255-(i-bbb3)/(4801-bbb3)*(104-128)/255,49/255-(i-bbb3)/(4801-bbb3)*(49-67)/255],'LineWidth',2.5);
   hold on
   end   
%    figure,
%    plot(M1_mean2(1:601,1),C_T_mean2(1:601,1),'r-','LineWidth',2.5);
%    hold on,plot(M1_mean2(601,1),C_T_mean2(601,1),'r *-','LineWidth',10);
%    hold on,plot(M1_mean2(601:bbb3,1),C_T_mean2(601:bbb3,1),'r--','LineWidth',2.5);
%    hold on,plot(M1_mean2(bbb3:4801,1),C_T_mean2(bbb3:4801,1),'r-','LineWidth',2.5);
   xlim([0 1]),ylim([0 1]);
   xlabel('Average M1 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
   title('Resistant patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);  
   
   figure,
   for i=1:601
   plot(M2_mean1(i,1),C_T_mean1(i,1),'.-','color',[250/255-i/601*(250-137)/255,227/255-i/601*(227-157)/255,113/255-i/601*(113-192)/255],'LineWidth',2.5);
   hold on
   end
   plot(M2_mean1(601,1),C_T_mean1(601,1),'*-','color',[137/255,157/255,192/255],'LineWidth',10);hold on,
   plot(M2_mean1(601:4801,1),C_T_mean1(601:4801,1),'--','color',[137/255,157/255,192/255],'LineWidth',2.5);hold on,  
%    figure,
%    plot(M2_mean1(1:601,1),C_T_mean1(1:601,1),'b-','LineWidth',2.5);
%    hold on,plot(M2_mean1(601,1),C_T_mean1(601,1),'b *-','LineWidth',10);
%    hold on,plot(M2_mean1(601:4801,1),C_T_mean1(601:4801,1),'b--','LineWidth',2.5);
   xlim([0 1]),ylim([0 1]);
   xlabel('Average M2 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
   title('Positive patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);
   figure,
   for i=1:601
   plot(M2_mean2(i,1),C_T_mean2(i,1),'.-','color',[137/255-i/601*(137-217)/255,157/255-i/601*(157-104)/255,192/255-i/601*(192-49)/255],'LineWidth',2.5);
   hold on
   end
   plot(M2_mean2(601,1),C_T_mean2(601,1),'*-','color',[217/255,104/255,49/255],'LineWidth',10);hold on,
   plot(M2_mean2(601:bbb3,1),C_T_mean2(601:bbb3,1),'--','color',[217/255,104/255,49/255],'LineWidth',2.5);hold on,
   for i=bbb3:4801
   plot(M2_mean2(i,1),C_T_mean2(i,1),'.-','color',[217/255-(i-bbb3)/(4801-bbb3)*(217-6)/255,104/255-(i-bbb3)/(4801-bbb3)*(104-128)/255,49/255-(i-bbb3)/(4801-bbb3)*(49-67)/255],'LineWidth',2.5);
   hold on
   end 
%    plot(M2_mean2(1:601,1),C_T_mean2(1:601,1),'r-','LineWidth',2.5);
%    hold on,plot(M2_mean2(601,1),C_T_mean2(601,1),'r *-','LineWidth',10);
%    hold on,plot(M2_mean2(601:bbb3,1),C_T_mean2(601:bbb3,1),'r--','LineWidth',2.5);
%    hold on,plot(M2_mean2(bbb3:4801,1),C_T_mean2(bbb3:4801,1),'r-','LineWidth',2.5);
   xlim([0 1]),ylim([0 1]);
   xlabel('Average M2 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
   title('Resistant patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);  
   
   
   
 %%%画图二  
   color_r=zeros(4801,1);
color_g=zeros(4801,1);
color_b=zeros(4801,1);
c1=[250 227 113];c2=[137 157 192];c3=[217 104 49];c4=[56 13 48];
color_r(1:601,1)=c1(1)/255:(c2(1)-c1(1))/255/600:c2(1)/255;
color_g(1:601,1)=c1(2)/255:(c2(2)-c1(2))/255/600:c2(2)/255;
color_b(1:601,1)=c1(3)/255:(c2(3)-c1(3))/255/600:c2(3)/255;
color_r(601:bbb3,1)=c2(1)/255:(c3(1)-c2(1))/255/(bbb3-601):c3(1)/255;
color_g(601:bbb3,1)=c2(2)/255:(c3(2)-c2(2))/255/(bbb3-601):c3(2)/255;
color_b(601:bbb3,1)=c2(3)/255:(c3(3)-c2(3))/255/(bbb3-601):c3(3)/255;
color_r(bbb3:4801,1)=c3(1)/255:(c4(1)-c3(1))/255/(4801-bbb3):c4(1)/255;
color_g(bbb3:4801,1)=c3(2)/255:(c4(2)-c3(2))/255/(4801-bbb3):c4(2)/255;
color_b(bbb3:4801,1)=c3(3)/255:(c4(3)-c3(3))/255/(4801-bbb3):c4(3)/255;
c=zeros(4801,3);
c(:,1)=color_r;
c(:,2)=color_g;
c(:,3)=color_b;
%敏感性病人
figure
% color_r=zeros(4801,1);
% color_g=zeros(4801,1);
% color_b=zeros(4801,1);
% c1=[250 227 113];c2=[137 157 192];c3=[217 104 49];c4=[56 13 48];
% color_r(1:601,1)=c1(1)/255:(c2(1)-c1(1))/255/600:c2(1)/255;
% color_g(1:601,1)=c1(2)/255:(c2(2)-c1(2))/255/600:c2(2)/255;
% color_b(1:601,1)=c1(3)/255:(c2(3)-c1(3))/255/600:c2(3)/255;
% color_r(601:aaa5,1)=c2(1)/255:(c3(1)-c2(1))/255/(aaa5-601):c3(1)/255;
% color_g(601:aaa5,1)=c2(2)/255:(c3(2)-c2(2))/255/(aaa5-601):c3(2)/255;
% color_b(601:aaa5,1)=c2(3)/255:(c3(3)-c2(3))/255/(aaa5-601):c3(3)/255;
% color_r(aaa5:4801,1)=c3(1)/255:(c4(1)-c3(1))/255/(4801-aaa5):c4(1)/255;
% color_g(aaa5:4801,1)=c3(2)/255:(c4(2)-c3(2))/255/(4801-aaa5):c4(2)/255;
% color_b(aaa5:4801,1)=c3(3)/255:(c4(3)-c3(3))/255/(4801-aaa5):c4(3)/255;
% c=zeros(4801,3);
% c(:,1)=color_r;
% c(:,2)=color_g;
% c(:,3)=color_b;
sz=25;
% map = colormap(c);
C=0:200/4800:200;
% map = colormap(C);
scatter(M1_mean1,C_T_mean1,sz,C,'filled');
% C=0:200/4800:200;
% scatter(M1_mean1,C_T_mean1,sz,C,'filled');
colorbar;
   caxis([0 200]);
   h=colorbar('xtick',[0,50,100,150,200]);
   xlim([0 1.2]),ylim([0 1]);
   set(gca,'xtick',0:0.2:1.2); set(gca,'xticklabel',0:0.2:1.2);
   xlabel('Average M1 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
%    title('Positive patient','FontWeight','Bold','FontSize',18);
   title('Sensitive patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);  
  
figure
% color_r=zeros(4801,1);
% color_g=zeros(4801,1);
% color_b=zeros(4801,1);
% c1=[250 227 113];c2=[137 157 192];c3=[217 104 49];c4=[56 13 48];
% color_r(1:601,1)=c1(1)/255:(c2(1)-c1(1))/255/600:c2(1)/255;
% color_g(1:601,1)=c1(2)/255:(c2(2)-c1(2))/255/600:c2(2)/255;
% color_b(1:601,1)=c1(3)/255:(c2(3)-c1(3))/255/600:c2(3)/255;
% color_r(601:aaa5,1)=c2(1)/255:(c3(1)-c2(1))/255/(aaa5-601):c3(1)/255;
% color_g(601:aaa5,1)=c2(2)/255:(c3(2)-c2(2))/255/(aaa5-601):c3(2)/255;
% color_b(601:aaa5,1)=c2(3)/255:(c3(3)-c2(3))/255/(aaa5-601):c3(3)/255;
% color_r(aaa5:4801,1)=c3(1)/255:(c4(1)-c3(1))/255/(4801-aaa5):c4(1)/255;
% color_g(aaa5:4801,1)=c3(2)/255:(c4(2)-c3(2))/255/(4801-aaa5):c4(2)/255;
% color_b(aaa5:4801,1)=c3(3)/255:(c4(3)-c3(3))/255/(4801-aaa5):c4(3)/255;
% c=zeros(4801,3);
% c(:,1)=color_r;
% c(:,2)=color_g;
% c(:,3)=color_b;
% C=1:4801;
sz=25;
% map = colormap(c);
% scatter(M2_mean1,C_T_mean1,sz,c,'filled');
scatter(M2_mean1,C_T_mean1,sz,C,'filled');
% pcolor(M1_mean1,C_T_mean1,C)
colorbar
   caxis([0 200])
   h=colorbar('xtick',[0,50,100,150,200])
   xlim([0 1.2]),ylim([0 1]);
   set(gca,'xtick',0:0.2:1.2); set(gca,'xticklabel',0:0.2:1.2);
   xlabel('Average M2 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
%    title('Positive patient','FontWeight','Bold','FontSize',18);
   title('Sensitive patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);    
 %耐药病人
   figure
% color_r=zeros(4801,1);
% color_g=zeros(4801,1);
% color_b=zeros(4801,1);
% c1=[250 227 113];c2=[137 157 192];c3=[217 104 49];c4=[56 13 48];
% color_r(1:601,1)=c1(1)/255:(c2(1)-c1(1))/255/600:c2(1)/255;
% color_g(1:601,1)=c1(2)/255:(c2(2)-c1(2))/255/600:c2(2)/255;
% color_b(1:601,1)=c1(3)/255:(c2(3)-c1(3))/255/600:c2(3)/255;
% color_r(601:bbb3,1)=c2(1)/255:(c3(1)-c2(1))/255/(bbb3-601):c3(1)/255;
% color_g(601:bbb3,1)=c2(2)/255:(c3(2)-c2(2))/255/(bbb3-601):c3(2)/255;
% color_b(601:bbb3,1)=c2(3)/255:(c3(3)-c2(3))/255/(bbb3-601):c3(3)/255;
% color_r(bbb3:4801,1)=c3(1)/255:(c4(1)-c3(1))/255/(4801-bbb3):c4(1)/255;
% color_g(bbb3:4801,1)=c3(2)/255:(c4(2)-c3(2))/255/(4801-bbb3):c4(2)/255;
% color_b(bbb3:4801,1)=c3(3)/255:(c4(3)-c3(3))/255/(4801-bbb3):c4(3)/255;
% c=zeros(4801,3);
% c(:,1)=color_r;
% c(:,2)=color_g;
% c(:,3)=color_b;
% C=1:4801;
sz=25;
% scatter(M1_mean2,C_T_mean2,sz,c,'filled');
scatter(M1_mean2,C_T_mean2,sz,C,'filled');
% pcolor(M1_mean1,C_T_mean1,C)
colorbar
% map = colormap(c);
   caxis([0 200])
   h=colorbar('xtick',[0,50,100,150,200])
   xlim([0 1.2]),ylim([0 1]);
   set(gca,'xtick',0:0.2:1.2); set(gca,'xticklabel',0:0.2:1.2);
   xlabel('Average M1 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
%    title('Positive patient','FontWeight','Bold','FontSize',18);
   title('Resistant patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);  
  
figure
% color_r=zeros(4801,1);
% color_g=zeros(4801,1);
% color_b=zeros(4801,1);
% c1=[250 227 113];c2=[137 157 192];c3=[217 104 49];c4=[56 13 48];
% color_r(1:601,1)=c1(1)/255:(c2(1)-c1(1))/255/600:c2(1)/255;
% color_g(1:601,1)=c1(2)/255:(c2(2)-c1(2))/255/600:c2(2)/255;
% color_b(1:601,1)=c1(3)/255:(c2(3)-c1(3))/255/600:c2(3)/255;
% color_r(601:bbb3,1)=c2(1)/255:(c3(1)-c2(1))/255/(bbb3-601):c3(1)/255;
% color_g(601:bbb3,1)=c2(2)/255:(c3(2)-c2(2))/255/(bbb3-601):c3(2)/255;
% color_b(601:bbb3,1)=c2(3)/255:(c3(3)-c2(3))/255/(bbb3-601):c3(3)/255;
% color_r(bbb3:4801,1)=c3(1)/255:(c4(1)-c3(1))/255/(4801-bbb3):c4(1)/255;
% color_g(bbb3:4801,1)=c3(2)/255:(c4(2)-c3(2))/255/(4801-bbb3):c4(2)/255;
% color_b(bbb3:4801,1)=c3(3)/255:(c4(3)-c3(3))/255/(4801-bbb3):c4(3)/255;
% c=zeros(4801,3);
% c(:,1)=color_r;
% c(:,2)=color_g;
% c(:,3)=color_b;
% C=1:4801;
sz=25;
% scatter(M2_mean2,C_T_mean2,sz,c,'filled');
scatter(M2_mean2,C_T_mean2,sz,C,'filled');
% pcolor(M1_mean1,C_T_mean1,C)
colorbar
% map = colormap(c);
   caxis([0 200])
   h=colorbar('xtick',[0,50,100,150,200])
   xlim([0 1.2]),ylim([0 1]);
   set(gca,'xtick',0:0.2:1.2); set(gca,'xticklabel',0:0.2:1.2);
   xlabel('Average M2 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
%    title('Positive patient','FontWeight','Bold','FontSize',18);
   title('Resistant patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);
   
   
   
   
   %%%画图三  
   c21=[137,157,192]/255;c22=[250,227,113]/255;c23=[0,90,171]/255;
   figure,
   plot(M1_mean1(1:601,1),C_T_mean1(1:601,1),'--','color',c21,'LineWidth',2.5);hold on,
   plot(M1_mean1(601,1),C_T_mean1(601,1),'*-','color',c21,'LineWidth',10);hold on,
   for i=601:aaa5
   plot(M1_mean1(i,1),C_T_mean1(i,1),'.-','color',c21-(i-601)/(aaa5-601)*(c21-c22),'LineWidth',2.5);
   hold on
   end
   for i=aaa5:4801
   plot(M1_mean1(i,1),C_T_mean1(i,1),'.-','color',c22-(i-aaa5)/(4801-aaa5)*(c22-c23),'LineWidth',2.5);
   hold on
   end
   xlim([0 1.1]),ylim([0 1]);
   xlabel('Average M1 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
   title('Positive patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);
   
   figure,
   plot(M2_mean1(1:601,1),C_T_mean1(1:601,1),'--','color',c21,'LineWidth',2.5);hold on,
   plot(M2_mean1(601,1),C_T_mean1(601,1),'*-','color',c21,'LineWidth',10);hold on,
   for i=601:aaa5
   plot(M2_mean1(i,1),C_T_mean1(i,1),'.-','color',c21-(i-601)/(aaa5-601)*(c21-c22),'LineWidth',2.5);
   hold on
   end
   for i=aaa5:4801
   plot(M2_mean1(i,1),C_T_mean1(i,1),'.-','color',c22-(i-aaa5)/(4801-aaa5)*(c22-c23),'LineWidth',2.5);
   hold on
   end
   xlim([0 1.1]),ylim([0 1]);
   xlabel('Average M2 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
   title('Positive patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);
   
   
   %耐药性
   c31=[137,157,192]/255;c32=[250,227,113]/255;c33=[90,13,67]/255;
   figure,
   plot(M1_mean2(1:601,1),C_T_mean2(1:601,1),'--','color',c31,'LineWidth',2.5);hold on,
   plot(M1_mean2(601,1),C_T_mean2(601,1),'*-','color',c31,'LineWidth',10);hold on,
   for i=601:bbb3
   plot(M1_mean2(i,1),C_T_mean2(i,1),'.-','color',c31-(i-601)/(bbb3-601)*(c31-c32),'LineWidth',2.5);
   hold on
   end
   for i=bbb3:4801
   plot(M1_mean2(i,1),C_T_mean2(i,1),'.-','color',c32-(i-bbb3)/(4801-bbb3)*(c32-c33),'LineWidth',2.5);
   hold on
   end
   xlim([0 1.1]),ylim([0 1]);
   xlabel('Average M1 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
   title('Resistant patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);
   
   figure,
   plot(M2_mean2(1:601,1),C_T_mean2(1:601,1),'--','color',c31,'LineWidth',2.5);hold on,
   plot(M2_mean2(601,1),C_T_mean2(601,1),'*-','color',c31,'LineWidth',10);hold on,
   for i=601:bbb3
   plot(M2_mean2(i,1),C_T_mean2(i,1),'.-','color',c31-(i-601)/(bbb3-601)*(c31-c32),'LineWidth',2.5);
   hold on
   end
   for i=bbb3:4801
   plot(M2_mean2(i,1),C_T_mean2(i,1),'.-','color',c32-(i-bbb3)/(4801-bbb3)*(c32-c33),'LineWidth',2.5);
   hold on
   end
   xlim([0 1.1]),ylim([0 1]);
   xlabel('Average M2 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
   title('Resistant patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);
   
   
   
   %%%画图四  
%    c21=[137,157,192]/255;c22=[250,227,113]/255;c23=[0,90,171]/255;
c21=[1,1,0];c22=[0.1 0.7 0.7];c23=[0.2,0.2,0.7];
   figure,
   plot(M1_mean1(1:601,1),C_T_mean1(1:601,1),'--','color',c21,'LineWidth',2.5);hold on,
   plot(M1_mean1(601,1),C_T_mean1(601,1),'*-','color',c21,'LineWidth',10);hold on,
   for i=601:aaa5
   plot(M1_mean1(i,1),C_T_mean1(i,1),'.-','color',c21-(i-601)/(aaa5-601)*(c21-c22),'LineWidth',2.5);
   hold on
   end
   for i=aaa5:4801
   plot(M1_mean1(i,1),C_T_mean1(i,1),'.-','color',c22-(i-aaa5)/(4801-aaa5)*(c22-c23),'LineWidth',2.5);
   hold on
   end
   xlim([0 1.1]),ylim([0 1]);
   xlabel('Average M1 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
   title('Positive patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);
   
   figure,
   plot(M2_mean1(1:601,1),C_T_mean1(1:601,1),'--','color',c21,'LineWidth',2.5);hold on,
   plot(M2_mean1(601,1),C_T_mean1(601,1),'*-','color',c21,'LineWidth',10);hold on,
   for i=601:aaa5
   plot(M2_mean1(i,1),C_T_mean1(i,1),'.-','color',c21-(i-601)/(aaa5-601)*(c21-c22),'LineWidth',2.5);
   hold on
   end
   for i=aaa5:4801
   plot(M2_mean1(i,1),C_T_mean1(i,1),'.-','color',c22-(i-aaa5)/(4801-aaa5)*(c22-c23),'LineWidth',2.5);
   hold on
   end
   xlim([0 1.1]),ylim([0 1]);
   xlabel('Average M2 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
   title('Positive patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);
   
   
   %耐药性
   c31=[137,157,192]/255;c32=[250,227,113]/255;c33=[90,13,67]/255;
   figure,
   plot(M1_mean2(1:601,1),C_T_mean2(1:601,1),'--','color',c21,'LineWidth',2.5);hold on,
   plot(M1_mean2(601,1),C_T_mean2(601,1),'*-','color',c21,'LineWidth',10);hold on,
   for i=601:aaa5
   plot(M1_mean2(i,1),C_T_mean2(i,1),'.-','color',c21-(i-601)/(aaa5-601)*(c21-c22),'LineWidth',2.5);
   hold on
   end
   for i=aaa5:4801
   plot(M1_mean2(i,1),C_T_mean2(i,1),'.-','color',c22-(i-aaa5)/(4801-aaa5)*(c22-c23),'LineWidth',2.5);
   hold on
   end
   xlim([0 1.1]),ylim([0 1]);
   xlabel('Average M1 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
   title('Resistant patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);
   
   figure,
   plot(M2_mean2(1:601,1),C_T_mean2(1:601,1),'--','color',c21,'LineWidth',2.5);hold on,
   plot(M2_mean2(601,1),C_T_mean2(601,1),'*-','color',c21,'LineWidth',10);hold on,
   for i=601:aaa5
   plot(M2_mean2(i,1),C_T_mean2(i,1),'.-','color',c21-(i-601)/(aaa5-601)*(c21-c22),'LineWidth',2.5);
   hold on
   end
   for i=aaa5:4801
   plot(M2_mean2(i,1),C_T_mean2(i,1),'.-','color',c22-(i-aaa5)/(4801-aaa5)*(c22-c23),'LineWidth',2.5);
   hold on
   end
   xlim([0 1.1]),ylim([0 1]);
   xlabel('Average M2 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
   title('Resistant patient','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);
   
   
   
   
   
   %%%画图五   
color_r=zeros(4801,1);
color_g=zeros(4801,1);
color_b=zeros(4801,1);
% c1=[250 227 113];c2=[137 157 192];c3=[217 104 49];c4=[56 13 48];
c1=[0.1 0.7 0.1]*255;c2=[137 157 192];c3=[0.201 0.5 1]*255;c4=[0.2 0.2 0.7]*255;
% c1=[1 1 0.1]*255;c2=[0.8 0.8 0.2]*255;c3=[0 0.7 0.7]*255;c4=[0.2 0.2 0.71]*255;
color_r(1:601,1)=c1(1)/255:(c2(1)-c1(1))/255/600:c2(1)/255;
color_g(1:601,1)=c1(2)/255:(c2(2)-c1(2))/255/600:c2(2)/255;
color_b(1:601,1)=c1(3)/255:(c2(3)-c1(3))/255/600:c2(3)/255;
color_r(601:bbb3,1)=c2(1)/255:(c3(1)-c2(1))/255/(bbb3-601):c3(1)/255;
color_g(601:bbb3,1)=c2(2)/255:(c3(2)-c2(2))/255/(bbb3-601):c3(2)/255;
color_b(601:bbb3,1)=c2(3)/255:(c3(3)-c2(3))/255/(bbb3-601):c3(3)/255;
color_r(bbb3:4801,1)=c3(1)/255:(c4(1)-c3(1))/255/(4801-bbb3):c4(1)/255;
color_g(bbb3:4801,1)=c3(2)/255:(c4(2)-c3(2))/255/(4801-bbb3):c4(2)/255;
color_b(bbb3:4801,1)=c3(3)/255:(c4(3)-c3(3))/255/(4801-bbb3):c4(3)/255;
c=zeros(4801,3);
c(:,1)=color_r;
c(:,2)=color_g;
c(:,3)=color_b;
%敏感性病人
figure
% color_r=zeros(4801,1);
% color_g=zeros(4801,1);
% color_b=zeros(4801,1);
% c1=[250 227 113];c2=[137 157 192];c3=[217 104 49];c4=[56 13 48];
% color_r(1:601,1)=c1(1)/255:(c2(1)-c1(1))/255/600:c2(1)/255;
% color_g(1:601,1)=c1(2)/255:(c2(2)-c1(2))/255/600:c2(2)/255;
% color_b(1:601,1)=c1(3)/255:(c2(3)-c1(3))/255/600:c2(3)/255;
% color_r(601:aaa5,1)=c2(1)/255:(c3(1)-c2(1))/255/(aaa5-601):c3(1)/255;
% color_g(601:aaa5,1)=c2(2)/255:(c3(2)-c2(2))/255/(aaa5-601):c3(2)/255;
% color_b(601:aaa5,1)=c2(3)/255:(c3(3)-c2(3))/255/(aaa5-601):c3(3)/255;
% color_r(aaa5:4801,1)=c3(1)/255:(c4(1)-c3(1))/255/(4801-aaa5):c4(1)/255;
% color_g(aaa5:4801,1)=c3(2)/255:(c4(2)-c3(2))/255/(4801-aaa5):c4(2)/255;
% color_b(aaa5:4801,1)=c3(3)/255:(c4(3)-c3(3))/255/(4801-aaa5):c4(3)/255;
% c=zeros(4801,3);
% c(:,1)=color_r;
% c(:,2)=color_g;
% c(:,3)=color_b;
sz=25;
% map = colormap(c);
% C=0:200/4800:200;
% map = colormap(C);
figure
% scatter(M1_mean1(1:601,1),C_T_mean1(1:601,1),sz,c(1:601,:));hold on,
% scatter(M1_mean1(601:4801,1),C_T_mean1(601:4801,1),sz,c(601:4801,:),'filled');
% C=0:200/4800:200;
scatter(M1_mean1,C_T_mean1,sz,c,'filled');
map = colormap(c);
colorbar;
   caxis([0 200]);
   h=colorbar('xtick',[0,50,100,150,200]);
   xlim([0 1.2]),ylim([0 1]);
   set(gca,'xtick',0:0.2:1.2); set(gca,'xticklabel',0:0.2:1.2);
   xlabel('Average M1 macrophage density','FontWeight','Bold','FontSize',18,'FontName','Arial');
   ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18,'FontName','Arial');
%    title('Positive patient','FontWeight','Bold','FontSize',18);
   title('Responder','FontWeight','Bold','FontSize',18,'FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',18,'FontName','Arial');  
  
figure
% color_r=zeros(4801,1);
% color_g=zeros(4801,1);
% color_b=zeros(4801,1);
% c1=[250 227 113];c2=[137 157 192];c3=[217 104 49];c4=[56 13 48];
% color_r(1:601,1)=c1(1)/255:(c2(1)-c1(1))/255/600:c2(1)/255;
% color_g(1:601,1)=c1(2)/255:(c2(2)-c1(2))/255/600:c2(2)/255;
% color_b(1:601,1)=c1(3)/255:(c2(3)-c1(3))/255/600:c2(3)/255;
% color_r(601:aaa5,1)=c2(1)/255:(c3(1)-c2(1))/255/(aaa5-601):c3(1)/255;
% color_g(601:aaa5,1)=c2(2)/255:(c3(2)-c2(2))/255/(aaa5-601):c3(2)/255;
% color_b(601:aaa5,1)=c2(3)/255:(c3(3)-c2(3))/255/(aaa5-601):c3(3)/255;
% color_r(aaa5:4801,1)=c3(1)/255:(c4(1)-c3(1))/255/(4801-aaa5):c4(1)/255;
% color_g(aaa5:4801,1)=c3(2)/255:(c4(2)-c3(2))/255/(4801-aaa5):c4(2)/255;
% color_b(aaa5:4801,1)=c3(3)/255:(c4(3)-c3(3))/255/(4801-aaa5):c4(3)/255;
% c=zeros(4801,3);
% c(:,1)=color_r;
% c(:,2)=color_g;
% c(:,3)=color_b;
% C=1:4801;
sz=25;
% map = colormap(c);
scatter(M2_mean1,C_T_mean1,sz,c,'filled');
% scatter(M2_mean1,C_T_mean1,sz,C,'filled');
% pcolor(M1_mean1,C_T_mean1,C)
colorbar
map = colormap(c);
   caxis([0 200])
   h=colorbar('xtick',[0,50,100,150,200])
   xlim([0 1.2]),ylim([0 1]);
   set(gca,'xtick',0:0.2:1.2); set(gca,'xticklabel',0:0.2:1.2);
   xlabel('Average M2 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
%    title('Positive patient','FontWeight','Bold','FontSize',18);
   title('Responder','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);    
 %耐药病人
   figure
% color_r=zeros(4801,1);
% color_g=zeros(4801,1);
% color_b=zeros(4801,1);
% c1=[250 227 113];c2=[137 157 192];c3=[217 104 49];c4=[56 13 48];
% color_r(1:601,1)=c1(1)/255:(c2(1)-c1(1))/255/600:c2(1)/255;
% color_g(1:601,1)=c1(2)/255:(c2(2)-c1(2))/255/600:c2(2)/255;
% color_b(1:601,1)=c1(3)/255:(c2(3)-c1(3))/255/600:c2(3)/255;
% color_r(601:bbb3,1)=c2(1)/255:(c3(1)-c2(1))/255/(bbb3-601):c3(1)/255;
% color_g(601:bbb3,1)=c2(2)/255:(c3(2)-c2(2))/255/(bbb3-601):c3(2)/255;
% color_b(601:bbb3,1)=c2(3)/255:(c3(3)-c2(3))/255/(bbb3-601):c3(3)/255;
% color_r(bbb3:4801,1)=c3(1)/255:(c4(1)-c3(1))/255/(4801-bbb3):c4(1)/255;
% color_g(bbb3:4801,1)=c3(2)/255:(c4(2)-c3(2))/255/(4801-bbb3):c4(2)/255;
% color_b(bbb3:4801,1)=c3(3)/255:(c4(3)-c3(3))/255/(4801-bbb3):c4(3)/255;
% c=zeros(4801,3);
% c(:,1)=color_r;
% c(:,2)=color_g;
% c(:,3)=color_b;
% C=1:4801;
sz=25;
scatter(M1_mean2,C_T_mean2,sz,c,'filled');
% scatter(M1_mean2,C_T_mean2,sz,C,'filled');
% pcolor(M1_mean1,C_T_mean1,C)
colorbar
   caxis([0 200])
   h=colorbar('xtick',[0,50,100,150,200])
   map = colormap(c);
   xlim([0 1.2]),ylim([0 1]);
   set(gca,'xtick',0:0.2:1.2); set(gca,'xticklabel',0:0.2:1.2);
   xlabel('Average M1 macrophage density','FontWeight','Bold','FontSize',18,'FontName','Arial');
   ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18,'FontName','Arial');
%    title('Positive patient','FontWeight','Bold','FontSize',18);
   title('Nonresponder','FontWeight','Bold','FontSize',18,'FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',18,'FontName','Arial');  
  
figure
% color_r=zeros(4801,1);
% color_g=zeros(4801,1);
% color_b=zeros(4801,1);
% c1=[250 227 113];c2=[137 157 192];c3=[217 104 49];c4=[56 13 48];
% color_r(1:601,1)=c1(1)/255:(c2(1)-c1(1))/255/600:c2(1)/255;
% color_g(1:601,1)=c1(2)/255:(c2(2)-c1(2))/255/600:c2(2)/255;
% color_b(1:601,1)=c1(3)/255:(c2(3)-c1(3))/255/600:c2(3)/255;
% color_r(601:bbb3,1)=c2(1)/255:(c3(1)-c2(1))/255/(bbb3-601):c3(1)/255;
% color_g(601:bbb3,1)=c2(2)/255:(c3(2)-c2(2))/255/(bbb3-601):c3(2)/255;
% color_b(601:bbb3,1)=c2(3)/255:(c3(3)-c2(3))/255/(bbb3-601):c3(3)/255;
% color_r(bbb3:4801,1)=c3(1)/255:(c4(1)-c3(1))/255/(4801-bbb3):c4(1)/255;
% color_g(bbb3:4801,1)=c3(2)/255:(c4(2)-c3(2))/255/(4801-bbb3):c4(2)/255;
% color_b(bbb3:4801,1)=c3(3)/255:(c4(3)-c3(3))/255/(4801-bbb3):c4(3)/255;
% c=zeros(4801,3);
% c(:,1)=color_r;
% c(:,2)=color_g;
% c(:,3)=color_b;
% C=1:4801;
sz=25;
scatter(M2_mean2,C_T_mean2,sz,c,'filled');
% scatter(M2_mean2,C_T_mean2,sz,C,'filled');
% pcolor(M1_mean1,C_T_mean1,C)
colorbar
map = colormap(c);
   caxis([0 200])
   h=colorbar('xtick',[0,50,100,150,200])
   xlim([0 1.2]),ylim([0 1]);
   set(gca,'xtick',0:0.2:1.2); set(gca,'xticklabel',0:0.2:1.2);
   xlabel('Average M2 macrophage density','FontWeight','Bold','FontSize',18,'FontName','Arial');
   ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18,'FontName','Arial');
%    title('Positive patient','FontWeight','Bold','FontSize',18);
   title('Nonresponder','FontWeight','Bold','FontSize',18,'FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',18,'FontName','Arial');
% Create a pseudo grid mapping with a fine mesh
clear
close all
clc
addpath ./Functions
addpath ./clcdalpha
addpath ./BT_Results
addpath ./BT1_Data
Global
up_BC = 'Slip';
initial_profile = 'Constant';
corr_block = 1.04; %or 1.08
Tip_Speed_Ratio = [12];
Angle_of_attack = figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(Tip_Speed_Ratio)
    TSR_theor = Tip_Speed_Ratio(i); TSR = TSR_theor;
%     folder = strcat(['./BT_Results/BT1_Dataset_final/RANS',num2str(TSR_theor)]);%,num2str(TSR_theor),'/',num2str(date),'/',up_BC,'/',out_BC,'/']);
    folder = strcat(['./BT_Results/TSR_',num2str(TSR_theor),'/',num2str(date),'/No_stall/',up_BC,'/',out_BC,'/']);
    load(strcat([folder,'Res.mat']));
%     fxt(i) = 0.02*TSR+0.4; 
    fxt(i) =.02*TSR_theor + .38;

    ftt(i) = 0.49;
    Lr = 0.5;

    Nr = length(Res.Ux(:,1));
    Pre_processing

    rGLC = Matrices.RmGLC(:,1);
    RANS = Res;
    Nt = length(TSR);
    D = .894; %m
    R = 0.5*D;
    U_infty = 10; %m/s
    Omega = TSR_theor*U_infty/R;
    r=flip(0:0.005:.55)';

% Psuedo blade characteristics
    Aerodyn_P=BlindTest1(r); % r chord twist Airfoil_class
    chord_P=Aerodyn_P(:,2);
    twist_P=Aerodyn_P(:,3);
    airfoil_P=Aerodyn_P(:,4);

    [~,R_ind_P]=min(abs(r-.5));
    [~,temp]=min(abs(r-nac));
    nac_ind_P=temp;

% Find relative incoming velocity and thus angle attack for Pseudo grid based on actual domain velocities
for k=1:Nt
    Index_C = Ind;
    Index = min(find(Matrices.XmGLC(1,:)-(x_T(k))<=0));
    [Vrel,alfa] = Angle_of_Attack([RANS.Ux(Ind_R:end-1,Index_C); RANS.Ux(end,Ind_hub)],RANS.Ut(Ind_R:end,Index_C));
    Vrel = [zeros(Nr-Nrd,1); Vrel]; alfa = [zeros(Nr-Nrd,1); alfa];
    Vrel_P(:,k)=interp1(rGLC,Vrel,r,'linear','extrap');
    Ux_int(:,k) = interp1(rGLC,RANS.Ux(:,Index_C),r);
    Ut_int(:,k) = interp1(rGLC,RANS.Ut(:,Index_C),r);
    alfa_P(:,k)=interp1(rGLC,alfa(:,k),r,'linear','extrap');
    
 % Calculate thrust and tangential forces for pseudo grid 
    AerodynForce_P=LiftDragBT1(r,alfa_P(:,k));
    CL_P(:,k)=AerodynForce_P(:,3);
    CD_P(:,k)=AerodynForce_P(:,4);
    CThrust_P(:,k)= CD_P(:,k).*sind(alfa_P(:,k)+twist_P) + CL_P(:,k).*cosd(alfa_P(:,k)+twist_P);
    CTang_P(:,k)= -CD_P(:,k).*cosd(alfa_P(:,k)+twist_P) + CL_P(:,k).*sind(alfa_P(:,k)+twist_P);
    CTang_P(:,k)=((r>r(nac_ind_P)).*(CTang_P(:,k)));
    CTang_P(:,k)=(CTang_P(:,k)>0).*CTang_P(:,k);
    CThrust_P(:,k)=((r>r(nac_ind_P)).*CThrust_P(:,k));
    CThrust_P(:,k)=((CThrust_P(:,k)>0).*CThrust_P(:,k));

%     g(i) = exp(-.2409*(3*TSR(k)-15.91))+.046;
    g(i)=exp(-.2431*(3*TSR_theor-16.18))+.03778;

    ghub = 1;
    
    % tip and hub corrections
    fxh(i) = nac; fth(i) = nac;
    
    f_tip_x_P=real(2/pi*acos(exp(-g(i)*3*(fxt(i)-r)./(2*r.*sind(alfa_P(:,k)+twist_P)))));
    f_hub_x_P=real(2/pi*acos(exp(-ghub*3*(r-fxh(i))./(2*r.*sind(alfa_P(:,k)+twist_P)))));
    f_x_P=f_tip_x_P(:,k).*f_hub_x_P(:,k);
    f_tip_t_P=real(2/pi*acos(exp(-g(i)*3*(ftt(i)-r)./(2*r.*sind(alfa_P(:,k)+twist_P)))));
    f_hub_t_P=real(2/pi*acos(exp(-ghub*3*(r-fth(i))./(2*r.*sind(alfa_P(:,k)+twist_P)))));
    f_t_P=f_tip_t_P(:,k).*f_hub_t_P(:,k);
    F_Thrust_P=((2*pi*r).^-1).*.5.*Vrel_P(:,k).^2.*chord_P.*f_x_P(:,k).*3.*CThrust_P(:,k);F_Thrust_P(end,k)=F_Thrust_P(end-1,k);
    F_Tang_P=((2*pi*r).^-1).*.5.*Vrel_P(:,k).^2.*chord_P.*f_t_P(:,k).*3.*(CTang_P(:,k));F_Tang_P(end,k)=0;
    Area_nac = sigma_r^2;
    F_nacelle_g_P=(.5.*Vrel_P(end,k).^2).*CDhub*exp(-0.5*(r/sigma_r).^2)*pi*nac^2/sigma_r;
    Fxtot=(F_nacelle_g_P > F_Thrust_P).*F_nacelle_g_P + (F_Thrust_P > F_nacelle_g_P).*F_Thrust_P;
end
%% Calculate Power and coefficients for this Pseudo grid
    Nr_P=length(r);
    for j=1:Nx
        Utinterp = interp1(rGLC,RANS.Ut(:,j),r);
        Uxinterp = interp1(rGLC,RANS.Ux(:,j),r);
        Ux_averaged(j) = 2*pi*trapz(flipud(r(R_ind_P:end)),...
        flipud(r(R_ind_P:end).*Uxinterp(R_ind_P:end)))/(0.25*pi);
        Ut_averaged(j) = 2*pi*trapz(flipud(r(R_ind_P:end)),...
        flipud(r(R_ind_P:end).*Utinterp(R_ind_P:end)))/(0.25*pi);
    end

    arg_P = chord_P.*Vrel_P(:,k).^2.*f_t_P(:,k).*CTang_P(:,k);
    Power_total_P = 3*rho*U_infty^3*TSR*D^2*trapz(flipud(r(R_ind_P:end)),flipud(r(R_ind_P:end)).*arg_P(R_ind_P:end));

    C_P(i) = Power_total_P/(0.5*rho*U_infty^3*pi*R^2);    

    arg_TH = chord_P.*Vrel_P(:,k).^2.*f_x_P(:,k).*CThrust_P(:,k);% + F_nacelle_g_P(:,k);
    arg_nac = exp(-0.5*(r/(sigma_r)).^2);
    F_nac_total_P = 0.5*rho*U_infty^2*Vrel_P(end,k)^2*pi*D^2*nac^2*CDhub;
    F_Thrust_total_P = 3/2*rho*U_infty^2*D^2*trapz(flipud(r(R_ind_P:end)),arg_TH(R_ind_P:end));
    C_TH(i) = (F_Thrust_total_P + F_nac_total_P)/(0.5*rho*U_infty^2*pi*R^2);

    Ux = Ux_averaged(:,Index); %Ux at turb. loc.
    Ut = Ut_averaged(:,Index); %Ut at turb. loc.

    a(i) = 1-Ux_averaged(Index);
        
    Res.C_P=C_P(i);
    Res.C_TH=C_TH(i);      
    [Fx1d,Ft1d] = Forcing([RANS.Ux(Ind_R:end-1,Ind);RANS.Ux(end,Ind_hub)],RANS.Ut(Ind_R:end,Ind));
%     clear chord_P Vrel_P CTang_P f_t_P
    Ux_external(i) = RANS.Ux(1,Ind);
    Ind_1=min(find(Matrices.XmGLC(1,:)-(x_T+1)<=0));
    Ux_external_1D(i) = RANS.Ux(1,Ind_1);
    
    hold on
    plot(rGLC,alfa,'linewidth',2);axis tight;Plot_options('r/D','\alpha','Relative incidence');
    drawnow
end
% savefig(Force_x,strcat(['./BT_Results/BT1_Dataset_final/Axial_force','.fig']));
% saveas(Force_x,strcat(['./BT_Results/BT1_Dataset_final/Axial_force','.png']));
%% figure on CPw
openfig('./BT1_DATA/Cp_BT1_Data.fig');
CP_fig = gcf;set(CP_fig,'units','normalized','outerposition',[0 0 0.5 1]);
hold on
plot(Tip_Speed_Ratio,C_P,'ks','markersize',10,'markerfacecolor','b');axis tight
Plot_options('TSR','C_P',[]);
leg = legend('Experiment','RANS global');set(leg,'location','south','fontsize',15,'fontname','arial');
% saveas(CP_fig,strcat(['./BT_Results/','CP','.png']));
%% figure on CTh
openfig('./BT1_DATA/Ct_BT1_Data.fig');
CT_fig = gcf;set(CT_fig,'units','normalized','outerposition',[0 0 0.5 1]);
hold on
plot(Tip_Speed_Ratio,C_TH,'ks','markersize',10);axis tight
Plot_options('TSR','C_T',[]);
leg = legend('Experiment','RANS global');set(leg,'location','south','fontsize',15,'fontname','arial');
% saveas(CT_fig,strcat(['./BT_Results/','C_T','.png']));
%% figure on axial forcing
Force_x = figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(3,1,1)
% plot(rGLC,Fx1d,'*-','linewidth',2);
% Plot_options('r/D','F_x',[]);axis tight
% subplot(3,1,2)
% plot(rGLC,Ft1d*2*pi.*rGLC,'*-','linewidth',2);
% Plot_options('r/D','F_{\theta}',[]);axis tight
% subplot(3,1,3)
% plot(rGLC,rGLC.^2.*Ft1d,'*-','linewidth',2);
% Plot_options('r/D','Power',[]);axis tight
% savefig(Force_x,strcat([folder,'/Force_x_TSR_',num2str(TSR_theor),'.fig']))
% saveas(Force_x,strcat([folder,'/Force_x_TSR_',num2str(TSR_theor),'.png']))
%% figure on f
% f_plot = figure('units','normalized','outerposition',[0 0 0.8 0.8]);
% plot(Tip_Speed_Ratio,fxt,'linewidth',5);
% hold on
% plot(Tip_Speed_Ratio,ftt,'linewidth',5);
% Plot_options('TSR','g',[]);
% leg = legend('fxt','ftt');
% set(leg,'location','southeast','fontsize',15);
% axis tight
% savefig(f_plot,'./BT_Results/ft.fig');
% saveas(f_plot,'./BT_Results/ft.png');
%% figure on blockage correction
block_plot = figure('units','normalized','outerposition',[0 0 0.8 0.8]);
plot(Tip_Speed_Ratio,Ux_external,'ko-','markersize',8,'linewidth',2,'markerfacecolor','k');
hold on
plot(Tip_Speed_Ratio,Ux_external_1D,'bo-','markersize',8,'linewidth',2,'markerfacecolor','b');
Plot_options('TSR','Blockage correction',[]);axis tight
leg = legend('x_T','x_T + 1');set(leg,'location','best','fontsize',15);
savefig(block_plot,strcat(['./BT_Results/BT1_Dataset_final/Blockage_correction','.fig']));
saveas(block_plot,strcat(['./BT_Results/BT1_Dataset_final/Blockage_correction','.png']));
%% Velocity plots
% folder = strcat(['./Results/Immagini/TSR_',num2str(TSR),'/xmap_',num2str(x_Map),'/']);
% mkdir(folder);

VEL = figure('units','normalized','outerposition',[0 0 0.5 1]);set(gcf,'Color','White');
subplot(4,1,1)
pcolor(Matrices.XmGLC,Matrices.RmGLC,RANS.Ur);
axis equal; axis ([xmin xmax 0 rmax]); colormap coolwarm; shading interp ; colorbar; 
hold on
plot_turbine
Plot_options('x/D','r/D','U_r');title(strcat(['TSR=',num2str(TSR_theor)]));
subplot(4,1,2)
pcolor(Matrices.XmGLC,Matrices.RmGLC,RANS.Ut);axis equal; axis ([xmin xmax 0 rmax]); colormap coolwarm; shading interp;colorbar;  
hold on
plot_turbine
Plot_options('x/D','r/D','U_{\theta}')
subplot(4,1,3)
pcolor(Matrices.XmGLC,Matrices.RmGLC,RANS.Ux); axis equal;axis ([xmin xmax 0 rmax]); colormap coolwarm; shading interp ;colorbar; caxis ([0 1.1]);
hold on
plot_turbine
Plot_options('x/D','r/D','U_x')
subplot(4,1,4)
pcolor(Matrices.XmGLC,Matrices.RmGLC,RANS.p); axis equal;axis ([xmin xmax 0 rmax]); colormap coolwarm; shading interp ;colorbar; 
hold on
plot_turbine
Plot_options('x/D','r/D','p')
savefig(VEL,strcat([folder,'Vel']));
saveas(VEL,strcat([folder,'Vel','.png']));
%% EV/ML
if isequal(Turbulence_model,'ML')==1 || isequal(Turbulence_model,'TKE')==1
    ML = ML_fun_MT(ML_coeff);
    mixing_length=figure('units','normalized','outerposition',[0 0 0.5 0.5]);set(gcf,'Color','White');
    plot(Matrices.XmGLC(1,:),ML(1,:),'o-');
    Plot_options('x/D','l_m(x)',[]);
%     savefig(mixing_length,strcat([folder,'/ml']));
%     saveas(mixing_length,strcat([folder,'/ml','.png']));
    if isequal(Turbulence_model,'ML')==1
        [NuT,DNuT_dUr,DNuT_dUt,DNuT_dUx,S] = EV_closure_All(0,ML,RANS.qc,NuT_smooth_par);
        NuT = reshape(NuT,Nx,Nr).';
    end
    if isequal(Turbulence_model,'TKE')==1
        NuT = 0.55*(RANS.tke).^0.5.*ML;
        TKE = figure('units','normalized','outerposition',[0 0 0.5 0.3]);set(gcf,'Color','White');
        pcolor(Matrices.XmGLC,Matrices.RmGLC,RANS.tke);axis equal; axis ([xmin xmax 0 rmax]); colormap coolwarm; shading interp ; colorbar
        hold on
        plot_turbine
        Plot_options('x/D','r/D','TKE');
        savefig(TKE,strcat([folder,'TKE']));
        saveas(TKE,strcat([folder,'TKE','.png']));
    end
    EV=figure('units','normalized','outerposition',[0 0 0.5 0.3]);set(gcf,'Color','White');
    pcolor(Matrices.XmGLC,Matrices.RmGLC,NuT);axis equal; axis ([xmin xmax 0 rmax]); colormap coolwarm; shading interp ; colorbar
    hold on
    plot_turbine
    Plot_options('x/D','r/D','\nu_T');
    savefig(EV,strcat([folder,'EV']));
    saveas(EV,strcat([folder,'EV','.png']));
end
%% Axial and radial velocity trends
% Ux_ax=figure('units','normalized','outerposition',[0 0 0.5 0.5]);set(gcf,'Color','White');
% plot(Matrices.XmGLC(1,:),Ux_averaged,'*-');
% Plot_options('x/D','U_x/U_{\infty}',strcat(['TSR=',num2str(TSR)]));
% savefig(Ux_ax,strcat([folder,'Ux_ax']));
% saveas(Ux_ax,strcat([folder,'Ux_ax','.png']));
%%
% Ux_rad=figure('units','normalized','outerposition',[0 0 0.5 0.5]);set(gcf,'Color','White');
% plot(Matrices.RmGLC(:,1),RANS.Ux(:,1),'*-');
% Plot_options('r/D','U_x/U_{\infty}',[]);
% savefig(Ux_rad,strcat([folder,'Ux_rad']));
% saveas(Ux_rad,strcat([folder,'Ux_rad','.png']));
%%
% Map=figure('units','normalized','outerposition',[0 0 0.5 0.5]);set(gcf,'Color','White');
% plot(Matrices.XmGLC(1,:),Matrices.xGLC,'*-');
% hold on
% plot(Matrices.XmGLC(1,:),Matrices.x1GLC,'o-');
% hold on
% plot(Matrices.XmGLC(1,:),Matrices.x2GLC,'+-');
% Plot_options('x/D','\eta',[]);
% leg = legend('\eta','d\eta/dx','d^2\eta/dx^2');set(leg,'fontsize',15,'location','southeast');
% savefig(Map,strcat([folder,'x_Map']));
% saveas(Map,strcat([folder,'x_Map','.png']));
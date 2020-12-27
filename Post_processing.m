%% Velocity plots
VEL = figure('units','normalized','outerposition',[0 0 0.5 1]);set(gcf,'Color','White');
subplot(4,1,1)
pcolor(Matrices.XmGLC,Matrices.RmGLC,RANS.Ur);
axis equal; axis ([xmin xmax 0 rmax]); colormap coolwarm; shading interp ; colorbar; 
hold on
plot_turbine
Plot_options('','$r/D$','$U_r$');%title(strcat(['TSR=',num2str(TSR_theor)]));
subplot(4,1,2)
pcolor(Matrices.XmGLC,Matrices.RmGLC,RANS.Ut);axis equal; axis ([xmin xmax 0 rmax]); colormap coolwarm; shading interp;colorbar;  
hold on
plot_turbine
Plot_options('','$r/D$','$U_{\theta}$')
subplot(4,1,3)
pcolor(Matrices.XmGLC,Matrices.RmGLC,RANS.Ux); axis equal;axis ([xmin xmax 0 rmax]); colormap coolwarm; shading interp ;colorbar; caxis ([0 1.1]);
hold on
plot_turbine
Plot_options('','$r/D$','$U_x$')
subplot(4,1,4)
pcolor(Matrices.XmGLC,Matrices.RmGLC,RANS.p); axis equal;axis ([xmin xmax 0 rmax]); colormap coolwarm; shading interp ;colorbar; 
hold on
plot_turbine
Plot_options('$x/D$','$r/D$','$p$')
savefig(VEL,strcat([folder,'Vel']));
saveas(VEL,strcat([folder,'Vel','.png']));
%% EV/ML
if isequal(Turbulence_model,'ML')==1 || isequal(Turbulence_model,'TKE')==1
    ML = ML_fun(ML_coeff);
    close all

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
        Plot_options('$x/D$','$r/D$','TKE');
        savefig(TKE,strcat([folder,'TKE']));
        saveas(TKE,strcat([folder,'TKE','.png']));
    end
    EV=figure('units','normalized','outerposition',[0 0 0.5 0.3]);set(gcf,'Color','White');
    pcolor(Matrices.XmGLC,Matrices.RmGLC,NuT);axis equal; axis ([xmin xmax 0 rmax]); colormap coolwarm; shading interp ; colorbar
    hold on
    plot_turbine
    Plot_options('$x/D$','$r/D$','$\nu_T$');
    savefig(EV,strcat([folder,'EV']));
    save2png(strcat([folder,'EV','.pdf']),EV,150);
end
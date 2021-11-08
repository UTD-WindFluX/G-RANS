clear all
close all
restoredefaultpath
addpath('./Functions')

%% Inputs
source_input='./Data/Input_file1.xlsx';

%% Main
GRANS_solver_v4

%% Plots
mkfig('max');
subplot(3,1,1)
pcolor(AXI.x,AXI.r,AXI.Ux);shading flat;axis equal;xlabel('$x$');ylabel('$r$');grid on;box on;colormap coolwarm
plot_turbine_fcn(0,0,0,0,1);axis([AXI.xmin AXI.xmax 0 AXI.rmax]);c=colorbar;c.Label.String='$u_x$';caxis([0.4 1.05]);
title(['\verb|',(AXI.source_input),'|: $c_t=',num2str(AXI.Ct),', \nu_T=',num2str(AXI.NuT),'$']);

subplot(3,1,2)
pcolor(AXI.x,AXI.r,AXI.Ur);shading flat;axis equal;xlabel('$x$');ylabel('$r$');grid on;box on;ax=gca;colormap(ax,'posneg');
plot_turbine_fcn(0,0,0,0,1);axis([AXI.xmin AXI.xmax 0 AXI.rmax]);c=colorbar;c.Label.String='$u_r$';caxis([-0.1 0.1]);

subplot(3,1,3)
pcolor(AXI.x,AXI.r,AXI.p);shading flat;axis equal;xlabel('$x$');ylabel('$r$');grid on;box on;ax=gca;colormap(ax,'posneg');
plot_turbine_fcn(0,0,0,0,1);axis([AXI.xmin AXI.xmax 0 AXI.rmax]);c=colorbar;c.Label.String='$p$';caxis([-0.1 0.1]);
TNR(12)
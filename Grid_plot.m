clear
close all
Global
% Nr = 40;
up_BC = 'Slip';
addpath ./Functions
initial_profile = 'Variable';
Pre_processing
%%
figure('units','normalized','outerposition',[0 0 1 0.5]);
plot(Matrices.XmGLC,Matrices.RmGLC,'k.','markersize',8);
hold on
plot_turbine
Plot_options('x/D','r/D',[]);axis equal
axis ([xmin xmax 0 rmax]);
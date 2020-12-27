clear
close all
clc
addpath ./Functions
addpath ./BT1_DATA
Tip_Speed_Ratio = 1:12;
Global
CL_CD %Evaluation of CL_alpha and CD_alpha
Buffer_zone = 0; %set 1 if EV must have a buffer zone on the last diameter; 0 otherwise
out_BC = 'Neumann'; %Choice of outlet boundary conditions
up_BC = 'Slip';
Under_relaxation_factor_tke = 1; %set 1 if tke production must be under-relaxed, 0 otherwise
Pre_ml_closure = 1;
initial_profile = 'Constant';
Turbulence_model = 'TKE';
corr_block = 1.04;
%% Direct solver
% for i=7:7
    TSR_theor = [3  3]; TSR = TSR_theor/corr_block; %
    Nr = 40;
    folder = strcat(['./BT_Results/TSR_',num2str(TSR_theor),'/',num2str(date),'/',up_BC,'/',out_BC,'/']);
    mkdir(folder);
    NuT_smooth_par = 0.0;
    Pre_processing %Here the input parameters are settled
    display(strcat(['TSR=',num2str(TSR_theor)]));
    BL_approx = 1;
    ML_coeff = [0 0.03 0.3142 x_T]; %ML modal coefficients
    DoFs = length(ML_coeff);
    tic
    if isequal(Turbulence_model,'TKE')==1 
        RANS = AxiSymDirect_TKE(2e-04,ML_coeff);end %Evaluate vel. field with TKE closure
    if isequal(Turbulence_model,'ML')==1  
        RANS = AxiSymDirect_ML(2e-04,ML_coeff);end %Evaluate vel. field with ML closure
    if isequal(Turbulence_model,'EV')==1  
        RANS = AxiSymDirect_EV([1e-03 5e-04]); end %Evaluate vel. field with EV closure
    time = toc;
%%
Post_processing
Res.Ux = RANS.Ux; Res.Ur = RANS.Ur; Res.Ut = RANS.Ut;Res.p = RANS.p;

if isequal(Turbulence_model,'TKE')==1 
    Res.tke = RANS.tke;end
Res.X = Matrices.XmGLC; Res.R = Matrices.RmGLC;
save(strcat([folder,'/Res.mat']),'Res');

if TSR_theor(1)==3 || TSR_theor(1)==6 || TSR_theor(1)==10
    Velocity_plot       %Plot axial velocity profiles over r @ 1,3 and 5 D past the rotor
    Velocity_plot_theta %Plot azimuthal velocity profiles over r @ 1,3 and 5 D past the rotor
end
close all
clc
% end
%% Validations
Mass_balance %Integral form of the mass conservation
Mom_balance_loc %Integral form of the axial momentum conservation
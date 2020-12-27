clear
close all
clc
addpath ./Functions                                         %Add the folder with functions to the main path
%% User-defined input
Global                                                      %Recall global variables
CL_CD                                                       %Evaluation of Cl,Cd, CL_alpha and CD_alpha
Pre_processing                                              %Here the input parameters are settled
%% Direct solver
folder = strcat(['./Results/TSR_',num2str(TSR),'/',Turbulence_model,'/',num2str(date),'/',up_BC,'/',out_BC,'/']);
mkdir(folder);                                              %Make folder to save results inside

display(strcat(['TSR=',num2str(TSR)]));               %Display the beginning of the iteration
BL_approx = 1;                                              %In case of Mixing Length, evaluate 2S:S through the boundary layer approximation
ML_coeff = [0 0.03 0.3142 x_T];                             %ML modal coefficients
DoFs = length(ML_coeff);                                    %Number of degrees of freedom

tic
if isequal(Turbulence_model,'TKE')                          %Evaluate vel. field with TKE closure
    RANS        = AxiSymDirect_TKE(2e-04,ML_coeff);
end             
    
if isequal(Turbulence_model,'ML')                           %Evaluate vel. field with ML closure
    RANS        = AxiSymDirect_ML(2e-04,ML_coeff);
end

if isequal(Turbulence_model,'EV')
    NuTCoeff    = [1e-03 1e-04];                            %Set constant value and slope of the EV model
    RANS        = AxiSymDirect_EV(NuTCoeff);                %AxiSymDirect_EV_mod evaluates EV closure from parabolic
end                 
time = toc;
%% Post-processing
Post_processing                                             %Plot of time-averaged flow quantities
Res.Ux  = RANS.Ux; 
Res.Ur  = RANS.Ur; 
Res.Ut  = RANS.Ut;
Res.p   = RANS.p;

if isequal(Turbulence_model,'TKE')
    Res.tke = RANS.tke;
end                                  %Add tke to the data structure if that is the chosen turbulence model
Res.X   = Matrices.XmGLC; 
Res.R   = Matrices.RmGLC;
save(strcat([folder,'/Res.mat']),'Res');                    %Save results
close all
clc
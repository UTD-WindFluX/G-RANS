clear
close all
clc
addpath ./Functions
Global
Nr = 40;Lr=0.5;
Tip_Speed_Ratio = 0;
corr_block = 1.08;
TSR_theor = Tip_Speed_Ratio; TSR = TSR_theor/corr_block;
initial_profile = 'Variable';
Pre_processing %Here the input parameters are settled
CL_CD
Turbulence_model = 'TKE';
Under_relaxation_factor_tke = 1; %set 1 if tke production must be under-relaxed, 0 otherwise
ML_coeff = [0 0.03 0.3142 x_T];
DoFs = length(ML_coeff);
BL_approx=1;

lm = ML_fun(ML_coeff);
out_BC = 'Neumann';
up_BC = 'Slip';
delta = 1e-09;
load('./BT_Results/Res_test.mat');

tke_fun = @(x,r) 1e-03*r+1e-03;
tke_ex = tke_fun(Matrices.XmGLC,Matrices.RmGLC);
% tke_ex = rand(Nr,Nx);
% tke_ex = tke_in*ones(Nr,Nx);
Initial_field = [reshape(Res.Ur.',Nx*Nr,1); reshape(Res.Ut.',Nx*Nr,1); reshape(Res.Ux.',Nx*Nr,1); zeros(NxGC*NrGC,1);reshape(tke_ex.',Nx*Nr,1)];

Ur = reshape(Initial_field(1:Nx*Nr),Nx,Nr).'; Ur2 = Ur; Ur1 = Ur;
Ux = reshape(Initial_field(2*Nx*Nr+1:3*Nx*Nr),Nx,Nr).'; Ux2 = Ux; Ux1 = Ux;
Ut = reshape(Initial_field(Nx*Nr+1:2*Nx*Nr),Nx,Nr).'; Ut2 = Ut; Ut1 = Ut;
%%
J_old.Ft_Ux=0;
J_old.Ft_Ut=0;
J_old.Fx_Ut=0;
J_old.Fx_Ux=0;
J_old.DP_DUr = zeros(Nx*Nr); J_old.DP_DUt = zeros(Nx*Nr);
J_old.DP_DUx = zeros(Nx*Nr); J_old.DP_Dtke = zeros(Nx*Nr);
[JA,~] = JacobianSponge_AD_TKE2(Initial_field,reshape(lm.',Nx*Nr,1),J_old,zeros(4*Nx*Nr+NxGC*NrGC,1),zeros(Nx*Nr,1),0);
%%
u_r = 30; u_x = 10; i_glob = (u_r-1)*Nx + u_x; %global residual index
v_r = 30; v_x = 10; j_glob = (v_r-1)*Nx + v_x; %Global velocity index
Ux2(v_r,v_x) = Ux(v_r,v_x)+delta; Ux1(v_r,v_x) = Ux(v_r,v_x)-delta;

qx2 = [Initial_field(1:Nx*Nr); reshape(Ut.',Nx*Nr,1); reshape(Ux2.',Nx*Nr,1); zeros(NxGC*NrGC,1);reshape(tke_ex.',Nx*Nr,1)];
qx1 = [Initial_field(1:Nx*Nr); reshape(Ut.',Nx*Nr,1); reshape(Ux1.',Nx*Nr,1); zeros(NxGC*NrGC,1);reshape(tke_ex.',Nx*Nr,1)];

Ut2(v_r,v_x) = Ut(v_r,v_x)+delta; Ut1(v_r,v_x) = Ut(v_r,v_x)-delta;
qt2 = [Initial_field(1:Nx*Nr); reshape(Ut2.',Nx*Nr,1); reshape(Ux.',Nx*Nr,1); zeros(NxGC*NrGC,1);reshape(tke_ex.',Nx*Nr,1)];
qt1 = [Initial_field(1:Nx*Nr); reshape(Ut1.',Nx*Nr,1); reshape(Ux.',Nx*Nr,1); zeros(NxGC*NrGC,1);reshape(tke_ex.',Nx*Nr,1)];

Ur2(v_r,v_x) = Ur(v_r,v_x)+delta; Ur1(v_r,v_x) = Ur(v_r,v_x)-delta;
qr2 = [reshape(Ur2.',Nx*Nr,1); reshape(Ut.',Nx*Nr,1); reshape(Ux.',Nx*Nr,1); zeros(NxGC*NrGC,1); reshape(tke_ex.',Nx*Nr,1)];
qr1 = [reshape(Ur1.',Nx*Nr,1); reshape(Ut.',Nx*Nr,1); reshape(Ux.',Nx*Nr,1); zeros(NxGC*NrGC,1); reshape(tke_ex.',Nx*Nr,1)];

tke_ex2 = tke_ex; 
tke_ex2(v_r,v_x) = tke_ex(v_r,v_x)+delta; 
tke_ex1 = tke_ex;tke_ex1(v_r,v_x) = tke_ex(v_r,v_x)-delta;
qtke2 =[reshape(Res.Ur.',Nx*Nr,1); reshape(Res.Ut.',Nx*Nr,1); reshape(Res.Ux.',Nx*Nr,1); zeros(NxGC*NrGC,1); reshape(tke_ex2.',Nx*Nr,1)];
qtke1 =[reshape(Res.Ur.',Nx*Nr,1); reshape(Res.Ut.',Nx*Nr,1); reshape(Res.Ux.',Nx*Nr,1); zeros(NxGC*NrGC,1); reshape(tke_ex1.',Nx*Nr,1)];
%%
% [J22,J23,J32,J33,~] = Jacobian_Forcing(Initial_field);
% [~,~,~,~,Fx2] = Jacobian_Forcing(qx2); [~,~,~,~,Fx1] = Jacobian_Forcing(qx1);
% [~,~,~,~,F2] = Jacobian_Forcing(qx2); 
% [~,~,~,~,F1] = Jacobian_Forcing(qx1); 
% [NuT,DNuTdUr,DNuTdUt,DNuT_dUx,S] = EV_closure_All(0,lm,Initial_field,NuT_smooth_par);
%%
% [NuT,DNuTdUr,DNuTdUt,DNuT_dUx,S] = EV_closure_All(0,lm,qt2,NuT_smooth_par);
tic
[~,F2] = JacobianSponge_AD_TKE2(qtke2,reshape(lm.',Nx*Nr,1),J_old,zeros(4*Nx*Nr+NxGC*NrGC,1),zeros(Nx*Nr,1),0); 

% [NuT,DNuTdUr,DNuTdUt,DNuT_dUx,S] = EV_closure_All(0,lm,qt1,NuT_smooth_par);
[~,F1] = JacobianSponge_AD_TKE2(qtke1,reshape(lm.',Nx*Nr,1),J_old,zeros(4*Nx*Nr+NxGC*NrGC,1),zeros(Nx*Nr,1),0);
time = toc;
%%
Dertke = (F2(3*Nx*Nr+NxGC*NrGC+1:end)-F1(3*Nx*Nr+NxGC*NrGC+1:end))/(2*delta);
%
J15 = JA(1:Nx*Nr,3*Nx*Nr+NxGC*NrGC+1:end);
J12 = JA(1:Nx*Nr,Nx*Nr+1:2*Nx*Nr);
J13 = JA(1:Nx*Nr,2*Nx*Nr+1:3*Nx*Nr);
J21 = JA(Nx*Nr+1:2*Nx*Nr,1:Nx*Nr);
J25 = JA(Nx*Nr+1:2*Nx*Nr,3*Nx*Nr+NxGC*NrGC+1:end);
J35 = JA(2*Nx*Nr+1:3*Nx*Nr,3*Nx*Nr+NxGC*NrGC+1:end);

J43 = JA(3*Nx*Nr+1:3*Nx*Nr+NxGC*NrGC,2*Nx*Nr+1:3*Nx*Nr);

J51 = JA(3*Nx*Nr+NxGC*NrGC+1:end,1:Nx*Nr);
J22 = JA(Nx*Nr+1:2*Nx*Nr,Nx*Nr+1:2*Nx*Nr);
J23 = JA(Nx*Nr+1:2*Nx*Nr,2*Nx*Nr+1:3*Nx*Nr);
J33 = JA(2*Nx*Nr+1:3*Nx*Nr,2*Nx*Nr+1:3*Nx*Nr);
J52 = JA(3*Nx*Nr+NxGC*NrGC+1:end,Nx*Nr+1:2*Nx*Nr);
J53 = JA(3*Nx*Nr+NxGC*NrGC+1:end,2*Nx*Nr+1:3*Nx*Nr);
J55 = JA(3*Nx*Nr+NxGC*NrGC+1:end,3*Nx*Nr+NxGC*NrGC+1:end);
% DFx_DUx_app = (Fx2(2*Nx*Nr+1:3*Nx*Nr)-Fx1(2*Nx*Nr+1:3*Nx*Nr))/(2*delta);
% Appr_x = (Fx2(i_glob) - Fx1(i_glob))/(2*delta);
% J33(i_glob,j_glob);

% DFx_DUx = reshape(J33((Ind+(i_r-1)*Nx),:),Nx,Nr).';
% Jx = reshape(JA(2*Nx*Nr+1:3*Nx*Nr,j_glob),Nx,Nr).';
% Dert_x = reshape(Dert(2*Nx*Nr+1:3*Nx*Nr),Nx,Nr).';
% DFt_DUt_app = reshape((Ft2(Nx*Nr+1:2*Nx*Nr)-Ft1(Nx*Nr+1:2*Nx*Nr))/(2*delta),Nx,Nr).';
% DFt_DUt = reshape(J22((Ind+(i_r-1)*Nx),:),Nx,Nr).';

% diffx = DFx_DUx_app-DFx_DUx; difft = DFt_DUt_app-DFt_DUt;
% figure
% plot(J33(:,j_glob),'*-')
% hold on
% plot(Dertke,'o-');grid on
difft=(J55(:,j_glob)-Dertke);
max(difft)
min(difft)
%
% figure
% plot(difft,'*-');
% 
% find(diffx==max(diffx))
% diff2d = reshape(diffx,Nx,Nr).';
% JA2d = reshape(J22(:,j_glob),Nx,Nr).';
% figure
% surf(Matrices.XmGLC,Matrices.RmGLC,diff2d);
clear
close all
addpath ./Functions
Tip_Speed_Ratio = 6;
Global
CL_CD %Evaluation of CL_alpha and CD_alpha

up_BC = 'Slip';
TSR = Tip_Speed_Ratio;
Nr = 40;Lr = 0.5;
BL_approx = 1;
initial_profile = 'Constant';
Turbulence_model = 'ML';
corr_block = 1.08;
folder = strcat(['./BT_Results/TSR_',num2str(TSR_theor),'/EV_parabolic/','20-Oct-2017','/',up_BC,'/',out_BC,'/']);
load(strcat([folder,'Res_',num2str(Tip_Speed_Ratio),'_',up_BC,'.mat'])); 
RANS = Res;
%%
Pre_processing %Here the input parameters are settled
up_ind = min(find(Matrices.RmGLC(:,1)<=1));
in = Nx; out = 1:Nx; %axial indeces of inlet and outlet boundaries
R = Matrices.RmGLC(:,1);
RANS.qc = [reshape(RANS.Ur.',Nx*Nr,1); reshape(RANS.Ut.',Nx*Nr,1); reshape(RANS.Ux.',Nr*Nx,1); reshape(RANS.p.',Nx*Nr,1)];
[~,~,~,~,F] = Jacobian_Forcing(RANS.qc); %Axial Forcing
Fx = reshape(F(2*Nx*Nr+1:3*Nx*Nr),Nx,Nr).';
S_S=Reynolds_stress(RANS.qc);
ML_coeff = [0 0.03 0.3142 x_T]; %ML modal coefficients
DoFs = length(ML_coeff);ML=ML_fun(ML_coeff);
[NuT,~,~,~,~] = EV_closure_All(0,ML,RANS.qc,NuT_smooth_par);
NuT = reshape(NuT,Nx,Nr).';
%%
for i=1:length(out)
    X=Matrices.XmGLC(1,out(i):in);
    Ux_inlet = RANS.Ux(up_ind:Nr,in); Ux_outlet = RANS.Ux(up_ind:Nr,out(i));
    Ux_up = RANS.Ux(up_ind,out(i):in); Ur_up = RANS.Ur(up_ind,out(i):in);
    p_in = RANS.p(up_ind:Nr,in); p_out = RANS.p(up_ind:Nr,out(i));
    
    Nu_in=1/Reynolds+NuT(up_ind:Nr,in); Nu_up=1/Reynolds+NuT(1,out(i):in); Nu_out=1/Reynolds+NuT(up_ind:Nr,out(i));
    S_in=S_S{3,3}(up_ind:Nr,in); S_up=S_S{3,1}(1,out(i):in); S_out=S_S{3,3}(up_ind:Nr,out(i));
    tau_in=(Nu_in.*S_in)'*diag(IWr(up_ind:Nr,up_ind:Nr)); 
    tau_up=max(R)*(Nu_up.*S_up)*diag(IWx(out(i):in,out(i):in)); 
    tau_out=(Nu_out.*S_out)'*diag(IWr(up_ind:Nr,up_ind:Nr));
    tau(i) = tau_up+tau_in+tau_out;
    
    IW2=kron(IWr(up_ind:Nr,up_ind:Nr),IWx(out(i):in,out(i):in));
    Force = Fx(up_ind:Nr,out(i):in);
    arg = reshape((Force).',(Nr-up_ind+1)*(Nx-i+1),1);
    Fx_AD_loc(i) = sqrt(arg).'*IW2*sqrt(arg);
    Qx_outlet_loc(i) = sqrt(Ux_outlet.^2).'*IWr(up_ind:Nr,up_ind:Nr)*sqrt(Ux_outlet.^2);
    Qx_up_loc(i) = R(up_ind)*sqrt(Ur_up.*Ux_up)*IWx(out(i):in,out(i):in)*sqrt(Ur_up.*Ux_up).';
    Qxout_loc(i) = Qx_up_loc(i) + Qx_outlet_loc(i);

    Qxin_loc(i) = sqrt(Ux_inlet.^2).'*IWr(up_ind:Nr,up_ind:Nr)*sqrt(Ux_inlet.^2);
    P(i) = (p_out-p_in).'*diag(IWr(up_ind:Nr,up_ind:Nr));
    Bal(i) = Qxout_loc(i)+Fx_AD_loc(i)+P(i)-tau(i);
end
%%
Balance = figure('units','normalized','outerposition',[0 0 0.5 1]);set(gcf,'Color','White');
plot(Matrices.XmGLC(1,:),Qxin_loc,'k*-');
hold on
plot(Matrices.XmGLC(1,:),Qx_outlet_loc,'k+-')
hold on
plot(Matrices.XmGLC(1,:),P,'kv-')
hold on
plot(Matrices.XmGLC(1,:),Fx_AD_loc,'ko-')
hold on
plot(Matrices.XmGLC(1,:),Bal,'k^-')
set(gco,'markersize',20);
leg = legend('Q_{in}','Q_{out}','Pressure','Forcing','Total');
set(leg,'location','best','fontsize',14);
Plot_options('x/D','Q',[]);
%
% savefig(Balance,strcat([folder,'Bal_x_new']));
% saveas(Balance,strcat([folder,'Bal_x_new','.png']));
%%
err_Balance = figure('units','normalized','outerposition',[0 0 0.5 1]);set(gcf,'Color','White');
plot(Matrices.XmGLC(1,:),(Bal-Qxin_loc)./Qxin_loc,'k+-');
Plot_options('x/D','Rel. error',[]);axis tight
%%
save(strcat([folder,'Forcing_TSR_6_Slip.mat']),'Fx');
save(strcat([folder,'Averaged_Forcing_TSR_6_Slip.mat']),'Fx_AD_loc');
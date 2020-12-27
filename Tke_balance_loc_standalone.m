clear
close all
addpath ./Functions
Tip_Speed_Ratio = 5;
Global
CL_CD %Evaluation of CL_alpha and CD_alpha

up_BC = 'Slip';
out_BC = 'Neumann'; %Choice of outlet boundary conditions

initial_profile = 'Constant';
Turbulence_model = 'TKE';
corr_block=1.08;
TSR = Tip_Speed_Ratio/corr_block;
Nr = 40;Lr = 0.5;
BL_approx = 1;
folder = strcat(['./BT_Results/TSR_',num2str(TSR_theor),'/',num2str(date),'/',up_BC,'/',out_BC,'/']);
load(strcat([folder,'Res','.mat'])); 
RANS = Res;
%%
Pre_processing %Here the input parameters are settled
up_ind = min(find(Matrices.RmGLC(:,1)<=0.7));
up_ind = 1;
in = Nx; out = 1:Nx; %axial indeces of inlet and outlet boundaries
R = Matrices.RmGLC(:,1);
RANS.qc = [reshape(RANS.Ur.',Nx*Nr,1); reshape(RANS.Ut.',Nx*Nr,1); reshape(RANS.Ux.',Nr*Nx,1); Matrices.MfrL*reshape(RANS.p.',Nx*Nr,1);reshape(RANS.tke.',Nr*Nx,1)];
ML_coeff = [0 0.03 0.3142 x_T]; %ML modal coefficients
DoFs = length(ML_coeff);
lm = ML_fun(ML_coeff);

[P,~,~,~,~] = Tke_production(reshape(lm.',Nx*Nr,1),RANS.qc); P2d = reshape(P,Nx,Nr).';%Production of tke
[T,~] = Tke_transportation(RANS.qc(3*Nx*Nr+NxGC*NrGC+1:end),reshape(lm.',Nx*Nr,1)); Transport2d = reshape(T,Nx,Nr).';%Transportation of tke
[eps_tke,~] = Tke_dissipation(RANS.qc(3*Nx*Nr+NxGC*NrGC+1:end)); Diss2d = reshape(eps_tke,Nx,Nr).';
%%
Tke_flow_in = zeros(1,Nx);Tke_flow_out = Tke_flow_in; Tke_flow_up=Tke_flow_in;Prod = Tke_flow_in;
Transp = Tke_flow_in;Diss = Transp; Bal=Transp;
for i=1:length(out)
    tke_inlet = RANS.tke(up_ind:Nr,in); tke_outlet = RANS.tke(up_ind:Nr,out(i)); tke_up = RANS.tke(up_ind,out(i):in);
    
    Ux_outlet = RANS.Ux(up_ind:Nr,out(i));
    Ur_up = RANS.Ur(up_ind,out(i):in);
    Ux_inlet = RANS.Ux(up_ind:Nr,in);
    
    Tke_flow_in(i) = sqrt(Ux_inlet.*tke_inlet).'*IWr(up_ind:Nr,up_ind:Nr)*sqrt(Ux_inlet.*tke_inlet);
    Tke_flow_out(i) = sqrt(Ux_outlet.*tke_outlet).'*IWr(up_ind:Nr,up_ind:Nr)*sqrt(Ux_outlet.*tke_outlet);
    Tke_flow_up(i) = R(up_ind)*sqrt(Ur_up.*tke_up)*IWx(out(i):in,out(i):in)*sqrt(Ur_up.*tke_up).';

    IW2=kron(IWr(up_ind:Nr,up_ind:Nr),IWx(out(i):in,out(i):in));

    Production = P2d(up_ind:Nr,out(i):in);
    arg_p = reshape((Production).',(Nr-up_ind+1)*(Nx-i+1),1);
    Prod(i) = sqrt(arg_p).'*IW2*sqrt(arg_p);
    
    Transport = Transport2d(up_ind:Nr,out(i):in);
    arg_t = reshape((Transport).',(Nr-up_ind+1)*(Nx-i+1),1);
    Transp(i) = sqrt(arg_t).'*IW2*sqrt(arg_t);
    
    Dissipation = Diss2d(up_ind:Nr,out(i):in);
    arg_d = reshape((Dissipation).',(Nr-up_ind+1)*(Nx-i+1),1);
    Diss(i) = sqrt(arg_d).'*IW2*sqrt(arg_d);
    
    Bal(i) =  Tke_flow_out(i) + Tke_flow_up(i) -  Prod(i) + Diss(i) - Transp(i);
end
%%
Balance = figure('units','normalized','outerposition',[0 0 0.5 1]);set(gcf,'Color','White');
plot(Matrices.XmGLC(1,:),Tke_flow_in,'k*-');
hold on
plot(Matrices.XmGLC(1,:),Tke_flow_out,'k+-')
hold on
plot(Matrices.XmGLC(1,:),Prod,'kv-')
hold on
plot(Matrices.XmGLC(1,:),Transp,'ko-')
hold on
plot(Matrices.XmGLC(1,:),Diss,'k.-','markersize',15)
hold on
plot(Matrices.XmGLC(1,:),Bal,'k^-')
set(gco,'markersize',20);
leg = legend('Tke_{in}','Tke_{out}','Production','Transport','Dissipation','Total');
set(leg,'location','best','fontsize',14);
Plot_options('x/D','Tke',[]);
savefig(Balance,strcat([folder,'Bal_tke']));
saveas(Balance,strcat([folder,'Bal_tke','.png']));
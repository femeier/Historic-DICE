clear all
%% Notation
%  - Damage functions: Nordhaus (N), Howard & Sterner (HS), Weizman (W)
%  - Optimization: fix (p) abatement optimized (o), abatement and investment optimized (oo)
%  - Columns 1-6: Column 1: DICE, Columns 2-5: SSPs
%  - Overview:
%       1. Import Data
%       2. Assign Data (to the 6 scenarios)
%       3. Parameters
%       4. Calculations: Accumulated Investments
%       5. Calculations: Accumulated Disinvestments
%       6. Calculations: Accumulated Net Investments
%       7. Plots
%% 1. IMPORT DATA
cd C:\Users\rickels\Dropbox\Historic_DICE_Dropbox\Optimization\SCC_results_1950_2018_kappa_FH
% Consumption 
Consumption_CHS=readmatrix('CHS.xlsx'); 
Consumption_CHSo=readmatrix('CHSo.xlsx'); 
Consumption_CHSIo=readmatrix('CHSoI.xlsx'); 
Consumption_CHSoo=readmatrix('CHSoo.xlsx'); 
Consumption_CN=readmatrix('CN.xlsx'); 
Consumption_CNo=readmatrix('CNo.xlsx'); 
Consumption_CNIo=readmatrix('CNoI.xlsx'); 
Consumption_CNoo=readmatrix('CNoo.xlsx'); 
Consumption_CW=readmatrix('CW.xlsx'); 
Consumption_CWo=readmatrix('CWo.xlsx'); 
Consumption_CWIo=readmatrix('CWoI.xlsx'); 
Consumption_CWoo=readmatrix('CWoo.xlsx'); 

% Emissions 
Emissions_EHS=readmatrix('EHS.xlsx'); 
Emissions_EHSo=readmatrix('EHSo.xlsx'); 
Emissions_EHSIo=readmatrix('EHSoI.xlsx'); 
Emissions_EHSoo=readmatrix('EHSoo.xlsx'); 
Emissions_EN=readmatrix('EN.xlsx'); 
Emissions_ENo=readmatrix('ENo.xlsx'); 
Emissions_ENIo=readmatrix('ENoI.xlsx'); 
Emissions_ENoo=readmatrix('ENoo.xlsx'); 
Emissions_EW=readmatrix('EW.xlsx'); 
Emissions_EWo=readmatrix('EWo.xlsx'); 
Emissions_EWIo=readmatrix('EWoI.xlsx'); 
Emissions_EWoo=readmatrix('EWoo.xlsx'); 
% Investments 

Investment_InvHS=readmatrix('InvHS.xlsx');
Investment_InvHSo=readmatrix('InvHSo.xlsx');
Investment_InvHSIo=readmatrix('InvHSoI.xlsx');
Investment_InvHSoo=readmatrix('InvHSoo.xlsx');
Investment_InvN=readmatrix('InvN.xlsx');
Investment_InvNo=readmatrix('InvNo.xlsx');
Investment_InvNIo=readmatrix('InvNoI.xlsx');
Investment_InvNoo=readmatrix('InvNoo.xlsx');
Investment_InvW=readmatrix('InvW.xlsx');
Investment_InvWo=readmatrix('InvWo.xlsx');
Investment_InvWIo=readmatrix('InvWoI.xlsx');
Investment_InvWoo=readmatrix('InvWoo.xlsx');

% SCC 
SCC_HS=readmatrix('sccHSp.xlsx'); 
SCC_HSo=readmatrix('sccHSpo.xlsx'); 
SCC_HSIo=readmatrix('sccHSpoI.xlsx'); 
SCC_HSoo=readmatrix('sccHSpoo.xlsx'); 
SCC_N=readmatrix('sccNp.xlsx'); 
SCC_No=readmatrix('sccNpo.xlsx'); 
SCC_NIo=readmatrix('sccNpoI.xlsx'); 
SCC_Noo=readmatrix('sccNpoo.xlsx'); 
SCC_W=readmatrix('sccWp.xlsx'); 
SCC_Wo=readmatrix('sccWpo.xlsx'); 
SCC_WIo=readmatrix('sccWpoI.xlsx'); 
SCC_Woo=readmatrix('sccWpoo.xlsx'); 

% Population external
Population_external=readmatrix('L_external.csv'); 

%% 2. ASSIGN DATA (to the 6 scenarios)
a1=2; %1950
a2=70; %2019
% Consumption 
for i=1:6
CHSp(:,i)=Consumption_CHS(a1:a2,i+1);
CHSo(:,i)=Consumption_CHSo(a1:a2,i+1);
CHSIo(:,i)=Consumption_CHSIo(a1:a2,i+1);
CHSoo(:,i)=Consumption_CHSoo(a1:a2,i+1);
CNp(:,i)=Consumption_CN(a1:a2,i+1);
CNo(:,i)=Consumption_CNo(a1:a2,i+1);
CNIo(:,i)=Consumption_CNIo(a1:a2,i+1);
CNoo(:,i)=Consumption_CNoo(a1:a2,i+1);
CWp(:,i)=Consumption_CW(a1:a2,i+1);
CWo(:,i)=Consumption_CWo(a1:a2,i+1);
CWIo(:,i)=Consumption_CWIo(a1:a2,i+1);
CWoo(:,i)=Consumption_CWoo(a1:a2,i+1);
end
% Emissions 
for i=1:6
EHSp(:,i)=Emissions_EHS(a1:a2,i+1);
EHSo(:,i)=Emissions_EHSo(a1:a2,i+1);
EHSIo(:,i)=Emissions_EHSIo(a1:a2,i+1);
EHSoo(:,i)=Emissions_EHSoo(a1:a2,i+1);
ENp(:,i)=Emissions_EN(a1:a2,i+1);
ENo(:,i)=Emissions_ENo(a1:a2,i+1);
ENIo(:,i)=Emissions_ENIo(a1:a2,i+1);
ENoo(:,i)=Emissions_ENoo(a1:a2,i+1);
EWp(:,i)=Emissions_EW(a1:a2,i+1);
EWo(:,i)=Emissions_EWo(a1:a2,i+1);
EWIo(:,i)=Emissions_EWIo(a1:a2,i+1);
EWoo(:,i)=Emissions_EWoo(a1:a2,i+1);
end
% Investments 
for i=1:6
InvHSp(:,i)=Investment_InvHS(a1:a2,i+1);
InvHSo(:,i)=Investment_InvHSo(a1:a2,i+1);
InvHSIo(:,i)=Investment_InvHSIo(a1:a2,i+1);
InvHSoo(:,i)=Investment_InvHSoo(a1:a2,i+1);
InvNp(:,i)=Investment_InvN(a1:a2,i+1);
InvNo(:,i)=Investment_InvNo(a1:a2,i+1);
InvNIo(:,i)=Investment_InvNIo(a1:a2,i+1);
InvNoo(:,i)=Investment_InvNoo(a1:a2,i+1);
InvWp(:,i)=Investment_InvW(a1:a2,i+1);
InvWo(:,i)=Investment_InvWo(a1:a2,i+1);
InvWIo(:,i)=Investment_InvWIo(a1:a2,i+1);
InvWoo(:,i)=Investment_InvWoo(a1:a2,i+1);
end
% SCC 
for i=1:6
SCCHSp(:,i)=SCC_HS(a1:a2,i+1);
SCCHSo(:,i)=SCC_HSo(a1:a2,i+1);
SCCHSIo(:,i)=SCC_HSIo(a1:a2,i+1);
SCCHSoo(:,i)=SCC_HSoo(a1:a2,i+1);
SCCNp(:,i)=SCC_N(a1:a2,i+1);
SCCNo(:,i)=SCC_No(a1:a2,i+1);
SCCNIo(:,i)=SCC_NIo(a1:a2,i+1);
SCCNoo(:,i)=SCC_Noo(a1:a2,i+1);
SCCWp(:,i)=SCC_W(a1:a2,i+1);
SCCWo(:,i)=SCC_Wo(a1:a2,i+1);
SCCWIo(:,i)=SCC_WIo(a1:a2,i+1);
SCCWoo(:,i)=SCC_Woo(a1:a2,i+1);
end
% Population
for i=1:6
Pop(:,i)=Population_external(a1:a2,i+1);
end
%% 3. PARAMETERS
rho_RMQ=0.01950802;
eta_RMQ=1.45;
for j=a1-1:a2-1
R_RMQ(j)=1/(1+rho_RMQ)^(j-1);
end
%% 4. CALCULATIONS: Accumulated Investment
for t=1:a2-1 %1950-2018
    for i=1:6 
% Howard & Sterner (p)
Factor_current_HSp(t,i)=(1000^(1-eta_RMQ)*(CHSp(t,i)/Pop(t,i))^(-eta_RMQ)*R_RMQ(t))/(1000^(1-eta_RMQ)*(CHSp(a2-1,i)/Pop(a2-1,i))^(-eta_RMQ)*R_RMQ(a2-1));
Current_inv_HSp(t,i)=Factor_current_HSp(t,i)*InvHSp(t,i);
Current_inv_acc_HSp=sum(Current_inv_HSp);
% Howard & Sterner (o)
Factor_current_HSo(t,i)=(1000^(1-eta_RMQ)*(CHSo(t,i)/Pop(t,i))^(-eta_RMQ)*R_RMQ(t))/(1000^(1-eta_RMQ)*(CHSo(a2-1,i)/Pop(a2-1,i))^(-eta_RMQ)*R_RMQ(a2-1));
Current_inv_HSo(t,i)=Factor_current_HSo(t,i)*InvHSo(t,i);
Current_inv_acc_HSo=sum(Current_inv_HSo);
% Howard & Sterner (Io)
Factor_current_HSIo(t,i)=(1000^(1-eta_RMQ)*(CHSIo(t,i)/Pop(t,i))^(-eta_RMQ)*R_RMQ(t))/(1000^(1-eta_RMQ)*(CHSIo(a2-1,i)/Pop(a2-1,i))^(-eta_RMQ)*R_RMQ(a2-1));
Current_inv_HSIo(t,i)=Factor_current_HSIo(t,i)*InvHSIo(t,i);
Current_inv_acc_HSIo=sum(Current_inv_HSIo);
% Howard & Sterner (oo)
Factor_current_HSoo(t,i)=(1000^(1-eta_RMQ)*(CHSoo(t,i)/Pop(t,i))^(-eta_RMQ)*R_RMQ(t))/(1000^(1-eta_RMQ)*(CHSoo(a2-1,i)/Pop(a2-1,i))^(-eta_RMQ)*R_RMQ(a2-1));
Current_inv_HSoo(t,i)=Factor_current_HSoo(t,i)*InvHSoo(t,i);
Current_inv_acc_HSoo=sum(Current_inv_HSoo);
% Nordhaus (p)
Factor_current_Np(t,i)=(1000^(1-eta_RMQ)*(CNp(t,i)/Pop(t,i))^(-eta_RMQ)*R_RMQ(t))/(1000^(1-eta_RMQ)*(CNp(a2-1,i)/Pop(a2-1,i))^(-eta_RMQ)*R_RMQ(a2-1));
Current_inv_Np(t,i)=Factor_current_Np(t,i)*InvNp(t,i);
Current_inv_acc_Np=sum(Current_inv_Np);
% Nordhaus (o)
Factor_current_No(t,i)=(1000^(1-eta_RMQ)*(CNo(t,i)/Pop(t,i))^(-eta_RMQ)*R_RMQ(t))/(1000^(1-eta_RMQ)*(CNo(a2-1,i)/Pop(a2-1,i))^(-eta_RMQ)*R_RMQ(a2-1));
Current_inv_No(t,i)=Factor_current_No(t,i)*InvNo(t,i);
Current_inv_acc_No=sum(Current_inv_No);
% Nordhaus (Io)
Factor_current_NIo(t,i)=(1000^(1-eta_RMQ)*(CNIo(t,i)/Pop(t,i))^(-eta_RMQ)*R_RMQ(t))/(1000^(1-eta_RMQ)*(CNIo(a2-1,i)/Pop(a2-1,i))^(-eta_RMQ)*R_RMQ(a2-1));
Current_inv_NIo(t,i)=Factor_current_NIo(t,i)*InvNIo(t,i);
Current_inv_acc_NIo=sum(Current_inv_NIo);
% Nordhaus (oo)
Factor_current_Noo(t,i)=(1000^(1-eta_RMQ)*(CNoo(t,i)/Pop(t,i))^(-eta_RMQ)*R_RMQ(t))/(1000^(1-eta_RMQ)*(CNoo(a2-1,i)/Pop(a2-1,i))^(-eta_RMQ)*R_RMQ(a2-1));
Current_inv_Noo(t,i)=Factor_current_Noo(t,i)*InvNoo(t,i);
Current_inv_acc_Noo=sum(Current_inv_Noo);
% Weitzman (p)
Factor_current_Wp(t,i)=(1000^(1-eta_RMQ)*(CWp(t,i)/Pop(t,i))^(-eta_RMQ)*R_RMQ(t))/(1000^(1-eta_RMQ)*(CWp(a2-1,i)/Pop(a2-1,i))^(-eta_RMQ)*R_RMQ(a2-1));
Current_inv_Wp(t,i)=Factor_current_Wp(t,i)*InvWp(t,i);
Current_inv_acc_Wp=sum(Current_inv_Wp);
% Weitzman (o)
Factor_current_Wo(t,i)=(1000^(1-eta_RMQ)*(CWo(t,i)/Pop(t,i))^(-eta_RMQ)*R_RMQ(t))/(1000^(1-eta_RMQ)*(CWo(a2-1,i)/Pop(a2-1,i))^(-eta_RMQ)*R_RMQ(a2-1));
Current_inv_Wo(t,i)=Factor_current_Wo(t,i)*InvWo(t,i);
Current_inv_acc_Wo=sum(Current_inv_Wo);
% Weitzman (Io)
Factor_current_WIo(t,i)=(1000^(1-eta_RMQ)*(CWIo(t,i)/Pop(t,i))^(-eta_RMQ)*R_RMQ(t))/(1000^(1-eta_RMQ)*(CWIo(a2-1,i)/Pop(a2-1,i))^(-eta_RMQ)*R_RMQ(a2-1));
Current_inv_WIo(t,i)=Factor_current_WIo(t,i)*InvWIo(t,i);
Current_inv_acc_WIo=sum(Current_inv_WIo);
% Weitzman (oo)
Factor_current_Woo(t,i)=(1000^(1-eta_RMQ)*(CWoo(t,i)/Pop(t,i))^(-eta_RMQ)*R_RMQ(t))/(1000^(1-eta_RMQ)*(CWoo(a2-1,i)/Pop(a2-1,i))^(-eta_RMQ)*R_RMQ(a2-1));
Current_inv_Woo(t,i)=Factor_current_Woo(t,i)*InvWoo(t,i);
Current_inv_acc_Woo=sum(Current_inv_Woo);
    end
end
% Mean and SD
%fix
Current_inv_acc_fix=[Current_inv_acc_HSp, Current_inv_acc_Np, Current_inv_acc_Wp];
Current_inv_acc_fix_mean=mean(Current_inv_acc_fix);
Current_inv_acc_fix_SD=std(Current_inv_acc_fix);
%opt
Current_inv_acc_opt=[Current_inv_acc_HSo, Current_inv_acc_No, Current_inv_acc_Wo];
Current_inv_acc_opt_mean=mean(Current_inv_acc_opt);
Current_inv_acc_opt_SD=std(Current_inv_acc_opt);
%Iopt
Current_inv_acc_Iopt=[Current_inv_acc_HSIo, Current_inv_acc_NIo, Current_inv_acc_WIo];
Current_inv_acc_Iopt_mean=mean(Current_inv_acc_Iopt);
Current_inv_acc_Iopt_SD=std(Current_inv_acc_Iopt);
%optopt
Current_inv_acc_optopt=[Current_inv_acc_HSoo, Current_inv_acc_Noo, Current_inv_acc_Woo];
Current_inv_acc_optopt_mean=mean(Current_inv_acc_optopt);
Current_inv_acc_optopt_SD=std(Current_inv_acc_optopt);

%% 5. CALCULATIONS: Accumulated Disinvestment
for t=1:a2-1 %1950-2018
    for i=1:6 
% Howard & Sterner (p)
Disinvestment_HSp(t,i)=EHSp(t)*3.664*SCCHSp(t,i)/1000;
Disinvestment_acc_HSp=sum(Disinvestment_HSp);
% Howard & Sterner (o)
Disinvestment_HSo(t,i)=EHSo(t)*3.664*SCCHSo(t,i)/1000;
Disinvestment_acc_HSo=sum(Disinvestment_HSo);
% Howard & Sterner (Io)
Disinvestment_HSIo(t,i)=EHSIo(t)*3.664*SCCHSIo(t,i)/1000;
Disinvestment_acc_HSIo=sum(Disinvestment_HSIo);
% Howard & Sterner (oo)
Disinvestment_HSoo(t,i)=EHSoo(t)*3.664*SCCHSoo(t,i)/1000;
Disinvestment_acc_HSoo=sum(Disinvestment_HSoo);
% Nordhaus (p)
Disinvestment_Np(t,i)=ENp(t)*3.664*SCCNp(t,i)/1000;
Disinvestment_acc_Np=sum(Disinvestment_Np);
% Nordhaus (o)
Disinvestment_No(t,i)=ENo(t)*3.664*SCCNo(t,i)/1000;
Disinvestment_acc_No=sum(Disinvestment_No);
% Nordhaus (Io)
Disinvestment_NIo(t,i)=ENIo(t)*3.664*SCCNIo(t,i)/1000;
Disinvestment_acc_NIo=sum(Disinvestment_NIo);
% Nordhaus (oo)
Disinvestment_Noo(t,i)=ENoo(t)*3.664*SCCNoo(t,i)/1000;
Disinvestment_acc_Noo=sum(Disinvestment_Noo);
% Weitzman (p)
Disinvestment_Wp(t,i)=EWp(t)*3.664*SCCWp(t,i)/1000;
Disinvestment_acc_Wp=sum(Disinvestment_Wp);
% Weitzman (o)
Disinvestment_Wo(t,i)=EWo(t)*3.664*SCCWo(t,i)/1000;
Disinvestment_acc_Wo=sum(Disinvestment_Wo);
% Weitzman (o)
Disinvestment_WIo(t,i)=EWIo(t)*3.664*SCCWIo(t,i)/1000;
Disinvestment_acc_WIo=sum(Disinvestment_WIo);
% Weitzman (oo)
Disinvestment_Woo(t,i)=EWoo(t)*3.664*SCCWoo(t,i)/1000;
Disinvestment_acc_Woo=sum(Disinvestment_Woo);
    end
end
% Mean and SD
%fix
Disinvestment_acc_fix=[Disinvestment_acc_HSp, Disinvestment_acc_Np, Disinvestment_acc_Wp];
Disinvestment_acc_fix_mean=mean(Disinvestment_acc_fix);
Disinvestment_acc_fix_SD=std(Disinvestment_acc_fix);
%opt
Disinvestment_acc_opt=[Disinvestment_acc_HSo, Disinvestment_acc_No, Disinvestment_acc_Wo];
Disinvestment_acc_opt_mean=mean(Disinvestment_acc_opt);
Disinvestment_acc_opt_SD=std(Disinvestment_acc_opt);
%Iopt
Disinvestment_acc_Iopt=[Disinvestment_acc_HSIo, Disinvestment_acc_NIo, Disinvestment_acc_WIo];
Disinvestment_acc_Iopt_mean=mean(Disinvestment_acc_Iopt);
Disinvestment_acc_Iopt_SD=std(Disinvestment_acc_Iopt);
%Iopt
Disinvestment_acc_optopt=[Disinvestment_acc_HSoo, Disinvestment_acc_Noo, Disinvestment_acc_Woo];
Disinvestment_acc_optopt_mean=mean(Disinvestment_acc_optopt);
Disinvestment_acc_optopt_SD=std(Disinvestment_acc_optopt);

%% 5. CALCULATIONS: Share of Disinvestment
sharefix=Disinvestment_acc_fix_mean/Current_inv_acc_fix_mean;
shareopt=Disinvestment_acc_opt_mean/Current_inv_acc_opt_mean;
shareIopt=Disinvestment_acc_Iopt_mean/Current_inv_acc_Iopt_mean;
shareoptopt=Disinvestment_acc_optopt_mean/Current_inv_acc_optopt_mean;
%for t=1:a2-1 %1950-2018
 %   for i=1:6 
% Howard & Sterner (p)
Disinvestment_share_HSp=Disinvestment_acc_HSp/Current_inv_acc_HSp;
% Howard & Sterner (o)
%Disinvestment_share_HSo=Disinvestment_acc_HSo/Current_inv_acc_HSo;
% Howard & Sterner (oo)
%Disinvestment_share_HSoo=Disinvestment_acc_HSoo/Current_inv_acc_HSoo;
% Nordhaus (p)
Disinvestment_share_Np=Disinvestment_acc_Np/Current_inv_acc_Np;
% Nordhaus (o)
%Disinvestment_share_No=Disinvestment_acc_No/Current_inv_acc_No;
% Nordhaus (oo)
%Disinvestment_share_Noo=Disinvestment_acc_Noo/Current_inv_acc_Noo;
% Weitzman (p)
Disinvestment_share_Wp=Disinvestment_acc_Wp/Current_inv_acc_Wp;
% Weitzman (o)
%Disinvestment_share_Wo=Disinvestment_acc_Wo/Current_inv_acc_Wo;
% Weitzman (oo)
%Disinvestment_share_Woo=Disinvestment_acc_Woo/Current_inv_acc_Woo;
  %end
  %end

% Mean and SD
%Disinvestment_share = [Disinvestment_share_HSp Disinvestment_share_HSo Disinvestment_share_HSoo Disinvestment_share_Np Disinvestment_share_No Disinvestment_share_Noo Disinvestment_share_Wp Disinvestment_share_Wo Disinvestment_share_Woo];
%Disinvestment_share_mean = mean(Disinvestment_share);
%Disinvestment_share_SD = std(Disinvestment_share);
% Mean and SD Only fix
Disinvestment_share_fix_table= Disinvestment_acc_fix./Current_inv_acc_fix;
Disinvestment_share_fix = [Disinvestment_share_HSp Disinvestment_share_Np Disinvestment_share_Wp];
Disinvestment_share_fix_mean = mean(Disinvestment_share_fix);
Disinvestment_share_fix_SD = std(Disinvestment_share_fix);
%writematrix(Disinvestment_share_fix_table,'Disinvestment_share.csv')
%% 7. PLOTS
%% 8. Timepaths
for n=1:69 
    Current_inv_mean(n)=mean(horzcat(Current_inv_HSp(n,:), Current_inv_Np(n,:), Current_inv_Wp(n,:))); 
    Current_inv_std(n)=std(horzcat(Current_inv_HSp(n,:), Current_inv_Np(n,:), Current_inv_Wp(n,:)));
    Current_inv_var(n)=var(horzcat(Current_inv_HSp(n,:), Current_inv_Np(n,:), Current_inv_Wp(n,:)));
    Current_disinv_mean(n)=mean(horzcat(Disinvestment_HSp(n,:), Disinvestment_Np(n,:), Disinvestment_Wp(n,:)));
    Current_disinv_std(n)=std(horzcat(Disinvestment_HSp(n,:), Disinvestment_Np(n,:), Disinvestment_Wp(n,:)));
    Current_disinv_var(n)=var(horzcat(Disinvestment_HSp(n,:), Disinvestment_Np(n,:), Disinvestment_Wp(n,:)));
end

for n=1:69 
    Current_invo_mean(n)=mean(horzcat(Current_inv_HSo(n,:), Current_inv_No(n,:), Current_inv_Wo(n,:))); 
    Current_invo_std(n)=std(horzcat(Current_inv_HSo(n,:), Current_inv_No(n,:), Current_inv_Wo(n,:)));
    Current_invo_var(n)=var(horzcat(Current_inv_HSo(n,:), Current_inv_No(n,:), Current_inv_Wo(n,:)));
    Current_disinvo_mean(n)=mean(horzcat(Disinvestment_HSo(n,:), Disinvestment_No(n,:), Disinvestment_Wo(n,:)));
    Current_disinvo_std(n)=std(horzcat(Disinvestment_HSo(n,:), Disinvestment_No(n,:), Disinvestment_Wo(n,:)));
    Current_disinvo_var(n)=var(horzcat(Disinvestment_HSo(n,:), Disinvestment_No(n,:), Disinvestment_Wo(n,:)));
end

for n=1:69
    current_inv_cum(n)=sum(Current_inv_mean(n:69));
    current_inv_cum_std(n)=sqrt(sum(Current_inv_var(n:69)));
    current_disinv_cum(n)=sum(Current_disinv_mean(n:69));
    current_disinv_cum_std(n)=sqrt(sum(Current_disinv_var(n:69)));
end

for n=1:69
    current_invo_cum(n)=sum(Current_invo_mean(n:69));
    current_invo_cum_std(n)=sqrt(sum(Current_invo_var(n:69)));
    current_disinvo_cum(n)=sum(Current_disinvo_mean(n:69));
    current_disinvo_cum_std(n)=sqrt(sum(Current_disinvo_var(n:69)));
end

cum_share_time=100*(current_disinv_cum./current_inv_cum);
cumo_share_time=100*current_disinvo_cum./current_invo_cum;
temp1=(current_disinv_cum_std./current_disinv_cum).^2;
temp1o=(current_disinvo_cum_std./current_disinvo_cum).^2;
temp2=(current_inv_cum_std./current_inv_cum).^2;
temp2o=(current_invo_cum_std./current_invo_cum).^2;
cum_share_time_std=sqrt(temp1+temp2).*cum_share_time;
cumo_share_time_std=sqrt(temp1o+temp2o).*cumo_share_time;
annual_share_time=100*Current_disinv_mean./Current_inv_mean;
temp3=(Current_disinv_std./Current_disinv_mean).^2;
temp4=(Current_inv_std./Current_inv_mean).^2;
annual_share_time_std=sqrt(temp3+temp4).*annual_share_time;
percapburden=10^6*current_disinv_cum./mean(Population_external(70,2:7));
percapburden_std=10^6*current_disinv_cum_std./mean(Population_external(70,2:7));


%cum_share_time=100*current_disinv_cum./current_inv_cum;
%temp1=(current_disinv_cum_std./current_disinv_cum).^2;
%temp2=(current_inv_cum_std./current_inv_cum).^2;
%cum_share_time_std=sqrt(temp1+temp2).*cum_share_time;
%annual_share_time=100*Current_disinv_mean./Current_inv_mean;
%temp3=(Current_disinv_std./Current_disinv_mean).^2;
%temp4=(Current_inv_std./Current_inv_mean).^2;
%annual_share_time_std=sqrt(temp3+temp4).*annual_share_time;
%percapburden=10^6*current_disinv_cum./mean(Population_external(70,2:7));
%percapburden_std=10^6*current_disinv_cum_std./mean(Population_external(70,2:7));

%cum_oshare_time=100*current_disinvo_cum./current_invo_cum;
%temp1o=(current_disinvo_cum_std./current_disinvo_cum).^2;
%temp2o=(current_invo_cum_std./current_invo_cum).^2;
%cum_oshare_time_std=sqrt(temp1o+temp2o).*cum_oshare_time;

cd C:\Users\rickels\Dropbox\Historic_DICE_Dropbox\IW_calculations\Global

output1=[cum_share_time; cum_share_time_std; annual_share_time; annual_share_time_std; percapburden; percapburden_std];
output2=[cum_share_time; cum_share_time_std; cumo_share_time; cumo_share_time_std;annual_share_time; annual_share_time_std; percapburden; percapburden_std];
writematrix(output1,'IW_output_global.csv')
colNames="y_"+(1950:2018);
rowNames={'cumshare','cumshare_std','cumoptshare','cumoptshare_std','annualshare','annualshare_std','percap','percap_std'};
rfix = array2table(output2,'RowNames',rowNames,'VariableNames',colNames);
writetable(rfix,'IW_table_global.csv','WriteRowNames',true)
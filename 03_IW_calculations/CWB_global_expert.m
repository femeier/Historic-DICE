clear all
%% 1. IMPORT DATA
cd C:\..\SCC_results_2022_constant_Drupp % Choose the file path 

% Capital 
Capital_CHS=readmatrix('KHS.xlsx'); 
Capital_CN=readmatrix('KN.xlsx'); 
Capital_CW=readmatrix('KW.xlsx'); 


% Consumption 
Consumption_CHS=readmatrix('CHS.xlsx'); 
Consumption_CN=readmatrix('CN.xlsx'); 
Consumption_CW=readmatrix('CW.xlsx'); 

% Emissions 
Emissions_EHS=readmatrix('EHS.xlsx'); 
Emissions_EN=readmatrix('EN.xlsx'); 
Emissions_EW=readmatrix('EW.xlsx'); 

% Investments 
Investment_InvHS=readmatrix('InvHS.xlsx');
Investment_InvN=readmatrix('InvN.xlsx');
Investment_InvW=readmatrix('InvW.xlsx');

% SCC 
SCC_HS=readmatrix('sccHSp.xlsx'); 
SCC_N=readmatrix('sccNp.xlsx'); 
SCC_W=readmatrix('sccWp.xlsx'); 

% Population external
Population_external=readmatrix('L_external.csv'); 

%% 2. ASSIGN DATA (to the 6 scenarios)
a1=2; %1950
a2=70; %2018
% NetInvestment
for i=1:6
NinvHSp(:,i)=Capital_CHS(a1+1:a2+1,i+1)-Capital_CHS(a1:a2,i+1);
NinvNp(:,i)=Capital_CN(a1+1:a2+1,i+1)-Capital_CN(a1:a2,i+1);
NinvWp(:,i)=Capital_CW(a1+1:a2+1,i+1)-Capital_CW(a1:a2,i+1);
end

% Consumption 
for i=1:6
CHSp(:,i)=Consumption_CHS(a1:a2,i+1);
CNp(:,i)=Consumption_CN(a1:a2,i+1);
CWp(:,i)=Consumption_CW(a1:a2,i+1);
end
% Emissions 
for i=1:6
EHSp(:,i)=Emissions_EHS(a1:a2,i+1);
ENp(:,i)=Emissions_EN(a1:a2,i+1);
EWp(:,i)=Emissions_EW(a1:a2,i+1);
end
% Investments 
for i=1:6
InvHSp(:,i)=Investment_InvHS(a1:a2,i+1);
InvNp(:,i)=Investment_InvN(a1:a2,i+1);
InvWp(:,i)=Investment_InvW(a1:a2,i+1);
end
% SCC 
for i=1:6
SCCHSp(:,i)=SCC_HS(a1:a2,i+1);
SCCNp(:,i)=SCC_N(a1:a2,i+1);
SCCWp(:,i)=SCC_W(a1:a2,i+1);
end

SCC=horzcat(SCCHSp,SCCNp,SCCWp);
SCCm=mean(SCC,2);
SCCstd=transpose(std(transpose(SCC)));

% Population
for i=1:6
Pop(:,i)=Population_external(a1:a2,i+1);
Popmean=mean(Pop,2)
end
%% 3. PARAMETERS
rho_RMQ=1.1/100; %0.01950802;
eta_RMQ=1.35; %1.45;
for j=a1-1:a2-1
R_RMQ(j)=1/(1+rho_RMQ)^(j-1);
end
%% 4. CALCULATIONS: Accumulated Investment
j1=1; %1950
j2=25; %1974
j3=26; %1975
j4=50; %1999
j5=51; %2000
j6=69; %2018

%Note that there is no statistical variation between the properties below,
%taking means is just to eliminate numerical variation; in the following calculations, only the means will be used. 

Inv=horzcat(InvHSp,InvNp,InvWp);
Inv_mean=mean(Inv,2);
Inv_std=std(transpose(Inv));

NInv=horzcat(NinvHSp,NinvNp,NinvWp);
NInv_mean=mean(NInv,2);
NInv_std=std(transpose(NInv));

Con=horzcat(CHSp,CNp,CWp);
Con_mean=mean(Con,2);
Con_std=std(transpose(Inv));

CO2=horzcat(EHSp,ENp,EWp);
CO2_mean=mean(CO2*3.664,2);
CO2_std=std(transpose(CO2*3.664));

for t=1:a2-1
Factor(t)=(1000^(1-eta_RMQ)*(Con_mean(t)/Popmean(t))^(-eta_RMQ)*R_RMQ(t))/(1000^(1-eta_RMQ)*(Con_mean(69)/Popmean(69))^(-eta_RMQ)*R_RMQ(69));
end

Present_inv_mean=Inv_mean.*transpose(Factor);
Present_Ninv_mean=NInv_mean.*transpose(Factor);
Present_inv_acc=sum(Present_inv_mean);
Present_inv_epo1=sum(Present_inv_mean(j1:j2));
Present_inv_epo2=sum(Present_inv_mean(j3:j4));
Present_inv_epo3=sum(Present_inv_mean(j5:j6));

Present_Ninv_acc=sum(Present_Ninv_mean);
Present_Ninv_epo1=sum(Present_Ninv_mean(j1:j2));
Present_Ninv_epo2=sum(Present_Ninv_mean(j3:j4));
Present_Ninv_epo3=sum(Present_Ninv_mean(j5:j6));

% Calculate CWB
CCO2_mean=SCCm.*CO2_mean/1000;
CCO2_std=SCCstd.*CO2_mean/1000;
CCO2_acc_mean=sum(CCO2_mean);
CCO2_epo1_mean=sum(CCO2_mean(j1:j2));
CCO2_epo2_mean=sum(CCO2_mean(j3:j4));
CCO2_epo3_mean=sum(CCO2_mean(j5:j6));

CCO2_acc_std=(sum(CCO2_std.^2)).^0.5;
CCO2_epo1_std=(sum(CCO2_mean(j1:j2).^2)).^0.5;
CCO2_epo2_std=(sum(CCO2_mean(j3:j4).^2)).^0.5;
CCO2_epo3_std=(sum(CCO2_mean(j5:j6).^2)).^0.5;


%Calculate CWB Share
share_meanCWBCO2gross=CCO2_acc_mean/Present_inv_acc;
share_meanCWBCO2grossepo1=CCO2_epo1_mean/Present_inv_epo1;
share_meanCWBCO2grossepo2=CCO2_epo2_mean/Present_inv_epo2;
share_meanCWBCO2grossepo3=CCO2_epo3_mean/Present_inv_epo3;
share_stdCWBCO2gross=share_meanCWBCO2gross*((CCO2_acc_std/CCO2_acc_mean)^2)^0.5;
share_stdCWBCO2grossepo1=share_meanCWBCO2grossepo1*((CCO2_epo1_std/CCO2_epo1_mean)^2)^0.5;
share_stdCWBCO2grossepo2=share_meanCWBCO2grossepo2*((CCO2_epo2_std/CCO2_epo2_mean)^2)^0.5;
share_stdCWBCO2grossepo3=share_meanCWBCO2grossepo3*((CCO2_epo3_std/CCO2_epo3_mean)^2)^0.5;

share_meanCWBCO2net=CCO2_acc_mean/Present_Ninv_acc;
share_meanCWBCO2netepo1=CCO2_epo1_mean/Present_Ninv_epo1;
share_meanCWBCO2netepo2=CCO2_epo2_mean/Present_Ninv_epo2;
share_meanCWBCO2netepo3=CCO2_epo3_mean/Present_Ninv_epo3;

share_stdCWBCO2net=share_meanCWBCO2net*((CCO2_acc_std/CCO2_acc_mean)^2)^0.5;
share_stdCWBCO2netepo1=share_meanCWBCO2netepo1*((CCO2_epo1_std/CCO2_epo1_mean)^2)^0.5;
share_stdCWBCO2netepo2=share_meanCWBCO2netepo2*((CCO2_epo2_std/CCO2_epo2_mean)^2)^0.5;
share_stdCWBCO2netepo3=share_meanCWBCO2netepo3*((CCO2_epo3_std/CCO2_epo3_mean)^2)^0.5;

c18=2.11/sqrt(18);

CWB_global_expert_summary=[CCO2_acc_mean, CCO2_epo1_mean,CCO2_epo2_mean,CCO2_epo3_mean;
    CCO2_acc_std, CCO2_epo1_std,CCO2_epo2_std,CCO2_epo3_std; c18*[CCO2_acc_std, CCO2_epo1_std, CCO2_epo2_std, CCO2_epo3_std]];

CWB_global_expert_share=[share_meanCWBCO2gross, share_meanCWBCO2grossepo1,share_meanCWBCO2grossepo2,share_meanCWBCO2grossepo3;
    share_stdCWBCO2gross, share_stdCWBCO2grossepo1,share_stdCWBCO2grossepo2,share_stdCWBCO2grossepo3;c18*[share_stdCWBCO2gross, share_stdCWBCO2grossepo1,share_stdCWBCO2grossepo2,share_stdCWBCO2grossepo3];
    share_meanCWBCO2net, share_meanCWBCO2netepo1,share_meanCWBCO2netepo2,share_meanCWBCO2netepo3;share_stdCWBCO2net, share_stdCWBCO2netepo1,share_stdCWBCO2netepo2,share_stdCWBCO2netepo3;c18*[share_stdCWBCO2net, share_stdCWBCO2netepo1,share_stdCWBCO2netepo2,share_stdCWBCO2netepo3]];

CWB_global_expert_percent=CWB_global_expert_share.*100;

cd C:\..\IW_calculations % Choose the file path 
writematrix(CWB_global_expert_summary,'CWB_global_expert.csv')
writematrix(CWB_global_expert_percent,'CWB_global_share_expert.csv')
writematrix(Factor,'social_discount_factor_expert.csv')

colNames={'Total','Epo1','Epo2','Epo3'};
rowNames={'Mean_gross','Std_gross','CI95_gross','Mean_net','Std_net','CI95_net'};
rowNames2={'Mean','Std','CI95'};

CWB_global_expert_table= array2table(CWB_global_expert_summary,'RowNames',rowNames2,'VariableNames',colNames);
CWB_global_share_expert_table= array2table(CWB_global_expert_percent,'RowNames',rowNames,'VariableNames',colNames);

writetable(CWB_global_expert_table,'CWB_global_expert_table.csv','WriteRowNames',true,'WriteVariableNames',true)
writetable(CWB_global_share_expert_table,'CWB_global_share_expert_table.csv','WriteRowNames',true,'WriteVariableNames',true)



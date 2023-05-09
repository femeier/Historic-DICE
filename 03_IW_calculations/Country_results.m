clear all

number_countries=217; % number of countries
number_years=69; % number of years from 1950 to 2018
C=3.664; % carbon to CO2 conversion factor
ci=2.11/sqrt(18);

%% READ IN DATA (Emissions, Population, Investment, Net Investment)
  num_sheets = 4;
  Inputdata = cell(1,num_sheets);
  
  for i=1:num_sheets
      %Inputdata{i} = readtable('Input_data_raw.xlsx', 'Sheet', i);
      Inputdata{i} = readtable('Input_data_raw_sorted.xlsx', 'Sheet',i); % Countries sorted by CO2 emissions
  end
  
  warning('off','MATLAB:table:ModifiedAndSavedVarnames');
  
    for i=2:number_countries+1
     Emissions(i-1,1:number_years)=Inputdata{1,1}{i,3:number_years+2}; 
     Population(i-1,1:number_years)=Inputdata{1,2}{i,3:number_years+2};
     Investment(i-1,1:number_years)=Inputdata{1,3}{i,3:number_years+2};
     Investment_net(i-1,1:number_years)=Inputdata{1,4}{i,3:number_years+2};
    end
  
Emissions(isnan(Emissions))=0;
Eclean=Emissions;
Eclean(isnan(Eclean))=0;
Ecum=sum(Eclean,2);
countries_label=importdata("country_list.xlsx");

%% READ IN DATA (Global SCC for calibrated and expert-based SDR)

% Calibrated SDR

%cd /Users/cocco/Desktop/Matlab/calib
cd C:\..\SCC_results_2022_constant % Choose the file path 

SCC_HS=readmatrix('sccHSp.xlsx'); 
SCC_N=readmatrix('sccNp.xlsx'); 
SCC_W=readmatrix('sccWp.xlsx'); 

% Drupp et al. SDR

%cd /Users/cocco/Desktop/Matlab/expert
cd C:\..\SCC_results_2022_constant_Drupp % Choose the file path 


SCC_HS_Drupp=readmatrix('sccHSp.xlsx'); 
%SCC_HSo_Drupp=readmatrix('sccHSpo.xlsx'); 
%SCC_HSIo_Drupp=readmatrix('sccHSpoI.xlsx'); 
%SCC_HSoo_Drupp=readmatrix('sccHSpoo.xlsx'); 
SCC_N_Drupp=readmatrix('sccNp.xlsx'); 
%SCC_No_Drupp=readmatrix('sccNpo.xlsx'); 
%SCC_NIo_Drupp=readmatrix('sccNpoI.xlsx'); 
%SCC_Noo_Drupp=readmatrix('sccNpoo.xlsx'); 
SCC_W_Drupp=readmatrix('sccWp.xlsx'); 
%SCC_Wo_Drupp=readmatrix('sccWpo.xlsx'); 
%SCC_WIo_Drupp=readmatrix('sccWpoI.xlsx'); 
%SCC_Woo_Drupp=readmatrix('sccWpoo.xlsx'); 

%% READ IN DATA (social discount factor, calculated in global matlab file)

%cd /Users/cocco/Desktop/Matlab
cd C:\..\IW_calculations\ % Choose the file path 


% Calibrated SDR
Factor_mean=readmatrix('social_discount_factor_calibrated.csv');

% Expert-based SDR
Factor_mean_Drupp=readmatrix('social_discount_factor_expert.csv'); 
%% ASSIGN SCC DATA (to the 6 scenarios)

a1=2; %1950
a2=70; %2018

% Calibrated SDR

for i=1:6
SCCHSp(:,i)=SCC_HS(a1:a2,i+1);
%SCCHSo(:,i)=SCC_HSo(a1:a2,i+1);
%SCCHSIo(:,i)=SCC_HSIo(a1:a2,i+1);
%SCCHSoo(:,i)=SCC_HSoo(a1:a2,i+1);
SCCNp(:,i)=SCC_N(a1:a2,i+1);
%SCCNo(:,i)=SCC_No(a1:a2,i+1);
%SCCNIo(:,i)=SCC_NIo(a1:a2,i+1);
%SCCNoo(:,i)=SCC_Noo(a1:a2,i+1);
SCCWp(:,i)=SCC_W(a1:a2,i+1);
%SCCWo(:,i)=SCC_Wo(a1:a2,i+1);
%SCCWIo(:,i)=SCC_WIo(a1:a2,i+1);
%SCCWoo(:,i)=SCC_Woo(a1:a2,i+1);
end

% Drupp et al. SDR

for i=1:6
SCCHSp_Drupp(:,i)=SCC_HS_Drupp(a1:a2,i+1);
%SCCHSo_Drupp(:,i)=SCC_HSo_Drupp(a1:a2,i+1);
%SCCHSIo_Drupp(:,i)=SCC_HSIo_Drupp(a1:a2,i+1);
%SCCHSoo_Drupp(:,i)=SCC_HSoo_Drupp(a1:a2,i+1);
SCCNp_Drupp(:,i)=SCC_N_Drupp(a1:a2,i+1);
%SCCNo_Drupp(:,i)=SCC_No_Drupp(a1:a2,i+1);
%SCCNIo_Drupp(:,i)=SCC_NIo_Drupp(a1:a2,i+1);
%SCCNoo_Drupp(:,i)=SCC_Noo_Drupp(a1:a2,i+1);
SCCWp_Drupp(:,i)=SCC_W_Drupp(a1:a2,i+1);
%SCCWo_Drupp(:,i)=SCC_Wo_Drupp(a1:a2,i+1);
%SCCWIo_Drupp(:,i)=SCC_WIo_Drupp(a1:a2,i+1);
%SCCWoo_Drupp(:,i)=SCC_Woo_Drupp(a1:a2,i+1);
end

%% Mean SCC and Std for different damage functions

% Calibrated SDR 

% Howard and Sterner
SCCm_HS=mean(SCCHSp,2);
SCCstd_HS=transpose(std(transpose(SCCHSp)));

% Nordhaus
SCCm_N=mean(SCCNp,2);
SCCstd_N=transpose(std(transpose(SCCNp)));

% Weitzman
SCCm_W=mean(SCCWp,2);
SCCstd_W=transpose(std(transpose(SCCWp)));

% All 3 damage functions
SCC=horzcat(SCCHSp,SCCNp,SCCWp);
SCCm=mean(SCC,2);
SCCstd=transpose(std(transpose(SCC)));

% Drupp et al. SDR

% Howard and Sterner
SCCm_HS_Drupp=mean(SCCHSp_Drupp,2);
SCCstd_HS_Drupp=transpose(std(transpose(SCCHSp_Drupp)));

% Nordhaus
SCCm_N_Drupp=mean(SCCNp_Drupp,2);
SCCstd_N_Drupp=transpose(std(transpose(SCCNp_Drupp)));

% Weitzman
SCCm_W_Drupp=mean(SCCWp_Drupp,2);
SCCstd_W_Drupp=transpose(std(transpose(SCCWp_Drupp)));

% All 3 damage functions
SCC_Drupp=horzcat(SCCHSp_Drupp,SCCNp_Drupp,SCCWp_Drupp);
SCCm_Drupp=mean(SCC_Drupp,2);
SCCstd_Drupp=transpose(std(transpose(SCC_Drupp)));

%% Calculation: Climate Wealth Borrowing

j1=1; %1950
j2=25; %1974
j3=26; %1975
j4=50; %1999
j5=51; %2000
j6=69; %2018

% Calibrated SDR

% Howard and Sterner
for j=1:number_years
    for i=1:number_countries
     SCC_PV_HS(i,j)=Emissions(i,j).*SCCm_HS(j)*C;
     SCCstd_PV_HS(i,j)=Emissions(i,j).*SCCstd_HS(j).*C;
    end
end
% Climate Wealth Borrowing for different time periods
% From 1950 to 2018
    for i=1:number_countries
     CWB_HS(i,:)=sum(SCC_PV_HS(i,j1:j6));
     CWBstd_HS(i,:)=(sum(SCCstd_PV_HS(i,j1:j6).^2)).^0.5;
    end
CWBstd_HS_agg=sqrt(sum((SCCstd_HS.*transpose(sum(Emissions*C))).^2))/10^6;
CWBci_HS_agg=ci*CWBstd_HS_agg;
% From 1950 to 1974
    for i=1:number_countries
     CWB_HS_epo1(i,:)=sum(SCC_PV_HS(i,j1:j2));
     CWBstd_HS_epo1(i,:)=(sum(SCCstd_PV_HS(i,j1:j2).^2)).^0.5;
    end
 CWBstd_HS_epo1_agg=sqrt(sum((SCCstd_HS(j1:j2).*transpose(sum(Emissions(j1:j2)*C))).^2))/10^6;
 CWBci_HS_epo1_agg=ci*CWBstd_HS_epo1_agg;

    % From 1974 to 1999
    for i=1:number_countries
     CWB_HS_epo2(i,:)=sum(SCC_PV_HS(i,j3:j4));
     CWBstd_HS_epo2(i,:)=(sum(SCCstd_PV_HS(i,j3:j4).^2)).^0.5;
    end
CWBstd_HS_epo2_agg=sqrt(sum((SCCstd_HS(j3:j4).*transpose(sum(Emissions(j3:j4)*C))).^2))/10^6;
CWBci_HS_epo2_agg=ci*CWBstd_HS_epo2_agg;

% From 1999 to 2018
    for i=1:number_countries
     CWB_HS_epo3(i,:)=sum(SCC_PV_HS(i,j5:j6));
     CWBstd_HS_epo3(i,:)=(sum(SCCstd_PV_HS(i,j5:j6).^2)).^0.5;
    end
CWBstd_HS_epo3_agg=sqrt(sum((SCCstd_HS(j5:j6).*transpose(sum(Emissions(j5:j6)*C))).^2))/10^6;
CWBci_HS_epo3_agg=ci*CWBstd_HS_epo3_agg;

% Nordhaus
for j=1:number_years
    for i=1:number_countries
     SCC_PV_N(i,j)=Emissions(i,j).*SCCm_N(j)*C;
     SCCstd_PV_N(i,j)=Emissions(i,j).*SCCstd_N(j).*C;
    end
end
% Climate Wealth Borrowing for different time periods
% From 1950 to 2018
    for i=1:number_countries
     CWB_N(i,:)=sum(SCC_PV_N(i,j1:j6));
     CWBstd_N(i,:)=(sum(SCCstd_PV_N(i,j1:j6).^2)).^0.5;
    end
CWBstd_N_agg=sqrt(sum((SCCstd_N.*transpose(sum(Emissions*C))).^2))/10^6;
CWBci_N_agg=ci*CWBstd_N_agg;

% From 1950 to 1974
    for i=1:number_countries
     CWB_N_epo1(i,:)=sum(SCC_PV_N(i,j1:j2));
     CWBstd_N_epo1(i,:)=(sum(SCCstd_PV_N(i,j1:j2).^2)).^0.5;
    end
CWBstd_N_epo1_agg=sqrt(sum((SCCstd_N(j1:j2).*transpose(sum(Emissions(j1:j2)*C))).^2))/10^6;
CWBci_N_epo1_agg=ci*CWBstd_N_epo1_agg;

    % From 1974 to 1999
    for i=1:number_countries
     CWB_N_epo2(i,:)=sum(SCC_PV_N(i,j3:j4));
     CWBstd_N_epo2(i,:)=(sum(SCCstd_PV_N(i,j3:j4).^2)).^0.5;
    end
CWBstd_N_epo2_agg=sqrt(sum((SCCstd_N(j3:j4).*transpose(sum(Emissions(j3:j4)*C))).^2))/10^6;
CWBci_N_epo2_agg=ci*CWBstd_N_epo2_agg;

% From 1999 to 2018
    for i=1:number_countries
     CWB_N_epo3(i,:)=sum(SCC_PV_N(i,j5:j6));
     CWBstd_N_epo3(i,:)=(sum(SCCstd_PV_N(i,j5:j6).^2)).^0.5;
    end
CWBstd_N_epo3_agg=sqrt(sum((SCCstd_N(j5:j6).*transpose(sum(Emissions(j5:j6)*C))).^2))/10^6;
CWBci_N_epo3_agg=ci*CWBstd_N_epo3_agg;

% Weitzman
for j=1:number_years
    for i=1:number_countries
     SCC_PV_W(i,j)=Emissions(i,j).*SCCm_W(j)*C;
     SCCstd_PV_W(i,j)=Emissions(i,j).*SCCstd_W(j).*C;
    end
end
% Climate Wealth Borrowing for different time periods
% From 1950 to 2018
for i=1:number_countries
     CWB_W(i,:)=sum(SCC_PV_W(i,j1:j6));
     CWBstd_W(i,:)=(sum(SCCstd_PV_W(i,j1:j6).^2)).^0.5;
end
CWBstd_W_agg=sqrt(sum((SCCstd_W.*transpose(sum(Emissions*C))).^2))/10^6;
CWBci_W_agg=ci*CWBstd_W_agg;

% From 1950 to 1974
    for i=1:number_countries
     CWB_W_epo1(i,:)=sum(SCC_PV_W(i,j1:j2));
     CWBstd_W_epo1(i,:)=(sum(SCCstd_PV_W(i,j1:j2).^2)).^0.5;
    end
CWBstd_W_epo1_agg=sqrt(sum((SCCstd_W(j1:j2).*transpose(sum(Emissions(j1:j2)*C))).^2))/10^6;
CWBci_W_epo1_agg=ci*CWBstd_W_epo1_agg;

% From 1974 to 1999
    for i=1:number_countries
     CWB_W_epo2(i,:)=sum(SCC_PV_W(i,j3:j4));
     CWBstd_W_epo2(i,:)=(sum(SCCstd_PV_W(i,j3:j4).^2)).^0.5;
    end
CWBstd_W_epo2_agg=sqrt(sum((SCCstd_W(j3:j4).*transpose(sum(Emissions(j3:j4)*C))).^2))/10^6;
CWBci_W_epo2_agg=ci*CWBstd_W_epo2_agg;

    % From 1999 to 2018
    for i=1:number_countries
     CWB_W_epo3(i,:)=sum(SCC_PV_W(i,j5:j6));
     CWBstd_W_epo3(i,:)=(sum(SCCstd_PV_W(i,j5:j6).^2)).^0.5;
    end
CWBstd_W_epo3_agg=sqrt(sum((SCCstd_W(j5:j6).*transpose(sum(Emissions(j5:j6)*C))).^2))/10^6;
CWBci_W_epo3_agg=ci*CWBstd_W_epo3_agg;

% Mean (all 3 damage functions)
for j=1:number_years
    for i=1:number_countries
     SCC_PV(i,j)=Emissions(i,j).*SCCm(j)*C;
     SCCstd_PV(i,j)=Emissions(i,j).*SCCstd(j).*C;
    end
end
% Climate Wealth Borrowing for different time periods
% From 1950 to 2018
    for i=1:number_countries
     CWB(i,:)=sum(SCC_PV(i,j1:j6));
     CWBstd(i,:)=(sum(SCCstd_PV(i,j1:j6).^2)).^0.5;
    end
CWBstd_agg=sqrt(sum((SCCstd.*transpose(sum(Emissions*C))).^2))/10^6;
CWBci_agg=ci*sqrt(sum((SCCstd.*transpose(sum(Emissions*C))).^2))/10^6;


% From 1950 to 1974
    for i=1:number_countries
     CWB_epo1(i,:)=sum(SCC_PV(i,j1:j2));
     CWBstd_epo1(i,:)=(sum(SCCstd_PV(i,j1:j2).^2)).^0.5;
    end
 CWBstd_epo1_agg=sqrt(sum((SCCstd(j1:j2).*transpose(sum(Emissions(j1:j2)*C))).^2))/10^6;
 CWBci_epo1_agg=ci*CWBstd_epo1_agg;
% From 1974 to 1999
    for i=1:number_countries
     CWB_epo2(i,:)=sum(SCC_PV(i,j3:j4));
     CWBstd_epo2(i,:)=(sum(SCCstd_PV(i,j3:j4).^2)).^0.5;
    end
CWBstd_epo2_agg=sqrt(sum((SCCstd(j3:j4).*transpose(sum(Emissions(j3:j4)*C))).^2))/10^6;
CWBci_epo2_agg=ci*CWBstd_epo2_agg;
% From 1999 to 2018
    for i=1:number_countries
     CWB_epo3(i,:)=sum(SCC_PV(i,j5:j6));
     CWBstd_epo3(i,:)=(sum(SCCstd_PV(i,j5:j6).^2)).^0.5;
    end
CWBstd_epo3_agg=sqrt(sum((SCCstd(j5:j6).*transpose(sum(Emissions(j5:j6)*C))).^2))/10^6;
CWBci_epo3_agg=ci*CWBstd_epo3_agg;

% Drupp et al. SDR

% Howard and Sterner
for j=1:number_years
    for i=1:number_countries
     SCC_PV_HS_Drupp(i,j)=Emissions(i,j).*SCCm_HS_Drupp(j)*C;
     SCCstd_PV_HS_Drupp(i,j)=Emissions(i,j).*SCCstd_HS_Drupp(j).*C;
    end
end
% Climate Wealth Borrowing for different time periods
% From 1950 to 2018
    for i=1:number_countries
     CWB_HS_Drupp(i,:)=sum(SCC_PV_HS_Drupp(i,j1:j6));
     CWBstd_HS_Drupp(i,:)=(sum(SCCstd_PV_HS_Drupp(i,j1:j6).^2)).^0.5;
    end
    CWBstd_HS_agg_Drupp=sqrt(sum((SCCstd_HS_Drupp.*transpose(sum(Emissions*C))).^2))/10^6;
     CWBci_HS_agg_Drupp=ci*CWBstd_HS_agg_Drupp;
% From 1950 to 1974
    for i=1:number_countries
     CWB_HS_epo1_Drupp(i,:)=sum(SCC_PV_HS_Drupp(i,j1:j2));
     CWBstd_HS_epo1_Drupp(i,:)=(sum(SCCstd_PV_HS_Drupp(i,j1:j2).^2)).^0.5;
    end
    CWBstd_HS_epo1_agg_Drupp=sqrt(sum((SCCstd_HS_Drupp(j1:j2).*transpose(sum(Emissions(j1:j2)*C))).^2))/10^6;
     CWBci_HS_epo1_agg_Drupp=ci*CWBstd_HS_epo1_agg_Drupp;
% From 1974 to 1999
    for i=1:number_countries
     CWB_HS_epo2_Drupp(i,:)=sum(SCC_PV_HS_Drupp(i,j3:j4));
     CWBstd_HS_epo2_Drupp(i,:)=(sum(SCCstd_PV_HS_Drupp(i,j3:j4).^2)).^0.5;
    end
    CWBstd_HS_epo2_agg_Drupp=sqrt(sum((SCCstd_HS_Drupp(j3:j4).*transpose(sum(Emissions(j3:j4)*C))).^2))/10^6;
     CWBci_HS_epo2_agg_Drupp=ci*CWBstd_HS_epo2_agg_Drupp;
% From 1999 to 2018
    for i=1:number_countries
     CWB_HS_epo3_Drupp(i,:)=sum(SCC_PV_HS_Drupp(i,j5:j6));
     CWBstd_HS_epo3_Drupp(i,:)=(sum(SCCstd_PV_HS_Drupp(i,j5:j6).^2)).^0.5;
    end
CWBstd_HS_epo3_agg_Drupp=sqrt(sum((SCCstd_HS_Drupp(j5:j6).*transpose(sum(Emissions(j5:j6)*C))).^2))/10^6;
CWBci_HS_epo3_agg_Drupp=ci*CWBstd_HS_epo3_agg_Drupp;

% Nordhaus
for j=1:number_years
    for i=1:number_countries
     SCC_PV_N_Drupp(i,j)=Emissions(i,j).*SCCm_N_Drupp(j)*C;
     SCCstd_PV_N_Drupp(i,j)=Emissions(i,j).*SCCstd_N_Drupp(j).*C;
    end
end
% Climate Wealth Borrowing for different time periods
% From 1950 to 2018
    for i=1:number_countries
     CWB_N_Drupp(i,:)=sum(SCC_PV_N_Drupp(i,j1:j6));
     CWBstd_N_Drupp(i,:)=(sum(SCCstd_PV_N_Drupp(i,j1:j6).^2)).^0.5;
    end
    CWBstd_N_agg_Drupp=sqrt(sum((SCCstd_N_Drupp.*transpose(sum(Emissions*C))).^2))/10^6;
     CWBci_N_agg_Drupp=ci*CWBstd_N_agg_Drupp;
% From 1950 to 1974
    for i=1:number_countries
     CWB_N_epo1_Drupp(i,:)=sum(SCC_PV_N_Drupp(i,j1:j2));
     CWBstd_N_epo1_Drupp(i,:)=(sum(SCCstd_PV_N_Drupp(i,j1:j2).^2)).^0.5;
    end
     CWBstd_N_epo1_agg_Drupp=sqrt(sum((SCCstd_N_Drupp(j1:j2).*transpose(sum(Emissions(j1:j2)*C))).^2))/10^6;
     CWBci_N_epo1_agg_Drupp=ci*CWBstd_N_epo1_agg_Drupp;
% From 1974 to 1999
    for i=1:number_countries
     CWB_N_epo2_Drupp(i,:)=sum(SCC_PV_N_Drupp(i,j3:j4));
     CWBstd_N_epo2_Drupp(i,:)=(sum(SCCstd_PV_N_Drupp(i,j3:j4).^2)).^0.5;
    end
    CWBstd_N_epo2_agg_Drupp=sqrt(sum((SCCstd_N_Drupp(j3:j4).*transpose(sum(Emissions(j3:j4)*C))).^2))/10^6;
    CWBci_N_epo2_agg_Drupp=ci*CWBstd_N_epo2_agg_Drupp;
% From 1999 to 2018
    for i=1:number_countries
     CWB_N_epo3_Drupp(i,:)=sum(SCC_PV_N_Drupp(i,j5:j6));
     CWBstd_N_epo3_Drupp(i,:)=(sum(SCCstd_PV_N_Drupp(i,j5:j6).^2)).^0.5;
    end
    CWBstd_N_epo3_agg_Drupp=sqrt(sum((SCCstd_N_Drupp(j5:j6).*transpose(sum(Emissions(j5:j6)*C))).^2))/10^6;
    CWBci_N_epo3_agg_Drupp=ci*CWBstd_N_epo3_agg_Drupp;

% Weitzman
for j=1:number_years
    for i=1:number_countries
     SCC_PV_W_Drupp(i,j)=Emissions(i,j).*SCCm_W_Drupp(j)*C;
     SCCstd_PV_W_Drupp(i,j)=Emissions(i,j).*SCCstd_W_Drupp(j).*C;
    end
end
% Climate Wealth Borrowing for different time periods
% From 1950 to 2018
for i=1:number_countries
     CWB_W_Drupp(i,:)=sum(SCC_PV_W_Drupp(i,j1:j6));
     CWBstd_W_Drupp(i,:)=(sum(SCCstd_PV_W_Drupp(i,j1:j6).^2)).^0.5;
end
CWBstd_W_agg_Drupp=sqrt(sum((SCCstd_W_Drupp.*transpose(sum(Emissions*C))).^2))/10^6;
     CWBci_W_agg_Drupp=ci*CWBstd_W_agg_Drupp;

% From 1950 to 1974
    for i=1:number_countries
     CWB_W_epo1_Drupp(i,:)=sum(SCC_PV_W_Drupp(i,j1:j2));
     CWBstd_W_epo1_Drupp(i,:)=(sum(SCCstd_PV_W_Drupp(i,j1:j2).^2)).^0.5;
    end
    CWBstd_W_epo1_agg_Drupp=sqrt(sum((SCCstd_W_Drupp(j1:j2).*transpose(sum(Emissions(j1:j2)*C))).^2))/10^6;
     CWBci_W_epo1_agg_Drupp=ci*CWBstd_W_epo1_agg_Drupp;
% From 1974 to 1999
    for i=1:number_countries
     CWB_W_epo2_Drupp(i,:)=sum(SCC_PV_W_Drupp(i,j3:j4));
     CWBstd_W_epo2_Drupp(i,:)=(sum(SCCstd_PV_W_Drupp(i,j3:j4).^2)).^0.5;
    end
    CWBstd_W_epo2_agg_Drupp=sqrt(sum((SCCstd_W_Drupp(j3:j4).*transpose(sum(Emissions(j3:j4)*C))).^2))/10^6;
     CWBci_W_epo2_agg_Drupp=ci*CWBstd_W_epo2_agg_Drupp;
% From 1999 to 2018
    for i=1:number_countries
     CWB_W_epo3_Drupp(i,:)=sum(SCC_PV_W_Drupp(i,j5:j6));
     CWBstd_W_epo3_Drupp(i,:)=(sum(SCCstd_PV_W_Drupp(i,j5:j6).^2)).^0.5;
    end
    CWBstd_W_epo3_agg_Drupp=sqrt(sum((SCCstd_W_Drupp(j5:j6).*transpose(sum(Emissions(j5:j6)*C))).^2))/10^6;
     CWBci_W_epo3_agg_Drupp=ci*CWBstd_W_epo3_agg_Drupp;


% Mean (all 3 damage functions)
for j=1:number_years
    for i=1:number_countries
     SCC_PV_Drupp(i,j)=Emissions(i,j).*SCCm_Drupp(j)*C;
     SCCstd_PV_Drupp(i,j)=Emissions(i,j).*SCCstd_Drupp(j).*C;
    end
end
% Climate Wealth Borrowing for different time periods
% From 1950 to 2018
    for i=1:number_countries
     CWB_Drupp(i,:)=sum(SCC_PV_Drupp(i,j1:j6));
     CWBstd_Drupp(i,:)=(sum(SCCstd_PV_Drupp(i,j1:j6).^2)).^0.5;
    end
   CWBstd_agg_Drupp=sqrt(sum((SCCstd_Drupp.*transpose(sum(Emissions*C))).^2))/10^6;
   CWBci_agg_Drupp=ci*sqrt(sum((SCCstd_Drupp.*transpose(sum(Emissions*C))).^2))/10^6;
% From 1950 to 1974
    for i=1:number_countries
     CWB_epo1_Drupp(i,:)=sum(SCC_PV_Drupp(i,j1:j2));
     CWBstd_epo1_Drupp(i,:)=(sum(SCCstd_PV_Drupp(i,j1:j2).^2)).^0.5;
    end
     CWBstd_epo1_agg_Drupp=sqrt(sum((SCCstd_Drupp(j1:j2).*transpose(sum(Emissions(j1:j2)*C))).^2))/10^6;
     CWBci_epo1_agg_Drupp=ci*CWBstd_epo1_agg_Drupp;
% From 1974 to 1999
    for i=1:number_countries
     CWB_epo2_Drupp(i,:)=sum(SCC_PV_Drupp(i,j3:j4));
     CWBstd_epo2_Drupp(i,:)=(sum(SCCstd_PV_Drupp(i,j3:j4).^2)).^0.5;
    end
    CWBstd_epo2_agg_Drupp=sqrt(sum((SCCstd_Drupp(j3:j4).*transpose(sum(Emissions(j3:j4)*C))).^2))/10^6;
     CWBci_epo2_agg_Drupp=ci*CWBstd_epo2_agg_Drupp;
% From 1999 to 2018
    for i=1:number_countries
     CWB_epo3_Drupp(i,:)=sum(SCC_PV_Drupp(i,j5:j6));
     CWBstd_epo3_Drupp(i,:)=(sum(SCCstd_PV_Drupp(i,j5:j6).^2)).^0.5;
     end
CWBstd_epo3_agg_Drupp=sqrt(sum((SCCstd_Drupp(j5:j6).*transpose(sum(Emissions(j5:j6)*C))).^2))/10^6;
CWBci_epo3_agg_Drupp=ci*CWBstd_epo3_agg_Drupp;
%% Calculation: Accumulated Investment

% Calibrated SDR 

% Present value investment (mean)
for j=j1:j6
    for i=1:number_countries
     Investment_PV(i,j)=Factor_mean(j).*Investment(i,j);
     Investment_net_PV(i,j)=Factor_mean(j).*Investment_net(i,j);
    end
end

% Accumulated investment for different time periods
% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     Investment_Acc(i,j)=sum(Investment_PV(i,j:j6));
     Investment_net_Acc(i,j)=sum(Investment_net_PV(i,j:j6));
    end
end
% From 1950 to 1974
    for i=1:number_countries
     Investment_Acc_epo1(i,:)=sum(Investment_PV(i,j1:j2));
     Investment_net_Acc_epo1(i,:)=sum(Investment_net_PV(i,j1:j2));
    end
% From 1975 to 1999
    for i=1:number_countries
     Investment_Acc_epo2(i,:)=sum(Investment_PV(i,j3:j4));
     Investment_net_Acc_epo2(i,:)=sum(Investment_net_PV(i,j3:j4));
    end
% From 1999 to 2018
    for i=1:number_countries
     Investment_Acc_epo3(i,:)=sum(Investment_PV(i,j5:j6));
     Investment_net_Acc_epo3(i,:)=sum(Investment_net_PV(i,j5:j6));
    end

% Drupp et al. SDR

% Present value investment (mean)
for j=j1:j6
    for i=1:number_countries
     Investment_PV_Drupp(i,j)=Factor_mean_Drupp(j).*Investment(i,j);
     Investment_net_PV_Drupp(i,j)=Factor_mean_Drupp(j).*Investment_net(i,j);
    end
end

% Accumulated investment for different time periods
% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     Investment_Acc_Drupp(i,j)=sum(Investment_PV_Drupp(i,j:j6));
     Investment_net_Acc_Drupp(i,j)=sum(Investment_net_PV_Drupp(i,j:j6));
    end
end
% From 1950 to 1974
    for i=1:number_countries
     Investment_Acc_epo1_Drupp(i,:)=sum(Investment_PV_Drupp(i,j1:j2));
     Investment_net_Acc_epo1_Drupp(i,:)=sum(Investment_net_PV_Drupp(i,j1:j2));
    end
% From 1975 to 1999
    for i=1:number_countries
     Investment_Acc_epo2_Drupp(i,:)=sum(Investment_PV_Drupp(i,j3:j4));
     Investment_net_Acc_epo2_Drupp(i,:)=sum(Investment_net_PV_Drupp(i,j3:j4));
    end
% From 1999 to 2018
    for i=1:number_countries
     Investment_Acc_epo3_Drupp(i,:)=sum(Investment_PV_Drupp(i,j5:j6));
     Investment_net_Acc_epo3_Drupp(i,:)=sum(Investment_net_PV_Drupp(i,j5:j6));
    end

%% CALCULATION: DISINVESTMENT SHARE

% Indicator assures that we accumulate only over years for which we have both emissions and investment data
Indicator=Investment;  
Indicator(Indicator>0)=1;
Indicator(Indicator<0)=1;
Indicator(isnan(Indicator))=0;
Emissions_dis=Indicator.*Emissions;

Indicator_net=Investment_net;  
Indicator_net(Indicator_net>0)=1;
Indicator_net(Indicator_net<0)=1;
Indicator_net(isnan(Indicator_net))=0;
Emissions_dis_net=Indicator_net.*Emissions;

% Calibrated SDR 

% Howard and Sterner

for j=1:number_years
    for i=1:number_countries
     SCC_PV_HS_dis(i,j)=Emissions_dis(i,j).*SCCm_HS(j)*C;
     SCC_PV_HS_dis_net(i,j)=Emissions_dis_net(i,j).*SCCm_HS(j)*C;
     SCCstd_PV_HS_dis(i,j)=Emissions_dis(i,j).*SCCstd_HS(j)*C;
     SCCstd_PV_HS_dis_net(i,j)=Emissions_dis_net(i,j).*SCCstd_HS(j)*C;
         end
end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     CWB_HS_dis(i,j)=sum(SCC_PV_HS_dis(i,j:j6));
     CWB_HS_dis_net(i,j)=sum(SCC_PV_HS_dis_net(i,j:j6));
     CWBstd_HS_dis(i,j)=(sum(SCCstd_PV_HS_dis(i,j:j6).^2)).^0.5;
     CWBstd_HS_dis_net(i,j)=(sum(SCCstd_PV_HS_dis_net(i,j:j6).^2)).^0.5;
    end
end
% From 1950 to 1974
    for i=1:number_countries
     CWB_HS_dis_epo1(i,:)=sum(SCC_PV_HS_dis(i,j1:j2));
     CWB_HS_dis_epo1_net(i,:)=sum(SCC_PV_HS_dis_net(i,j1:j2));
     CWBstd_HS_dis_epo1(i,:)=(sum(SCCstd_PV_HS_dis(i,j1:j2).^2)).^0.5;
     CWBstd_HS_dis_epo1_net(i,:)=(sum(SCCstd_PV_HS_dis_net(i,j1:j2).^2)).^0.5;
    end
% From 1974 to 1999
    for i=1:number_countries
     CWB_HS_dis_epo2(i,:)=sum(SCC_PV_HS_dis(i,j3:j4));
     CWB_HS_dis_epo2_net(i,:)=sum(SCC_PV_HS_dis_net(i,j3:j4));
     CWBstd_HS_dis_epo2(i,:)=(sum(SCCstd_PV_HS_dis(i,j3:j4).^2)).^0.5;
     CWBstd_HS_dis_epo2_net(i,:)=(sum(SCCstd_PV_HS_dis_net(i,j3:j4).^2)).^0.5;
    end
% From 1999 to 2018
    for i=1:number_countries
     CWB_HS_dis_epo3(i,:)=sum(SCC_PV_HS_dis(i,j5:j6));
     CWB_HS_dis_epo3_net(i,:)=sum(SCC_PV_HS_dis_net(i,j5:j6));
     CWBstd_HS_dis_epo3(i,:)=(sum(SCCstd_PV_HS_dis(i,j5:j6).^2)).^0.5;
     CWBstd_HS_dis_epo3_net(i,:)=(sum(SCCstd_PV_HS_dis_net(i,j5:j6).^2)).^0.5;

    end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     Dis_inv_share_HS(i,j)=CWB_HS_dis(i,j)./Investment_Acc(i,j);
     Dis_inv_share_HS_net(i,j)=CWB_HS_dis_net(i,j)./Investment_net_Acc(i,j);
     Dis_inv_share_std_HS(i,j)=CWBstd_HS_dis(i,j)./Investment_Acc(i,j);
     Dis_inv_share_std_HS_net(i,j)=CWBstd_HS_dis_net(i,j)./Investment_net_Acc(i,j);
     end
end
% From 1950 to 1974
    for i=1:number_countries
     Dis_inv_share_HS_epo1(i,:)=CWB_HS_dis_epo1(i,:)./Investment_Acc_epo1(i,:);
     Dis_inv_share_HS_epo1_net(i,:)=CWB_HS_dis_epo1_net(i,:)./Investment_net_Acc_epo1(i,:);
     Dis_inv_share_std_HS_epo1(i,:)=CWBstd_HS_dis_epo1(i,:)./Investment_Acc_epo1(i,:);
     Dis_inv_share_std_HS_epo1_net(i,:)=CWBstd_HS_dis_epo1_net(i,:)./Investment_net_Acc_epo1(i,:);
    end
% From 1974 to 1999
    for i=1:number_countries
     Dis_inv_share_HS_epo2(i,:)=CWB_HS_dis_epo2(i,:)./Investment_Acc_epo2(i,:);
     Dis_inv_share_HS_epo2_net(i,:)=CWB_HS_dis_epo2_net(i,:)./Investment_net_Acc_epo2(i,:);
     Dis_inv_share_std_HS_epo2(i,:)=CWBstd_HS_dis_epo2(i,:)./Investment_Acc_epo2(i,:);
     Dis_inv_share_std_HS_epo2_net(i,:)=CWBstd_HS_dis_epo2_net(i,:)./Investment_net_Acc_epo2(i,:);
    end
% From 1999 to 2018
    for i=1:number_countries
     Dis_inv_share_HS_epo3(i,:)=CWB_HS_dis_epo3(i,:)./Investment_Acc_epo3(i,:);
     Dis_inv_share_HS_epo3_net(i,:)=CWB_HS_dis_epo3_net(i,:)./Investment_net_Acc_epo3(i,:);
     Dis_inv_share_std_HS_epo3(i,:)=CWBstd_HS_dis_epo3(i,:)./Investment_Acc_epo3(i,:);
     Dis_inv_share_std_HS_epo3_net(i,:)=CWBstd_HS_dis_epo3_net(i,:)./Investment_net_Acc_epo3(i,:);

    end

% Nordhaus

for j=1:number_years
    for i=1:number_countries
     SCC_PV_N_dis(i,j)=Emissions_dis(i,j).*SCCm_N(j)*C;
     SCC_PV_N_dis_net(i,j)=Emissions_dis_net(i,j).*SCCm_N(j)*C;
     SCCstd_PV_N_dis(i,j)=Emissions_dis(i,j).*SCCstd_N(j)*C;
     SCCstd_PV_N_dis_net(i,j)=Emissions_dis_net(i,j).*SCCstd_N(j)*C;
    end
end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     CWB_N_dis(i,j)=sum(SCC_PV_N_dis(i,j:j6));
     CWB_N_dis_net(i,j)=sum(SCC_PV_N_dis_net(i,j:j6));
     CWBstd_N_dis(i,j)=(sum(SCCstd_PV_N_dis(i,j:j6).^2)).^0.5;
     CWBstd_N_dis_net(i,j)=(sum(SCCstd_PV_N_dis_net(i,j:j6).^2)).^0.5;
    end
end
% From 1950 to 1974
    for i=1:number_countries
     CWB_N_dis_epo1(i,:)=sum(SCC_PV_N_dis(i,j1:j2));
     CWB_N_dis_epo1_net(i,:)=sum(SCC_PV_N_dis_net(i,j1:j2));
     CWBstd_N_dis_epo1(i,:)=(sum(SCCstd_PV_N_dis(i,j1:j2).^2)).^0.5;
     CWBstd_N_dis_epo1_net(i,:)=(sum(SCCstd_PV_N_dis_net(i,j1:j2).^2)).^0.5;
    end
% From 1974 to 1999
    for i=1:number_countries
     CWB_N_dis_epo2(i,:)=sum(SCC_PV_N_dis(i,j3:j4));
     CWB_N_dis_epo2_net(i,:)=sum(SCC_PV_N_dis_net(i,j3:j4));
     CWBstd_N_dis_epo2(i,:)=(sum(SCCstd_PV_N_dis(i,j3:j4).^2)).^0.5;
     CWBstd_N_dis_epo2_net(i,:)=(sum(SCCstd_PV_N_dis_net(i,j3:j4).^2)).^0.5;
    end
% From 1999 to 2018
    for i=1:number_countries
     CWB_N_dis_epo3(i,:)=sum(SCC_PV_N_dis(i,j5:j6));
     CWB_N_dis_epo3_net(i,:)=sum(SCC_PV_N_dis_net(i,j5:j6));
     CWBstd_N_dis_epo3(i,:)=(sum(SCCstd_PV_N_dis(i,j5:j6).^2)).^0.5;
     CWBstd_N_dis_epo3_net(i,:)=(sum(SCCstd_PV_N_dis_net(i,j5:j6).^2)).^0.5;
    end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     Dis_inv_share_N(i,j)=CWB_N_dis(i,j)./Investment_Acc(i,j);
     Dis_inv_share_N_net(i,j)=CWB_N_dis_net(i,j)./Investment_net_Acc(i,j);
     Dis_inv_share_std_N(i,j)=CWBstd_N_dis(i,j)./Investment_Acc(i,j);
     Dis_inv_share_std_N_net(i,j)=CWBstd_N_dis_net(i,j)./Investment_net_Acc(i,j);

    end
end
% From 1950 to 1974
    for i=1:number_countries
     Dis_inv_share_N_epo1(i,:)=CWB_N_dis_epo1(i,:)./Investment_Acc_epo1(i,:);
     Dis_inv_share_N_epo1_net(i,:)=CWB_N_dis_epo1_net(i,:)./Investment_net_Acc_epo1(i,:);
     Dis_inv_share_std_N_epo1(i,:)=CWBstd_N_dis_epo1(i,:)./Investment_Acc_epo1(i,:);
     Dis_inv_share_std_N_epo1_net(i,:)=CWBstd_N_dis_epo1_net(i,:)./Investment_net_Acc_epo1(i,:);
    end
% From 1974 to 1999
    for i=1:number_countries
     Dis_inv_share_N_epo2(i,:)=CWB_N_dis_epo2(i,:)./Investment_Acc_epo2(i,:);
     Dis_inv_share_N_epo2_net(i,:)=CWB_N_dis_epo2_net(i,:)./Investment_net_Acc_epo2(i,:);
     Dis_inv_share_std_N_epo2(i,:)=CWBstd_N_dis_epo2(i,:)./Investment_Acc_epo2(i,:);
     Dis_inv_share_std_N_epo2_net(i,:)=CWBstd_N_dis_epo2_net(i,:)./Investment_net_Acc_epo2(i,:);
    end
% From 1999 to 2018
    for i=1:number_countries
     Dis_inv_share_N_epo3(i,:)=CWB_N_dis_epo3(i,:)./Investment_Acc_epo3(i,:);
     Dis_inv_share_N_epo3_net(i,:)=CWB_N_dis_epo3_net(i,:)./Investment_net_Acc_epo3(i,:);
     Dis_inv_share_std_N_epo3(i,:)=CWBstd_N_dis_epo3(i,:)./Investment_Acc_epo3(i,:);
     Dis_inv_share_std_N_epo3_net(i,:)=CWBstd_N_dis_epo3_net(i,:)./Investment_net_Acc_epo3(i,:);
    end


% Weitzman

for j=1:number_years
    for i=1:number_countries
     SCC_PV_W_dis(i,j)=Emissions_dis(i,j).*SCCm_W(j)*C;
     SCCstd_PV_W_dis(i,j)=Emissions_dis(i,j).*SCCstd_W(j).*C;
     SCC_PV_W_dis_net(i,j)=Emissions_dis_net(i,j).*SCCm_W(j)*C;
     SCCstd_PV_W_dis_net(i,j)=Emissions_dis_net(i,j).*SCCstd_W(j)*C;

    end
end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     CWB_W_dis(i,j)=sum(SCC_PV_W_dis(i,j:j6));
     CWBstd_W_dis(i,j)=(sum(SCCstd_PV_W_dis(i,j:j6).^2)).^0.5;
     CWB_W_dis_net(i,j)=sum(SCC_PV_W_dis_net(i,j:j6));
     CWBstd_W_dis_net(i,j)=(sum(SCCstd_PV_W_dis_net(i,j:j6).^2)).^0.5;
    end
end
% From 1950 to 1974
    for i=1:number_countries
     CWB_W_dis_epo1(i,:)=sum(SCC_PV_W_dis(i,j1:j2));
     CWBstd_W_dis_epo1(i,:)=(sum(SCCstd_PV_W_dis(i,j1:j2).^2)).^0.5;
     CWB_W_dis_epo1_net(i,:)=sum(SCC_PV_W_dis_net(i,j1:j2));
     CWBstd_W_dis_epo1_net(i,:)=(sum(SCCstd_PV_W_dis_net(i,j1:j2).^2)).^0.5;
    end
% From 1974 to 1999
    for i=1:number_countries
     CWB_W_dis_epo2(i,:)=sum(SCC_PV_W_dis(i,j3:j4));
     CWBstd_W_dis_epo2(i,:)=(sum(SCCstd_PV_W_dis(i,j3:j4).^2)).^0.5;
     CWB_W_dis_epo2_net(i,:)=sum(SCC_PV_W_dis_net(i,j3:j4));
     CWBstd_W_dis_epo2_net(i,:)=(sum(SCCstd_PV_W_dis_net(i,j3:j4).^2)).^0.5;
    end
% From 1999 to 2018
    for i=1:number_countries
     CWB_W_dis_epo3(i,:)=sum(SCC_PV_W_dis(i,j5:j6));
     CWBstd_W_dis_epo3(i,:)=(sum(SCCstd_PV_W_dis(i,j5:j6).^2)).^0.5;
     CWB_W_dis_epo3_net(i,:)=sum(SCC_PV_W_dis_net(i,j5:j6));
     CWBstd_W_dis_epo3_net(i,:)=(sum(SCCstd_PV_W_dis_net(i,j5:j6).^2)).^0.5;
    end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     Dis_inv_share_W(i,j)=CWB_W_dis(i,j)./Investment_Acc(i,j);
     Dis_inv_share_W_net(i,j)=CWB_W_dis_net(i,j)./Investment_net_Acc(i,j);
     Dis_inv_share_std_W(i,j)=CWBstd_W_dis(i,j)./Investment_Acc(i,j);
     Dis_inv_share_std_W_net(i,j)=CWBstd_W_dis_net(i,j)./Investment_net_Acc(i,j);
    end
end
% From 1950 to 1974
    for i=1:number_countries
     Dis_inv_share_W_epo1(i,:)=CWB_W_dis_epo1(i,:)./Investment_Acc_epo1(i,:);
     Dis_inv_share_W_epo1_net(i,:)=CWB_W_dis_epo1_net(i,:)./Investment_net_Acc_epo1(i,:);
     Dis_inv_share_std_W_epo1(i,:)=CWBstd_W_dis_epo1(i,:)./Investment_Acc_epo1(i,:);
     Dis_inv_share_std_W_epo1_net(i,:)=CWBstd_W_dis_epo1_net(i,:)./Investment_net_Acc_epo1(i,:);
    end
% From 1974 to 1999
    for i=1:number_countries
     Dis_inv_share_W_epo2(i,:)=CWB_W_dis_epo2(i,:)./Investment_Acc_epo2(i,:);
     Dis_inv_share_W_epo2_net(i,:)=CWB_W_dis_epo2_net(i,:)./Investment_net_Acc_epo2(i,:);
     Dis_inv_share_std_W_epo2(i,:)=CWBstd_W_dis_epo2(i,:)./Investment_Acc_epo2(i,:);
     Dis_inv_share_std_W_epo2_net(i,:)=CWBstd_W_dis_epo2_net(i,:)./Investment_net_Acc_epo2(i,:);
    end
% From 1999 to 2018
    for i=1:number_countries
     Dis_inv_share_W_epo3(i,:)=CWB_W_dis_epo3(i,:)./Investment_Acc_epo3(i,:);
     Dis_inv_share_W_epo3_net(i,:)=CWB_W_dis_epo3_net(i,:)./Investment_net_Acc_epo3(i,:);
     Dis_inv_share_std_W_epo3(i,:)=CWBstd_W_dis_epo3(i,:)./Investment_Acc_epo3(i,:);
     Dis_inv_share_std_W_epo3_net(i,:)=CWBstd_W_dis_epo3_net(i,:)./Investment_net_Acc_epo3(i,:);
    end

% Mean (all 3 damage functions)

for j=1:number_years
    for i=1:number_countries
     SCC_PV_dis(i,j)=Emissions_dis(i,j).*SCCm(j)*C;
     SCCstd_PV_dis(i,j)=Emissions_dis(i,j).*SCCstd(j).*C;
     SCC_PV_dis_net(i,j)=Emissions_dis_net(i,j).*SCCm(j)*C;
     SCCstd_PV_dis_net(i,j)=Emissions_dis_net(i,j).*SCCstd(j)*C;
    end
end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     CWB_dis(i,j)=sum(SCC_PV_dis(i,j:j6));
     CWBstd_dis(i,j)=(sum(SCCstd_PV_dis(i,j:j6).^2)).^0.5;
     CWB_dis_net(i,j)=sum(SCC_PV_dis_net(i,j:j6));
     CWBstd_dis_net(i,j)=(sum(SCCstd_PV_dis_net(i,j:j6).^2)).^0.5;
    end
end
% From 1950 to 1974
    for i=1:number_countries
     CWB_dis_epo1(i,:)=sum(SCC_PV_dis(i,j1:j2));
     CWBstd_dis_epo1(i,:)=(sum(SCCstd_PV_dis(i,j1:j2).^2)).^0.5;
     CWB_dis_epo1_net(i,:)=sum(SCC_PV_dis_net(i,j1:j2));
     CWBstd_dis_epo1_net(i,:)=(sum(SCCstd_PV_dis_net(i,j1:j2).^2)).^0.5;
    end
% From 1974 to 1999
    for i=1:number_countries
     CWB_dis_epo2(i,:)=sum(SCC_PV_dis(i,j3:j4));
     CWBstd_dis_epo2(i,:)=(sum(SCCstd_PV_dis(i,j3:j4).^2)).^0.5;
     CWB_dis_epo2_net(i,:)=sum(SCC_PV_dis_net(i,j3:j4));
     CWBstd_dis_epo2_net(i,:)=(sum(SCCstd_PV_dis_net(i,j3:j4).^2)).^0.5;
    end
% From 1999 to 2018
    for i=1:number_countries
     CWB_dis_epo3(i,:)=sum(SCC_PV_dis(i,j5:j6));
     CWBstd_dis_epo3(i,:)=(sum(SCCstd_PV_dis(i,j5:j6).^2)).^0.5;
     CWB_dis_epo3_net(i,:)=sum(SCC_PV_dis_net(i,j5:j6));
     CWBstd_dis_epo3_net(i,:)=(sum(SCCstd_PV_dis_net(i,j5:j6).^2)).^0.5;
    end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     Dis_inv_share(i,j)=CWB_dis(i,j)./Investment_Acc(i,j);
     Dis_inv_share_std(i,j)=Dis_inv_share(i,j).*((CWBstd_dis(i,j)./CWB_dis(i,j))^2)^0.5;
     Dis_inv_share_net(i,j)=CWB_dis_net(i,j)./Investment_net_Acc(i,j);
     Dis_inv_share_std_net(i,j)=CWBstd_dis_net(i,j)./Investment_net_Acc(i,j);  
    end
end
% From 1950 to 1974
    for i=1:number_countries
     Dis_inv_share_epo1(i,:)=CWB_dis_epo1(i,:)./Investment_Acc_epo1(i,:);
     Dis_inv_share_std_epo1(i,:)=Dis_inv_share_epo1(i,:).*((CWBstd_dis_epo1(i,:)./CWB_dis_epo1(i,:)).^2).^0.5;
     Dis_inv_share_epo1_net(i,:)=CWB_dis_epo1_net(i,:)./Investment_net_Acc_epo1(i,:);
     Dis_inv_share_std_epo1_net(i,:)=CWBstd_dis_epo1_net(i,:)./Investment_net_Acc_epo1(i,:);
    end
% From 1974 to 1999
    for i=1:number_countries
     Dis_inv_share_epo2(i,:)=CWB_dis_epo2(i,:)./Investment_Acc_epo2(i,:);
     Dis_inv_share_std_epo2(i,:)=Dis_inv_share_epo2(i,:).*((CWBstd_dis_epo2(i,:)./CWB_dis_epo2(i,:)).^2).^0.5;
     Dis_inv_share_epo2_net(i,:)=CWB_dis_epo2_net(i,:)./Investment_net_Acc_epo2(i,:);
     Dis_inv_share_std_epo2_net(i,:)=CWBstd_dis_epo2_net(i,:)./Investment_net_Acc_epo2(i,:);
    end
% From 1999 to 2018
    for i=1:number_countries
     Dis_inv_share_epo3(i,:)=CWB_dis_epo3(i,:)./Investment_Acc_epo3(i,:);
     Dis_inv_share_std_epo3(i,:)=Dis_inv_share_epo3(i,:).*((CWBstd_dis_epo3(i,:)./CWB_dis_epo3(i,:)).^2).^0.5;
     Dis_inv_share_epo3_net(i,:)=CWB_dis_epo3_net(i,:)./Investment_net_Acc_epo3(i,:);
     Dis_inv_share_std_epo3_net(i,:)=CWBstd_dis_epo3_net(i,:)./Investment_net_Acc_epo3(i,:);
    end

% Drupp et al. SDR 

% Howard and Sterner

for j=1:number_years
    for i=1:number_countries
     SCC_PV_HS_dis_Drupp(i,j)=Emissions_dis(i,j).*SCCm_HS_Drupp(j)*C;
     SCC_PV_HS_dis_net_Drupp(i,j)=Emissions_dis_net(i,j).*SCCm_HS_Drupp(j)*C;
     SCCstd_PV_HS_dis_Drupp(i,j)=Emissions_dis(i,j).*SCCstd_HS_Drupp(j)*C;
     SCCstd_PV_HS_dis_net_Drupp(i,j)=Emissions_dis_net(i,j).*SCCstd_HS_Drupp(j)*C;
    end
end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     CWB_HS_dis_Drupp(i,j)=sum(SCC_PV_HS_dis_Drupp(i,j:j6));
     CWB_HS_dis_net_Drupp(i,j)=sum(SCC_PV_HS_dis_net_Drupp(i,j:j6));
     CWBstd_HS_dis_Drupp(i,j)=(sum(SCCstd_PV_HS_dis(i,j:j6).^2)).^0.5;
     CWBstd_HS_dis_net_Drupp(i,j)=(sum(SCCstd_PV_HS_dis_net(i,j:j6).^2)).^0.5;
    end
end
% From 1950 to 1974
    for i=1:number_countries
     CWB_HS_dis_epo1_Drupp(i,:)=sum(SCC_PV_HS_dis_Drupp(i,j1:j2));
     CWB_HS_dis_epo1_net_Drupp(i,:)=sum(SCC_PV_HS_dis_net_Drupp(i,j1:j2));
     CWBstd_HS_dis_epo1_Drupp(i,:)=(sum(SCCstd_PV_HS_dis(i,j1:j2).^2)).^0.5;
     CWBstd_HS_dis_epo1_net_Drupp(i,:)=(sum(SCCstd_PV_HS_dis_net(i,j1:j2).^2)).^0.5;
    end
% From 1974 to 1999
    for i=1:number_countries
     CWB_HS_dis_epo2_Drupp(i,:)=sum(SCC_PV_HS_dis_Drupp(i,j3:j4));
     CWB_HS_dis_epo2_net_Drupp(i,:)=sum(SCC_PV_HS_dis_net_Drupp(i,j3:j4));
     CWBstd_HS_dis_epo2_Drupp(i,:)=(sum(SCCstd_PV_HS_dis(i,j3:j4).^2)).^0.5;
     CWBstd_HS_dis_epo2_net_Drupp(i,:)=(sum(SCCstd_PV_HS_dis_net(i,j3:j4).^2)).^0.5;
    end
% From 1999 to 2018
    for i=1:number_countries
     CWB_HS_dis_epo3_Drupp(i,:)=sum(SCC_PV_HS_dis_Drupp(i,j5:j6));
     CWB_HS_dis_epo3_net_Drupp(i,:)=sum(SCC_PV_HS_dis_net_Drupp(i,j5:j6));
     CWBstd_HS_dis_epo3_Drupp(i,:)=(sum(SCCstd_PV_HS_dis(i,j5:j6).^2)).^0.5;
     CWBstd_HS_dis_epo3_net_Drupp(i,:)=(sum(SCCstd_PV_HS_dis_net(i,j5:j6).^2)).^0.5;
    end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     Dis_inv_share_HS_Drupp(i,j)=CWB_HS_dis_Drupp(i,j)./Investment_Acc_Drupp(i,j);
     Dis_inv_share_HS_net_Drupp(i,j)=CWB_HS_dis_net_Drupp(i,j)./Investment_net_Acc_Drupp(i,j);
     Dis_inv_share_std_HS_Drupp(i,j)=CWBstd_HS_dis_Drupp(i,j)./Investment_Acc_Drupp(i,j);
     Dis_inv_share_std_HS_net_Drupp(i,j)=CWBstd_HS_dis_net_Drupp(i,j)./Investment_net_Acc_Drupp(i,j);
    end
end
% From 1950 to 1974
    for i=1:number_countries
     Dis_inv_share_HS_epo1_Drupp(i,:)=CWB_HS_dis_epo1_Drupp(i,:)./Investment_Acc_epo1_Drupp(i,:);
     Dis_inv_share_HS_epo1_net_Drupp(i,:)=CWB_HS_dis_epo1_net_Drupp(i,:)./Investment_net_Acc_epo1_Drupp(i,:);
     Dis_inv_share_std_HS_epo1_Drupp(i,:)=CWBstd_HS_dis_epo1_Drupp(i,:)./Investment_Acc_epo1_Drupp(i,:);
     Dis_inv_share_std_HS_epo1_net_Drupp(i,:)=CWBstd_HS_dis_epo1_net_Drupp(i,:)./Investment_net_Acc_epo1_Drupp(i,:);
    end
% From 1974 to 1999
    for i=1:number_countries
     Dis_inv_share_HS_epo2_Drupp(i,:)=CWB_HS_dis_epo2_Drupp(i,:)./Investment_Acc_epo2_Drupp(i,:);
     Dis_inv_share_HS_epo2_net_Drupp(i,:)=CWB_HS_dis_epo2_net_Drupp(i,:)./Investment_net_Acc_epo2_Drupp(i,:);
     Dis_inv_share_std_HS_epo2_Drupp(i,:)=CWBstd_HS_dis_epo2_Drupp(i,:)./Investment_Acc_epo2_Drupp(i,:);
     Dis_inv_share_std_HS_epo2_net_Drupp(i,:)=CWBstd_HS_dis_epo2_net_Drupp(i,:)./Investment_net_Acc_epo2_Drupp(i,:);
    end
% From 1999 to 2018
    for i=1:number_countries
     Dis_inv_share_HS_epo3_Drupp(i,:)=CWB_HS_dis_epo3_Drupp(i,:)./Investment_Acc_epo3_Drupp(i,:);
     Dis_inv_share_HS_epo3_net_Drupp(i,:)=CWB_HS_dis_epo3_net_Drupp(i,:)./Investment_net_Acc_epo3_Drupp(i,:);
     Dis_inv_share_std_HS_epo3_Drupp(i,:)=CWBstd_HS_dis_epo3_Drupp(i,:)./Investment_Acc_epo3_Drupp(i,:);
     Dis_inv_share_std_HS_epo3_net_Drupp(i,:)=CWBstd_HS_dis_epo3_net_Drupp(i,:)./Investment_net_Acc_epo3_Drupp(i,:);
    end

% Nordhaus

for j=1:number_years
    for i=1:number_countries
     SCC_PV_N_dis_Drupp(i,j)=Emissions_dis(i,j).*SCCm_N_Drupp(j)*C;
     SCC_PV_N_dis_net_Drupp(i,j)=Emissions_dis_net(i,j).*SCCm_N_Drupp(j)*C;
     SCCstd_PV_N_dis_Drupp(i,j)=Emissions_dis(i,j).*SCCstd_N_Drupp(j)*C;
     SCCstd_PV_N_dis_net_Drupp(i,j)=Emissions_dis_net(i,j).*SCCstd_N_Drupp(j)*C;
    end
end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     CWB_N_dis_Drupp(i,j)=sum(SCC_PV_N_dis_Drupp(i,j:j6));
     CWB_N_dis_net_Drupp(i,j)=sum(SCC_PV_N_dis_net_Drupp(i,j:j6));
     CWBstd_N_dis_Drupp(i,j)=(sum(SCCstd_PV_N_dis_Drupp(i,j:j6).^2)).^0.5;
     CWBstd_N_dis_net_Drupp(i,j)=(sum(SCCstd_PV_N_dis_net_Drupp(i,j:j6).^2)).^0.5;
    end
end
% From 1950 to 1974
    for i=1:number_countries
     CWB_N_dis_epo1_Drupp(i,:)=sum(SCC_PV_N_dis_Drupp(i,j1:j2));
     CWB_N_dis_epo1_net_Drupp(i,:)=sum(SCC_PV_N_dis_net_Drupp(i,j1:j2));
     CWBstd_N_dis_epo1_Drupp(i,:)=(sum(SCCstd_PV_N_dis(i,j1:j2).^2)).^0.5;
     CWBstd_N_dis_epo1_net_Drupp(i,:)=(sum(SCCstd_PV_N_dis_net(i,j1:j2).^2)).^0.5;
    end
% From 1974 to 1999
    for i=1:number_countries
     CWB_N_dis_epo2_Drupp(i,:)=sum(SCC_PV_N_dis_Drupp(i,j3:j4));
     CWB_N_dis_epo2_net_Drupp(i,:)=sum(SCC_PV_N_dis_net_Drupp(i,j3:j4));
     CWBstd_N_dis_epo2_Drupp(i,:)=(sum(SCCstd_PV_N_dis(i,j3:j4).^2)).^0.5;
     CWBstd_N_dis_epo2_net_Drupp(i,:)=(sum(SCCstd_PV_N_dis_net(i,j3:j4).^2)).^0.5;
    end
% From 1999 to 2018
    for i=1:number_countries
     CWB_N_dis_epo3_Drupp(i,:)=sum(SCC_PV_N_dis_Drupp(i,j5:j6));
     CWB_N_dis_epo3_net_Drupp(i,:)=sum(SCC_PV_N_dis_net_Drupp(i,j5:j6));
     CWBstd_N_dis_epo3_Drupp(i,:)=(sum(SCCstd_PV_N_dis(i,j5:j6).^2)).^0.5;
     CWBstd_N_dis_epo3_net_Drupp(i,:)=(sum(SCCstd_PV_N_dis_net(i,j5:j6).^2)).^0.5;
    end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     Dis_inv_share_N_Drupp(i,j)=CWB_N_dis_Drupp(i,j)./Investment_Acc_Drupp(i,j);
     Dis_inv_share_N_net_Drupp(i,j)=CWB_N_dis_net_Drupp(i,j)./Investment_net_Acc_Drupp(i,j);
     Dis_inv_share_std_N_Drupp(i,j)=CWBstd_N_dis_Drupp(i,j)./Investment_Acc_Drupp(i,j);
     Dis_inv_share_std_N_net_Drupp(i,j)=CWBstd_N_dis_net_Drupp(i,j)./Investment_net_Acc_Drupp(i,j);

    end
end
% From 1950 to 1974
    for i=1:number_countries
     Dis_inv_share_N_epo1_Drupp(i,:)=CWB_N_dis_epo1_Drupp(i,:)./Investment_Acc_epo1_Drupp(i,:);
     Dis_inv_share_N_epo1_net_Drupp(i,:)=CWB_N_dis_epo1_net_Drupp(i,:)./Investment_net_Acc_epo1_Drupp(i,:);
     Dis_inv_share_std_N_epo1_Drupp(i,:)=CWBstd_N_dis_epo1_Drupp(i,:)./Investment_Acc_epo1_Drupp(i,:);
     Dis_inv_share_std_N_epo1_net_Drupp(i,:)=CWBstd_N_dis_epo1_net_Drupp(i,:)./Investment_net_Acc_epo1_Drupp(i,:);
    end
% From 1974 to 1999
    for i=1:number_countries
     Dis_inv_share_N_epo2_Drupp(i,:)=CWB_N_dis_epo2_Drupp(i,:)./Investment_Acc_epo2_Drupp(i,:);
     Dis_inv_share_N_epo2_net_Drupp(i,:)=CWB_N_dis_epo2_net_Drupp(i,:)./Investment_net_Acc_epo2_Drupp(i,:);
     Dis_inv_share_std_N_epo2_Drupp(i,:)=CWBstd_N_dis_epo2_Drupp(i,:)./Investment_Acc_epo2_Drupp(i,:);
     Dis_inv_share_std_N_epo2_net_Drupp(i,:)=CWBstd_N_dis_epo2_net_Drupp(i,:)./Investment_net_Acc_epo2_Drupp(i,:);
    end
% From 1999 to 2018
    for i=1:number_countries
     Dis_inv_share_N_epo3_Drupp(i,:)=CWB_N_dis_epo3_Drupp(i,:)./Investment_Acc_epo3_Drupp(i,:);
     Dis_inv_share_N_epo3_net_Drupp(i,:)=CWB_N_dis_epo3_net_Drupp(i,:)./Investment_net_Acc_epo3_Drupp(i,:);
     Dis_inv_share_std_N_epo3_Drupp(i,:)=CWBstd_N_dis_epo3_Drupp(i,:)./Investment_Acc_epo3_Drupp(i,:);
     Dis_inv_share_std_N_epo3_net_Drupp(i,:)=CWBstd_N_dis_epo3_net_Drupp(i,:)./Investment_net_Acc_epo3_Drupp(i,:);
    end

% Weitzman

for j=1:number_years
    for i=1:number_countries
     SCC_PV_W_dis_Drupp(i,j)=Emissions_dis(i,j).*SCCm_W_Drupp(j)*C;
     SCC_PV_W_dis_net_Drupp(i,j)=Emissions_dis_net(i,j).*SCCm_W_Drupp(j)*C;
     SCCstd_PV_W_dis_Drupp(i,j)=Emissions_dis(i,j).*SCCstd_W_Drupp(j)*C;
     SCCstd_PV_W_dis_net_Drupp(i,j)=Emissions_dis_net(i,j).*SCCstd_W_Drupp(j)*C;
    end
end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     CWB_W_dis_Drupp(i,j)=sum(SCC_PV_W_dis_Drupp(i,j:j6));
     CWB_W_dis_net_Drupp(i,j)=sum(SCC_PV_W_dis_net_Drupp(i,j:j6));
     CWBstd_W_dis_Drupp(i,j)=(sum(SCCstd_PV_W_dis_Drupp(i,j:j6).^2)).^0.5;
     CWBstd_W_dis_net_Drupp(i,j)=(sum(SCCstd_PV_W_dis_net_Drupp(i,j:j6).^2)).^0.5;

    end
end
% From 1950 to 1974
    for i=1:number_countries
     CWB_W_dis_epo1_Drupp(i,:)=sum(SCC_PV_W_dis_Drupp(i,j1:j2));
     CWB_W_dis_epo1_net_Drupp(i,:)=sum(SCC_PV_W_dis_net_Drupp(i,j1:j2));
     CWBstd_W_dis_epo1_Drupp(i,:)=(sum(SCCstd_PV_W_dis_Drupp(i,j1:j2).^2)).^0.5;
     CWBstd_W_dis_epo1_net_Drupp(i,:)=(sum(SCCstd_PV_W_dis_net_Drupp(i,j1:j2).^2)).^0.5;
    end
% From 1974 to 1999
    for i=1:number_countries
     CWB_W_dis_epo2_Drupp(i,:)=sum(SCC_PV_W_dis_Drupp(i,j3:j4));
     CWB_W_dis_epo2_net_Drupp(i,:)=sum(SCC_PV_W_dis_net_Drupp(i,j3:j4));
     CWBstd_W_dis_epo2_Drupp(i,:)=(sum(SCCstd_PV_W_dis_Drupp(i,j3:j4).^2)).^0.5;
     CWBstd_W_dis_epo2_net_Drupp(i,:)=(sum(SCCstd_PV_W_dis_net_Drupp(i,j3:j4).^2)).^0.5;
    end
% From 1999 to 2018
    for i=1:number_countries
     CWB_W_dis_epo3_Drupp(i,:)=sum(SCC_PV_W_dis_Drupp(i,j5:j6));
     CWB_W_dis_epo3_net_Drupp(i,:)=sum(SCC_PV_W_dis_net_Drupp(i,j5:j6));
     CWBstd_W_dis_epo3_Drupp(i,:)=(sum(SCCstd_PV_W_dis_Drupp(i,j5:j6).^2)).^0.5;
     CWBstd_W_dis_epo3_net_Drupp(i,:)=(sum(SCCstd_PV_W_dis_net_Drupp(i,j5:j6).^2)).^0.5;
    end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     Dis_inv_share_W_Drupp(i,j)=CWB_W_dis_Drupp(i,j)./Investment_Acc_Drupp(i,j);
     Dis_inv_share_W_net_Drupp(i,j)=CWB_W_dis_net_Drupp(i,j)./Investment_net_Acc_Drupp(i,j);
     Dis_inv_share_std_W_Drupp(i,j)=CWBstd_W_dis_Drupp(i,j)./Investment_Acc_Drupp(i,j);
     Dis_inv_share_std_W_net_Drupp(i,j)=CWBstd_W_dis_net_Drupp(i,j)./Investment_net_Acc_Drupp(i,j);
    end
end
% From 1950 to 1974
    for i=1:number_countries
     Dis_inv_share_W_epo1_Drupp(i,:)=CWB_W_dis_epo1_Drupp(i,:)./Investment_Acc_epo1_Drupp(i,:);
     Dis_inv_share_W_epo1_net_Drupp(i,:)=CWB_W_dis_epo1_net_Drupp(i,:)./Investment_net_Acc_epo1_Drupp(i,:);
     Dis_inv_share_std_W_epo1_Drupp(i,:)=CWBstd_W_dis_epo1_Drupp(i,:)./Investment_Acc_epo1_Drupp(i,:);
     Dis_inv_share_std_W_epo1_net_Drupp(i,:)=CWBstd_W_dis_epo1_net_Drupp(i,:)./Investment_net_Acc_epo1_Drupp(i,:);
    end
% From 1974 to 1999
    for i=1:number_countries
     Dis_inv_share_W_epo2_Drupp(i,:)=CWB_W_dis_epo2_Drupp(i,:)./Investment_Acc_epo2_Drupp(i,:);
     Dis_inv_share_W_epo2_net_Drupp(i,:)=CWB_W_dis_epo2_net_Drupp(i,:)./Investment_net_Acc_epo2_Drupp(i,:);
     Dis_inv_share_std_W_epo2_Drupp(i,:)=CWBstd_W_dis_epo2_Drupp(i,:)./Investment_Acc_epo2_Drupp(i,:);
     Dis_inv_share_std_W_epo2_net_Drupp(i,:)=CWBstd_W_dis_epo2_net_Drupp(i,:)./Investment_net_Acc_epo2_Drupp(i,:);
    end
% From 1999 to 2018
    for i=1:number_countries
     Dis_inv_share_W_epo3_Drupp(i,:)=CWB_W_dis_epo3_Drupp(i,:)./Investment_Acc_epo3_Drupp(i,:);
     Dis_inv_share_W_epo3_net_Drupp(i,:)=CWB_W_dis_epo3_net_Drupp(i,:)./Investment_net_Acc_epo3_Drupp(i,:);
     Dis_inv_share_std_W_epo3_Drupp(i,:)=CWBstd_W_dis_epo3_Drupp(i,:)./Investment_Acc_epo3_Drupp(i,:);
     Dis_inv_share_std_W_epo3_net_Drupp(i,:)=CWBstd_W_dis_epo3_net_Drupp(i,:)./Investment_net_Acc_epo3_Drupp(i,:);
    end

% Mean (all 3 damage functions)

for j=1:number_years
    for i=1:number_countries
     SCC_PV_dis_Drupp(i,j)=Emissions_dis(i,j).*SCCm_Drupp(j)*C;
     SCCstd_PV_dis_Drupp(i,j)=Emissions_dis(i,j).*SCCstd_Drupp(j).*C;
     SCC_PV_dis_net_Drupp(i,j)=Emissions_dis_net(i,j).*SCCm_Drupp(j)*C;
     SCCstd_PV_dis_net_Drupp(i,j)=Emissions_dis_net(i,j).*SCCstd_Drupp(j)*C;
    end
end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     CWB_dis_Drupp(i,j)=sum(SCC_PV_dis_Drupp(i,j:j6));
     CWBstd_dis_Drupp(i,j)=(sum(SCCstd_PV_dis_Drupp(i,j:j6).^2)).^0.5;
     CWB_dis_net_Drupp(i,j)=sum(SCC_PV_dis_net_Drupp(i,j:j6));
     CWBstd_dis_net_Drupp(i,j)=(sum(SCCstd_PV_dis_net_Drupp(i,j:j6).^2)).^0.5;
    end
end
% From 1950 to 1974
    for i=1:number_countries
     CWB_dis_epo1_Drupp(i,:)=sum(SCC_PV_dis_Drupp(i,j1:j2));
     CWBstd_dis_epo1_Drupp(i,:)=(sum(SCCstd_PV_dis_Drupp(i,j1:j2).^2)).^0.5;
     CWB_dis_epo1_net_Drupp(i,:)=sum(SCC_PV_dis_net_Drupp(i,j1:j2));
     CWBstd_dis_epo1_net_Drupp(i,:)=(sum(SCCstd_PV_dis_net_Drupp(i,j1:j2).^2)).^0.5;
    end
% From 1974 to 1999
    for i=1:number_countries
     CWB_dis_epo2_Drupp(i,:)=sum(SCC_PV_dis_Drupp(i,j3:j4));
     CWBstd_dis_epo2_Drupp(i,:)=(sum(SCCstd_PV_dis_Drupp(i,j3:j4).^2)).^0.5;
     CWB_dis_epo2_net_Drupp(i,:)=sum(SCC_PV_dis_net_Drupp(i,j3:j4));
     CWBstd_dis_epo2_net_Drupp(i,:)=(sum(SCCstd_PV_dis_net_Drupp(i,j3:j4).^2)).^0.5;
    end
% From 1999 to 2018
    for i=1:number_countries
     CWB_dis_epo3_Drupp(i,:)=sum(SCC_PV_dis_Drupp(i,j5:j6));
     CWBstd_dis_epo3_Drupp(i,:)=(sum(SCCstd_PV_dis_Drupp(i,j5:j6).^2)).^0.5;
     CWB_dis_epo3_net_Drupp(i,:)=sum(SCC_PV_dis_net_Drupp(i,j5:j6));
     CWBstd_dis_epo3_net_Drupp(i,:)=(sum(SCCstd_PV_dis_net_Drupp(i,j5:j6).^2)).^0.5;
    end

% From 1950 to 2018
for j=j1:j6
    for i=1:number_countries
     Dis_inv_share_Drupp(i,j)=CWB_dis_Drupp(i,j)./Investment_Acc_Drupp(i,j);
     Dis_inv_share_std_Drupp(i,j)=Dis_inv_share_Drupp(i,j).*((CWBstd_dis_Drupp(i,j)./CWB_dis_Drupp(i,j))^2)^0.5;
     Dis_inv_share_net_Drupp(i,j)=CWB_dis_net_Drupp(i,j)./Investment_net_Acc_Drupp(i,j);
     Dis_inv_share_std_net_Drupp(i,j)=CWBstd_dis_net_Drupp(i,j)./Investment_net_Acc_Drupp(i,j);

    end
end
% From 1950 to 1974
    for i=1:number_countries
     Dis_inv_share_epo1_Drupp(i,:)=CWB_dis_epo1_Drupp(i,:)./Investment_Acc_epo1_Drupp(i,:);
     Dis_inv_share_std_epo1_Drupp(i,:)=Dis_inv_share_epo1_Drupp(i,:).*((CWBstd_dis_epo1_Drupp(i,:)./CWB_dis_epo1_Drupp(i,:)).^2).^0.5;
     Dis_inv_share_epo1_net_Drupp(i,:)=CWB_dis_epo1_net_Drupp(i,:)./Investment_net_Acc_epo1_Drupp(i,:);
     Dis_inv_share_std_epo1_net_Drupp(i,:)=CWBstd_dis_epo1_net_Drupp(i,:)./Investment_net_Acc_epo1_Drupp(i,:);

    end
% From 1974 to 1999
    for i=1:number_countries
     Dis_inv_share_epo2_Drupp(i,:)=CWB_dis_epo2_Drupp(i,:)./Investment_Acc_epo2_Drupp(i,:);
     Dis_inv_share_std_epo2_Drupp(i,:)=Dis_inv_share_epo2_Drupp(i,:).*((CWBstd_dis_epo2_Drupp(i,:)./CWB_dis_epo2_Drupp(i,:)).^2).^0.5;
     Dis_inv_share_epo2_net_Drupp(i,:)=CWB_dis_epo2_net_Drupp(i,:)./Investment_net_Acc_epo2_Drupp(i,:);
     Dis_inv_share_std_epo2_net_Drupp(i,:)=CWBstd_dis_epo2_net_Drupp(i,:)./Investment_net_Acc_epo2_Drupp(i,:);
    end
% From 1999 to 2018
    for i=1:number_countries
     Dis_inv_share_epo3_Drupp(i,:)=CWB_dis_epo3_Drupp(i,:)./Investment_Acc_epo3_Drupp(i,:);
     Dis_inv_share_std_epo3_Drupp(i,:)=Dis_inv_share_epo3_Drupp(i,:).*((CWBstd_dis_epo3_Drupp(i,:)./CWB_dis_epo3_Drupp(i,:)).^2).^0.5;
     Dis_inv_share_epo3_net_Drupp(i,:)=CWB_dis_epo3_net_Drupp(i,:)./Investment_net_Acc_epo3_Drupp(i,:);
     Dis_inv_share_std_epo3_net_Drupp(i,:)=CWBstd_dis_epo3_net_Drupp(i,:)./Investment_net_Acc_epo3_Drupp(i,:);
    end
    
 %% Climate Wealth Borrowing per capita (only for mean)
 
 % Calibrated SDR
 
 CWB_per_capita=CWB./Population(:,69);
 CWBstd_per_capita=CWBstd./Population(:,69);
 
 % Drupp et al. SDR 
 
  %for j=j1:j6
  % for i=1:number_countries
     CWB_per_capita_Drupp=CWB_Drupp./Population(:,69);
     CWBstd_per_capita_Drupp=CWBstd_Drupp./Population(:,69);
  %  end
  %end
 
  CWBci_per_capita=CWBstd_per_capita*2.11/sqrt(18);
  CWBci_per_capita_Drupp=CWBstd_per_capita_Drupp*2.11/sqrt(18);
  %% Results
  
%cd /Users/cocco/Desktop/Matlab
cd C:\..\IW_calculations % Choose the file path 
  
CWBresults=[Ecum, CWB(:,1), CWBstd(:,1)*ci, CWB_N(:,1), CWBstd_N(:,1)*ci, CWB_W(:,1), CWBstd_W(:,1)*ci, CWB_HS(:,1), CWBstd_HS(:,1)*ci, CWB_Drupp(:,1), CWBstd_Drupp(:,1)*ci, CWB_N_Drupp(:,1), CWBstd_N_Drupp(:,1)*ci, CWB_W_Drupp(:,1), CWBstd_W_Drupp(:,1)*ci, CWB_HS_Drupp(:,1), CWBstd_HS_Drupp(:,1)*ci];
%writematrix(CWBresults,'CWBresults.csv');

CWBcapita=[CWB_per_capita,CWBci_per_capita,CWB_per_capita_Drupp,CWBci_per_capita_Drupp];

colNames1={'All_impact_calib_m','All_impact_calib_ci','N_calib_m','N_calib_ci','W_calib_m','W_calib_ci','HS_calib_m','HS_calib_ci','All_impact_expert_m','All_impact_expert_ci','N_expert_m','N_expert_ci','W_expert_m','W_expert_ci','HS_expert_m','HS_expert_ci'};
colNames2={'CO2','All_impact_calib_m','All_impact_calib_ci','N_calib_m','N_calib_ci','W_calib_m','W_calib_ci','HS_calib_m','HS_calib_ci','All_impact_expert_m','All_impact_expert_ci','N_expert_m','N_expert_ci','W_expert_m','W_expert_ci','HS_expert_m','HS_expert_ci'};
colNames3={'Mean_calibrated','CI_calib','Mean_expert','CI_expert'};
rowNames=countries_label(:,1);
CWBresults_table = array2table(CWBresults,'RowNames',rowNames,'VariableNames',colNames2);
writetable(CWBresults_table,'CWBcountry_results_table.csv','WriteRowNames',true,'WriteVariableNames',true);
CWBcapita_table = array2table(CWBcapita,'RowNames',rowNames,'VariableNames',colNames3);
writetable(CWBcapita_table,'CWBcountry_capita_table.csv','WriteRowNames',true,'WriteVariableNames',true);

%CWBresults2=[Ecum, CWB(:,1), Dis_inv_share(:,1), Dis_inv_share_epo1(:,1), Dis_inv_share_epo2(:,1), Dis_inv_share_epo3(:,1), Dis_inv_share_epo1_net(:,1), Dis_inv_share_epo2_net(:,1), Dis_inv_share_epo3_net(:,1), CWB_Drupp(:,1), Dis_inv_share_Drupp(:,1), Dis_inv_share_epo1_Drupp(:,1), Dis_inv_share_epo2_Drupp(:,1), Dis_inv_share_epo3_Drupp(:,1), Dis_inv_share_epo1_net_Drupp(:,1), Dis_inv_share_epo2_net_Drupp(:,1), Dis_inv_share_epo3_net_Drupp(:,1)];
%writematrix(CWBresults2,'CWBresults2.csv')

%CWB_share=[Dis_inv_share(:,1),Dis_inv_share_std(:,1),Dis_inv_share_N(:,1), Dis_inv_share_std_N(:,1),Dis_inv_share_W(:,1),Dis_inv_share_std_W(:,1),Dis_inv_share_HS(:,1),Dis_inv_share_std_HS(:,1),Dis_inv_share_Drupp(:,1),Dis_inv_share_std_Drupp(:,1),Dis_inv_share_N_Drupp(:,1), Dis_inv_share_std_N_Drupp(:,1),Dis_inv_share_W_Drupp(:,1),Dis_inv_share_std_W_Drupp(:,1),Dis_inv_share_HS_Drupp(:,1),Dis_inv_share_std_HS(:,1)_Drupp];
%CWB_netshare=[Dis_inv_share_net(:,1),Dis_inv_share_std_net(:,1),Dis_inv_share_N_net(:,1),Dis_inv_share_std_N_net(:,1),Dis_inv_share_W_net(:,1),Dis_inv_share_std_W_net(:,1),Dis_inv_share_HS_net(:,1),Dis_inv_share_std_HS_net,Dis_inv_share_net_Drupp(:,1),Dis_inv_share_std_net_Drupp(:,1),Dis_inv_share_N_net_Drupp(:,1),Dis_inv_share_std_N_net_Drupp(:,1),Dis_inv_share_W_net_Drupp(:,1),Dis_inv_share_std_W_net_Drupp(:,1),Dis_inv_share_HS_net_Drupp(:,1),Dis_inv_share_std_HS_net_Drupp(:,1)]
%CWB_share_epo1=[Dis_inv_share_epo1(:,1), Dis_inv_share_std_epo1(:,1),Dis_inv_share_N_epo1(:,1), Dis_inv_share_std_N_epo1(:,1),Dis_inv_share_W_epo1(:,1), Dis_inv_share_std_W_epo1(:,1),Dis_inv_share_HS_epo1(:,1),Dis_inv_share_std_HS_epo1(:,1),Dis_inv_share_epo1_Drupp(:,1), Dis_inv_share_std_epo1_Drupp(:,1),Dis_inv_share_N_epo1_Drupp(:,1), Dis_inv_share_std_N_epo1_Drupp(:,1),Dis_inv_share_W_epo1_Drupp(:,1), Dis_inv_share_std_W_epo1_Drupp(:,1),Dis_inv_share_HS_epo1_Drupp(:,1),Dis_inv_share_std_HS_epo1_Drupp(:,1)];
%CWB_netshare_epo1=[Dis_inv_share_epo1_net(:,1),Dis_inv_share_std_epo1_net(:,1),Dis_inv_share_N_epo1_net(:,1), Dis_inv_share_std_N_epo1_net(:,1),Dis_inv_share_W_epo1_net(:,1), Dis_inv_share_std_W_epo1_net(:,1), Dis_inv_share_HS_epo1_net(:,1), Dis_inv_share_HS_std_epo1_net(:,1),Dis_inv_share_epo1_net_Drupp(:,1),Dis_inv_share_std_epo1_net_Drupp(:,1),Dis_inv_share_N_epo1_net_Drupp(:,1), Dis_inv_share_std_N_epo1_net_Drupp(:,1),Dis_inv_share_W_epo1_net_Drupp(:,1), Dis_inv_share_std_W_epo1_net_Drupp(:,1),Dis_inv_share_HS_epo1_net_Drupp(:,1), Dis_inv_share_HS_std_epo1_net_Drupp(:,1)];
%CWB_share_epo2=[Dis_inv_share_epo2(:,1),Dis_inv_share_std_epo2(:,1),Dis_inv_share_N_epo2(:,1), Dis_inv_share_std_N_epo2(:,1),Dis_inv_share_W_epo2(:,1), Dis_inv_share_std_W_epo2(:,1) ,Dis_inv_share_HS_epo2(:,1), Dis_inv_share_std_HS_epo2(:,1),Dis_inv_share_epo2_Drupp(:,1),Dis_inv_share_std_epo2_Drupp(:,1),Dis_inv_share_N_epo2_Drupp(:,1),Dis_inv_share_std_N_epo2_Drupp(:,1),Dis_inv_share_W_epo2_Drupp(:,1),Dis_inv_share_std_W_epo2_Drupp(:,1),Dis_inv_share_HS_epo2_Drupp(:,1), Dis_inv_share_std_HS_epo2_Drupp(:,1)];
%CWB_netshare_epo2=[Dis_inv_share_epo2_net(:,1),Dis_inv_share_std_epo2_net(:,1),Dis_inv_share_N_epo2_net(:,1),Dis_inv_share_std_N_epo2_net(:,1),Dis_inv_share_W_epo2_net(:,1), Dis_inv_share_std_W_epo2_net(:,1),Dis_inv_share_HS_epo2_net(:,1), Dis_inv_share_HS_std_epo2_net(:,1),Dis_inv_share_epo2_net_Drupp(:,1),Dis_inv_share_std_epo2_net_Drupp(:,1),Dis_inv_share_N_epo2_net_Drupp(:,1),Dis_inv_share_std_N_epo2_net_Drupp(:,1),Dis_inv_share_W_epo2_net_Drupp(:,1), Dis_inv_share_std_W_epo2_net_Drupp(:,1),Dis_inv_share_HS_epo2_net_Drupp(:,1),Dis_inv_share_HS_std_epo2_net(:,1)_Drupp];
%CWB_share_epo3=[Dis_inv_share_epo3(:,1),Dis_inv_share_std_epo3(:,1),Dis_inv_share_N_epo3(:,1), Dis_inv_share_std_N_epo3(:,1),Dis_inv_share_W_epo3(:,1),Dis_inv_share_std_W_epo3(:,1),Dis_inv_share_HS_epo3(:,1), Dis_inv_share_std_HS_epo3(:,1),Dis_inv_share_epo3_Drupp(:,1),Dis_inv_share_std_epo3_Drupp(:,1),Dis_inv_share_N_epo3_Drupp(:,1),Dis_inv_share_std_N_epo3_Drupp(:,1),Dis_inv_share_W_epo3_Drupp(:,1),Dis_inv_share_std_W_epo3_Drupp(:,1),Dis_inv_share_HS_epo3_Drupp(:,1), Dis_inv_share_std_HS_epo3_Drupp(:,1)];
%CWB_netshare_epo3=[Dis_inv_share_epo3_net(:,1),Dis_inv_share_std_epo3_net(:,1),Dis_inv_share_N_epo3_net(:,1), Dis_inv_share_std_N_epo3_net(:,1),Dis_inv_share_W_epo3_net(:,1), Dis_inv_share_std_W_epo3_net(:,1), Dis_inv_share_HS_epo3_net(:,1),Dis_inv_share_std_HS_epo3_net(:,1),Dis_inv_share_epo3_net_Drupp(:,1),Dis_inv_share_std_epo3_net_Drupp(:,1),Dis_inv_share_N_epo3_net_Drupp(:,1), Dis_inv_share_std_N_epo3_net_Drupp(:,1),Dis_inv_share_W_epo3_net_Drupp(:,1), Dis_inv_share_std_W_epo3_net_Drupp(:,1),Dis_inv_share_HS_epo3_net_Drupp(:,1),Dis_inv_share_std_HS_epo3_net_Drupp(:,1)];

CWB_share_ci=[Dis_inv_share(:,1),Dis_inv_share_std(:,1)*ci,Dis_inv_share_N(:,1),Dis_inv_share_std_N(:,1)*ci,Dis_inv_share_W(:,1),Dis_inv_share_std_W(:,1)*ci,Dis_inv_share_HS(:,1),Dis_inv_share_std_HS(:,1),Dis_inv_share_Drupp(:,1),Dis_inv_share_std_Drupp(:,1),Dis_inv_share_N_Drupp(:,1),Dis_inv_share_std_N_Drupp(:,1),Dis_inv_share_W_Drupp(:,1),Dis_inv_share_std_W_Drupp(:,1),Dis_inv_share_HS_Drupp(:,1),Dis_inv_share_std_HS_Drupp(:,1)];
CWB_netshare_ci=[Dis_inv_share_net(:,1),Dis_inv_share_std_net(:,1),Dis_inv_share_N_net(:,1),Dis_inv_share_std_N_net(:,1),Dis_inv_share_W_net(:,1),Dis_inv_share_std_W_net(:,1),Dis_inv_share_HS_net(:,1),Dis_inv_share_std_HS_net,Dis_inv_share_net_Drupp(:,1),Dis_inv_share_std_net_Drupp(:,1),Dis_inv_share_N_net_Drupp(:,1),Dis_inv_share_std_N_net_Drupp(:,1),Dis_inv_share_W_net_Drupp(:,1),Dis_inv_share_std_W_net_Drupp(:,1),Dis_inv_share_HS_net_Drupp(:,1),Dis_inv_share_std_HS_net_Drupp(:,1)]
CWB_share_epo1_ci=[Dis_inv_share_epo1(:,1), Dis_inv_share_std_epo1(:,1),Dis_inv_share_N_epo1(:,1), Dis_inv_share_std_N_epo1(:,1),Dis_inv_share_W_epo1(:,1), Dis_inv_share_std_W_epo1(:,1),Dis_inv_share_HS_epo1(:,1),Dis_inv_share_std_HS_epo1(:,1),Dis_inv_share_epo1_Drupp(:,1), Dis_inv_share_std_epo1_Drupp(:,1),Dis_inv_share_N_epo1_Drupp(:,1), Dis_inv_share_std_N_epo1_Drupp(:,1),Dis_inv_share_W_epo1_Drupp(:,1), Dis_inv_share_std_W_epo1_Drupp(:,1),Dis_inv_share_HS_epo1_Drupp(:,1),Dis_inv_share_std_HS_epo1_Drupp(:,1)];
CWB_netshare_epo1_ci=[Dis_inv_share_epo1_net(:,1),Dis_inv_share_std_epo1_net(:,1),Dis_inv_share_N_epo1_net(:,1), Dis_inv_share_std_N_epo1_net(:,1),Dis_inv_share_W_epo1_net(:,1), Dis_inv_share_std_W_epo1_net(:,1), Dis_inv_share_HS_epo1_net(:,1), Dis_inv_share_std_HS_epo1_net(:,1),Dis_inv_share_epo1_net_Drupp(:,1),Dis_inv_share_std_epo1_net_Drupp(:,1),Dis_inv_share_N_epo1_net_Drupp(:,1), Dis_inv_share_std_N_epo1_net_Drupp(:,1),Dis_inv_share_W_epo1_net_Drupp(:,1), Dis_inv_share_std_W_epo1_net_Drupp(:,1),Dis_inv_share_HS_epo1_net_Drupp(:,1), Dis_inv_share_std_HS_epo1_net_Drupp(:,1)];
CWB_share_epo2_ci=[Dis_inv_share_epo2(:,1),Dis_inv_share_std_epo2(:,1),Dis_inv_share_N_epo2(:,1), Dis_inv_share_std_N_epo2(:,1),Dis_inv_share_W_epo2(:,1), Dis_inv_share_std_W_epo2(:,1) ,Dis_inv_share_HS_epo2(:,1), Dis_inv_share_std_HS_epo2(:,1),Dis_inv_share_epo2_Drupp(:,1),Dis_inv_share_std_epo2_Drupp(:,1),Dis_inv_share_N_epo2_Drupp(:,1),Dis_inv_share_std_N_epo2_Drupp(:,1),Dis_inv_share_W_epo2_Drupp(:,1),Dis_inv_share_std_W_epo2_Drupp(:,1),Dis_inv_share_HS_epo2_Drupp(:,1), Dis_inv_share_std_HS_epo2_Drupp(:,1)];
CWB_netshare_epo2_ci=[Dis_inv_share_epo2_net(:,1),Dis_inv_share_std_epo2_net(:,1),Dis_inv_share_N_epo2_net(:,1),Dis_inv_share_std_N_epo2_net(:,1),Dis_inv_share_W_epo2_net(:,1), Dis_inv_share_std_W_epo2_net(:,1),Dis_inv_share_HS_epo2_net(:,1), Dis_inv_share_std_HS_epo2_net(:,1),Dis_inv_share_epo2_net_Drupp(:,1),Dis_inv_share_std_epo2_net_Drupp(:,1),Dis_inv_share_N_epo2_net_Drupp(:,1),Dis_inv_share_std_N_epo2_net_Drupp(:,1),Dis_inv_share_W_epo2_net_Drupp(:,1), Dis_inv_share_std_W_epo2_net_Drupp(:,1),Dis_inv_share_HS_epo2_net_Drupp(:,1),Dis_inv_share_std_HS_epo2_net_Drupp(:,1)];
CWB_share_epo3_ci=[Dis_inv_share_epo3(:,1),Dis_inv_share_std_epo3(:,1),Dis_inv_share_N_epo3(:,1), Dis_inv_share_std_N_epo3(:,1),Dis_inv_share_W_epo3(:,1),Dis_inv_share_std_W_epo3(:,1),Dis_inv_share_HS_epo3(:,1), Dis_inv_share_std_HS_epo3(:,1),Dis_inv_share_epo3_Drupp(:,1),Dis_inv_share_std_epo3_Drupp(:,1),Dis_inv_share_N_epo3_Drupp(:,1),Dis_inv_share_std_N_epo3_Drupp(:,1),Dis_inv_share_W_epo3_Drupp(:,1),Dis_inv_share_std_W_epo3_Drupp(:,1),Dis_inv_share_HS_epo3_Drupp(:,1), Dis_inv_share_std_HS_epo3_Drupp(:,1)];
CWB_netshare_epo3_ci=[Dis_inv_share_epo3_net(:,1),Dis_inv_share_std_epo3_net(:,1),Dis_inv_share_N_epo3_net(:,1), Dis_inv_share_std_N_epo3_net(:,1),Dis_inv_share_W_epo3_net(:,1), Dis_inv_share_std_W_epo3_net(:,1), Dis_inv_share_HS_epo3_net(:,1),Dis_inv_share_std_HS_epo3_net(:,1),Dis_inv_share_epo3_net_Drupp(:,1),Dis_inv_share_std_epo3_net_Drupp(:,1),Dis_inv_share_N_epo3_net_Drupp(:,1), Dis_inv_share_std_N_epo3_net_Drupp(:,1),Dis_inv_share_W_epo3_net_Drupp(:,1), Dis_inv_share_std_W_epo3_net_Drupp(:,1),Dis_inv_share_HS_epo3_net_Drupp(:,1),Dis_inv_share_std_HS_epo3_net_Drupp(:,1)];

CWB_share_ci=100*[Dis_inv_share(:,1),Dis_inv_share_std(:,1)*ci,Dis_inv_share_N(:,1), Dis_inv_share_std_N(:,1)*ci,Dis_inv_share_W(:,1),Dis_inv_share_std_W(:,1)*ci,Dis_inv_share_HS(:,1),Dis_inv_share_std_HS(:,1)*ci,Dis_inv_share_Drupp(:,1),Dis_inv_share_std_Drupp(:,1)*ci,Dis_inv_share_N_Drupp(:,1),Dis_inv_share_std_N_Drupp(:,1)*ci,Dis_inv_share_W_Drupp(:,1),Dis_inv_share_std_W_Drupp(:,1)*ci,Dis_inv_share_HS_Drupp(:,1),Dis_inv_share_std_HS_Drupp(:,1)*ci];
CWB_netshare_ci=100*[Dis_inv_share_net(:,1),Dis_inv_share_std_net(:,1)*ci,Dis_inv_share_N_net(:,1),Dis_inv_share_std_N_net(:,1)*ci,Dis_inv_share_W_net(:,1)*ci,Dis_inv_share_std_W_net(:,1)*ci,Dis_inv_share_HS_net(:,1),Dis_inv_share_std_HS_net(:,1)*ci,Dis_inv_share_net_Drupp(:,1),Dis_inv_share_std_net_Drupp(:,1)*ci,Dis_inv_share_N_net_Drupp(:,1),Dis_inv_share_std_N_net_Drupp(:,1)*ci,Dis_inv_share_W_net_Drupp(:,1),Dis_inv_share_std_W_net_Drupp(:,1)*ci,Dis_inv_share_HS_net_Drupp(:,1),Dis_inv_share_std_HS_net_Drupp(:,1)*ci];
CWB_share_epo1_ci=100*[Dis_inv_share_epo1(:,1),Dis_inv_share_std_epo1(:,1)*ci,Dis_inv_share_N_epo1(:,1), Dis_inv_share_std_N_epo1(:,1)*ci,Dis_inv_share_W_epo1(:,1),Dis_inv_share_std_W_epo1(:,1)*ci,Dis_inv_share_HS_epo1(:,1),Dis_inv_share_std_HS_epo1(:,1)*ci,Dis_inv_share_epo1_Drupp(:,1), Dis_inv_share_std_epo1_Drupp(:,1)*ci,Dis_inv_share_N_epo1_Drupp(:,1),Dis_inv_share_std_N_epo1_Drupp(:,1)*ci,Dis_inv_share_W_epo1_Drupp(:,1), Dis_inv_share_std_W_epo1_Drupp(:,1)*ci,Dis_inv_share_HS_epo1_Drupp(:,1),Dis_inv_share_std_HS_epo1_Drupp(:,1)*ci];
CWB_netshare_epo1_ci=100*[Dis_inv_share_epo1_net(:,1),Dis_inv_share_std_epo1_net(:,1)*ci,Dis_inv_share_N_epo1_net(:,1), Dis_inv_share_std_N_epo1_net(:,1)*ci,Dis_inv_share_W_epo1_net(:,1),Dis_inv_share_std_W_epo1_net(:,1)*ci, Dis_inv_share_HS_epo1_net(:,1), Dis_inv_share_std_HS_epo1_net(:,1)*ci,Dis_inv_share_epo1_net_Drupp(:,1),Dis_inv_share_std_epo1_net_Drupp(:,1)*ci,Dis_inv_share_N_epo1_net_Drupp(:,1), Dis_inv_share_std_N_epo1_net_Drupp(:,1)*ci,Dis_inv_share_W_epo1_net_Drupp(:,1), Dis_inv_share_std_W_epo1_net_Drupp(:,1)*ci,Dis_inv_share_HS_epo1_net_Drupp(:,1),Dis_inv_share_std_HS_epo1_net_Drupp(:,1)*ci];
CWB_share_epo2_ci=100*[Dis_inv_share_epo2(:,1),Dis_inv_share_std_epo2(:,1)*ci,Dis_inv_share_N_epo2(:,1),Dis_inv_share_std_N_epo2(:,1)*ci,Dis_inv_share_W_epo2(:,1),Dis_inv_share_std_W_epo2(:,1)*ci,Dis_inv_share_HS_epo2(:,1), Dis_inv_share_std_HS_epo2(:,1)*ci,Dis_inv_share_epo2_Drupp(:,1),Dis_inv_share_std_epo2_Drupp(:,1)*ci,Dis_inv_share_N_epo2_Drupp(:,1),Dis_inv_share_std_N_epo2_Drupp(:,1)*ci,Dis_inv_share_W_epo2_Drupp(:,1),Dis_inv_share_std_W_epo2_Drupp(:,1)*ci,Dis_inv_share_HS_epo2_Drupp(:,1), Dis_inv_share_std_HS_epo2_Drupp(:,1)*ci];
CWB_netshare_epo2_ci=100*[Dis_inv_share_epo2_net(:,1),Dis_inv_share_std_epo2_net(:,1)*ci,Dis_inv_share_N_epo2_net(:,1),Dis_inv_share_std_N_epo2_net(:,1)*ci,Dis_inv_share_W_epo2_net(:,1),Dis_inv_share_std_W_epo2_net(:,1)*ci,Dis_inv_share_HS_epo2_net(:,1),Dis_inv_share_std_HS_epo2_net(:,1)*ci,Dis_inv_share_epo2_net_Drupp(:,1),Dis_inv_share_std_epo2_net_Drupp(:,1)*ci,Dis_inv_share_N_epo2_net_Drupp(:,1),Dis_inv_share_std_N_epo2_net_Drupp(:,1)*ci,Dis_inv_share_W_epo2_net_Drupp(:,1),Dis_inv_share_std_W_epo2_net_Drupp(:,1)*ci,Dis_inv_share_HS_epo2_net_Drupp(:,1),Dis_inv_share_std_HS_epo2_net_Drupp(:,1)*ci];
CWB_share_epo3_ci=100*[Dis_inv_share_epo3(:,1),Dis_inv_share_std_epo3(:,1)*ci,Dis_inv_share_N_epo3(:,1),Dis_inv_share_std_N_epo3(:,1)*ci,Dis_inv_share_W_epo3(:,1),Dis_inv_share_std_W_epo3(:,1)*ci,Dis_inv_share_HS_epo3(:,1), Dis_inv_share_std_HS_epo3(:,1)*ci,Dis_inv_share_epo3_Drupp(:,1),Dis_inv_share_std_epo3_Drupp(:,1)*ci,Dis_inv_share_N_epo3_Drupp(:,1),Dis_inv_share_std_N_epo3_Drupp(:,1)*ci,Dis_inv_share_W_epo3_Drupp(:,1),Dis_inv_share_std_W_epo3_Drupp(:,1)*ci,Dis_inv_share_HS_epo3_Drupp(:,1),Dis_inv_share_std_HS_epo3_Drupp(:,1)*ci];
CWB_netshare_epo3_ci=100*[Dis_inv_share_epo3_net(:,1),Dis_inv_share_std_epo3_net(:,1)*ci,Dis_inv_share_N_epo3_net(:,1),Dis_inv_share_std_N_epo3_net(:,1)*ci,Dis_inv_share_W_epo3_net(:,1),Dis_inv_share_std_W_epo3_net(:,1)*ci, Dis_inv_share_HS_epo3_net(:,1),Dis_inv_share_std_HS_epo3_net(:,1)*ci,Dis_inv_share_epo3_net_Drupp(:,1),Dis_inv_share_std_epo3_net_Drupp(:,1)*ci,Dis_inv_share_N_epo3_net_Drupp(:,1),Dis_inv_share_std_N_epo3_net_Drupp(:,1)*ci,Dis_inv_share_W_epo3_net_Drupp(:,1), Dis_inv_share_std_W_epo3_net_Drupp(:,1)*ci,Dis_inv_share_HS_epo3_net_Drupp(:,1),Dis_inv_share_std_HS_epo3_net_Drupp(:,1)*ci];

%colNames="y_"+(1950:2018);
colNames={'All_impact_calib_m','All_impact_calib_ci','N_calib_m','N_calib_ci','W_calib_m','W_calib_ci','HS_calib_m','HS_calib_ci','All_impact_expert_m','All_impact_expert_ci','N_expert_m','N_expert_ci','W_expert_m','W_expert_ci','HS_expert_m','HS_expert_ci'};
%rfix = array2table(output2,'RowNames',rowNames,'VariableNames',colNames);
%CWB_table = array2table(CWBresults,'RowNames',rowNames,'VariableNames',colNames);
CWBshare_table = array2table(CWB_share_ci,'RowNames',rowNames,'VariableNames',colNames);
CWBnetshare_table = array2table(CWB_netshare_ci,'RowNames',rowNames,'VariableNames',colNames);
CWBshareepo1_table = array2table(CWB_share_epo1_ci,'RowNames',rowNames,'VariableNames',colNames);
CWBnetshareepo1_table = array2table(CWB_netshare_epo1_ci,'RowNames',rowNames,'VariableNames',colNames);
CWBshareepo2_table = array2table(CWB_share_epo2_ci,'RowNames',rowNames,'VariableNames',colNames);
CWBnetshareepo2_table = array2table(CWB_netshare_epo2_ci,'RowNames',rowNames,'VariableNames',colNames);
CWBshareepo3_table = array2table(CWB_share_epo3_ci,'RowNames',rowNames,'VariableNames',colNames);
CWBnetshareepo3_table = array2table(CWB_netshare_epo3_ci,'RowNames',rowNames,'VariableNames',colNames);

%writetable(CWB_table,'CWB_results.csv','WriteRowNames',true,'WriteVariableNames',true)
writetable(CWBshare_table,'CWBcountryshare_table.csv','WriteRowNames',true,'WriteVariableNames',true)
writetable(CWBnetshare_table,'CWBcountrynetshare_table.csv','WriteRowNames',true,'WriteVariableNames',true)
writetable(CWBshareepo1_table,'CWBcountryshareepo1_table.csv','WriteRowNames',true,'WriteVariableNames',true)
writetable(CWBnetshareepo1_table,'CWBcountrynetshareepo1_table.csv','WriteRowNames',true,'WriteVariableNames',true)
writetable(CWBshareepo2_table,'CWBcountryshareepo2_table.csv','WriteRowNames',true,'WriteVariableNames',true)
writetable(CWBnetshareepo2_table,'CWBcountrynetshareepo2_table.csv','WriteRowNames',true,'WriteVariableNames',true)
writetable(CWBshareepo3_table,'CWBcountryshareepo3_table.csv','WriteRowNames',true,'WriteVariableNames',true)
writetable(CWBnetshareepo3_table,'CWBcountrynetshareepo3_table.csv','WriteRowNames',true,'WriteVariableNames',true)

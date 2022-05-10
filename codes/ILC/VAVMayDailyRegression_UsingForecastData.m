clear all
close all
swEPSfigure
load Data4RegressionDaily2022.mat
%load VAVdataMay2014.mat
%VAVdata = VAV201TO216301TO316May2014r;

%Number of VAVs
VSign = 32;

%Define the order of the state and control
m = 1;
p = 1;

%Left-hand-side of regression
xt = Trm;

%Right-hand-side of regression
%Time-delay states
xt1 = xt; %(n-1) to (n-m)
xt0 = []; % for order iteration
%Time-delay controls
u = Tdisch.*flow;
u1 = u; %(n-1) to (n-p)
u0 = [];  % for order iteration
%Disturbance - Forecast
w_forecast = Toa_forecast;
w = Toa;

for j = 1 : m

    xt0 = [xt1 xt0];
    if (j <= p)
        u0 = [u1 u0];
    end
% delete extra measurements to compenstate the time sequence difference
% among inputs and output
    NIndex = (NumDayPartition-(j-1)) * linspace(0, length(xt)/(NumDayPartition-(j-1))-1,length(xt)/(NumDayPartition-(j-1)));
    AHUdataDay(NIndex'+1,:)=[]; VAVdataDay(NIndex'+1,:)=[];
    xt(NIndex'+1,:)=[]; %[n] discard first sample everyday for measurement
    w(NIndex'+1,:)=[]; w_forecast(NIndex'+1,:)=[];   
    xt0(NIndex'+NumDayPartition-(j-1),:)=[]; %[n-j] discard last sample everyday for Y1
    xt1(NIndex'+1,:)=[]; %[n-j-1]
    u(NIndex'+1,:)=[]; %[n] discard first sample everyday for measurement
    u0(NIndex'+NumDayPartition-(j-1),:)=[]; %[n-j] discard last sample everyday for Y1
    if (j <= p)
        u1(NIndex'+1,:)=[]; %[n-j-1]
    end
end

tt = linspace(0, DayNum, length(xt));
N = length(tt);

RegMatrixT = zeros([N VSign+2 VSign]);
RegCoefT = zeros([VSign+2 VSign]);

for i = 1:VSign
    RegMatrixT(:,:,i) = [xt0 u0(:,i) w_forecast];
    RegMatrixW(:,:,i) = [xt0 u0(:,i) w];
    RegCoefT(:,i) = RegMatrixT(1:(NumDayPartition-m) * DayNum,:,i) ...
        \ xt(1:(NumDayPartition-m) * DayNum,i); %LS regression of coefficient of each X term
    RegCoefW(:,i) = RegMatrixW(1:(NumDayPartition-m) * DayNum,:,i) ...
        \ xt(1:(NumDayPartition-m) * DayNum,i); %LS regression of coefficient of each X term
end

for i = 1:VSign

    Predxt(:,i) = RegMatrixT(:,:,i) * RegCoefT(:,i);
    PredxtW(:,i) = RegMatrixW(:,:,i) * RegCoefT(:,i);
    PredxtPerfectModel(:,i) = RegMatrixW(:,:,i) * RegCoefW(:,i);
end

A = RegCoefT(1:VSign,:)';
B = RegCoefT(VSign+1,:)';
B = diag(B);
C = RegCoefT(VSign+2,:)';

% for i = 1:VSign
%     figure
%     plot(tt,k2f(xt(1:(NumDayPartition-m)*DayNum,i)),'b-',tt, k2f(PredxtW(1:(NumDayPartition-m)*DayNum,i)),'r--',...
%         tt, k2f(PredxtPerfectModel(1:(NumDayPartition-m)*DayNum,i)),'g-')
%     xlim([0 DayNum])
%     legend("Ground Truth", "Forecast Model", "Perfect Model")
%     xlabel('Days (May)')
%     ylabel('Temperature /F')
% end

Rtrue = zeros(VSign, 5);

for i = 1:VSign
    Rtrue(i, 1) = gfit2(xt(:,i), PredxtPerfectModel(:,i),'5'); %MAE
    Rtrue(i, 2) = gfit2(xt(:,i), PredxtPerfectModel(:,i),'1'); %MSE
    Rtrue(i, 3) = gfit2(xt(:,i), PredxtPerfectModel(:,i),'3'); %RMSE
    Rtrue(i, 4) = gfit2(xt(:,i), PredxtPerfectModel(:,i),'8'); %r^2
    Rtrue(i, 5) = gfit2(xt(:,i), PredxtPerfectModel(:,i),'10'); %MaxAE
end

Rforecast = zeros(VSign, 5);

for i = 1:VSign
    Rforecast(i, 1) = gfit2(xt(:,i), PredxtW(:,i),'5'); %MAE
    Rforecast(i, 2) = gfit2(xt(:,i), PredxtW(:,i),'1'); %MSE
    Rforecast(i, 3) = gfit2(xt(:,i), PredxtW(:,i),'3'); %RMSE
    Rforecast(i, 4) = gfit2(xt(:,i), PredxtW(:,i),'8'); %r^2
    Rforecast(i, 5) = gfit2(xt(:,i), PredxtW(:,i),'10'); %MaxAE
end

save VAVMPCDailyForecast.mat
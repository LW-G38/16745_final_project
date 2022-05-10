close all
clear all
swEPSfigure

%loading the building parameters
BuildingParameter = load ('SEBuildingParameterSuspendedRoofWithoutV358R.txt');
Surf = BuildingParameter (1:32,2);
Wd  = BuildingParameter (1:32,3);
Wa  = BuildingParameter (1:32,4);
Vol  = BuildingParameter (1:32,5);
r = Vol/sum(Vol);

load VAVdataMay2014.mat
load forecast.mat
VAVdata = VAV201TO216301TO316May2014r;
FlowDemand2f = csvread('AHU9VAV2fSPr.csv',0,1); %flow rate set point
FlowDemand3f = csvread('AHU9VAV3fSPr.csv',0,1);
VSign = 32;
for i = 1: VSign
    if i<=16 %
        flow(:,i) = VAVdata(:,(i-1)*8 + 1); %actual flow rate
        Trm(:,i) = VAVdata(:,(i-1)*8 + 2); %actual room temp
        Tdisch(:,i) = VAVdata(:,(i-1)*8 + 3);
        Hw(:,i) = VAVdata(:,(i-1)*8 + 4);
        Clg(:,i) = VAVdata(:,(i-1)*8 + 5);
        TspC(:,i) = VAVdata(:,(i-1)*8 + 6);
        Htg(:,i) = VAVdata(:,(i-1)*8 + 7);
        TspH(:,i) = VAVdata(:,i*8);
    else
        flow(:,i) = VAVdata(:,(i-1)*8 + 1);
        Trm(:,i) = VAVdata(:,(i-1)*8 + 2);
        Tdisch(:,i) = VAVdata(:,(i-1)*8 + 3);
        Clg(:,i) = VAVdata(:,(i-1)*8 + 4);
        Htg(:,i) = VAVdata(:,(i-1)*8 + 5);
        Hw(:,i) = VAVdata(:,(i-1)*8 + 6);
        TspC(:,i) = VAVdata(:,(i-1)*8 + 7);
        TspH(:,i) = VAVdata(:,i*8);        
    end
end
VAVdata = [Trm Tdisch (TspC+TspH)/2 flow FlowDemand2f FlowDemand3f Htg]; %FlowDemand2f + FlowDemand3f: 32 columns
TrmOri = Trm;
% clear VAV201TO216301TO316May2014r flow Trm Tdisch Hw Clg TspC TspH Htg FlowDemand2f FlowDemand3f

%load AHU data
[date oat OA_RH wind] = textread('OATMayr.csv', '%s %f %f %f');
Temp = csvread('AHU9Tr.csv',0,1); %Temp info
AHUdata = csvread('AHU9Fr.csv',0,1); %flow rate & damper position info
SPdata = csvread('AHU9Pr.csv',0,1); %static pressure & fan info

Wsf = SPdata(:,6);   %Fan power: Supply Fan
Wrf = SPdata(:,7);   %Fan power: Return Fan

Temp = [oat Temp oat_forecast]; %Toa Tsa Tra Tma
AHUflow = AHUdata(:,1:3); %OA flow RA flow SA flow
AHUDamper = AHUdata(:,4:6); %EA damper OA damper RA damper

AHUdata = [Temp AHUflow AHUDamper Wsf Wrf];

clear Temp SPdata date oat OA_RH wind AHUflow AHUDamper Wsf Wrf oat_forecast

%Sampling Info
day = 1440;
sample = 15;
NumDailySamples = day / sample;
DayNum = length(AHUdata) / NumDailySamples - 1;
NumDayPartition = 18 * 60 / sample;
NumNightPartition = NumDailySamples - NumDayPartition;
NumSampleBeforeNightPartitoin = 1 * 60 / sample;

% discard the samples from 12am to 1am in the first day
VAVdata (1:NumSampleBeforeNightPartitoin,:) = [];
AHUdata (1:NumSampleBeforeNightPartitoin,:) = [];
% discard the samples from 1am to 12am in the last day
VAVdata(length(VAVdata)-(NumDailySamples-NumSampleBeforeNightPartitoin)+1:length(VAVdata),:) = [];
AHUdata(length(AHUdata)-(NumDailySamples-NumSampleBeforeNightPartitoin)+1:length(AHUdata),:) = [];
% seperate the samples into Day Partition and Night Partition
for i = 1 : length(VAVdata) / NumDailySamples
    VAVdataNight((i-1)*NumNightPartition+1 : i*NumNightPartition ,:) = VAVdata((i-1)*NumDailySamples+1 : (i-1)*NumDailySamples+NumNightPartition ,:);
    VAVdataDay((i-1)*NumDayPartition+1 : i*NumDayPartition ,:) = VAVdata((i-1)*NumDailySamples+NumNightPartition+1 : i*NumDailySamples ,:);
    AHUdataNight((i-1)*NumNightPartition+1 : i*NumNightPartition ,:) = AHUdata((i-1)*NumDailySamples+1 : (i-1)*NumDailySamples+NumNightPartition ,:);
    AHUdataDay((i-1)*NumDayPartition+1 : i*NumDayPartition ,:) = AHUdata((i-1)*NumDailySamples+NumNightPartition+1 : i*NumDailySamples ,:);
end

AData = AHUdataDay;
VData = VAVdataDay;

%Coefficients list
% conversion coefficient from feet^3/min to m^3/sec
CFMCoeff = 4.719333 * 10 ^ ( - 4 );
a = 0.0254;
% volumetric heat capacity of air at 20 C [KJ * m^(-3) * K^(-1)]
Cv = 1.211;
% density of air at 20 C [kg * m^(-3)]
%density = 1.205; %// this should be a function for higher accuracy
% one standard atm in [Pa]
atm = 101325;
% [WC] to [Pa]
WC2Pa = 249.088908333;

%Data conversion: AHU
Toa = f2k(AData(:,1));
Tsa = f2k(AData(:,2));
Tra = f2k(AData(:,3));
Tma = f2k(AData(:,4));
Toa_forecast = f2k(AData(:,5));

OAflow = CFMCoeff * AData(:,5);
RAflow = CFMCoeff * AData(:,6);
SAflow = CFMCoeff * AData(:,7);

RAflowOrign = RAflow;

%Data conversion: VAV
Trm = f2k(VData(:,1:32));
Tdisch = f2k(VData(:,(1:32) + 32));
Tsp = f2k(VData(:,(1:32) + 2*32));
flow = CFMCoeff * VData(:,(1:32) + 3*32);
flowSP = CFMCoeff * VData(:,(1:32) + 4*32);
reheat = VData(:,(1:32) + 5*32);

time = linspace(0,26,length(SAflow));

OAflow = sgolayfilt(OAflow,3,5);
RAflow = sgolayfilt(RAflow,3,5);
SAflow = sgolayfilt(SAflow,3,5);

EAD = AData(:,8);
OAD = AData(:,9);
RAD = AData(:,10);

betaD = RAD ./ (EAD+RAD);
OAflow = SAflow - betaD .* RAflow;

%Data conversion: VAV
Trm = f2k(VData(:,1:VSign));
Tin = Trm * r;
SCoeffT = linfit( Tin, Tra);
TrmReg = SCoeffT(1) * Tin + SCoeffT(2);
FlowDemand = CFMCoeff * VData(:,VSign * 4 + 1 : VSign * 5);
SAflowIn = sum(FlowDemand,2);
SCoeffS  = linfit( SAflowIn, SAflow );
SAflowIn =SCoeffS(1) * SAflowIn+ SCoeffS(2);
Tdd = f2k(VData(:,1 + VSign:2 * VSign));
Tddd = Tdd * r;
SCoeffTdd = linfit(Tddd, Tsa);
TddReg = SCoeffTdd(1) * Tddd + SCoeffTdd(2);

% for i = 1 : VSign
%     figure
%     plot(time,k2f(Trm(:,i)),time,k2f(Tsp(:,i)),'--r')
%     xlim([DayNum - 9 DayNum])
%     xlabel('Days (May)')
%     ylabel('Temperature (F)')
%     legend('$T_{RM}$','$T_{SP}$')
% end

Dsign = 7;
figure
plot(6+18*time(1:72),k2f(Toa((Dsign-1)*72+1:Dsign*72)),6+18*time(1:72),k2f(Toa_forecast((Dsign-1)*72+1:Dsign*72)),'--k')
xlim([6 24])
%xlim([DayNum - 9 DayNum])
xlabel('Hour')
ylabel('Temperature ($^{\circ}$F)')
legend('Groundtruth $T_{OA}$','Forecast $T_{OA}$')
title(['Real-time $T_{OA}$ vs Forecast $T_{OA}$ - Day', num2str(Dsign)])
save Data4RegressionDaily2022.mat

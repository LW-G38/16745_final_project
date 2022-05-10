clear all
close all

load VAVMPCDailyForecast.mat

%Parameter Indentifiction: Forecasting Model
A = RegCoefT(1:VSign,:);
B = diag(RegCoefT(VSign+1,:));
C = RegCoefT(VSign+2,:);

Atrue = RegCoefW(1:VSign,:);
Btrue = diag(RegCoefW(VSign+1,:));
Ctrue = RegCoefW(VSign+2,:);

% back to full length data
xt = Trm;
u = Tdisch.*flow;
w_forecast = Toa_forecast;
w = Toa;
tt = linspace(0, DayNum, length(xt));
N = length(tt);
RegMatrixT = zeros([N VSign+2 VSign]);
RegMatrixW = RegMatrixT;
Predxt = zeros(size(xt));
PredxtW = Predxt;
PredxtPerfectModel = Predxt;

for i = 1:VSign
    RegMatrixT(:,:,i) = [xt u(:,i) w_forecast];
    RegMatrixW(:,:,i) = [xt u(:,i) w];
    Predxt(:,i) = RegMatrixT(:,:,i) * RegCoefT(:,i);
    PredxtW(:,i) = RegMatrixW(:,:,i) * RegCoefT(:,i);
    PredxtPerfectModel(:,i) = RegMatrixW(:,:,i) * RegCoefW(:,i);
end

% for i = 1:VSign
%     figure
%     plot(tt,k2f(xt(1:NumDayPartition*DayNum,i)),'b-',tt, k2f(PredxtW(1:NumDayPartition*DayNum,i)),'r--',...
%         tt, k2f(PredxtPerfectModel(1:NumDayPartition*DayNum,i)),'g-')
%     xlim([0 DayNum])
%     legend("Ground Truth", "Forecast Model with true OA", "Perfect Model")
%     xlabel('Days (May)')
%     ylabel('Temperature /F')
% end
Dsign = 15;
for i = 1:VSign
    figure
    plot(6+18*time(1:72),k2f(xt((Dsign-1)*72+1:Dsign*72,i)),'b-',6+18*time(1:72), k2f(PredxtW((Dsign-1)*72+1:Dsign*72,i)),'r--',...
        6+18*time(1:72), k2f(PredxtPerfectModel((Dsign-1)*72+1:Dsign*72,i)),'m-')
    xlim([6 24])
    legend("Ground Truth", "Forecast Model with true OA", "Ideal Model","Location","best")
    xlabel('Hour')
    ylabel('Temperature ($^{\circ}$F)')
    title(['Room' num2str(i) ' Model Comparison - Day' num2str(Dsign)])
end

x_rollout = zeros(size(xt));
x_rollout_forecast = x_rollout;
for j = 1:DayNum
    x_rollout(NumDayPartition*(j-1)+1,:) = xt(NumDayPartition*(j-1)+1,:);
    x_rollout_forecast(NumDayPartition*(j-1)+1,:) = xt(NumDayPartition*(j-1)+1,:);   
    for k = 2:NumDayPartition
        x_rollout(NumDayPartition*(j-1)+k,:) = x_rollout(NumDayPartition*(j-1)+k-1,:)*Atrue + ...
            u(NumDayPartition*(j-1)+k,:)*Btrue + w(NumDayPartition*(j-1)+k-1)*Ctrue;
        x_rollout_forecast(NumDayPartition*(j-1)+k,:) = x_rollout_forecast(NumDayPartition*(j-1)+k-1,:)*A + ...
            u(NumDayPartition*(j-1)+k,:)*B + w(NumDayPartition*(j-1)+k-1)*C;
    end
end

% for i = 1:VSign
%     figure
%     plot(tt,k2f(xt(1:NumDayPartition*DayNum,i)),'b-',tt, k2f(x_rollout_forecast(1:NumDayPartition*DayNum,i)),'r--',...
%         tt, k2f(x_rollout(1:NumDayPartition*DayNum,i)),'g-')
%     xlim([0 DayNum])
%     legend("Ground Truth", "Daily Rollout - forecast model", "Daily Rollout - perfect model")
%     xlabel('Days (May)')
%     ylabel('Temperature /F')
% end

Rx_rollout = zeros(VSign, 5);

for i = 1:VSign
    Rx_rollout(i, 1) = gfit2(xt(:,i), x_rollout(:,i),'5'); %MAE
    Rx_rollout(i, 2) = gfit2(xt(:,i), x_rollout(:,i),'1'); %MSE
    Rx_rollout(i, 3) = gfit2(xt(:,i), x_rollout(:,i),'3'); %RMSE
    Rx_rollout(i, 4) = gfit2(xt(:,i), x_rollout(:,i),'8'); %r^2
    Rx_rollout(i, 5) = gfit2(xt(:,i), x_rollout(:,i),'10'); %MaxAE
end

Rx_forecast = zeros(VSign, 5);

for i = 1:VSign
    Rx_forecast(i, 1) = gfit2(xt(:,i), x_rollout_forecast(:,i),'5'); %MAE
    Rx_forecast(i, 2) = gfit2(xt(:,i), x_rollout_forecast(:,i),'1'); %MSE
    Rx_forecast(i, 3) = gfit2(xt(:,i), x_rollout_forecast(:,i),'3'); %RMSE
    Rx_forecast(i, 4) = gfit2(xt(:,i), x_rollout_forecast(:,i),'8'); %r^2
    Rx_forecast(i, 5) = gfit2(xt(:,i), x_rollout_forecast(:,i),'10'); %MaxAE
end
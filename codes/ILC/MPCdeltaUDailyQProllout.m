clear all
close all

load VAVMPCDailyForecast.mat
% back to full length data
x0 = Trm;
Ubas = Tdisch.*flow;
w_forecast = Toa_forecast;
w = Toa;
tt = linspace(0, DayNum, length(xt));
N = length(tt);

%% torque limits
% airflow limits
uf_2f_min = [200 190 160 110 120 130 110 110 120 130 110 75 60 110 110 510];
uf_2f_max = [800 750 630 450 470 530 450 450 470 530 450 300 230 440 430 2030];

uf_3f_min = [195 170 150 110 120 120 110 110 120 120 110 75 60 105 105 225];
uf_3f_max = [780 690 600 450 470 490 450 450 470 490 450 300 230 420 420 1020];

% discharge air temperature limits
ut_min = 55;
ut_max = 70;

% control limits

u_min = f2k(ut_min)*CFMCoeff*[uf_2f_min uf_3f_min];
u_max = f2k(ut_max)*CFMCoeff*[uf_2f_max uf_3f_max];

%Parameter Indentifiction: Forecasting Model
A = RegCoefT(1:VSign,:)';
B = RegCoefT(VSign+1,:)';
B = diag(B);
C = RegCoefT(VSign+2,:)';

%Real-time designed data before for energy consumption calculation
Toa_forecast = w_forecast;
Tsp = f2k(VData(:,(1:32) + 2*32));
%%
p = 3;
NumDayMPC = NumDayPartition - p;
%Quadratic Weight Factor
Q = 1000000*eye(32);
R = 0.01*eye(32);

%accessory
I = eye(32);
O = zeros(32,32);
Xcalc = zeros(size(xt0));
Xsim = zeros(size(xt0));
Usim = zeros(size(Ubas));

% Hqp = blkdiag(Q,R);
Hqp = blkdiag(Q,Q,Q,R,R,R);
for i = 1 : 1
    Dsign = i - 1; %Select the day for MPC
    Xsim(Dsign * NumDayPartition + 1,:) = xt0(Dsign * NumDayPartition + 1,:)';
    Xcalc(Dsign * NumDayPartition + 1,:) = xt0(Dsign * NumDayPartition + 1,:)';
    
    for j = 1 : NumDayMPC
        Hsign = j; %Select the hour after occupied mode turns on(7am). Ex: 2*4 means 9 am for 4 samples per hour
        %Define the disigned data in MPC steps
        Wf = Toa_forecast(Dsign * NumDayPartition + Hsign: Dsign * NumDayPartition + Hsign + 2);
        Wt = Toa(Dsign * NumDayPartition + Hsign: Dsign * NumDayPartition + Hsign + 2);
        Xprev = Xsim(Dsign * NumDayPartition + Hsign,:)';
        Xprev_calc = Xcalc(Dsign * NumDayPartition + Hsign,:)';
        Xref = Tsp(Dsign * NumDayPartition + Hsign + 1 : Dsign * NumDayPartition + Hsign + 3,:)';        
        %qp for MPC
        Fqp = [-Q'*Xref(:,1); -Q'*Xref(:,2); -Q'*Xref(:,3); zeros(32*3,1)];
        Aeq = [-I O O B O O; A -I O O B O; O A -I O O B];
        beq = -[A*Xprev + C*Wf(1); C*Wf(2); C*Wf(3)];
        lb = [-Inf*ones(32*3,1); u_min'; u_min'; u_min'];
        ub = [Inf*ones(32*3,1); u_max'; u_max'; u_max'];
        Sol = quadprog(Hqp,Fqp,[],[],Aeq,beq,lb,ub);
        % roll out on real outdoor temperature
        Ucurr = Sol((1:32) + 3*32);
        Usim(Dsign * NumDayPartition + Hsign,:) = Ucurr;
        Xsim(Dsign * NumDayPartition + Hsign + 1,:) = (A*Xprev + B*Ucurr + C*Wt(1))';
        Xcalc(Dsign * NumDayPartition + Hsign + 1,:) = (A*Xprev_calc + B*Ucurr + C*Wf(1))';
    end
    % roll out the remainder of the day
    for k = 1 : 2
        Xprev = Xsim(Dsign * NumDayPartition + Hsign + k,:)';
        Xprev_calc = Xcalc(Dsign * NumDayPartition + Hsign + k,:)';
        Wprev = Toa(Dsign * NumDayPartition + Hsign + k,:)';
        Wprev_f = Toa_forecast(Dsign * NumDayPartition + Hsign + k,:)';
        UCurr = Sol((1:32) + (3+k)*32);
        Usim(Dsign * NumDayPartition + Hsign + k,:) = Ucurr;
        Xsim(Dsign * NumDayPartition + Hsign + k + 1,:) = (A*Xprev + B*Ucurr + C*Wprev)';
        Xcalc(Dsign * NumDayPartition + Hsign + k + 1,:) = (A*Xprev_calc + B*Ucurr + C*Wprev_f)';
    end
end
deltaXref = Xsim(1: NumDayPartition,:) - Tsp(1: NumDayPartition,:);
% deltaXref = Xsim(1: NumDayPartition,:) - Xcalc(1: NumDayPartition,:);
deltaXref = reshape(deltaXref',1,[])';
deltaW = w_forecast(1:NumDayPartition-1) - w(1:NumDayPartition-1);
Uref = Usim(1: NumDayPartition-1,:);
Uref = reshape(Uref',1,[])';
Uref = zeros(32*(NumDayPartition-1),1);

Q_ilc = eye(NumDayPartition*32);
R_ilc = eye((NumDayPartition-1)*32);
Hqp_ilc = blkdiag(Q_ilc,R_ilc);
Fqp_ilc = [-Q_ilc*deltaXref; zeros(32*(NumDayPartition-1),1)];
Aeq_ilc = zeros(32*(NumDayPartition-1),32*(2*NumDayPartition-1));
beq_ilc = zeros(32*(NumDayPartition-1),1);
for k = 1 : NumDayPartition - 1
    Aeq_ilc((1:32)+(k-1)*32, (1:32)+(k-1)*32) = A;
    Aeq_ilc((1:32)+(k-1)*32, (1:32)+k*32) = -eye(32);
    Aeq_ilc((1:32)+(k-1)*32, (1:32)+(k-1)*32 + (NumDayPartition-1)*32) = B;
    beq_ilc((1:32)+(k-1)*32) = -C*deltaW(k);
end
lb_ilc = [-Inf*ones(32*NumDayPartition,1); repmat(u_min',NumDayPartition-1,1) - Uref];
ub_ilc = [Inf*ones(32*NumDayPartition,1); repmat(u_max',NumDayPartition-1,1) - Uref];
Sol_ilc = quadprog(Hqp_ilc,Fqp_ilc,[],[],Aeq_ilc,beq_ilc,lb_ilc,ub_ilc);
deltaU = Sol_ilc(NumDayPartition*32 + 1:end);
%%
deltaU = reshape(deltaU, [32 NumDayPartition-1]);
for i = 2 : 2
    Dsign = i - 1; %Select the day for MPC
    Xsim(Dsign * NumDayPartition + 1,:) = xt0(Dsign * NumDayPartition + 1,:);
    Xcalc(Dsign * NumDayPartition + 1,:) = xt0(Dsign * NumDayPartition + 1,:)';
    
    for j = 1 : NumDayMPC
        Hsign = j; %Select the hour after occupied mode turns on(7am). Ex: 2*4 means 9 am for 4 samples per hour
        %Define the disigned data in MPC steps
        Wf = Toa_forecast(Dsign * NumDayPartition + Hsign: Dsign * NumDayPartition + Hsign + 2);
        Wt = Toa(Dsign * NumDayPartition + Hsign: Dsign * NumDayPartition + Hsign + 2);
        Xprev = Xsim(Dsign * NumDayPartition + Hsign,:)';
        Xprev_calc = Xcalc(Dsign * NumDayPartition + Hsign,:)';
        Xref = Tsp(Dsign * NumDayPartition + Hsign + 1 : Dsign * NumDayPartition + Hsign + 3,:)';        
        %qp for MPC
        Fqp = [-Q'*Xref(:,1); -Q'*Xref(:,2); -Q'*Xref(:,3); zeros(32*3,1)];
        Aeq = [-I O O B O O; A -I O O B O; O A -I O O B];
        beq = -[A*Xprev + C*Wf(1); C*Wf(2); C*Wf(3)];
        lb = [-Inf*ones(32*3,1); u_min'; u_min'; u_min'];
        ub = [Inf*ones(32*3,1); u_max'; u_max'; u_max'];
        Sol = quadprog(Hqp,Fqp,[],[],Aeq,beq,lb,ub);
        % roll out on real outdoor temperature
        Ucurr = Sol((1:32) + 3*32) + deltaU(:,Hsign);
        Usim(Dsign * NumDayPartition + Hsign,:) = Ucurr;
        Xsim(Dsign * NumDayPartition + Hsign + 1,:) = (A*Xprev + B*Ucurr + C*Wt(1))';
        Xcalc(Dsign * NumDayPartition + Hsign + 1,:) = (A*Xprev_calc + B*Ucurr + C*Wf(1))';
    end
    % roll out the remainder of the day
    for k = 1 : 2
        Xprev = Xsim(Dsign * NumDayPartition + Hsign + k,:)';
        Xprev_calc = Xcalc(Dsign * NumDayPartition + Hsign + k,:)';
        Wprev = Toa(Dsign * NumDayPartition + Hsign + k,:)';
        Wprev_f = Toa_forecast(Dsign * NumDayPartition + Hsign + k,:)';
        UCurr = Sol((1:32) + (3+k)*32) + deltaU(:,Hsign) + k;
        Usim(Dsign * NumDayPartition + Hsign + k,:) = Ucurr;
        Xsim(Dsign * NumDayPartition + Hsign + k + 1,:) = (A*Xprev + B*Ucurr + C*Wprev)';
        Xcalc(Dsign * NumDayPartition + Hsign + k + 1,:) = (A*Xprev_calc + B*Ucurr + C*Wprev_f)';
    end
end
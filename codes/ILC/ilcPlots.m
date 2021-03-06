clear all
close all
load ILCresults.mat

%% temperature plots
for Rsign = 1:1
    styVec = ["-b","--r","--oc","--*g", "--dm","--^k","--xc","--|g"];
    figure
    for m = 1 : 8
        plot(6+18*time(1:72),k2f(Xsim((m-1)*NumDayPartition + (1:NumDayPartition), Rsign)),styVec(m))
        legappend(['day ',num2str(m)])
        hold on;
    end
    %yline(72.5,'r-','Day 1-6 SP');
    %yline(74,'r-','Day 7-8 SP');
    %yline(72.5,'r-','SP');
    yline(74,'r-','SP');
    xlim([6 24])
    xlabel('Hour')
    ylabel('Temperature ($^{\circ}$F)')
    title(['Temperature of Room' num2str(Rsign)])
end

%% control plots
%close all
for Rsign = 1:1
    styVec = ["-b","--r","--oc","--*g", "--dm","--^k","--xc","--|g"];
    figure
    for m = 1 : 8
        plot(6+18*time(1:72),Usim((m-1)*NumDayPartition + (1:NumDayPartition), Rsign),styVec(m))
        legappend(['day ',num2str(m)])
        hold on;
    end
    %yline(72.5,'r-','Day 1-6 SP');
    %yline(74,'r-','Day 7-8 SP');
    yline(u_max(Rsign),'r-','ub');
    yline(u_min(Rsign),'r-','lb');
    xlim([6 24])
    xlabel('Hour')
    ylabel('Control Signal')
    title(['Control inputs of Room' num2str(Rsign)])
end
% plotMultSimsTime_ters plots the channel base current and wavefront
% position as a function of time for multiple simulations
%
% Author: Caitano L. da Silva, NMT
% Created: Jun/14/2024
% Last modification: Jun/15/2025
%
clc; clear all


%% Input
simCases = {'testSim1.mat', ...
            'testSim2.mat', ...
            };
        

%% Initialize plot
Nsims = length(simCases);
Iunit = 1e-3; % Convert to kA
Vunit = 1e-6; % Convert to MV
tunit = 1e6;  % Convert to us
LFS   = 24;
FS    = 20;
LW    = 5;
MS    = 12;
lcolors = [1 0 0; 0 .6 0; 0 0 1;  0.6275 0.3216 0.1765; 0.55 0 0.55];
lstyl = {'-','-','--','-','-'};
lgray = 0.7*[1 1  1];


%% Plot
f2 = figure(2); clf
f2.Position = [100 100 1400 500];

for j=1:Nsims
    
    % Load simulation
    load(['./sims/' simCases{j}]);
    
    % Generate legend
    % this snippet needs to be adjusted so the proper 
    % input params appear in the legend
    simCasestr{j} = ['R_0 = ' num2str(drs.params.R0) ' \Omega/m,' blanks(2-length(num2str(drs.params.R0))) ' \tau_a = ' num2str(drs.params.tauA*1e6) ' \mus'];
    
    % Plot
    subplot(1,2,1)
    pIg = plot(drs.t*tunit,drs.Igt*Iunit,lstyl{j},'Color',lcolors(j,:),'LineWidth',((1+Nsims/10) - (j-1)*Nsims/30)*LW);
    if j==1
        hold on
        Ip = drs.params.Ip;        
    end    
    if j==Nsims
        hold off
        xlim([-0.05 1]*max(drs.t)*tunit)
        set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
        ylabel('Channel-base current (kA)','FontSize',LFS)
        xlabel('Time (\mus)','FontSize',LFS)               
    end
    
    subplot(1,2,2)
    pz(j) = plot(drs.t*tunit,drs.zIwavpos,lstyl{j},'Color',lcolors(j,:),'LineWidth',((1+Nsims/10) - (j-1)*Nsims/30)*LW);    
    if j==1
        hold on
    end
    if j==Nsims
       
        trange = drs.t(drs.t>=8e-6);
        pvt = plot(drs.t*tunit,drs.zfrntpos,':','Color',lgray,'LineWidth',LW);        
        hold off
        ylim([0 max(drs.zfrntpos)])
        xlim([-0.05 1]*max(drs.t)*tunit)
        legend([pz pvt],[simCasestr, 'v*t'],'Location','SouthEast')
        set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
        ylabel('Wave front position (m)','FontSize',LFS)
        xlabel('Time (\mus)','FontSize',LFS)
    end
        
end
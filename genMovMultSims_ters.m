% genMovMultSims_ters generates a movie including the results of one
% or more TERS simulations. The movie displays vertical profiles of
% current, potential, and resistance. This works best when tsave is defined
% in TERS like this, for example: Ntsave = 101; tsave = linspace(0,tmax,Ntsave);
% where Ntsave will be the number of frames in the movie
%
% Author: Caitano L. da Silva, NMT
% Created: Mar/26/2025
% Last modification: Jun/15/2025
clc; clear all


%% Input
simCases = {'testSim1.mat';
            'testSim2.mat';
            };

legName = {'Case 1','Case 2'}; % Legend names

iSaveMov = 1; % =1 to save movie


%% Load simulation that will define the boundaries
Nsims  = length(simCases);
load(['./sims/' simCases{1}]); % 1st one in this example
Ntsave = length(drs.tsave);
zmax   =  max(drs.zI);
v      = drs.params.v;
Vlims  = drs.Vlims;
Ilims  = drs.Ilims;
Rlims  = drs.Rlims;


%% Movie
fh1 = figure(1); clf
Iunit = 1e-3; % Convert to kA
Vunit = 1e-6; % Convert to MV
tunit = 1e6;  % Convert to us
LFS   = 22;
FS    = 18;
LW    = 2.2;
MS    = 12;

% Initialize file to save movie
if iSaveMov==1
    % Make directory
    if exist('./sims','dir')
        % do nothing else
    else
        mkdir 'sims';
    end
    % Initialize movie
    movPanels = VideoWriter('./sims/multSimsMov.mp4','MPEG-4');
    open(movPanels);
end

% Main loop
for ksave=1:Ntsave    
    
    for j=1:Nsims
        load(['./sims/' simCases{j}]);
        
        subplot(1,3,1) % Current
        pI  = plot(drs.Isave(ksave,:)*Iunit,drs.zI,'-','LineWidth',LW);
        
        if j==1 && Nsims>1
            hold on
        end
        if j==Nsims
            hold off
            if Nsims>1
                legend(legName,'Location','NorthEast')
            end
            xlim(Ilims*Iunit)
            ylim([0 zmax])
            set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
            ylabel('Height (m)','FontSize',LFS)
            xlabel('Current (kA)','FontSize',LFS)
        end
        
        subplot(1,3,2) % Potential
        pV = plot(drs.Vsave(ksave,:)*Vunit,drs.zV,'-','LineWidth',LW);
        if j==1 && Nsims>1
            hold on
        end
        if j==Nsims
            hold off
            xlim(Vlims*Vunit)
            ylim([0 zmax])
            set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
            ylabel('Height (m)','FontSize',LFS)
            xlabel('Electric potential (MV)','FontSize',LFS)
            title(['t = ' sprintf('%4.1f',drs.tsave(ksave)*tunit) ' \mus'],'FontWeight','Normal');
        end
        
        subplot(1,3,3) % Resistance
        plot(drs.Rsave(ksave,:),drs.zI,'-','LineWidth',LW);
        if j==1 && Nsims>1
            hold on
        end
        if j==Nsims
            line(Rlims,drs.tsave(ksave)*v*[1 1],'LineStyle','--','Color',.3*[1 1 1],'LineWidth',.9*LW)
            hold off
            xlim(Rlims)
            ylim([0 zmax])
            set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
            ylabel('Height (m)','FontSize',LFS)
            xlabel('Resistance (\Omega/m)','FontSize',LFS)
        end
        
    end    
    drawnow;
    
    % Save movie?
    if iSaveMov==1
        currFrame = getframe(fh1);
        writeVideo(movPanels,currFrame);
    else
        pause(.0001)
    end
    
end

% Close file with movie
if iSaveMov==1
    close(movPanels);
end
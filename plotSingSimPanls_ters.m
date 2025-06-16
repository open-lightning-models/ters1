% plotSingSimPanls_ters plots Current & Potential as a function of height
% for a few selected time instants. This works best when tsave is defined
% in TERS like this, for example: tsave = [2 6 10]*1e-6;
%
% Author: Caitano L. da Silva, NMT
% Created: May/27/2025
% Last modification: Jun/15/2025
%
clc; clear all

%% Input
simCase = 'testSim.mat';


%% Load sim
load(['./sims/' simCase]);
iEvAnlyt = drs.flags.iEvAnlyt;
iupdateR = drs.flags.iupdateR;
Ntsave   = length(drs.tsave);
Imin     = min(drs.Isave(:));
Imax     = max(drs.Isave(:));
dI       = 0.05*(Imax-Imin);
zmax     = max(drs.zI);
Rlims    = drs.Rlims;


%% Sanity check
if Ntsave>10
    error('These are too many curves to be plotted in a single graph.');
end


%% Plot
Iunit = 1e-3; % Convert to kA
Vunit = 1e-6; % Convert to MV
tunit = 1e6;  % Convert to us
LFS   = 24;
FS    = 20;
LW    = 3;
MS    = 14;
lcolors = [1 0 0; 0 .6 0; 0 0 1];
lstyl = {'-','-','--','-','--'};
lgray = 0.7*[1 1  1];

f2 = figure(2); clf
f2.Position = [100 100 1400 500];

for ksave=1:Ntsave
    
    tsavestr{ksave} = [num2str(drs.tsave(ksave)*tunit) ' \mus'];
    
    sbp1 = subplot(1,2+iupdateR,1);    
    pwp = plot(interp1(drs.zI,drs.Isave(ksave,:)*Iunit,drs.zIwavsav(ksave)),drs.zIwavsav(ksave),'^','MarkerFaceColor',lcolors(ksave,:),'MarkerEdgeColor','k','MarkerSize',MS);
    
    if ksave==1
        hold on
    end
    pI  = plot(drs.Isave(ksave,:)*Iunit,drs.zI,'-','Color',lcolors(ksave,:),'LineWidth',1.9*LW);
    if iEvAnlyt==1
        pIa = plot(drs.Iasave(ksave,:)*Iunit,drs.zI,'--','Color',lgray,'LineWidth',1.1*LW);
    end
    pfp = plot((Imin-0.5*dI)*Iunit,drs.zfrntsav(ksave),'o','MarkerFaceColor',lcolors(ksave,:),'MarkerEdgeColor','k','MarkerSize',MS);
        
    if ksave==Ntsave
        hold off
        xlim(drs.Ilims*Iunit)
        ylim([0 zmax])
        if iEvAnlyt==1
            legend([pIa, pI, pfp, pwp],'Analytical','Numerical','v*t','Wave front position','Location','Best')
        elseif iEvAnlyt==0
            legend([pI, pfp, pwp],'Calculation','v*t','Wave front position')
        end
        set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
        ylabel('Height (m)','FontSize',LFS)
        xlabel('Current (kA)','FontSize',LFS)
        pos_sbp1 = get(sbp1,'Position');
    end
            
    sbp2 = subplot(1,2+iupdateR,2);    
    pV(ksave)  = plot(drs.Vsave(ksave,:)*Vunit,drs.zV,'-','Color',lcolors(ksave,:),'LineWidth',1.9*LW);
    if ksave==1
        hold on
    end
    if iEvAnlyt==1
        plot(drs.Vasave(ksave,:)*Vunit,drs.zV,'--','Color',lgray,'LineWidth',1.1*LW);
    end
    if ksave==Ntsave        
        hold off
        xlim(drs.Vlims*Vunit)
        ylim([0 zmax])
        legend(pV,tsavestr,'Location','Best')
        set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on') % ,'YAxisLocation','right'
        ylabel('Height (m)','FontSize',LFS)
        xlabel('Potential (MV)','FontSize',LFS)
        pos_sbp2 = get(sbp2,'Position');
    end
    
    
    if iupdateR==1
        sbp3 = subplot(1,2+iupdateR,3);
        plot(drs.Rsave(ksave,:),drs.zI,'-','Color',lcolors(ksave,:),'LineWidth',1.9*LW);
        if ksave==1
            hold on
        end
        xlim(Rlims)

        if ksave==Ntsave
            ylim([0 zmax])
            set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on') % ,'YAxisLocation','right'
            ylabel('Height (m)','FontSize',LFS)
            xlabel('Resistance (\Omega/m)','FontSize',LFS)
            pos_sbp3 = get(sbp3,'Position');
            
            % Fix all positions
            pos_sbp1(3) = pos_sbp1(3)*1.25;
            pos_sbp2(3) = pos_sbp2(3)*1.25;            
            
            pos_sbp2(1) = pos_sbp2(1) + pos_sbp1(3)*0.2;
            pos_sbp3(1) = pos_sbp3(1) + pos_sbp1(3)*0.205 + pos_sbp2(3)*0.205;
            pos_sbp3(3) = pos_sbp3(3)*0.68;
            
            set(sbp1,'Position',pos_sbp1)
            set(sbp2,'Position',pos_sbp2)
            set(sbp3,'Position',pos_sbp3)
        end
    end        
    
end
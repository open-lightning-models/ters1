% Telegrapher's Equation Return Stroke (TERS) model
%
% Author: Caitano L. da Silva, NMT
% Created: circa Nov/2022
% Last modification: Jun/15/2025
%
clc; clear all


%% Universal constants
e0 = 8.85e-12;
mu0 = 4*pi*1e-7;


%% Input params

% Spatial
zmax = 2000;
dz = 1;

% Temporal
tmax = 1.1e-5;
dt = 1e-10;

% Trans. line params. (linear params)
rq = 50;       % Radius of charge sheath [m]
rI = 5e-3;     % Current-carrying radius [m]
R0 = 3;        % Linear resistance [Ohm/m]
Vi = -10e6;    % Initial channel potential [V]
tauA = 1e-7;   % Attachment time [s]

% Non-linear resistance params
iupdateR = 1;       % =0 R=const, =1 use nonlinear R model
Ess = 2.9e3;        % Steady-state field
tauRss = 1e-6;      % Time scale to achieve steady state resistance

% Save (which time instants should be saved?)
Ntsave = 101; tsave = linspace(0,tmax,Ntsave); % Defined via number of instants to be save
% tsave = [2 6 10]*1e-6; Ntsave = length(tsave);  % Defined via particular instants

% Flags
iEvAnlyt  = 0; % Set =1 to evaluate analytical solution
iUseAbsBC = 0; %>0 uses absorbing boundary which encompasses the top iUseAbsBC fraction of the domain (functionality under development)
iMakeMov  = 1; % =0 no movie, =1 make it, =2 make and save movie
iSavVars  = 1; % =0 don't save vars, =1 does save it
caseName  = 'testSim'; % file name with saved vars


%% Secondary constants: TL elecrical params
zmaxLC = 2000;               % I am putting this number here so that L & C do not change when I change zmax
C = 2*pi*e0/log(zmaxLC/rq);  % Linear capacitance [F/m]
L = mu0/2/pi*log(zmaxLC/rI); % Linear inductance [H/m]
v = 1/sqrt(L*C);             % Speed of EM waves guided by the TL [m/s]
Z = sqrt(L/C);               % Impedance [Ohm]
Ip = -Vi/Z;                  % Peak current [A]
tauR = 2*L/R0;               % Time scale for resistive wave decay [s]
b = dt/tauRss/2;
sIp = sign(Ip);
aVi = abs(Vi);


%% Check numerical resolution and issue warnings
if dt>dz/v/2
    warning('Insufficient temporal resolution.');
end

if iEvAnlyt==1
    if iupdateR>0 && tauRss>1e-7
        warning('The nonlinear resistance steady state time scale may be too long to yield a satisfactory comparison.');
    end
    
    if tauA>1e-6
        warning('The ground attachment time scale may be too long to yield a satisfactory comparison.');
    end
end


%% Independent vars
zI = 0:dz:zmax;
Nz = length(zI);
zV = zI(1:Nz-1) + dz/2;
t  = 0:dt:tmax;
Nt = length(t);


%% Boundary conditions: Voltage at ground level
Vgf = @(to) -Vi*(1 - exp(-to/tauA)).*(to>0);


%% Preallocation
% Model vars
V        = zeros(1,Nz-1);
I        = zeros(1,Nz);
R        = R0*ones(1,Nz);
Vo       = zeros(1,Nz-1);
Io       = zeros(1,Nz);
Ro       = zeros(1,Nz);
% Save vars
Vsave    = zeros(Ntsave,Nz-1);
Isave    = zeros(Ntsave,Nz);
Rsave    = zeros(Ntsave,Nz);
Imaxsave = zeros(1,Nt);
Igt      = zeros(1,Nt);
zfrntpos = NaN(1,Nt);
zIwavpos = NaN(1,Nt);
zfrntsav = NaN(1,Ntsave);
zIwavsav = NaN(1,Ntsave);
if iEvAnlyt==1
    Iasave = zeros(Ntsave,Nz);
    Vasave = zeros(Ntsave,Nz-1);
end


%% Implement absorbing boundary condition
% This is very simple! The resistance  just increases exponentially
% with height in the damping layer which occupies the top iUseAbsBC
% fraction of the simulation domain
if iUseAbsBC>0
    ldamp = iUseAbsBC*zmax/2/log(16*Z/iUseAbsBC/zmax/R0);
    [~, idamp] = min(abs(zI - (1-iUseAbsBC)*zmax));
    R(1,idamp:end) = R0*exp( (zI(idamp:end) - (1-iUseAbsBC)*zmax)/ldamp );
    idamp = min(idamp - 1,Nz);
    figure(3302); clf; semilogy(zI,R); ylim([0.1*R0 1.1*max(R)])
    
else
    idamp = Nz;
    
end


%% Define modified bessel funtion I1/rho
mymbesseli1 = @(ro) besseli(1,ro)./max(ro,1e-30).*(ro>0) + 0.5.*(ro==0);


%% Main loop
ksave = 1;
saveIntFlag = 1;
R0vec = R;
for k=1:Nt
    
    % Display time and front position, if interested
    %     fprintf('%1.1f\n',t(k)*v);
    %     fprintf('%1.1f\n',round(t(k)*v));
    
    % Max interation index to take care of Heaviside functions
    imax = min(ceil(t(k)*v/dz),Nz-1);
    
    % Store info
    Io = I;
    Vo = V;
    Ro = R;
    
    % Set BCs
    V(1) = Vgf(t(k));
    
    % First half step
    for i=2:imax
        a = dt*R(i)/2/L;
        I(i) = I(i)*(1-a)/(1+a) - dt/L*(V(i)-V(i-1))/dz/(1+a);
        V(i) = V(i) - dt/C*(I(i+1)-I(i))/dz;
    end
    I(1) = I(2);
    if iupdateR==1
        R(1:idamp) = ( (1-b)*R(1:idamp) + 2*b*min(R0vec(1:idamp),Ess./abs(I(1:idamp))) )/(1+b);
    end
    
    % Second half step
    for i=2:imax
        a = dt*(R(i) + Ro(i))/4/L;
        I(i) = Io(i)*(1-a)/(1+a) - dt/L*(V(i)-V(i-1))/dz/(1+a);
        V(i) = Vo(i) - dt/C*(I(i+1)-I(i))/dz;
    end
    I(1) = I(2);
    if iupdateR==1
        R(1:idamp) = ( (1-b)*Ro(1:idamp) + 2*b*min(R0vec(1:idamp),Ess./abs(I(1:idamp))) )/(1+b);
    end
    
    % Wave front position
    zfrntpos(k) = t(k)*v;
    if t(k)>tauA
        for iwp=imax:-1:1
            if sIp*I(iwp)>=0.1*max(sIp*I)
                break
            end
        end
        zIwavpos(k) = zI(iwp);
    end
    
    % Save variables as a function of both space and time to plot it
    if t(k)>=tsave(ksave) && saveIntFlag==1
                
        Isave(ksave,:) = I;
        Vsave(ksave,:) = V + Vi;
        Rsave(ksave,:) = R;
        
        % Update analytical sol.
        if iEvAnlyt==1
            for i=1:imax
                toI = t(k) - zI(i)/v;
                toV = t(k) - zV(i)/v;
                rhoo = sqrt(t(k)^2 - zI(i)^2/v^2)/tauR;
                
                Vasave(ksave,i) = Vgf(toV);    % Perfect conductor
                Iasave(ksave,i) = Vgf(toI)/Z;  % solution
                
                if iupdateR==0 && R0>0 % Constant resistance solution
                    expterm = exp(-t(k)/tauR);
                    besstermI = besseli(0,rhoo);
                    
                    % Current
                    Iasave(ksave,i) = Iasave(ksave,i)*expterm*besstermI;                
                    
                    % Potential
                    s = linspace(0,zV(i)-v*t(k),300);
                    ts = -s/v;
                    rho = sqrt(t(k)^2 - (zV(i) - s).^2/v^2)/tauR;                    
                    integ = (t(k).*mymbesseli1(rho)/tauR + besseli(0,rho)).*Vgf(ts);                                        
                    besstermV = -trapz(s,integ)/tauR/v*expterm;
                    
                    Vasave(ksave,i) = Vasave(ksave,i)*expterm + besstermV;
                    
                elseif iupdateR==1 % Particular case of nonlinear resistance solution
                    Vasave(ksave,i) = Vasave(ksave,i)*(1 - Ess*zV(i)/2/aVi).*(zV(i)<v*t(k));
                    Iasave(ksave,i) = Iasave(ksave,i)*(1 - Ess*t(k)*Z/L/2/aVi).*(zI(i)<v*t(k));
                    
                end
            end
        end
        
        % Save wavefront position
        zfrntsav(ksave) = zfrntpos(k);
        zIwavsav(ksave) = zIwavpos(k);
        
        % Update ksave
        ksave = ksave + 1;
        if ksave>Ntsave
            saveIntFlag = 0;
            ksave = Ntsave;
        end
    end
    
    % Save variables as a function of time only
    Imaxsave(k) = max(I);
    Igt(k) = I(1);
    
end


%% Rescale Vasave and Iasave to proper units
if iEvAnlyt==1
    Vasave = Vasave + Vi;
end


%% Set plot limits
Rmin = min(Rsave(:));
if Rmin==0
    Rmin = -1;
else
    Rmin = 0;
end
Rmax = max(Rsave(:));
if Rmax==0
    Rmax = 1;
else
    Rmax = 1.2*Rmax;
end
Rlims = [Rmin Rmax];

Vmin = min(Vsave(:));
Vmax = max(Vsave(:));
dV = 0.05*(Vmax-Vmin);
Vlims = [Vmin-dV Vmax+dV];

Imin = min(Isave(:));
Imax = max(Isave(:));
dI = 0.05*(Imax-Imin);
Ilims = [Imin-dI Imax+dI];


%% Save variables (drs is structure with important vars)
if iSavVars>0
    
    % Make directory (this can be changed by the user)
    if exist('./sims','dir')
        % do nothing else
    else
        mkdir 'sims';
    end
    
    % Save input parameters
    drs.params.dt   = dt;
    drs.params.dz   = dz;
    drs.params.rq   = rq;
    drs.params.rI   = rI;
    drs.params.R0   = R0;
    drs.params.L    = L;
    drs.params.C    = C;
    drs.params.Vi   = Vi;
    drs.params.Ip   = Ip; % Max theoretical peak current, not actual
    drs.params.tauA = tauA;
    drs.params.v    = v;
    drs.params.Z    = Z;
    drs.params.tauR = tauR;

    if iupdateR==1
        drs.params.tauRss = tauRss;
        drs.params.Ess = Ess;
    else
        drs.params.tauRss = Inf;
    end
    
    % Save flags
    drs.flags.iupdateR  = iupdateR;
    drs.flags.iEvAnlyt  = iEvAnlyt;
    drs.flags.iUseAbsBC = iUseAbsBC;
    drs.flags.iMakeMov  = iMakeMov;
    drs.flags.iSavVars  = iSavVars;
    drs.flags.caseName  = caseName;

    % Save variables
    drs.t = t;
    drs.tsave = tsave;
    drs.zI = zI;
    drs.zV = zV;
    drs.Vsave = Vsave;
    drs.Isave = Isave;
    drs.Rsave = Rsave;
    if iEvAnlyt==1 % Analytical solutions
        drs.Vasave = Vasave;
        drs.Iasave = Iasave;
    end
    drs.Igt = Igt;
    drs.zfrntpos = zfrntpos;
    drs.zIwavpos = zIwavpos;
    drs.zfrntsav = zfrntsav;
    drs.zIwavsav = zIwavsav;
    
    % Save var limits
    drs.Vlims = Vlims;
    drs.Ilims = Ilims;
    drs.Rlims = Rlims;
    
    save(['./sims/' caseName],'drs');
end


%% Initialize plots
Iunit = 1e-3; % Convert to kA
Vunit = 1e-6; % Convert to MV
tunit = 1e6;  % Convert to us
LFS   = 24;
FS    = 20;
LW    = 2.2;
MS    = 12;


%% Movie animation
if iMakeMov>0
    fh1 = figure(1); clf
    
    % Initialize file to save movie
    if iMakeMov>1
        % Make directory
        if exist('./sims','dir')
            % do nothing else
        else
            mkdir 'sims';
        end
        % Initialize movie
        movPanels = VideoWriter(['./sims/' caseName '_movie.mp4'],'MPEG-4');
        open(movPanels);
    end
    
    for ksave=1:Ntsave
        
        subplot(1,3,1)
        pI  = plot(Isave(ksave,:)*Iunit,zI,'-k','LineWidth',LW);
        if iEvAnlyt==1
            hold on
            pIa = plot(Iasave(ksave,:)*Iunit,zI,'--r','LineWidth',1.1*LW);
            hold off
        end
        xlim(Ilims*Iunit)
        ylim([0 zmax])
        if iEvAnlyt==1
            legend([pIa, pI],'Analytical','Numerical')
        end
        set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
        ylabel('Height (m)','FontSize',LFS)
        xlabel('Current (kA)','FontSize',LFS)
        
        subplot(1,3,2)
        pV = plot(Vsave(ksave,:)*Vunit,zV,'-k','LineWidth',LW);
        if iEvAnlyt==1
            hold on
            pVa = plot(Vasave(ksave,:)*Vunit,zV,'--r','LineWidth',1.1*LW);
            hold off
        end
        xlim(Vlims*Vunit)
        ylim([0 zmax])
        set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
        ylabel('Height (m)','FontSize',LFS)
        xlabel('Electric potential (MV)','FontSize',LFS)
        title(['t = ' sprintf('%4.1f',tsave(ksave)*tunit) ' \mus'],'FontWeight','Normal');
        
        subplot(1,3,3)
        plot(Rsave(ksave,:),zI,'-k','LineWidth',LW);
        hold on
        line(Rlims,tsave(ksave)*v*[1 1],'LineStyle','--','Color','b','LineWidth',.9*LW)
        hold off
        xlim(Rlims)
        ylim([0 zmax])
        set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
        ylabel('Height (m)','FontSize',LFS)
        xlabel('Resistance (\Omega/m)','FontSize',LFS)
        
        drawnow;
        
        if iMakeMov>1
            currFrame = getframe(fh1);
            writeVideo(movPanels,currFrame);
        else
            pause(.0001)
        end
        
    end
    
    % Close file with movie
    if iMakeMov>1
        close(movPanels);
    end
    
end


%% Static vertical plot
lcolors = [1 0 0; 0 .6 0; 0 0 1];
lgray = 0.7*[1 1  1];

if Ntsave<=10    
    f2 = figure(2); clf
    f2.Position = [100 100 1400 500];

    
    for ksave=1:Ntsave
        
        tsavestr{ksave} = [num2str(tsave(ksave)*tunit) ' \mus'];
        
        subplot(1,2,1)
        pwp = plot(interp1(zI,Isave(ksave,:)*Iunit,zIwavsav(ksave)),zIwavsav(ksave),'^','MarkerFaceColor',lcolors(ksave,:),'MarkerEdgeColor','k','MarkerSize',MS);
        
        if ksave==1
            hold on
        end
        pI  = plot(Isave(ksave,:)*Iunit,zI,'-','Color',lcolors(ksave,:),'LineWidth',1.9*LW);
        if iEvAnlyt==1
            pIa = plot(Iasave(ksave,:)*Iunit,zI,'--','Color',lgray,'LineWidth',1.1*LW);
        end
        pfp = plot((Imin-0.5*dI)*Iunit,zfrntsav(ksave),'o','MarkerFaceColor',lcolors(ksave,:),'MarkerEdgeColor','k','MarkerSize',MS);                
        
        if ksave==Ntsave
            hold off
            xlim(Ilims*Iunit)
            ylim([0 zmax])
            if iEvAnlyt==1
                legend([pIa, pI, pfp, pwp],'Analytical','Numerical','v*t','Wave front position')
            elseif iEvAnlyt==0
                legend([pI, pfp, pwp],'Calculation','v*t','Wave front position')
            end
            set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
            ylabel('Height (m)','FontSize',LFS)
            xlabel('Current (kA)','FontSize',LFS)
        end
        
        subplot(1,2,2)        
        pV(ksave)  = plot(Vsave(ksave,:)*Vunit,zV,'-','Color',lcolors(ksave,:),'LineWidth',1.9*LW);
        if ksave==1
            hold on
        end
        if iEvAnlyt==1
            plot(Vasave(ksave,:)*Vunit,zV,'--','Color',lgray,'LineWidth',1.1*LW);
        end
        if ksave==Ntsave
            hold off
            xlim(Vlims*Vunit)
            ylim([0 zmax])
            legend(pV,tsavestr)
            set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
            ylabel('Height (m)','FontSize',LFS)
            xlabel('Potential (MV)','FontSize',LFS)
        end
        
    end
end


%% Plot time series of current at ground
figure(3); clf
subplot(1,2,1)
plot(t*tunit,Igt*Iunit,'-k','LineWidth',LW)
% ylim([min(0,1.01*Ip) max(0,1.01*Ip)]*Iunit)
xlim([-0.05 1]*tmax*tunit)
set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
ylabel('Channel-base current (kA)','FontSize',LFS)
xlabel('Time (\mus)','FontSize',LFS)

subplot(1,2,2)
plot(t*tunit,zIwavpos,'-r','LineWidth',1.1*LW)
hold on
plot(t*tunit,zfrntpos,'--k','LineWidth',LW)
hold off
ylim([0 max(zfrntpos)])
xlim([-0.05 1]*tmax*tunit)
legend('Wave front','v*t','Location','SouthEast')
set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
ylabel('Position (m)','FontSize',LFS)
xlabel('Time (\mus)','FontSize',LFS)
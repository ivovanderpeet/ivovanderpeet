%% Solves: Steady, compressible convection-diffusion problems.
clear; close all; clc;

%% Global constants and variables
% Global contants
global NPI NPJ XMAX YMAX JBOT JMID JTOP HBOT HTOP HMID
global LARGE SMALL U_IN COFLOW Dy
global sigmak sigmaom gamma1 beta1 betastar Cmu kappa ERough Ti 

% Variables
global x y u v pc p T rho mu Gamma Cp 
global aP aE aW aN aS b d_u d_v SMAX SAVG relax_rho 
global k om
    
%% Configuration parameters
COFLOW = -1; % set to -1 for counterflow, or 1 for coflow

%% Constants
% Domain
NPI        = 100;       % number of grid cells in x-direction [-]
NPJ        = 40;        % number of grid cells in y-direction [-]
XMAX       = 0.30;      % width of the domain [m]
HBOT       = 0.01;
HTOP       = 0.01;
HMID       = 0.002;
YMAX       = HBOT + HMID + HTOP;  % height of the domain [m]

JBOT = 2:ceil((NPJ+1)/YMAX*HBOT);
JMID = ceil((NPJ+1)/YMAX*HBOT)+1:ceil((NPJ+1)/YMAX*(HBOT+HMID));
JTOP = ceil((NPJ+1)/YMAX*(HBOT+HMID))+1:NPJ+1;

% Iterations
MAX_ITER   = 1000;  % maximum number of outer iterations [-]
U_ITER     = 1;     % number of Newton iterations for u equation [-]
V_ITER     = 1;     % number of Newton iterations for u equation [-]
PC_ITER    = 200;   % number of Newton iterations for pc equation [-]
K_ITER     = 1;     % number of Newton iterations for k equation [-]
OM_ITER    = 1;     % number of Newton iterations for omega equation [-]
T_ITER     = 1;     % number of Newton iterations for T equation [-]

% Accuracy
SMAXneeded = 1E-8;      % maximum accepted error in mass balance [kg/s]
SAVGneeded = 1E-9;      % maximum accepted average error in mass balance [kg/s]
LARGE      = 1E30;      % arbitrary very large value [-]
SMALL      = 1E-30;     % arbitrary very small value [-]
NPRINT     = 1;         % number of iterations between printing output to screen

% Input constants
P_ATM      = 101000.;   % atmospheric pressure [Pa]
U_IN       = .2;      % in flow velocity [m/s]

%% Coefficients for wilcox k-omega model 
sigmak     = 2.;
sigmaom    = 2.;
gamma1     = 0.553;
beta1      = 0.075;
betastar   = 0.09;

Cmu        = 0.09;
kappa      = 0.4187;
ERough     = 9.793;
Ti         = 0.04;

%% Main calculations
init();  %call initialization function

iter = 1;
% outer iteration loop
while (iter <= MAX_ITER && SMAX > SMAXneeded && SAVG > SAVGneeded)
    
    %% u calculation
    ucoeff(); %call ucoeffe.m function to calculate the coefficients for u function
    for iter_u = 1:U_ITER
        u = solve(u, b, aE, aW, aN, aS, aP); %solve u function
    end
    
    %% v calculation
    vcoeff(); %call vcoeffe.m function to calculate the coefficients for v function
    for iter_v = 1:V_ITER
        v = solve(v, b, aE, aW, aN, aS, aP); %solve v function
    end
    
    %% Boundary conditions
    bound(); %apply boundary conditions 
    
    %% p calculation
    pccoeff(); %call pccoeffe.m function to calculate the coefficients for p function
    for iter_pc = 1:PC_ITER
        pc = solve(pc, b, aE, aW, aN, aS, aP); %solve p function
    end
    
    %% Pressure and velocity correction
    velcorr(); % Correct pressure and velocity
    
    %% k calculation
    kcoeff();
    for iter_k = 1:K_ITER
        k = solve(k, b, aE, aW, aN, aS, aP);
    end
    
    %% omega calculation
    omcoeff();
    for iter_om = 1:OM_ITER
        om = solve(om, b, aE, aW, aN, aS, aP);
    end

    %% T calculation
    Tcoeff(); %call Tcoeff.m function to calculate the coefficients for T function
    for iter_T = 1:T_ITER
        T = solve(T, b, aE, aW, aN, aS, aP); %solve T function
    end
    
    viscosity();
    bound();

    %% rho calculation
    for I = 1:NPI+2
        for J = 2:NPJ+1
            if (I == 1) % Since p(1, J) doesn't exist, we set:
                rho(I,J) = (1 - relax_rho)*rho(I,J) + relax_rho*(p(I+1,J) + P_ATM)/(287.*T(I,J));
            else
                rho(I,J) = (1 - relax_rho)*rho(I,J) + relax_rho*(p(I,J) + P_ATM)/(287.*T(I,J));
            end
        end
    end
    % end of density calculation======================================================================
    
    %% mu calculation
    mu(1:NPI+2,2:NPJ+1) = (0.0395*T(1:NPI+2,2:NPJ+1) + 6.58)*1.E-6;   
    
    %% Gamma calculation    
    for I = 1:NPI+2
        for J = 2:NPJ+1
            if max(J==JMID)
                Gamma(I,J) = 400;
            else
                Gamma(I,J) = (6.1E-5*T(I,J) + 8.4E-3)/Cp(I,J);
            end
            if (Gamma(I,J) < 0.)
                output();
                fprintf('Error: Gamma(%d,%d) = %e\n', I, J, Gamma(I,J));
                exit(1);
            end
        end
    end
    
    %% Print temporary results
    if iter == 1
        fprintf ('Iter.\t d_u/u\t\t d_v/v\t\t SMAX\t\t SAVG\n');
    end
    if mod(iter,NPRINT) == 0
        I = round((NPI+1)/2);
        J = round((NPJ+1)/4);
        du = d_u(I,J)*(pc(I-1,J) - pc(I,J));
        dv = d_v(I,J)*(pc(I,J-1) - pc(I,J));
        fprintf ('%3d\t%10.2e\t%10.2e\t%10.2e\t%10.2e\n', iter,du/u(I,J), dv/v(I,J), SMAX, SAVG);
    end

    iter = iter + 1; % increase interation number
end

%% visualize the velocity profile
[X,Y]=meshgrid(x,y);

figure(1)
quiver(X,Y, u', v',1.5); hold on; title("Velocity")
plot(x, (HBOT+HMID)*ones(1,length(x)), 'k');
plot(x, HBOT*ones(1,length(x)), 'k');
plot(x, YMAX*ones(1,length(x)), 'k');
xlim([0,XMAX])
% ylim([-XMAX/4,XMAX/4])
grid on

figure(2)
surf(X,Y,T'); colorbar; title("Temperature")
view(0,90); hold on
plot3(x, (HBOT+HMID)*ones(1,length(x)), max(max(T))*ones(1,length(x)), 'r');
plot3(x, HBOT*ones(1,length(x)), max(max(T))*ones(1,length(x)), 'r');
xlim([0,XMAX])
% ylim([-XMAX/4,XMAX/4])


figure(3)
surf(X,Y,p'); colorbar; title("Pressure")
view(0,90)
xlim([0,XMAX])
% ylim([-XMAX/4,XMAX/4])


% figure(4)
% surf(X,Y,Gamma'); colorbar; title("Gamma")
% view(0,90)
% 
% figure(5)
% surf(X,Y,Cp'); colorbar; title("Heat capacity")
% view(0,90)
% 
% figure(6)
% surf(X,Y,rho'); colorbar; title("Density")
% view(0,90)

%% Check balance
if COFLOW == 1
    m_in_bot = u(2,JBOT).*rho(2,JBOT)*Dy;
    m_out_bot = u(NPI+1,JBOT).*rho(NPI+1,JBOT)*Dy;

    Q_in_bot = sum(m_in_bot.*Cp(2,JBOT).*T(2,JBOT));
    Q_out_bot = sum(m_out_bot.*Cp(NPI+1,JBOT).*T(NPI+1,JBOT));
elseif COFLOW == -1
    m_in_bot = abs(u(NPI+1,JBOT).*rho(NPI+1,JBOT)*Dy);
    m_out_bot = abs(u(2,JBOT).*rho(2,JBOT)*Dy);

    Q_in_bot = sum(m_in_bot.*Cp(NPI+1,JBOT).*T(NPI+1,JBOT));
    Q_out_bot = sum(m_out_bot.*Cp(2,JBOT).*T(2,JBOT));
end

m_in_top = u(2,JTOP).*rho(2,JTOP)*Dy;
m_out_top = u(NPI+1,JTOP).*rho(NPI+1,JTOP)*Dy;

Q_in_top = sum(m_in_top.*Cp(2,JTOP).*T(2,JTOP));
Q_out_top = sum(m_out_top.*Cp(NPI+1,JTOP).*T(NPI+1,JTOP));

dQ_bot = Q_out_bot - Q_in_bot;
dQ_top = Q_out_top - Q_in_top;

Q_in = Q_in_top + Q_in_bot;
Q_out = Q_out_top + Q_out_bot;

fprintf("\ndQ_bot =%10.2e \t Q_in  =%10.2e \ndQ_top =%10.2e \t Q_out =%10.2e \n\n", dQ_bot, Q_in, dQ_top, Q_out)
%% Solves: Unsteady, compressible convection-diffusion problems.
% Description:
% This program solves unsteady convection-diffusion problems
% using the transient simple algorithm described in ch. 8.7.1 in "Computational
% Fluid Dynamics" by H.K. Versteeg and W. Malalasekera. Symbols and
% variables follow exactly the notations in this reference, and all
% equations cited are from this reference unless mentioned otherwise.

% Converted from C to Matlab by YTANG
% References: 1. Computational Fluid Dynamics, H.K. Versteeg and W. Malalasekera, Longman Group Ltd, 1995

clear
close all
clc
%% declare all variables and contants
% variables
global u v d_u d_v pc p T rho mu Gamma Cp
global b aP aE aW aN aS omega k
global u_old v_old pc_old T_old omega_old k_old 
global uplus yplus yplus1 yplus2  
% constants
global x x_u y y_v NPI NPJ XMAX YMAX JBOT JMID JTOP HBOT HMID HTOP Dy
global SMAX SAVG LARGE SMALL BIG
global sigmak sigmaom gamma1 beta1 betastar ERough Ti Cmu kappa 
global NPRINT COFLOW Dt U_IN

%% Configuration parameters
COFLOW = -1; % set to -1 for counterflow, or 1 for coflow

%% Constants
% Domain
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
MAX_ITER   = 1000;      % maximum number of outer iterations [-]
U_ITER     = 1;         % number of Newton iterations for u equation [-]
V_ITER     = 15;         % number of Newton iterations for v equation [-]
PC_ITER    = 30;        % number of Newton iterations for pc equation [-]
T_ITER     = 1;         % number of Newton iterations for T equation [-]
OMEGA_ITER = 1;         % number of Newton iterations for omega equation [-]
K_ITER     = 1;         % number of Newton iterations for K equation [-]
NPRINT     = 1;

% Accuracy
SMAXneeded = 1E-6;      % maximum accepted error in mass balance [kg/s]
SAVGneeded = 1E-7;      % maximum accepted average error in mass balance [kg/s]
LARGE      = 1E30;      % arbitrary very large value [-]
SMALL      = 1E-30;     % arbitrary very small value [-]
BIG        = 1E10;

% Input constants
P_ATM      = 101000.;   % athmospheric pressure [Pa]
U_IN       = 0.02;       % in flow velocity [m/s]

% k-epsilon
sigmak     = 2.;
sigmaom    = 2.;
gamma1     = 0.553;
beta1      = 0.075;
betastar   = 0.09;
ERough     = 9.793;
Ti         = 0.04;
Cmu        = 0.09;
kappa      = 0.4187;


Dt         = 0.05;
TOTAL_TIME = 100;

%% start main function here
init(); % initialization
bound(); % apply boundary conditions

f = waitbar(0,'Nou we gaan beginnen hoor');
TMID = zeros(length(x),round(TOTAL_TIME/Dt));

for time = Dt:Dt:TOTAL_TIME
    waitbar(time/TOTAL_TIME,f,'Even geduld pik');
    iter = 0;
    
    % outer iteration loop
    while iter < MAX_ITER && SMAX > SMAXneeded && SAVG > SAVGneeded
        
        derivatives();
        ucoeff();
        for iter_u = 1:U_ITER
            u = solve(u, b, aE, aW, aN, aS, aP);
        end
        
        vcoeff();
        for iter_v = 1:V_ITER
            v = solve(v, b, aE, aW, aN, aS, aP);
        end
        
        bound();
        
        pccoeff();
        for iter_pc = 1:PC_ITER
            pc = solve(pc, b, aE, aW, aN, aS, aP);
        end
        
        velcorr(); % Correct pressure and velocity
        
        kcoeff();
        for iter_k = 1:K_ITER
            k = solve(k, b, aE, aW, aN, aS, aP);
            if min(min(k)) < 0
                return;
            end
        end
        
        omcoeff();
        for iter_eps = 1:OMEGA_ITER
            omega = solve(omega, b, aE, aW, aN, aS, aP);
        end
        
        Tcoeff();
        for iter_T = 1:T_ITER
            T = solve(T, b, aE, aW, aN, aS, aP);
        end
        
        viscosity();
        bound();
        
        % begin:storeresults()=============================================
        % Store data at current time level in arrays for "old" data
        % To newly calculated variables are stored in the arrays ...
        % for old variables, which can be used in the next timestep       
        u_old(3:NPI+1,2:NPJ+1)   = u(3:NPI+1,2:NPJ+1);       
        v_old(2:NPI+1,3:NPJ+1)   = v(2:NPI+1,3:NPJ+1);        
        pc_old(2:NPI+1,2:NPJ+1)  = pc(2:NPI+1,2:NPJ+1);        
        T_old(2:NPI+1,2:NPJ+1)   = T(2:NPI+1,2:NPJ+1);        
        omega_old(2:NPI+1,2:NPJ+1) = omega(2:NPI+1,2:NPJ+1);
        k_old(2:NPI+1,2:NPJ+1)   = k(2:NPI+1,2:NPJ+1);
        
        % end: storeresults()==============================================
        
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
        
        % increase iteration number
        iter = iter +1;        
    end % end of while loop (outer interation)
    
    % begin: printConv(time,iter)=========================================
    % print convergence to the screen
    if time == Dt
        fprintf ('Iter\t Time\t u\t v\t T\t SMAX\t SAVG\n');
    end    
    fprintf ("%4d %10.3e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\n",iter,...
        time,u(ceil(3*(NPI+1)/10),ceil(2*(NPJ+1)/5)),v(ceil(3*(NPI+1)/10),ceil(2*(NPJ+1)/5)),...
        T(ceil(3*(NPI+1)/10),ceil(2*(NPJ+1)/5)), SMAX, SAVG);
%     end: printConv(time, iter)===========================================
    
    % reset SMAX and SAVG
    SMAX = LARGE;
    SAVG = LARGE;   

    TMID(:,round(time/Dt)) = mean(T(:,JMID)')';


end % end of calculation

close(f)

%% plot vector map
[X,Y]=meshgrid(x,y);
[Xu,Yv]=meshgrid(x_u,y_v);
quiver(Xu,Yv,u',v',0.05);

figure(1)
quiver(X,Y, u', v',1.5); hold on; title("Velocity")
plot(x, (HBOT+HMID)*ones(1,length(x)), 'k');
plot(x, HBOT*ones(1,length(x)), 'k');
plot(x, YMAX*ones(1,length(x)), 'k');
grid on

figure(2)
surf(X,Y,T'); colorbar; title("Temperature")
view(0,90); hold on
plot3(x, (HBOT+HMID)*ones(1,length(x)), max(max(T))*ones(1,length(x)), 'r');
plot3(x, HBOT*ones(1,length(x)), max(max(T))*ones(1,length(x)), 'r');
xlim([0,XMAX])

p(:,JMID) = 0;

figure(3)
surf(X,Y,p'); colorbar; title("Pressure")
hold on
plot3(x, (HBOT+HMID)*ones(1,length(x)), max(max(p))*ones(1,length(x)), 'r');
plot3(x, HBOT*ones(1,length(x)), max(max(p))*ones(1,length(x)), 'r');
view(0,90)
xlim([0,XMAX])

figure(4)
surf(X,Y,k'); colorbar; title("k")
view(0,90)
hold on
plot3(x, (HBOT+HMID)*ones(1,length(x)), max(max(k))*ones(1,length(x)), 'r');
plot3(x, HBOT*ones(1,length(x)), max(max(k))*ones(1,length(x)), 'r');
xlim([0,XMAX])
figure(5)
surf(X,Y,omega'); colorbar; title("omega")
view(0,90)
hold on
plot3(x, (HBOT+HMID)*ones(1,length(x)), max(max(omega))*ones(1,length(x)), 'r');
plot3(x, HBOT*ones(1,length(x)), max(max(omega))*ones(1,length(x)), 'r');

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

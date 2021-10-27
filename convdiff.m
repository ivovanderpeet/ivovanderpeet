%% Solves: Steady, compressible convection-diffusion problems.
clear; close all; clc;

%% Declare all constants and variables
% global contants
global NPI NPJ XMAX YMAX LARGE U_IN JBOT JMID JTOP Dy COFLOW

% variables
global x x_u y y_v u v pc p T rho mu Gamma Cp aP aE aW aN aS b d_u d_v  SMAX SAVG relax_rho 
    
% configuration parameters
COFLOW = -1; % set to -1 for counterflow, or 1 for coflow

% constants
NPI        = 96;        % number of grid cells in x-direction [-]
NPJ        = 48;        % number of grid cells in y-direction [-]
XMAX       = 0.24;      % width of the domain [m]
YMAX       = 0.12;      % height of the domain [m]
HBOT       = 0.05;
HTOP       = 0.05;
HMID       = YMAX - HBOT - HTOP;
MAX_ITER   = 1000;      % maximum number of outer iterations [-]
U_ITER     = 1;         % number of Newton iterations for u equation [-]
V_ITER     = 1;         % number of Newton iterations for u equation [-]
PC_ITER    = 200;       % number of Newton iterations for pc equation [-]
T_ITER     = 100;       % number of Newton iterations for T equation [-]
SMAXneeded = 1E-8;      % maximum accepted error in mass balance [kg/s]
SAVGneeded = 1E-9;     % maximum accepted average error in mass balance [kg/s]
LARGE      = 1E30;      % arbitrary very large value [-]
P_ATM      = 101000.;   % atmospheric pressure [Pa]
U_IN       = 0.02;      % in flow velocity [m/s]
NPRINT     = 1;         % number of iterations between printing output to screen

JBOT = 2:ceil((NPJ+1)/0.12*0.05);
JMID = ceil((NPJ+1)/0.12*0.05)+1:ceil((NPJ+1)/0.12*0.07);
JTOP = ceil((NPJ+1)/0.12*0.07)+1:NPJ+1;

%% main calculations
init();  %call initialization function

% This is part of the init function
for JJ = 1:NPJ+2
    if max(JJ == JBOT)
        u(:,JJ) = COFLOW*U_IN*(1.-(2.*(y(JJ)-HBOT/2.)/HBOT)^2); % inlet bot
    elseif max(JJ == JTOP)
        u(:,JJ) = U_IN*(1.-(2.*((y(JJ) - (HBOT+HMID))-HTOP/2.)/HTOP)^2); % inlet top
    else
        u(:,JJ) = zeros(size(u(:,JJ)));
    end
end

iter = 1;
% outer iteration loop
while (iter <= MAX_ITER && SMAX > SMAXneeded && SAVG > SAVGneeded)
    
%     bound(); %call boundary function
    
    ucoeff(); %call ucoeffe.m function to calculate the coefficients for u function
    for iter_u = 1:U_ITER
        u = solve(u, b, aE, aW, aN, aS, aP); %solve u function
    end
    
    vcoeff(); %call vcoeffe.m function to calculate the coefficients for v function
    for iter_v = 1:V_ITER
        v = solve(v, b, aE, aW, aN, aS, aP); %solve v function
    end
    
    bound(); %apply boundary conditions 
    % Setting the velocity at inlet
    % This is part of the bound script
    for J = 1:NPJ+2
        if max(J == JBOT)
            u(2,J) = COFLOW*U_IN*(1.-(2.*(y(J)-HBOT/2.)/HBOT)^2); % inlet bot
        elseif max(J == JTOP)
            u(2,J) = U_IN*(1.-(2.*((y(J) - (HBOT+HMID))-HTOP/2.)/HTOP)^2); % inlet top
        end
    end
    
    pccoeff(); %call pccoeffe.m function to calculate the coefficients for p function
    for iter_pc = 1:PC_ITER
        pc = solve(pc, b, aE, aW, aN, aS, aP); %solve p function
    end
    
    velcorr(); % Correct pressure and velocity
    
    Tcoeff(); %call Tcoeff.m function to calculate the coefficients for T function
    for iter_T = 1:T_ITER
        T = solve(T, b, aE, aW, aN, aS, aP); %solve T function
    end
    
    % begin: density()==============================================================================
    % Purpose: Calculate the density rho(I, J) in the fluid as a function of the ideal gas law. 
    % Note: rho at the walls are not needed in this case, and therefore not calculated.    
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
    
    % begin: viscosity()==============================================================================
    % Purpose: Calculate the viscosity in the fluid as a function of temperature. 
    % Max error in the actual temp. interval is 0.5%   
    mu(1:NPI+2,2:NPJ+1) = (0.0395*T(1:NPI+2,2:NPJ+1) + 6.58)*1.E-6;   
    % end of viscosity calculation======================================================================
    
    % begin: conductivity()===========================================================================
    % Purpose: Calculate the thermal conductivity in the fluid as a function of temperature.    
    for I = 1:NPI+2
        for J = 2:NPJ+1
            Gamma(I,J) = (6.1E-5*T(I,J) + 8.4E-3)/Cp(I,J);
            if (Gamma(I,J) < 0.)
                output();
                fprintf('Error: Gamma(%d,%d) = %e\n', I, J, Gamma(I,J));
                exit(1);
            end
        end
    end
    % end of thermal conductivity calculation========================================================
    
    % begin: printConv(iter)========================================================================
    % print temporary results
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
    % end of print temporaty results=================================================================
    
    % increase interation number
    iter = iter + 1;   
end
%% begin: output()
% print out results in files
% fp   = fopen('output.dat','w');
% str  = fopen('str.dat','w');
% velu = fopen('velu.dat','w');
% velv = fopen('velv.dat','w');
% 
% for I = 1:NPI+1
%     i = I;
%     for J = 2:NPJ+1
%         j = J;
%         ugrid = 0.5*(u(i,J)+u(i+1,J));
%         vgrid = 0.5*(v(I,j)+v(I,j+1));
%         fprintf(fp, '%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\n',...
%             x(I), y(J), ugrid, vgrid, p(I,J), T(I,J), rho(I,J), mu(I,J), Gamma(I,J));
%     end
%     fprintf(fp, '\n');
% end
% fclose(fp);
% 
% for I = 1:NPI+1
%     i = I;
%     for J = 2:NPJ+1
%         j = J;
%         stream = -(u(i,J+1)-u(i,J))/(y(J+1)-y(J))+(v(I+1,j)-v(I,j))/(x(I+1)-x(I));
%         fprintf(str, '%10.2e\t%10.2e\t%10.5e\n',x_u(i), y_v(j), stream);
%         fprintf(velu,'%10.2e\t%10.2e\t%10.5e\n',x_u(i), y(J)  , u(i,J));
%         fprintf(velv,'%10.2e\t%10.2e\t%10.5e\n',x(I)  , y_v(j), v(I,j));
%     end
%     fprintf(str, '\n');
%     fprintf(velu,'\n');
%     fprintf(velv,'\n');
% end
% 
% fclose(str);
% fclose(velu);
% fclose(velv);
%% visualize the velocity profile
[X,Y]=meshgrid(x,y);

figure(1)
quiver(X,Y, u', v',1.5); hold on; title("Velocity")
plot(x, 0.07*ones(1,length(x)), 'k');
plot(x, 0.05*ones(1,length(x)), 'k');
plot(x, 0.12*ones(1,length(x)), 'k');

figure(2)
surf(X,Y,T'); colorbar; title("Temperature")
view(0,90); hold on
plot3(x, 0.07*ones(1,length(x)), max(max(T))*ones(1,length(x)), 'r');
plot3(x, 0.05*ones(1,length(x)), max(max(T))*ones(1,length(x)), 'r');

figure(3)
surf(X,Y,p'); colorbar; title("Pressure")
view(0,90)

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
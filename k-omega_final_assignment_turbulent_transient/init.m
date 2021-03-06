function [] = init()
% Purpose: To initilise all parameters.

% constants
global NPI NPJ LARGE U_IN XMAX YMAX JBOT JMID JTOP COFLOW HBOT HMID HTOP Dy
% variables
global x x_u y y_v u v pc p T rho mu mut mueff Gamma Cp k omega delta E E2 yplus yplus1 ...
    yplus2 uplus tw b SP Su d_u d_v SMAX SAVG relax_u relax_v ...
    relax_pc relax_T aP aE aW aN aS F_u F_v u_old v_old pc_old T_old k_old ...
    omega_old dudx dudy dvdx dvdy gamma_k gamma_om kinematic_mu ...
    A_mu B_mu C_mu D_mu relax_rho U_IN_bot rho_ratio

% begin: memalloc()========================================================
% allocate memory for variables
x   = zeros(1,NPI+2);
x_u = zeros(1,NPI+2);
y   = zeros(1,NPJ+2);
y_v = zeros(1,NPJ+2);

u   = zeros(NPI+2,NPJ+2);
v   = zeros(NPI+2,NPJ+2);
pc  = zeros(NPI+2,NPJ+2);
p   = zeros(NPI+2,NPJ+2);
T   = zeros(NPI+2,NPJ+2);
rho = zeros(NPI+2,NPJ+2);
mu  = zeros(NPI+2,NPJ+2);
kinematic_mu = zeros(NPI+2,NPJ+2);
mut  = zeros(NPI+2,NPJ+2);
mueff  = zeros(NPI+2,NPJ+2);
Gamma = zeros(NPI+2,NPJ+2);
gamma_k = zeros(NPI+2,NPJ+2);
gamma_om = zeros(NPI+2,NPJ+2);
Cp  = zeros(NPI+2,NPJ+2);
k  = zeros(NPI+2,NPJ+2);
omega  = zeros(NPI+2,NPJ+2);
delta  = zeros(NPI+2,NPJ+2);
E  = zeros(NPI+2,NPJ+2);
E2  = zeros(NPI+2,NPJ+2);

u_old  = zeros(NPI+2,NPJ+2);
v_old  = zeros(NPI+2,NPJ+2);
pc_old = zeros(NPI+2,NPJ+2);
T_old  = zeros(NPI+2,NPJ+2);
k_old  = zeros(NPI+2,NPJ+2);
omega_old  = zeros(NPI+2,NPJ+2);

dudx   = zeros(NPI+2,NPJ+2);
dudy   = zeros(NPI+2,NPJ+2);
dvdx   = zeros(NPI+2,NPJ+2);
dvdy   = zeros(NPI+2,NPJ+2);

aP  = zeros(NPI+2,NPJ+2);
aE  = zeros(NPI+2,NPJ+2);
aW  = zeros(NPI+2,NPJ+2);
aN  = zeros(NPI+2,NPJ+2);
aS  = zeros(NPI+2,NPJ+2);
b   = zeros(NPI+2,NPJ+2);

SP  = zeros(NPI+2,NPJ+2);
Su  = zeros(NPI+2,NPJ+2);

F_u = zeros(NPI+2,NPJ+2);
F_v = zeros(NPI+2,NPJ+2);

d_u = zeros(NPI+2,NPJ+2);
d_v = zeros(NPI+2,NPJ+2);
% end of memory allocation=================================================

% begin: grid()===========================================================
% Purpose: Defining the geometrical variables See fig. 6.2-6.4 in ref. 1
% Length of volume element
Dx = XMAX/NPI;
Dy = YMAX/NPJ;

% Length variable for the scalar points in the x direction
x(1) = 0.;
x(2) = 0.5*Dx;
for I = 3:NPI+1
    x(I) = x(I-1) + Dx;
end
x(NPI+2) = x(NPI+1) + 0.5*Dx;

% Length variable for the scalar points T(i,j) in the y direction
y(1) = 0.;
y(2) = 0.5*Dy;
for J = 3:NPJ+1
    y(J) = y(J-1) + Dy;
end
y(NPJ+2) = y(NPJ+1) + 0.5*Dy;

% Length variable for the velocity components u(i,j) in the x direction
x_u(1) = 0.;
x_u(2) = 0.;
for i = 3:NPI+2
    x_u(i) = x_u(i-1) + Dx;
end

% Length variable for the velocity components v(i,j) in the y direction 
y_v(1) = 0.;
y_v(2) = 0.;
for j = 3:NPJ+2
    y_v(j) = y_v(j-1) + Dy;
end
% end of grid setting======================================================

% begin: init()===========================================================
% Initialising all other variables

% Initialize convergence parameters at large values
SMAX = LARGE;
SAVG = LARGE;

% m_in  = 1.;
% m_out = 1.;

% if COFLOW == 1
%     rho_ratio = (rho(2,JBOT)/JBOT)/(rho(2,JTOP)/JTOP);
% elseif COFLOW == -1
%     rho_ratio = rho(NPI+1,JBOT)./rho(NPI+1,JTOP);
% end

% for J = 1:NPJ+2
%     if max(J == JBOT)
%         u(:,J) = COFLOW*U_IN_bot*(1.-(2.*(y(J)-HBOT/2.)/HBOT)^2); % inlet bot
%     elseif max(J == JTOP)
%         u(:,J) = U_IN_bot*1*(1.-(2.*((y(J) - (HBOT+HMID))-HTOP/2.)/HTOP)^2); % inlet top
%     end
% end

v(:,:)     = 0.;       % Velocity in y-direction
p(:,:)     = 0.;       % Relative pressure
T(:,:)     = 273+10.;     % Temperature
T(:,JMID)  = 273+50.;
T(:,JTOP)  = 273+90.;
rho(:,:) = (999.83952+16.945176.*(T-273)-7.9870401*10^(-3)*(T-273).^2-46.170461*10^(-6)*(T-273).^3+105.56302*10^(-9)*(T-273).^4-280.54253*10^(-12)*(T-273).^5)./(1+16.897850*10^(-3)*(T-273));
rho(:,JMID) = 8960;
mu(:,:) = A_mu*exp((B_mu./T)+C_mu.*T+D_mu.*T.^2)*10^(-3);
kinematic_mu(:,:) = mu./rho;

Cp(:,:) = 4200;
Cp(:,JMID) = 385;
Gamma = 0.598./Cp;
Gamma(:,JMID) = 400./Cp(:,JMID);
k(:,:)     = 1e-3;     % k
omega(:,:) = 1e-4;     % omega

rho_ratio = rho(2,2)/rho(2,NPJ+1);
for J = 1:NPJ+2
    if max(J == JBOT)
        u(:,J) = COFLOW*U_IN_bot*(1.-(2.*(y(J)-HBOT/2.)/HBOT)^2); % inlet bot
    elseif max(J == JTOP)
        u(:,J) = U_IN_bot*rho_ratio*(1.-(2.*((y(J) - (HBOT+HMID))-HTOP/2.)/HTOP)^2); % inlet top
    end
end

u_old      = u;        % Velocity in x-direction old timestep
v_old      = v;        % Velocity in y-direction old timestep
pc_old     = pc;       % Pressure correction old timestep
T_old      = T;        % Temperature old timestep
omega_old  = omega;      % epsilon old timestep
k_old      = k;        % k old timestep

% Setting the relaxation parameters
relax_u   = 0.8;            % See eq. 6.36
relax_v   = relax_u;        % See eq. 6.37
relax_pc  = 1.1 - relax_u;  % See eq. 6.33
relax_T   = 1.0;            % Relaxation factor for temperature
relax_rho = 0.1;
% end of initilization=====================================================
end


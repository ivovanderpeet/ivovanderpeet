function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ JBOT JMID JTOP HBOT HTOP HMID
global U_IN COFLOW 
global Cmu Ti
% variables
global y u v T y_v y F_u
global m_in_TOP m_out_TOP m_in_BOT m_out_BOT
global k eps

%% Inflow properties
% Temperature
if COFLOW == -1
    T(NPI+2,JBOT) = 273;
elseif COFLOW == 1
    T(1,JBOT) = 273;
end
T(1,JTOP) = 573;

% Velocity profile
u(2,JBOT) = COFLOW*U_IN*(1.-(2.*(y(JBOT)-HBOT/2.)/HBOT).^2);
u(2,JTOP) = U_IN*(1.-(2.*((y(JTOP) - (HBOT+HMID))-HTOP/2.)/HTOP).^2);

% Set k and eps at the inlet
k(1,JBOT)     = 1.5*(U_IN*Ti)^2; % at inlet
k(1,JTOP)     = 1.5*(U_IN*Ti)^2;
% k(:,JMID) = 0;

eps(1,JBOT)   = Cmu^0.75 *k(1,JBOT).^1.5/(0.07*HBOT*0.5); % at inlet
eps(1,JTOP)   = Cmu^0.75 *k(1,JTOP).^1.5/(0.07*HTOP*0.5); % at inlet
% eps(:,JMID) = 0;

%% Outer wall boundary (adiabatic)
T(:, NPJ+2) = T(:, NPJ+1);
T(:,1) = T(:,2);

%% Continuity equations (mass balance)
convect();

m_in_TOP = 0.;
m_out_TOP = 0.;
m_in_BOT = 0.;
m_out_BOT = 0.;

for J = min(JTOP):max(JTOP)
    j = J;
    AREAw = y_v(j+1) - y_v(j); % See fig. 6.3
    m_in_TOP  = m_in_TOP  + F_u(2,J)*AREAw;
    m_out_TOP = m_out_TOP + F_u(NPI+1,J)*AREAw;
end
for J = min(JBOT):max(JBOT)
    j = J;
    AREAw = y_v(j+1) - y_v(j); % See fig. 6.3
    m_in_BOT  = m_in_BOT  + F_u(2,J)*AREAw;
    m_out_BOT = m_out_BOT + F_u(NPI+1,J)*AREAw;
end


%% Corection varibles: Correction factor m_in/m_out is used to satisfy global continuity
% TOP
u(NPI+2,JTOP) = u(NPI+1,JTOP)*m_in_TOP/m_out_TOP;
v(NPI+2,JTOP) = v(NPI+1,JTOP);
T(NPI+2,JTOP) = T(NPI+1,JTOP);
k(NPI+2,JTOP) = k(NPI+1,JTOP);
eps(NPI+2,JTOP) = eps(NPI+1,JTOP);

% BOT
u(NPI+2,JBOT) = u(NPI+1,JBOT)*m_in_BOT/m_out_BOT;
v(NPI+2,JBOT) = v(NPI+1,JBOT);
if COFLOW == -1
    T(1,JBOT) = T(2,JBOT);
    k(1,JBOT) = k(2,JBOT);
    eps(1,JBOT) = eps(2,JBOT);
elseif COFLOW == 1
    T(NPI+2,JBOT) = T(NPI+1,JBOT);
    k(NPI+2,JBOT) = k(NPI+1,JBOT);
    eps(NPI+2,JBOT) = eps(NPI+1,JBOT);
end

% MID (Otherwise both ends of the wall will always have the initial temperature)
T(NPI+2,JMID) = T(NPI+1,JMID);
T(1,JMID)     = T(2,JMID);
end

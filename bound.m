function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ JBOT JMID JTOP HBOT HTOP HMID
global U_IN COFLOW 
global Cmu Ti
% variables
global u v T y_v y F_u
global m_in_TOP m_out_TOP m_in_BOT m_out_BOT
global k om

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

% Turbulent model parameters
k(1,:)  = 2/3*(U_IN*Ti)^2; % Waar komt 1.5 vandaan, en deze equation in general
om(1,JBOT)   = Cmu^(3/4) * k(1,JBOT).^0.5 / (0.07*HBOT); % Page 77 of book
om(1,JTOP)   = Cmu^(3/4) * k(1,JTOP).^0.5 / (0.07*HTOP); % Page 77 of book

%% Outer wall boundary (adiabatic)
T(:, NPJ+2) = T(:, NPJ+1);
k(:,NPJ+2) = k(:,NPJ+1);
om(:,NPJ+2) = om(:,NPJ+1);

T(:,1) = T(:,2);
k(:,1) = k(:,2);
om(:,1) = om(:,2);

%% TEMPORARY EXPERIMENTAL Inner wall boundary (perfect conductor)
% T(:,max(JBOT)) = T(:,min(JTOP));
% T(:,max(JMID)) = T(:,min(JTOP));
% T(:,max(JMID)) = T(:,max(JBOT));

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
om(NPI+2,JTOP) = om(NPI+1,JTOP);

% BOT
u(NPI+2,JBOT) = u(NPI+1,JBOT)*m_in_BOT/m_out_BOT;
v(NPI+2,JBOT) = v(NPI+1,JBOT);
if COFLOW == -1
    T(1,JBOT) = T(2,JBOT);
    k(1,JBOT) = k(2,JBOT);
    om(1,JBOT) = om(2,JBOT);
elseif COFLOW == 1
    T(NPI+2,JBOT) = T(NPI+1,JBOT);
    k(NPI+2,JBOT) = k(NPI+1,JBOT);
    om(NPI+2,JBOT) = om(NPI+1,JBOT);
end

% MID (Otherwise both ends of the wall will always have the initial temperature)
T(NPI+2,JMID) = T(NPI+1,JMID);
T(1,JMID)     = T(2,JMID);

end


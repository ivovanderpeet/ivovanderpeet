function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ U_IN JBOT JTOP JMID
% variables
global  u v T y_v F_u v m_in_TOP m_out_TOP m_in_BOT m_out_BOT

% Fixed temperature in Kelvin of the incoming fluid
T(1,JBOT) = 273;
T(1,JTOP) = 573;

% Setting the velocity at inlet
u(2,JBOT) = U_IN;
u(2,JTOP) = U_IN;

% u(:,JMID) = 0;
% v(:,JMID) = 0;

% Zero temperature gradient at outer walls
T(:, NPJ+2) = T(:, NPJ+1);
T(:, 1) = T(:, 2);

% begin: globcont(): Velocity and temperature gradient at outlet = zero:
% Purpose: Calculate mass in and out of the calculation domain to correct for the continuity at outlet.
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
% end: globcont():

% corection varibles: Correction factor m_in/m_out is used to satisfy global continuity
u(NPI+2,JTOP) = u(NPI+1,JTOP)*m_in_TOP/m_out_TOP;
u(NPI+2,JBOT) = u(NPI+1,JBOT)*m_in_BOT/m_out_BOT;

v(NPI+2,JTOP) = v(NPI+1,JTOP);
v(NPI+2,JBOT) = v(NPI+1,JBOT);

T(NPI+2,JTOP) = T(NPI+1,JTOP);
T(NPI+2,JMID) = T(NPI+1,JMID);  % Otherwise both ends of the wall will always have the initial temperature
T(1,JMID)     = T(2,JMID);          % Otherwise both ends of the wall will always have the initial temperature
T(NPI+2,JBOT) = T(NPI+1,JBOT);
end


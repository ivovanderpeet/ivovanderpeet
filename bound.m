function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ U_IN JBOT JTOP JMID
% variables
global  u v T m_in m_out y_v F_u v

% Fixed temperature in Kelvin of the incoming fluid
T(1,JBOT) = 273;
T(1,JTOP) = 573;

% Setting the velocity at inlet
u(2,JBOT) = U_IN;
u(2,JTOP) = 2*U_IN;

% u(:,JMID) = 0;
% v(:,JMID) = 0;

% Zero temperature gradient at outer walls
T(:, NPJ+2) = T(:, NPJ+1);
T(:, 1) = T(:, 2);

% begin: globcont(): Velocity and temperature gradient at outlet = zero:
% Purpose: Calculate mass in and out of the calculation domain to correct for the continuity at outlet.
convect();

m_in = 0.;
m_out = 0.;

for J = 2:NPJ+1
    j = J;
    AREAw = y_v(j+1) - y_v(j); % See fig. 6.3
    m_in  = m_in  + F_u(2,J)*AREAw;
    m_out = m_out + F_u(NPI+1,J)*AREAw;
end
% end: globcont():

% corection varibles: Correction factor m_in/m_out is used to satisfy global continuity
u(NPI+2,2:NPJ+1) = u(NPI+1,2:NPJ+1)*m_in/m_out;
v(NPI+2,2:NPJ+1) = v(NPI+1,2:NPJ+1);
T(NPI+2,2:NPJ+1) = T(NPI+1,2:NPJ+1);
end


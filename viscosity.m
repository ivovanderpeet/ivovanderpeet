function [] = viscosity()
% Purpose: To calculate the viscosity in the fluid as a function of temperature.

% constants
global NPI NPJ SMALL
% variables
global rho k mu mut mueff om

for I = 2:NPI+1
    for J = 2:NPJ+1
        mut(I,J)   = rho(I,J)*k(I,J)/(om(I,J)+SMALL);   % Eqn 3.71 of page 90
        mueff(I,J)  = mu(I,J) + mut(I,J);
    end
end
end
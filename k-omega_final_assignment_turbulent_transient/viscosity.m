function [] = viscosity()
% Purpose: To calculate the viscosity in the fluid as a function of temperature.

% constants
global NPI NPJ SMALL
% variables
global rho k mu mut mueff omega

for I = 1:NPI+1
    for J = 2: NPJ+2
        mut(I,J)   = rho(I,J)*k(I,J)/(omega(I,J)+SMALL);
        mueff(I,J) = mu(I,J) + mut(I,J);
    end
end
end
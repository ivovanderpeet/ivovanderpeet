function [] = viscosity()
% Purpose: To calculate the viscosity in the fluid as a function of temperature.

% constants
global NPI NPJ SMALL sigmaom sigmak
% variables
global rho k mu mut mueff omega gamma_om gamma_k

for I = 1:NPI+1
    for J = 2: NPJ+2
        mut(I,J)   = rho(I,J)*k(I,J)/(omega(I,J)+SMALL);
        mueff(I,J) = mu(I,J) + mut(I,J);
        gamma_om(I,J) = mut(I,J)/sigmaom + mu(I,J);
        gamma_k(I,J) = mut(I,J)/sigmak + mu(I,J);
    end
end
end
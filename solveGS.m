function fi = solveGS(fi,b,aE,aW,aN,aS,aP)
% Purpose: To solve the algebraic equation 7.7. using Gauss-Seidel

% variables
global Istart Iend Jstart Jend

for I = Istart:Iend
    for J = Jstart:Jend
        fi(I,J) = (aE(I,J)*fi(I+1,J) + aW(I,J)*fi(I-1,J) ...
            + aN(I,J)*fi(I,J+1) + aS(I,J)*fi(I,J-1) + b(I,J))/aP(I,J);
    end
end
end

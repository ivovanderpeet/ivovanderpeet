function [] = kcoeff()
% Purpose: To calculate the coefficient for the k equation.

% constants
global NPI NPJ Cmu sigmak betastar

% variables
global x x_u y y_v SP Su F_u F_v mut rho u uplus tw Istart Iend ...
    Jstart Jend b aE aW aN aS aP k E2 mu om

Istart = 2;
Iend = NPI+1;
Jstart = 2;
Jend = NPJ+1;

convect();
viscosity();
calculateuplus();

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;        
        % Geometrical parameters: Areas of the cell faces
        AREAw = y_v(j+1) - y_v(j); % = A(i,J) See fig. 6.3
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i); % = A(I,j)
        AREAn = AREAs;
        
        % eq. 6.9a-6.9d - the convective mass flux defined in eq. 5.8a
        % note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition.
        Fw = F_u(i,J)*AREAw;
        Fe = F_u(i+1,J)*AREAe;
        Fs = F_v(I,j)*AREAs;
        Fn = F_v(I,j+1)*AREAn;
        
        % eq. 6.9e-6.9h - the transport by diffusion defined in eq. 5.8b
        % note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition
        % The conductivity, Gamma, at the interface is calculated with the use of a harmonic mean.
        Dw = ((mu(I,J)+(mut(I-1,J)*mut(I,J)))/sigmak)/(mut(I-1,J)*(x(I) - x_u(i)) + ...
            mut(I,J)*(x_u(i)-x(I-1)))*AREAw;
        De = ((mu(I,J)+(mut(I,J)*mut(I+1,J)))/sigmak)/(mut(I,J)*(x(I+1) - x_u(i+1)) + ...
            mut(I+1,J)*(x_u(i+1)-x(I)))*AREAe;
        Ds = ((mu(I,J)+(mut(I,J-1)*mut(I,J)))/sigmak)/(mut(I,J-1)*(y(J) - y_v(j)) + ...
            mut(I,J)*(y_v(j)-y(J-1)))*AREAs;
        Dn = ((mu(I,J)+(mut(I,J)*mut(I,J+1)))/sigmak)/(mut(I,J)*(y(J+1) - y_v(j+1)) + ...
            mut(I,J+1)*(y_v(j+1)-y(J)))*AREAn;
        
        % The source terms
        if J==2 || J==NPJ+1
            SP(I,J) = -rho(I,J)*Cmu^0.75*k(I,J)^0.5*uplus(I,J)/(0.5*AREAw)*AREAs*AREAw;
            Su(I,J) = tw(I,J)*0.5*(u(i,J) + u(i+1,J))/(0.5*AREAw)*AREAs*AREAw;
        else
            SP(I,J) = -rho(I,J)*betastar*k(I,J)*om(I,J);
            Su(I,J) = (2.0*mut(I,J)*E2(I,J));
        end
        
        Su(I,J) =  Su(I,J)*AREAw*AREAs;
        SP(I,J) =  SP(I,J)*AREAw*AREAs;
        
        % The coefficients (hybrid differencing scheme)
        aW(I,J) = max([ Fw, Dw + Fw/2, 0.]);
        aE(I,J) = max([-Fe, De - Fe/2, 0.]);
        if J==2
            aS(I,J) = 0;
        else
            aS(I,J) = max([ Fs, Ds + Fs/2, 0.]);
        end
        
        if J==NPJ+1
            aN(I,J) = 0;
        else
            aN(I,J) = max([-Fn, Dn - Fn/2, 0.]);
        end
                
        % eq. 8.31 with time dependent terms (see also eq. 5.14):
        aP(I,J) = aW(I,J) + aE(I,J) + aS(I,J) + aN(I,J) + Fe - Fw + Fn - Fs - SP(I,J);
        
        % setting the source term equal to b
        b(I,J) = Su(I,J);
        
        % now the TDMA algorithm can be called to solve the equation.
        % This is done in the next step of the main program.       
    end
end
end

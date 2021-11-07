function [] = omcoeff()
% Purpose: To calculate the coefficients for the epsilon.

% constants
global NPI NPJ JBOT JMID JTOP LARGE SMALL BIG
global Dt gamma1 beta1 Dy
% variables
global x x_u y y_v SP Su F_u F_v rho Istart Iend ...
    Jstart Jend b aE aW aN aS aP omega E2 omega_old delta gamma_om

Istart = 2;
Iend = NPI+1;
Jstart = 2;
Jend = NPJ+1;

convect();
viscosity();

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
        Dw = gamma_om(I-1,J)*gamma_om(I,J)/(gamma_om(I-1,J)*(x(I) - x_u(i)) + ...
            gamma_om(I,J)*(x_u(i)-x(I-1)))*AREAw;
        De = gamma_om(I,J)*gamma_om(I+1,J)/(gamma_om(I,J)*(x(I+1) - x_u(i+1)) + ...
            gamma_om(I+1,J)*(x_u(i+1)-x(I)))*AREAe;
        Ds = gamma_om(I,J-1)*gamma_om(I,J)/(gamma_om(I,J-1)*(y(J) - y_v(j)) + ...
            gamma_om(I,J)*(y_v(j)-y(J-1)))*AREAs;
        Dn = gamma_om(I,J)*gamma_om(I,J+1)/(gamma_om(I,J)*(y(J+1) - y_v(j+1)) + ...
            gamma_om(I,J+1)*(y_v(j+1)-y(J)))*AREAn;

        SP(I,J) = -beta1*rho(I,J)*omega(I,J);
        Su(I,J) = gamma1*(2.0*rho(I,J)*E2(I,J) - 2/3*rho(I,J)*omega(I,J)*delta(I,J));
        
        Su(I,J) =  Su(I,J)*AREAw*AREAs;
        SP(I,J) =  SP(I,J)*AREAw*AREAs;

        if max(J == JMID)
            SP(I,J) = -LARGE;
            Su(I,J) = LARGE*BIG;
        end

        if J == 2 || J == NPJ+1
            SP(I,J) = -LARGE;
            Su(I,J) = LARGE*6*1e-6/(beta1*(Dy/2)^2); % Kinematic viscosity v = 1e-6 (at 20C)

        elseif isequal(J,max(JBOT)) || isequal(J,min(JTOP))
            SP(I,J) = -LARGE;
            Su(I,J) = LARGE*6*1e-6/(beta1*Dy^2); % Kinematic viscosity v = 1e-6 (at 20C)

        end

        % The coefficients (hybrid differencing scheme)
        aW(I,J) = max([ Fw, Dw + Fw/2, 0.]);
        aE(I,J) = max([-Fe, De - Fe/2, 0.]);
        aS(I,J) = max([ Fs, Ds + Fs/2, 0.]);
        aN(I,J) = max([-Fn, Dn - Fn/2, 0.]);
        aPold   =  rho(I,J)*AREAe*AREAn/Dt;

        % eq. 8.31 with time dependent terms (see also eq. 5.14):
        aP(I,J) = aW(I,J) + aE(I,J) + aS(I,J) + aN(I,J) + Fe - Fw + Fn - Fs - SP(I,J) + aPold;

        % setting the source term equal to b
        b(I,J) = Su(I,J) + aPold*omega_old(I,J);

        % now the TDMA algorithm can be called to solve the equation.
        % This is done in the next step of the main program.
    end
end
end

function [] = init_u()
global u JBOT JTOP HBOT HMID COFLOW U_IN HTOP NPJ y

for JJ = 1:NPJ+2
    if max(JJ == JBOT)
        u(:,JJ) = COFLOW*U_IN*(1.-(2.*(y(JJ)-HBOT/2.)/HBOT)^2); % inlet bot
    elseif max(JJ == JTOP)
        u(:,JJ) = U_IN*(1.-(2.*((y(JJ) - (HBOT+HMID))-HTOP/2.)/HTOP)^2); % inlet top
    else
        u(:,JJ) = zeros(size(u(:,JJ)));
    end
end

end
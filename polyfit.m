p = polyfit(y,u,2);


y1 = 0:0.0001:0.01;
u1 = polyval(p,y1);

u2 = 0.1*(1.-(2.*(y1 - 0.01/2.)/0.01).^2);

figure
plot(y,u,'o')
hold on
plot(y1,u1)
plot(y1,u2)
hold off
grid on
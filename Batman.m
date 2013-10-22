%% Draw Batman
x = -7.5:.001:7.5;
y1=1.5.*sqrt((-abs(abs(x) - 1)) .* abs(3 - abs(x))./((abs(x) - 1).*(3 - abs(x)))) .* (1+abs(abs(x) - 3)./(abs(x)- 3)) .* sqrt(1 - (x./7).^2)+(4.5+0.75 .* (abs(x - 0.5)+abs(x+0.5)) - 2.75 .* (abs(x-0.75)+abs(x+0.75))) .* (1+abs(1 - abs(x))./(1 - abs(x)));
y2=(-3).*sqrt(1 -(x./7).^2) .* sqrt(abs(abs(x) - 4)./(abs(x)-4));
y3=abs(x./2) - 0.0913722 .* x.^2-3+sqrt(1 - (abs(abs(x) - 2) - 1).^2);
y4=(2.71052+1.5 - 0.5 .* abs(x) - 1.35526 .* sqrt(4 - (abs(x) - 1).^2)) .* sqrt(abs(abs(x) - 1)./(abs(x) - 1));
figure(1);
hold on;
axis equal
plot(x,y1);
plot(x,y2);
plot(x,y3);
plot(x,y4);
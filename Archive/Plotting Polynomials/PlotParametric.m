
function [] = PlotParametric()

t = linspace (0,100);
x = sin(t);
y = cos(t);
z = t;

figure()
hold on
plot3(x,y,z,'-s')
hold off

xt = @(t) sin(t);
yt = @(t) cos(t);
zt = @(t) t;
figure()
hold on
fplot3(xt,yt,zt)
hold off

end
%% E谱画图(风域X=30km，风向角0，U10=5：5：30）
clear
kmin = 1e-5;
kmax = 1e5;
k = logspace(log10(kmin),log10(kmax),1e5);
X= 30e3;
X0= 2.2e4;
g = 9.81;
figure(1);
for U10 = 5:5:30
    X_= X.*g./(U10^2);
    age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
    [S,kp] = Elfouhaily(k,U10,age);
    loglog(k,S);
    hold on;
    grid on;
end
% line([370,370],[1e-15,1e10],'Color','black','LineWidth',1.3);
axis([1e-5 1e5 1e-15 1e10]);
legend("U10 = 5 m/s","U10 = 10 m/s","U10 = 15 m/s","U10 = 20 m/s","U10 = 25 m/s","U10 = 30 m/s");
xlabel("log(k)");
ylabel("log(S(k))");
savefig(gcf,"result/1-Espectrum");
exportgraphics(gcf,'result/1-Espectrum.pdf');
exportgraphics(gcf,'result/1-Espectrum.png');
%% E谱画图(风速U10=5，风向角0，X=30：80km） 与论文结论相反 【要不要再说】【不要】
clear
kmin = 1e-5;
kmax = 1e5;
k = logspace(log10(kmin),log10(kmax),1e5);
U10 = 5;
X0= 2.2e4;
g = 9.81;
figure(2);
%  for age = 0.1:0.1:0.84
for X = 30e3:10e3:80e3
    X_= X.*g./(U10^2);
%     age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
    age = 0.84*(tanh((X_./X0).^(0.4))).^(-0.75);
    [S,kp] = Elfouhaily(k,U10,age);
    loglog(k,S);
    hold on;
    grid on;
end
% line([370,370],[1e-15,1e10],'Color','black','LineWidth',1.3);
axis([1e-2 1e4 1e-15 1e0]);
legend("X = 30 m/s","X = 40 m/s","X = 50 m/s","X = 60 m/s","X = 70 m/s","X = 80 m/s");
xlabel("log(k)");
ylabel("log(S(k))");
savefig(gcf,"result/1-Espectrum-wind");
exportgraphics(gcf,'result/1-Espectrum-wind.pdf');
exportgraphics(gcf,'result/1-Espectrum-wind.png');
%% E谱 曲率 画图（K^3s）(风域X=30km，风向角0，U10=5：5：30）
clear
kmin = 1e-5;
kmax = 1e5;
k = logspace(log10(kmin),log10(kmax),1e5);
X= 30e3;
X0= 2.2e4;
g = 9.81;
figure(1);
for U10 = 5:5:30
    X_= X.*g./(U10^2);
    age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
    [S,kp] = Elfouhaily(k,U10,age);
    loglog(k,S.*k.^3);
    hold on;
    grid on;
end
line([370,370],[1e-15,1e10],'Color','black','LineWidth',1.3);
axis([1e-5 1e5 1e-4 1e0]);
legend("U10 = 5 m/s","U10 = 10 m/s","U10 = 15 m/s","U10 = 20 m/s","U10 = 25 m/s","U10 = 30 m/s","k = 370 rad/s");
xlabel("log(k)");
ylabel("log(k^3S(k))");
savefig(gcf,"result/2-EspectrumCurv");
exportgraphics(gcf,'result/2-EspectrumCurv.pdf');
%% E谱角度分布函数(K=0.3,U10=5,X=30)
clear
k = 0.3;
X= 30e3;
U10 = 5;
X0= 2.2e4;
g = 9.81;
km = 370.0; % rad/s
cm = 0.23; %minimum phase speed at wavenumber km?
k0 = g/(U10^2);%?
X_= X.*g./(U10^2);
age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
kp = k0 * age^2; %wavenumber of the spectral peak?
cp = sqrt(g/kp); %phase speed at the spectral peak cp = U10/age?
c = sqrt((g./k).*(1 + (k/km).^2)); %wave phase speed 波浪相速度
Cd10N = 0.00144; %drag coefficient？？
ustar = sqrt(Cd10N)*U10;%friction velocity at the water surface，论文2-18下面uf，为啥？

% 计算
a0 = log(2)/4;%改动，原本为log2/2
ap = 4;
am = 0.13*ustar/cm;
Delk = tanh(a0 + ap*(c/cp).^(2.5) + am*(cm./c).^(2.5));
phi = linspace(0, 2*pi, 100);
PSI = 1/(2*pi).*(1 + Delk.*cos(2*phi));%论文2-24关系

figure(1);
polarplot(phi, PSI, 'LineWidth', 1.5);
hold on;
PSI = 1/(2*pi).*(1 + Delk.*cos(2*(phi-pi/2)));%论文2-24关系
polarplot(phi, PSI, 'LineWidth', 1.5);
legend("\phi_w = 0^°/180^°","\phi_w = 90^°");
savefig(gcf,"result/3-Eforward");
exportgraphics(gcf,'result/3-Eforward.pdf');
exportgraphics(gcf,'result/3-Eforward.png');
%% E谱二维图(phi = 0 and phi = pi/2）
clear
U10 = 5;
X0= 2.2e4;
X= 30e3;
g = 9.81;
X_= X.*g./(U10^2);
age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);

figure(1);
phi_w = 0;
kmin = 1e-5;
kmax = 1e5;
K = logspace(log10(kmin),log10(kmax),256);
PHI = linspace(0,2*pi,256);
[k, phi] = meshgrid(K, PHI);

PSI = Elfouhaily2D(k,phi,U10,age,phi_w);
contour(k, phi, PSI);
clabel(contour(k, phi, PSI), 'FontSize', 10);
axis([0 0.5 0 pi]);
xlabel('k $(m^{-1})$', 'Interpreter', 'latex');
ylabel("\phi (rad)");
savefig(gcf,"result/4-Espectrum2D");
exportgraphics(gcf,'result/4-Espectrum2D.pdf');
exportgraphics(gcf,'result/4-Espectrum2D.png');

figure(2);
phi_w = pi/2;
PSI = Elfouhaily2D(k,phi,U10,age,phi_w);
contour(k, phi, PSI);
clabel(contour(k, phi, PSI), 'FontSize', 10);
axis([0 0.5 0 pi]);
xlabel('k $(m^{-1})$', 'Interpreter', 'latex');
ylabel("\phi (rad)");
savefig(gcf,"result/4-Espectrum2D90");
exportgraphics(gcf,'result/4-Espectrum2D90.pdf');
exportgraphics(gcf,'result/4-Espectrum2D90.png');

figure(3);
phi_w = pi;
PSI = Elfouhaily2D(k,phi,U10,age,phi_w);
contour(k, phi, PSI);
clabel(contour(k, phi, PSI), 'FontSize', 10);
axis([0 0.5 0 pi]);
xlabel('k $(m^{-1})$', 'Interpreter', 'latex');
ylabel("\phi (rad)");
savefig(gcf,"result/4-Espectrum2D180");
exportgraphics(gcf,'result/4-Espectrum2D180.pdf');
exportgraphics(gcf,'result/4-Espectrum2D180.png');
%% 二维海面图风速【4个图】（海面尺度500m，N=2L，风向角=0 or pi/4，U10 = 5/10, X = 30KM)
clear
%设置海面参数
g=9.81;
X0= 2.2e4;

L=500;
N=2*L;
X=30e3;


%1
phi_w=0;
U10=3;
X_= X.*g./(U10^2);
age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
[h, k, S, V, kx, ky,x,y] = generateSeaSurface2D(L, N, U10, age, phi_w);
% plot
% figure(1);
% % mesh(x,y,h),ylabel('y(m)'),xlabel('x(m)');
% % axis([0 L 0 L -2 2])% 设置画图范围，可以去掉
% [X1, Y1] = meshgrid(x, y);
% scatter(X1(:), Y1(:), 10, h(:), 'filled');
% ylabel('y(m)'),xlabel('x(m)');
% colorbar;
% savefig(gcf,"result/5-sea-3-0");
% exportgraphics(gcf,'result/5-sea-3-0.pdf');
% exportgraphics(gcf,'result/5-sea-3-0.png');
% 
% %2
% phi_w=0;
% U10=5;
% X_= X.*g./(U10^2);
% age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
% [h, k, S, V, kx, ky,x,y] = generateSeaSurface2D(L, N, U10, age, phi_w);
% % plot
% figure(2);
% % mesh(x,y,h),ylabel('y(m)'),xlabel('x(m)');
% % axis([0 L 0 L -5 5])% 设置画图范围，可以去掉
% [X1, Y1] = meshgrid(x, y);
% scatter(X1(:), Y1(:), 10, h(:), 'filled');
% ylabel('y(m)'),xlabel('x(m)');
% colorbar;
% savefig(gcf,"result/5-sea-5-0");
% exportgraphics(gcf,'result/5-sea-5-0.pdf');
% exportgraphics(gcf,'result/5-sea-5-0.png');
% 
% %3
% phi_w=0;
% U10=7;
% X_= X.*g./(U10^2);
% age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
% [h, k, S, V, kx, ky,x,y] = generateSeaSurface2D(L, N, U10, age, phi_w);
% % plot
% figure(3);
% % mesh(x,y,h),ylabel('y(m)'),xlabel('x(m)');
% % axis([0 L 0 L -2 2])% 设置画图范围，可以去掉
% [X1, Y1] = meshgrid(x, y);
% scatter(X1(:), Y1(:), 10, h(:), 'filled');
% ylabel('y(m)'),xlabel('x(m)');
% colorbar;
% savefig(gcf,"result/5-sea-7-0");
% exportgraphics(gcf,'result/5-sea-7-0.pdf');
% exportgraphics(gcf,'result/5-sea-7-0.png');

%4
phi_w=0;
U10=9;
X_= X.*g./(U10^2);
age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
[h, k, S, V, kx, ky,x,y] = generateSeaSurface2D(L, N, U10, age, phi_w);
% plot
figure(4);
% mesh(x,y,h),ylabel('y(m)'),xlabel('x(m)');
% axis([0 L 0 L -5 5])% 设置画图范围，可以去掉
[X1, Y1] = meshgrid(x, y);
scatter(X1(:), Y1(:), 10, h(:), 'filled');
ylabel('y(m)'),xlabel('x(m)');
colorbar;
savefig(gcf,"result/5-sea-9-0");
exportgraphics(gcf,'result/5-sea-9-0.pdf');
exportgraphics(gcf,'result/5-sea-9-0.png');

phi_w=0;
U10=15;
X_= X.*g./(U10^2);
age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
[h, k, S, V, kx, ky,x,y] = generateSeaSurface2D(L, N, U10, age, phi_w);
% plot
figure(5);
% mesh(x,y,h),ylabel('y(m)'),xlabel('x(m)');
% axis([0 L 0 L -5 5])% 设置画图范围，可以去掉
[X1, Y1] = meshgrid(x, y);
scatter(X1(:), Y1(:), 10, h(:), 'filled');
ylabel('y(m)'),xlabel('x(m)');
colorbar;
savefig(gcf,"result/5-sea-15-0");
exportgraphics(gcf,'result/5-sea-15-0.pdf');
exportgraphics(gcf,'result/5-sea-15-0.png');
%% 二维海面图风向角【4个图】（海面尺度500m，N=2L，风向角=0 pi/4 90 180，U10 = 5, X = 30KM)
clear
%设置海面参数
g=9.81;
X0= 2.2e4;

L=500;
N=2*L;
X=30e3;


%1
phi_w=0;
U10=5;
X_= X.*g./(U10^2);
age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
[h, k, S, V, kx, ky,x,y] = generateSeaSurface2D(L, N, U10, age, phi_w);
% plot
figure(1);
% mesh(x,y,h),ylabel('y(m)'),xlabel('x(m)');
% axis([0 L 0 L -2 2])% 设置画图范围，可以去掉
[X1, Y1] = meshgrid(x, y);
scatter(X1(:), Y1(:), 10, h(:), 'filled');
ylabel('y(m)'),xlabel('x(m)');
colorbar;
savefig(gcf,"result/5-sea-5-0");
exportgraphics(gcf,'result/5-sea-5-0.pdf');
exportgraphics(gcf,'result/5-sea-5-0.png');

%2
phi_w=pi/4;
U10=5;
X_= X.*g./(U10^2);
age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
[h, k, S, V, kx, ky,x,y] = generateSeaSurface2D(L, N, U10, age, phi_w);
% plot
figure(2);
% mesh(x,y,h),ylabel('y(m)'),xlabel('x(m)');
% axis([0 L 0 L -5 5])% 设置画图范围，可以去掉
[X1, Y1] = meshgrid(x, y);
scatter(X1(:), Y1(:), 10, h(:), 'filled');
ylabel('y(m)'),xlabel('x(m)');
colorbar;
savefig(gcf,"result/5-sea-5-45");
exportgraphics(gcf,'result/5-sea-5-45.pdf');
exportgraphics(gcf,'result/5-sea-5-45.png');

%3
phi_w=pi/2;
U10=5;
X_= X.*g./(U10^2);
age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
[h, k, S, V, kx, ky,x,y] = generateSeaSurface2D(L, N, U10, age, phi_w);
% plot
figure(3);
% mesh(x,y,h),ylabel('y(m)'),xlabel('x(m)');
% axis([0 L 0 L -2 2])% 设置画图范围，可以去掉
[X1, Y1] = meshgrid(x, y);
scatter(X1(:), Y1(:), 10, h(:), 'filled');
ylabel('y(m)'),xlabel('x(m)');
colorbar;
savefig(gcf,"result/5-sea-5-90");
exportgraphics(gcf,'result/5-sea-5-90.pdf');
exportgraphics(gcf,'result/5-sea-5-90.png');

%4
phi_w=pi;
U10=5;
X_= X.*g./(U10^2);
age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
[h, k, S, V, kx, ky,x,y] = generateSeaSurface2D(L, N, U10, age, phi_w);
% plot
figure(4);
% mesh(x,y,h),ylabel('y(m)'),xlabel('x(m)');
% axis([0 L 0 L -5 5])% 设置画图范围，可以去掉
[X1, Y1] = meshgrid(x, y);
scatter(X1(:), Y1(:), 10, h(:), 'filled');
ylabel('y(m)'),xlabel('x(m)');
colorbar;
savefig(gcf,"result/5-sea-5-180");
exportgraphics(gcf,'result/5-sea-5-180.pdf');
exportgraphics(gcf,'result/5-sea-5-180.png');
%% 二维海面恢复功率谱（海面尺度500m，N=2L，风向角=0，U10 = 5, X = 30KM)【结果不好看】
clear
%设置海面参数
g=9.81;
X0= 2.2e4;

L=500;
N=2*L;
X=30e3;
U10=5;
X_= X.*g./(U10^2);
age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);

% y=0时候x的切片
phi_w=0;
[h, k, S, V, kx, ky,x,y] = generateSeaSurface2D(L, N, U10, age, phi_w);
% plot
figure(1);
h_slice = h(:,1);
plot(x,h_slice);
title("Sea Surface xSlice");
xlabel("x(m)");
ylabel("h(m)");

figure(2)
[S,kp] = Elfouhaily(k(:,1),U10,age);
S(1,1)=0.0;
loglog(k(:,1),S);
axis([1e-2 1e1 1e-10 1e1]);
hold on;
F = fft(h_slice');
recover = abs(F).^2 * 2*(L/N)/N;
scatter(k(:,1),recover);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
title("Spectrum comparison");
xlabel("x(m)");
ylabel("S_{x}(m^{4}/rad)");

% x=0时候y的切片
phi_w=0;
[h, k, S, V, kx, ky,x,y] = generateSeaSurface2D(L, N, U10, age, phi_w);
% plot
figure(3);
h_slice = h(1,:);
plot(y,h_slice);
title("Sea Surface ySlice");
xlabel("y(m)");
ylabel("h(m)");

figure(4)
[S,kp] = Elfouhaily(k(1,:),U10,age);
S(1,1)=0.0;
loglog(k(1,:),S);
axis([1e-2 1e1 1e-10 1e1]);
hold on;
F = fft(h_slice');
recover = abs(F).^2 * 2*(L/N)/N;
scatter(k(1,:),recover);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
title("Spectrum comparison");
xlabel("y(m)");
ylabel("S_{y}(m^{4}/rad)");


%% 时变海浪 【w=sqrt(gk)】（海面尺度500m，N=2L，风向角=0 or pi/4，U10 = 5/10, X = 30KM）
clear
%设置海面参数
g=9.81;
X0= 2.2e4;

L=500;
N=2*L;
X=30e3;
phi_w=0;
U10=5;
X_= X.*g./(U10^2);
age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
% w = sqrt(g*2*pi/L);
dt = 1;

[h, k, S, V, kx, ky,x,y] = generateSeaSurface2D_time_new(L, N, U10, age, phi_w, dt);




function [h, he, he_t, k, S, Se, V, kx, ky,x,y] = generateSeaSurface2D_time_boat_paper(L, N, U10, age, phi_w, dt, time, vs, beta)
%% parameter
% clear;
% L=4;
% N=256; % need to be even
% % U10=2,age=1会出bug
% U10=5;
% age=0.84;
% phi_w=0;
% vs = 5;
% beta = pi/4;
% time = 50;
%% Frequency Mesh
dk = 2*pi/L; %frequency grid spacing (1/m)
%build up the matrix of wave numbers
kx = (-N/2: N/2 - 1) * dk; 
ky = kx';
[kxx,kyy] = meshgrid(kx,ky);
%shift so that element (1,1) is k = 0
kxx = ifftshift(kxx);
kyy = ifftshift(kyy);

%convert to polar coordinates
[phi,k] = cart2pol(kxx,kyy);

g = 9.81;
w = sqrt(g.*k);
we = abs(w - k.* vs .* cos(phi-beta));%论文7-3
he_t = zeros(1,time+1);
counter = 1;
%% Spectral Representation
%get the elfouhaily spectrum
S = Elfouhaily2D(k,phi,U10,age,phi_w);  
S(1,1) = 0.0;

%we now have a 2-way spectrum (positive and negative frequencies) and need to scale the total power by 2
S = 2*S;

%now we get encounter frequency
Se = S./abs(1-2.*w./g.*vs.*cos(phi-beta));

%build up the random representation in frequency space 构建高斯变量
V = zeros(N);
Ve = zeros(N);

%need to generate two sub matrices Va and Vb
Va = zeros(N/2-1);
Vb = zeros(N/2-1);
Vae = zeros(N/2-1);
Vbe = zeros(N/2-1);

%matrices of random variables for Va and Vb (0~1 random num，using gauss)
W1 = randn(N/2-1);
W2 = randn(N/2-1);
W3 = randn(N/2-1);
W4 = randn(N/2-1);

%now need to handle the cases: kx = 0, and ky = 0 起始点
w1 = randn(1,N/2); % 1 * N/2的随机数
u1 = randn(1,N/2);
w2 = randn(1,N/2);
u2 = randn(1,N/2);

Vx0 = [];
Vy0 = [];
Vx0e = [];
Vy0e = [];

%build the zero frequency lines
Vx0(1) = 0;
Vy0(1) = 0;
Vx0e(1) = 0;
Vy0e(1) = 0;


%now need to handle the cases: kx = N/2, and ky = N/2
w3 = randn(1,N/2);
u3 = randn(1,N/2);
w4 = randn(1,N/2);
u4 = randn(1,N/2);


%引入时间
% figure;
%初始化视频
% writerObj = VideoWriter('wave_new_en.avi'); % 指定视频文件名
% writerObj.FrameRate = 5; % 设置帧率
% open(writerObj);

for t=0:dt:time
%loop and populate Va and Vb
for u = 1:N/2-1
    for v = 1:N/2-1
        ua = u+1;
        va = v+1;
        ub = u + 1;
        vb = v + N/2 + 1;
        Va(u,v) = pi*sqrt(S(ua,va)*L^2)*(W1(u,v) + 1i*W2(u,v))*exp(1i*w(ua,va)*t);
        Vb(u,v) = pi*sqrt(S(ub,vb)*L^2)*(W3(u,v) + 1i*W4(u,v))*exp(1i*w(ub,vb)*t);
        Vae(u,v) = pi*sqrt(Se(ua,va)*L^2)*(W1(u,v) + 1i*W2(u,v))*exp(1i*we(ua,va)*t);
        Vbe(u,v) = pi*sqrt(Se(ub,vb)*L^2)*(W3(u,v) + 1i*W4(u,v))*exp(1i*we(ub,vb)*t);
    end
end

%now place the submatrices and their Hermitian conjugates in place
V(2:N/2,2:N/2) = Va; % （1，1）=0
V(N/2+2:N,N/2+2:N) = conj(flipud(fliplr(Va))); % 对角线翻转，取共轭=厄密共轭
V(2:N/2,N/2+2:N) = Vb;
V(N/2+2:N,2:N/2) = conj(flipud(fliplr(Vb)));

Ve(2:N/2,2:N/2) = Vae; % （1，1）=0
Ve(N/2+2:N,N/2+2:N) = conj(flipud(fliplr(Vae))); % 对角线翻转，取共轭=厄密共轭
Ve(2:N/2,N/2+2:N) = Vbe;
Ve(N/2+2:N,2:N/2) = conj(flipud(fliplr(Vbe)));
%now is [0000000
%        0Va Vb0
%        0Vb~Va~0]
% operation after is to fill the 0 blank


for j = 2:N/2
    Vx0(j) = pi*sqrt(S(1,j)*L^2)*(w1(j) + 1i*u1(j))*exp(1i*w(1,j)*t);
    Vy0(j) = pi*sqrt(S(j,1)*L^2)*(w2(j) + 1i*u2(j))*exp(1i*w(j,1)*t);
    
    Vx0e(j) = pi*sqrt(Se(1,j)*L^2)*(w1(j) + 1i*u1(j))*exp(1i*we(1,j)*t);
    Vy0e(j) = pi*sqrt(Se(j,1)*L^2)*(w2(j) + 1i*u2(j))*exp(1i*we(j,1)*t);
end

Vx0(N/2+1) = 2*pi*sqrt(S(1,N/2+1)*L^2)*u1(1)*exp(1i*w(1,N/2+1)*t);
Vy0(N/2+1) = 2*pi*sqrt(S(N/2+1,1)*L^2)*u2(1)*exp(1i*w(N/2+1,1)*t);

Vx0e(N/2+1) = 2*pi*sqrt(Se(1,N/2+1)*L^2)*u1(1)*exp(1i*we(1,N/2+1)*t);
Vy0e(N/2+1) = 2*pi*sqrt(Se(N/2+1,1)*L^2)*u2(1)*exp(1i*we(N/2+1,1)*t);

for j = N/2+2:N
Vx0(j) = conj(Vx0(N-j + 2));
Vy0(j) = conj(Vy0(N-j + 2));

Vx0e(j) = conj(Vx0e(N-j + 2));
Vy0e(j) = conj(Vy0e(N-j + 2));
end

%place the zero frequency lines in V - rows are y, columns are x (need to
%transpose to a column vector)
V(1,:) = Vy0 ;
V(:,1) = Vx0';
Ve(1,:) = Vy0e ;
Ve(:,1) = Vx0e';


%build the N/2 frequency lines
Vx2(1) = 2*pi*sqrt(S(N/2+1,1)*L^2)*w3(1)*exp(1i*w(N/2+1,1)*t);
Vy2(1) = 2*pi*sqrt(S(1,N/2+1)*L^2)*w4(1)*exp(1i*w(1,N/2+1)*t);
Vx2e(1) = 2*pi*sqrt(Se(N/2+1,1)*L^2)*w3(1)*exp(1i*we(N/2+1,1)*t);
Vy2e(1) = 2*pi*sqrt(Se(1,N/2+1)*L^2)*w4(1)*exp(1i*we(1,N/2+1)*t);
for j = 2:N/2
    Vx2(j) = pi*sqrt(S(N/2+1,j)*L^2)*(w3(j) + 1i*u3(j))*exp(1i*w(N/2+1,j)*t);
    Vy2(j) = pi*sqrt(S(j,N/2+1)*L^2)*(w4(j) + 1i*u4(j))*exp(1i*w(j,N/2+1)*t);
    Vx2e(j) = pi*sqrt(Se(N/2+1,j)*L^2)*(w3(j) + 1i*u3(j))*exp(1i*we(N/2+1,j)*t);
    Vy2e(j) = pi*sqrt(Se(j,N/2+1)*L^2)*(w4(j) + 1i*u4(j))*exp(1i*we(j,N/2+1)*t);
end

Vx2(N/2+1) = 2*pi*sqrt(S(N/2+1,N/2+1)*L^2)*u4(1)*exp(1i*w(N/2+1,N/2+1)*t);
Vy2(N/2+1) = 2*pi*sqrt(S(N/2+1,N/2+1)*L^2)*u4(1)*exp(1i*w(N/2+1,N/2+1)*t);

Vx2e(N/2+1) = 2*pi*sqrt(Se(N/2+1,N/2+1)*L^2)*u4(1)*exp(1i*we(N/2+1,N/2+1)*t);
Vy2e(N/2+1) = 2*pi*sqrt(Se(N/2+1,N/2+1)*L^2)*u4(1)*exp(1i*we(N/2+1,N/2+1)*t);

for j = N/2+2:N
Vx2(j) = conj(Vx2(N-j + 2));
Vy2(j) = conj(Vy2(N-j + 2));

Vx2e(j) = conj(Vx2e(N-j + 2));
Vy2e(j) = conj(Vy2e(N-j + 2));
end

%now place the N/2 frequency lines - rows are y, columns are x (need to
%transpose to a column vector)
V(N/2+1,:) = Vy2;
V(:,N/2+1) = Vx2';

Ve(N/2+1,:) = Vy2e;
Ve(:,N/2+1) = Vx2e';

% compute
h = real(ifft2(V)*length(V)^2)/L^2;
x = (0:N-1)*L/N;
y = (0:N-1)*L/N;
he = real(ifft2(Ve)*length(Ve)^2)/L^2;

% 做y（N/2-1）的平面，求交线
intersLine = he(N/2-1,:);
% plot(x,intersLine);
% hold on;
% deltah = (intersLine(N/2)-intersLine(N/2-2))
% deltax = (x(N/2)-x(N/2-2))
slope = (intersLine(N/2)-intersLine(N/2-2))/(x(N/2)-x(N/2-2));
he_t(counter)= atan(slope);

% he_t(counter) = he(N/2-1, N/2-1);
counter = counter + 1;
%% plot
% clf;
%whole plot
% [X, Y] = meshgrid(x, y);
% clf
% scatter(X(:), Y(:), 10, h(:), 'filled');
% colorbar;

% [X, Y] = meshgrid(x, y);
% clf
% scatter(X(:), Y(:), 10, he(:), 'filled');
% colorbar;

% %xslice plot
% h_slice = h(:,1);
% plot(x,h_slice);
% title("Sea Surface xSlice");
% xlabel("x(m)");
% ylabel("h(m)");
% axis([0 L -0.6 0.6]);

%yslice plot
% h_slice = h(1,:);
% plot(y,h_slice);
% title("Sea Surface ySlice");
% xlabel("y(m)");
% ylabel("h(m)");

% frame = getframe(gcf);
% writeVideo(writerObj, frame);
end
% close(writerObj);
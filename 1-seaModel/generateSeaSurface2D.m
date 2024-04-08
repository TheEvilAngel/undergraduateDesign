function [h, k, S, V, kx, ky,x,y] = generateSeaSurface2D(L, N, U10, age, phi_w)
%% parameter
% L=4;
% N=256; % need to be even
% % U10=2,age=1会出bug
% U10=5;
% age=0.84;
% phi_w=0;
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

%% Spectral Representation
%get the elfouhaily spectrum
S = Elfouhaily2D(k,phi,U10,age,phi_w);  
S(1,1) = 0.0;

%we now have a 2-way spectrum (positive and negative frequencies) and need to scale the total power by 2
S = 2*S;

%build up the random representation in frequency space 构建高斯变量
V = zeros(N);

%need to generate two sub matrices Va and Vb
Va = zeros(N/2-1);
Vb = zeros(N/2-1);

%matrices of random variables for Va and Vb (0~1 random num，using gauss)
W1 = randn(N/2-1);
W2 = randn(N/2-1);
W3 = randn(N/2-1);
W4 = randn(N/2-1);

%loop and populate Va and Vb
for u = 1:N/2-1
    for v = 1:N/2-1
        ua = u+1;
        va = v+1;
        ub = u + 1;
        vb = v + N/2 + 1;
        Va(u,v) = 1/2*sqrt(S(ua,va)*dk^2)*(W1(u,v) + 1i*W2(u,v));
        Vb(u,v) = 1/2*sqrt(S(ub,vb)*dk^2)*(W3(u,v) + 1i*W4(u,v));
    end
end

%now place the submatrices and their Hermitian conjugates in place
V(2:N/2,2:N/2) = Va; % （1，1）=0
V(N/2+2:N,N/2+2:N) = conj(flipud(fliplr(Va))); % 对角线翻转，取共轭=厄密共轭
V(2:N/2,N/2+2:N) = Vb;
V(N/2+2:N,2:N/2) = conj(flipud(fliplr(Vb)));

%now is [0000000
%        0Va Vb0
%        0Vb~Va~0]
% operation after is to fill the 0 blank

%now need to handle the cases: kx = 0, and ky = 0 起始点
w1 = randn(1,N/2); % 1 * N/2的随机数
u1 = randn(1,N/2);
w2 = randn(1,N/2);
u2 = randn(1,N/2);

Vx0 = [];
Vy0 = [];

%build the zero frequency lines
Vx0(1) = 0;
Vy0(1) = 0;
for j = 2:N/2
    Vx0(j) = 1/2*sqrt(S(1,j)*dk^2)*(w1(j) + 1i*u1(j));
    Vy0(j) = 1/2*sqrt(S(j,1)*dk^2)*(w2(j) + 1i*u2(j));
end
Vx0(N/2+1) = sqrt(S(1,N/2+1)*dk^2)*u1(1);
Vy0(N/2+1) = sqrt(S(N/2+1,1)*dk^2)*u2(1);

for j = N/2+2:N
Vx0(j) = conj(Vx0(N-j + 2));
Vy0(j) = conj(Vy0(N-j + 2));
end

%place the zero frequency lines in V - rows are y, columns are x (need to
%transpose to a column vector)
V(1,:) = Vy0 ;
V(:,1) = Vx0';

%now need to handle the cases: kx = N/2, and ky = N/2
w3 = randn(1,N/2);
u3 = randn(1,N/2);
w4 = randn(1,N/2);
u4 = randn(1,N/2);

%build the N/2 frequency lines
Vx2(1) = sqrt(S(N/2+1,1)*dk^2)*w3(1);
Vy2(1) = sqrt(S(1,N/2+1)*dk^2)*w4(1);

for j = 2:N/2
    Vx2(j) = 1/2*sqrt(S(N/2+1,j)*dk^2)*(w3(j) + 1i*u3(j));
    Vy2(j) = 1/2*sqrt(S(j,N/2+1)*dk^2)*(w4(j) + 1i*u4(j));
end

Vx2(N/2+1) = sqrt(S(N/2+1,N/2+1)*dk^2)*u4(1);
Vy2(N/2+1) = sqrt(S(N/2+1,N/2+1)*dk^2)*u4(1);

for j = N/2+2:N
Vx2(j) = conj(Vx2(N-j + 2));
Vy2(j) = conj(Vy2(N-j + 2));
end

%now place the N/2 frequency lines - rows are y, columns are x (need to
%transpose to a column vector)
V(N/2+1,:) = Vy2;
V(:,N/2+1) = Vx2';

% generate the surface and create the x and y vectors for return
h = real(ifft2(V)*length(V)^2);
x = (0:N-1)*L/N;
y = (0:N-1)*L/N;

% % plot
% mesh(x,y,h),ylabel('y(m)'),xlabel('x(m)');
% axis([0 L 0 L -1 1])% 设置画图范围，可以去掉
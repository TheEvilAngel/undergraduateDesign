%% 波浪波倾角切线版本
clear
clc
%设置海面参数
g=9.81;
X0= 2.2e4;

%Sea parameter
L=100;
N=2*L;%精度得满足采样率30足够
X=30e3;
phi_w=0;
U10=5;
X_= X.*g./(U10^2);
age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
dt = 0.05;
%boat parameter
vs = 0;
beta = pi/2;
time = 1; % total 300s
[h, he, he_t, k, S, Se, V, kx, ky,x,y] = generateSeaSurface2D_time_boat(L, N, U10, age, phi_w, dt,time, vs, beta);

t = 0:dt:time;
boat_depth = 1;%船吃水深度
Cw = 0.8072;%水线面系数
B = 4.8;%船宽度
K1 = exp(-1*2*pi/L*boat_depth/2);
K2 = 1-sqrt(Cw*(B/L)^2);
he_t_deg = rad2deg(K1*K2*he_t);
if vs == 0
    figure;
    plot(t,he_t_deg);
    xlabel('t(s)');ylabel('\alpha(°)');
    savefig(gcf,"result/7-waveAngle");
    exportgraphics(gcf,'result/7-waveAngle.pdf');
    exportgraphics(gcf,'result/7-waveAngle.png');
else 
    figure;
    plot(t,he_t_deg);
    xlabel('t(s)');ylabel('\alpha(°)');
    savefig(gcf,"result/7-waveAngle1");
    exportgraphics(gcf,'result/7-waveAngle1.pdf');
    exportgraphics(gcf,'result/7-waveAngle1.png');
end

%% 横摇数据模型（画幅频特性曲线）
D = 135.3;%排水量
h_boat = 0.73;%横稳心高
L_ship = 38.53;%船长
d = 2.48;%船高
B = 4.8;%船宽度
% V = 1;%方形系数面积
Cb = 0.65;%方形系数（文件给出）
% Aw = 1;%水线面面积
Cw = 0.8072;%水线面系数
H = 1.92;%型深
Ix = D/9.81*(B^2*Cw^2/11.4/Cb+H^2/12);%论文赵晔7-14
delta_Ix = Ix * (-0.186+1.179*Cb-0.615*Cb^2);%课本164页
w_phi = sqrt(D*h_boat/(Ix+delta_Ix));%课本5-84
% % miu = 1;%无因次衰减系数（matlab中n_phi,经验公式）
% miu = 0.06*L_ship*B^4/(D*(B^2+H^2))*0.6;%尼古拉耶夫公式
kesi_phi = 0.687;


A = [0, 1; -(w_phi^2), -2*kesi_phi*w_phi];
B = [0, w_phi^2]';
C = [1, 0];
D = 0;
G_phi = ss(A, B, C, D);
tt = 0:dt:time;
u = he_t_deg;
if vs==0
[ship_phi,t]=lsim(G_phi, u, tt);
figure;
plot(t,ship_phi);
hold on
plot(t,u,'Color',[0.5 0.5 0.5]);
xlabel('t(s)');ylabel('deg(°)');
legend("波倾角","横摇角");
savefig(gcf,"result/8-shakeAngle");
exportgraphics(gcf,'result/8-shakeAngle.pdf');
exportgraphics(gcf,'result/8-shakeAngle.png');
save('result/vs0.mat');
else
[ship_phi,t]=lsim(G_phi, u, tt);
figure;
plot(t,ship_phi);
hold on
plot(t,u,'Color',[0.5 0.5 0.5]);
xlabel('t(s)');ylabel('deg(°)');
legend("波倾角","横摇角");
savefig(gcf,"result/8-shakeAngle1");
exportgraphics(gcf,'result/8-shakeAngle1.pdf');
exportgraphics(gcf,'result/8-shakeAngle1.png');
save('result/vs1.mat')
end

%% 只导出海面模型
[X, Y] = meshgrid(x, y);

%另一种更好的方法
figure
sea = surf(X,Y,h);
surf2stl(sea, '1_sea_origin');


%% 修改船只模型（自己写的）【把库文件添加到路径】
[Vscaled, VFaces, ~, ~] = stlRead('serie_60_una_pieza.stl');
% 船只缩放倍数
Vscaled = Vscaled*0.3;

% % 船只修改倾角（空间几何变换）
% phi_rad = deg2rad(ship_phi(time+1));
rot_matrix = roty(0);
% rot_matrix = roty(30);% 大角度测试
Vscaled = rot_matrix * Vscaled';
Vscaled = Vscaled';

% 船只移动到海面中心
x_offset = mean(Vscaled(:,1)) - mean(x);
Vscaled(:,1) = Vscaled(:,1) - x_offset;
y_offset = mean(Vscaled(:,2)) - mean(y);
Vscaled(:,2) = Vscaled(:,2) - y_offset;
z_offset = mean(Vscaled(:,3)) - (d/2) + boat_depth;%需要知道船只吃水深度代替z_offset，把海浪平均高度作为吃水深度
Vscaled(:,3) = Vscaled(:,3) - z_offset;


% 删除船只内海面的点
[X, Y] = meshgrid(x, y);%建立海面网格
figure;
Vscaled_circle = convhull(Vscaled(:,1),Vscaled(:,2));
plot(Vscaled(Vscaled_circle,1),Vscaled(Vscaled_circle,2));
[in, on] = inpolygon(X, Y, Vscaled(Vscaled_circle,1),Vscaled(Vscaled_circle,2));%找到海面在船里面的点
% 删除有交线的点
% 遍历网格中的每个线
for i = 1:size(x,2)
    if ~ismember(1, in(:,i))%如果这列没有点在船内（固定x，线垂直x轴）
        xi = [x(i) x(i)];
        yi = [y(1) y(size(y,2))];
        [xi,yi] = polyxpoly(xi,yi,Vscaled(Vscaled_circle,1),Vscaled(Vscaled_circle,2));
        if ~isempty(yi)%有交点
            scale = L/N;
            nearest_scale = round(yi / scale) * scale;  
            idy = [0 0];
            eps = 1e-5;
            idy(1) = find(abs(y - nearest_scale(1)) < eps);
            idy(2) = find(abs(y - nearest_scale(2)) < eps);
            in(idy(1),i) = 1;%in 中行为y坐标，列为x坐标
            in(idy(2),i) = 1;
        end
    end
end
for j = 1:size(y,2)
    if ~ismember(1, in(j,:))%如果这行没有点在船内
        xj = [x(1) x(size(x,2))];
        yj = [y(j) y(j)];
        [xj,yj] = polyxpoly(xj,yj,Vscaled(Vscaled_circle,1),Vscaled(Vscaled_circle,2));
        if ~isempty(xj)%有交点
            scale = L/N;
            nearest_scale = round(xj / scale) * scale; 
            idx = [0 0];
            eps = 1e-5;
            idx(1) = find(abs(x - nearest_scale(1)) < eps);
            idx(2) = find(abs(x - nearest_scale(2)) < eps);
            in(j,idx(1)) = 1; %in 中行为y坐标，列为x坐标
            in(j,idx(2)) = 1;
        end
    end
end
h_filt = h;
h_filt(in) = NaN;%海面在船内的点删除
h_filt(on) = NaN;

hold on;
scatter(X(in),Y(in),'r+');

% 输出生成的海面和船只
figure;
sea_filt = surf(X,Y,h_filt);
surf2stl(sea_filt, '1_sea_filt');
stlwrite('2_boat_filt.stl', VFaces, Vscaled);


% 重新读入并合并
[ship_new_verticles, ship_new_faces, ~, ~] = stlRead('2_boat_filt.stl');
[sea_new_verticles, sea_new_faces, ~, ~] = stlRead('1_sea_filt.stl');
F_whole = [ship_new_faces; sea_new_faces+length(ship_new_verticles)];
V_whole = [ship_new_verticles; sea_new_verticles];
stlwrite('3_whole.stl', F_whole,V_whole);


% % 画图
% % figure;
% % scatter3(x_new,y_new,h_new),ylabel('y(m)'),xlabel('x(m)');
% % [X, Y] = meshgrid(x_new, y_new);
% % H = reshape(h_new, size(X));
% % surf(X, Y, H);
% % xlabel('x');
% % ylabel('y');
% % zlabel('h');



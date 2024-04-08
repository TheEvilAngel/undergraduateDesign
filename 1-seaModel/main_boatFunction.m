% %% 海浪波倾角（海面尺度50m，N=2L，风向角=0，U10 = 5, X = 30KM，vs=5，航向角beta=30°）
% clear
% %设置海面参数
% g=9.81;
% X0= 2.2e4;
% 
% %Sea parameter
% L=50;
% N=2*L;%精度得满足采样率30足够
% X=30e3;
% phi_w=0;
% U10=5;
% X_= X.*g./(U10^2);
% age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
% dt = 0.01;
% 
% %boat parameter
% vs = 5;
% beta = 0;
% time = 50; % total 300s
% [h, he, he_t, k, S, Se, V, kx, ky,x,y] = generateSeaSurface2D_time_boat(L, N, U10, age, phi_w, dt,time, vs, beta);
% 
% t = 0:dt:time;
% K1 = 0.8564;
% K2 = 0.7281;
% % he_t = he_t .* 2*pi/L .* sin(beta);
% he_t = he_t .* 2*pi/L;
% % he_t = he_t .* k(N/2-1,N/2-1);
% he_t_deg = rad2deg(he_t);
% figure;
% plot(t,he_t_deg);
% xlabel('t(s)');ylabel('\alpha(°)');
% % savefig(gcf,"result/7-waveAngle");
% % exportgraphics(gcf,'result/7-waveAngle.pdf');
% % exportgraphics(gcf,'result/7-waveAngle.png');
% figure;
% mesh(x,y,h),ylabel('y(m)'),xlabel('x(m)');
% % savefig(gcf,"result/7-waveSea");
% % exportgraphics(gcf,'result/7-waveSea.pdf');
% % exportgraphics(gcf,'result/7-waveSea.png');
% figure;
% mesh(x,y,he),ylabel('y(m)'),xlabel('x(m)');
% % savefig(gcf,"result/7-waveSeaEncounter");
% % exportgraphics(gcf,'result/7-waveSeaEncounter.pdf');
% % exportgraphics(gcf,'result/7-waveSeaEncounter.png');

%% 波浪波倾角切线版本
clear
%设置海面参数
g=9.81;
X0= 2.2e4;

%Sea parameter
L=50;
N=2*L;%精度得满足采样率30足够
X=30e3;
phi_w=0;
U10=5;
X_= X.*g./(U10^2);
age = 0.84*power(tanh(power(X_./X0,0.4)),0.75);
dt = 0.01;

%boat parameter
vs = 5;
beta = 0;
time = 5; % total 300s
[h, he, he_t, k, S, Se, V, kx, ky,x,y] = generateSeaSurface2D_time_boat(L, N, U10, age, phi_w, dt,time, vs, beta);

t = 0:dt:time;
K1 = 0.8564;
K2 = 0.7281;
% he_t is the height of the point with time
% he_t = he_t .* 2*pi/L .* sin(beta);


% he_t = he_t .* 2*pi/L;
% he_t = he_t .* k(N/2-1,N/2-1);
he_t_deg = rad2deg(he_t);
figure;
plot(t,he_t_deg);
xlabel('t(s)');ylabel('\alpha(°)');
% savefig(gcf,"result/7-waveAngle");
% exportgraphics(gcf,'result/7-waveAngle.pdf');
% exportgraphics(gcf,'result/7-waveAngle.png');
figure;
mesh(x,y,h),ylabel('y(m)'),xlabel('x(m)');
% savefig(gcf,"result/7-waveSea");
% exportgraphics(gcf,'result/7-waveSea.pdf');
% exportgraphics(gcf,'result/7-waveSea.png');
figure;
mesh(x,y,he),ylabel('y(m)'),xlabel('x(m)');
% savefig(gcf,"result/7-waveSeaEncounter");
% exportgraphics(gcf,'result/7-waveSeaEncounter.pdf');
% exportgraphics(gcf,'result/7-waveSeaEncounter.png');
%% 横摇数据模型（画幅频特性曲线）[w_phi = 0.3]
kesi_phi = 0.1;%无因次衰减系数，经验公式求解
w_phi = 0.3;%用书上公式5-84，Ix由赵烨7-47，deltaIx由书中公式求解
A = [0, 1; -(w_phi^2), -2*kesi_phi*w_phi];
B = [0, w_phi^2]';
C = [1, 0];
D = 0;
G_phi1 = ss(A, B, C, D);

kesi_phi = 0.15;%无因次衰减系数，经验公式求解
w_phi = 0.3;%用书上公式5-84，Ix由赵烨7-47，deltaIx由书中公式求解
A = [0, 1; -(w_phi^2), -2*kesi_phi*w_phi];
B = [0, w_phi^2]';
C = [1, 0];
D = 0;
G_phi2 = ss(A, B, C, D);

kesi_phi = 0.3;%无因次衰减系数，经验公式求解
w_phi = 0.3;%用书上公式5-84，Ix由赵烨7-47，deltaIx由书中公式求解
A = [0, 1; -(w_phi^2), -2*kesi_phi*w_phi];
B = [0, w_phi^2]';
C = [1, 0];
D = 0;
G_phi3 = ss(A, B, C, D);

kesi_phi = 0.6;%无因次衰减系数，经验公式求解
w_phi = 0.3;%用书上公式5-84，Ix由赵烨7-47，deltaIx由书中公式求解
A = [0, 1; -(w_phi^2), -2*kesi_phi*w_phi];
B = [0, w_phi^2]';
C = [1, 0];
D = 0;
G_phi4 = ss(A, B, C, D);

kesi_phi = 1;%无因次衰减系数，经验公式求解
w_phi = 0.3;%用书上公式5-84，Ix由赵烨7-47，deltaIx由书中公式求解
A = [0, 1; -(w_phi^2), -2*kesi_phi*w_phi];
B = [0, w_phi^2]';
C = [1, 0];
D = 0;
G_phi5 = ss(A, B, C, D);

kesi_phi = 2;%无因次衰减系数，经验公式求解
w_phi = 0.3;%用书上公式5-84，Ix由赵烨7-47，deltaIx由书中公式求解
A = [0, 1; -(w_phi^2), -2*kesi_phi*w_phi];
B = [0, w_phi^2]';
C = [1, 0];
D = 0;
G_phi6 = ss(A, B, C, D);

bode(G_phi1,G_phi2,G_phi3,G_phi4,G_phi5,G_phi6);

legend('$\xi=0.1$','$\xi=0.15$','$\xi=0.3$','$\xi=0.6$','$\xi=1.0$','$\xi=2.0$','Interpreter','latex');
xlim([0.01 10]);
grid on;

%% 横摇数据模型（画幅频特性曲线）
% 计算相关参数
% D = 1;%排水量
% h = 1;%横稳心高
% B = 1;%船宽
% L_ship = 1;%船长
% d = 1;%船高
% V = 1;%方形系数面积
% Cb = V/(L_ship*B*d);%方形系数（文件给出）
% Aw = 1;%水线面面积
% Cw = Aw/(L_ship*B);%水线面系数
% H = 1;%型深
% Ix = D/9.81*(B^2*Cw^2/11.4/Cb+H^2/12);
% delta_Ix = Ix * (-0.186+1.179*Cb-0.615*Cb^2);
% w_phi = sqrt(D*h/(Ix+delta_Ix));
% % miu = 1;%无因次衰减系数（matlab中n_phi,经验公式）
% miu = 0.058*L_ship*B^4/(D*(B^2+H^2))*0.55;
% kesi_phi = miu/w_phi;


kesi_phi = 0.1;%无因次衰减系数，经验公式求解
w_phi = 0.3;%用书上公式5-84，Ix由赵烨7-47，deltaIx由书中公式求解
A = [0, 1; -(w_phi^2), -2*kesi_phi*w_phi];
B = [0, w_phi^2]';
C = [1, 0];
D = 0;
G_phi = ss(A, B, C, D);

% dt = 1;
% time = 400;
tt = 0:dt:time;
u = he_t_deg;
[ship_phi,t]=lsim(G_phi, u, tt);
plot(t,ship_phi,'-');


%% 只导出海面模型
[X, Y] = meshgrid(x, y);

%可行的一种方法
% k_data = convhulln(data);
% trisurf2stl('surf.stl',k_data,data);

%另一种更好的方法
sea = surf(X,Y,h);
surf2stl(sea, '1_sea_origin');

% x_new = data(:,1);
% y_new = data(:,2);
% h_new = data(:,3);
% DT = delaunayTriangulation(x_new,y_new,h_new);
% tetramesh(DT,'FaceAlpha',0.3);
% DT_tri = triangulation(DT.ConnectivityList,DT.Points);
% stlwrite(DT_tri,'tritext.stl','binary');
% patch(x_new,y_new,h_new, 'EdgeColor', 'k');
% axis([0 L 0 L -1.5 1.5]);
% axis equal;


%% 修改船只模型（自己写的）【把库文件别添加到路径】
ship = stlread('serie_60_una_pieza.stl');
Vscaled = ship.Points; %先重新新建数组（原只为可读文件）

% 船只缩放倍数
Vscaled = Vscaled*0.1;

% % 船只修改倾角（空间几何变换）
% phi_rad = deg2rad(ship_phi(time+1));
rot_matrix = roty(ship_phi(time+1));
% rot_matrix = roty(30);% 大角度测试
Vscaled = rot_matrix * Vscaled';
Vscaled = Vscaled';

% 船只移动到海面中心
x_offset = mean(Vscaled(:,1)) - mean(x);
Vscaled(:,1) = Vscaled(:,1) - x_offset;
y_offset = mean(Vscaled(:,2)) - mean(y);
Vscaled(:,2) = Vscaled(:,2) - y_offset;
sea_depth = mean(h(:));
z_offset = mean(Vscaled(:,3)) - sea_depth;%需要知道船只吃水深度代替z_offset，把海浪平均高度作为吃水深度
Vscaled(:,3) = Vscaled(:,3) - z_offset;


% 删除船只内海面的点
[X, Y] = meshgrid(x, y);%建立海面网格
% indice = find(Vscaled(:,3) < sea_depth);%找到吃水线下面坐标索引
% ship_polygon = [Vscaled(indice,1), Vscaled(indice,2)];%绘制多边形(船只范围）【要和水高度一样】
% ship_polygon = [Vscaled(:,1), Vscaled(:,1)];%绘制多边形(船只范围）【最长的范围】
figure;
% plot(Vscaled(:,1),Vscaled(:,2));
Vscaled_circle = convhull(Vscaled(:,1),Vscaled(:,2));
% plot(Vscaled(Vscaled_circle,1),Vscaled(Vscaled_circle,2));
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
            idy(1) = find(y == nearest_scale(1));
            idy(2) = find(y == nearest_scale(2));
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
            idx(1) = find(x == nearest_scale(1));
            idx(2) = find(x == nearest_scale(2));
            in(j,idx(1)) = 1; %in 中行为y坐标，列为x坐标
            in(j,idx(2)) = 1;
        end
    end
end
h_filt = h;
h_filt(in) = NaN;%海面在船内的点删除
h_filt(on) = NaN;
% hold on;
% scatter(X(in),Y(in),'r+');

% 输出生成的海面和船只
sea_filt = surf(X,Y,h_filt);
surf2stl(sea_filt, '1_sea_filt');
stlwrite('2_boat_filt.stl', ship.ConnectivityList, Vscaled);


% 重新读入并合并
ship_new = stlread('2_boat_filt.stl');
sea_new = stlread('1_sea_filt.stl');
F_whole = [ship_new.ConnectivityList; sea_new.ConnectivityList+length(ship_new.Points)];
V_whole = [ship_new.Points; sea_new.Points];
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


%% 修改船只模型gptool版本【失败】
[V_shipOrigin,F_shipOrigin,N_shipOrigin] = readSTL('serie_60_una_pieza.stl');
[V_seaOrigin,F_seaOrigin,N_seaOrigin] = readSTL('1_sea_origin.stl');

V_shipFilt = V_shipOrigin*0.1;

% % 船只修改倾角（空间几何变换）
% phi_rad = deg2rad(ship_phi(time+1));
rot_matrix = roty(ship_phi(time+1));
% rot_matrix = roty(30);% 大角度测试
V_shipFilt = rot_matrix * V_shipFilt';
V_shipFilt = V_shipFilt';

% 船只移动到海面中心
x_offset = mean(V_shipFilt(:,1)) - mean(x);
V_shipFilt(:,1) = V_shipFilt(:,1) - x_offset;
y_offset = mean(V_shipFilt(:,2)) - mean(y);
V_shipFilt(:,2) = V_shipFilt(:,2) - y_offset;
sea_depth = mean(h(:));
z_offset = mean(V_shipFilt(:,3)) - sea_depth;%需要知道船只吃水深度代替z_offset，把海浪平均高度作为吃水深度
V_shipFilt(:,3) = V_shipFilt(:,3) - z_offset;

% % Inputs:
% %   V  #V by 3 list of vertex positions of first mesh
% %   F  #F by 3 list of triangle indices into V
% %   U  #U by 3 list of vertex positions of second mesh
% %   G  #G by 3 list of triangle indices into U
% %   operation  followed by operation to perform as a string, one of: 'union',
% %     'intersect', 'minus', 'xor', or 'resolve'
% %     Optional:
% %       'BooleanLib' followed by boolean library back-end to use, one of:
% %         {'libigl'}  uses CGAL's exact arithmetic kernel and is believed to be
% %                     correct.
% %         'cork'  is faster but may give incorrect results. 
% %         'libigl-try-cork-resolve'  libigl boolean extraction but tries to use
% %                                    cork's fast resolve, if intersections
% %                                    persist, then resolves remaining with
% %                                    libigl's resolve. This adds a "layer of
% %                                    robustness" on top of cork, but since it's
% %                                    not understood _how_ cork is failing, it
% %                                    is unknown whether this will lead to
% %                                    correct results.
% % Outputs:
% %   W  #W by 3 list of vertex positions of boolean result mesh
% %   H  #H by 3 list of triangle indices into W
% %   J  #H list of indices into [FA;FB] of facet birth parents
[V_seaFilt,F_seaFilt,N_seaFilt] = mesh_boolean_winding_number(V_seaOrigin,F_seaOrigin,V_shipFilt,F_shipOrigin,'minus');

%% 修改船只模型iso2mesh版本
[V_shipOrigin,F_shipOrigin,N_shipOrigin] = readSTL('serie_60_una_pieza.stl');
[V_seaOrigin,F_seaOrigin,N_seaOrigin] = readSTL('1_sea_origin.stl');

V_shipFilt = V_shipOrigin*0.1;

% % 船只修改倾角（空间几何变换）
% phi_rad = deg2rad(ship_phi(time+1));
rot_matrix = roty(ship_phi(time+1));
% rot_matrix = roty(30);% 大角度测试
V_shipFilt = rot_matrix * V_shipFilt';
V_shipFilt = V_shipFilt';

% 船只移动到海面中心
x_offset = mean(V_shipFilt(:,1)) - mean(x);
V_shipFilt(:,1) = V_shipFilt(:,1) - x_offset;
y_offset = mean(V_shipFilt(:,2)) - mean(y);
V_shipFilt(:,2) = V_shipFilt(:,2) - y_offset;
sea_depth = mean(h(:));
z_offset = mean(V_shipFilt(:,3)) - sea_depth;%需要知道船只吃水深度代替z_offset，把海浪平均高度作为吃水深度
V_shipFilt(:,3) = V_shipFilt(:,3) - z_offset;

[V_seaFilt,F_seaFilt]=surfboolean(V_seaOrigin,F_seaOrigin,'-',V_shipFilt,F_shipOrigin);
plotmesh(V_seaFilt,F_seaFilt);
%% 海水介电常数计算（密度为1.05）
salt = 0.35;
temper = 17;
eps1 = (87.134-0.1949*temper-0.01276*temper^2+0.0002491*temper^3)*(1+1.63e-5*temper*salt-0.003656*salt+3.21e-5*salt^2-4.232e-7*salt^3);
tau = (1.768e-11-6.086e-13*temper+1.104e-14*temper^2-8.111e-17*temper^3)*(1+2.282e-5*temper*salt-7.638e-4*salt-7.76e-6*salt^2+1.105e-8*salt^3);
fr = (1/tau)*2*pi;

%% 修改船只模型（alphashape版本）【包括gptool库函数】【失败】
[V_shipOrigin,F_shipOrigin,~] = readSTL('serie_60_una_pieza.stl');

V_shipFilt = V_shipOrigin*0.1;

% % 船只修改倾角（空间几何变换）
% phi_rad = deg2rad(ship_phi(time+1));
rot_matrix = roty(ship_phi(time+1));
% rot_matrix = roty(30);% 大角度测试
V_shipFilt = rot_matrix * V_shipFilt';
V_shipFilt = V_shipFilt';

% 船只移动到海面中心
x_offset = mean(V_shipFilt(:,1)) - mean(x);
V_shipFilt(:,1) = V_shipFilt(:,1) - x_offset;
y_offset = mean(V_shipFilt(:,2)) - mean(y);
V_shipFilt(:,2) = V_shipFilt(:,2) - y_offset;
sea_depth = mean(h(:));
z_offset = mean(V_shipFilt(:,3)) - sea_depth;%需要知道船只吃水深度代替z_offset，把海浪平均高度作为吃水深度
V_shipFilt(:,3) = V_shipFilt(:,3) - z_offset;

shp1=alphaShape(V_shipFilt(:,1),V_shipFilt(:,2),V_shipFilt(:,3));
[X, Y] = meshgrid(x, y);%建立海面网格
V_seaOrigin = reshape([X(:), Y(:), h(:)], [], 3);
shp2=alphaShape(V_seaOrigin(:,1),V_seaOrigin(:,2),V_seaOrigin(:,3));

figure;
plot(shp1);
figure;
plot(shp2);

id2=inShape(shp1,V_seaOrigin(:,1),V_seaOrigin(:,2),V_seaOrigin(:,3));
h_filt = h;
h_filt(id2)=NaN;


sea_filt = surf(X,Y,h_filt);
surf2stl(sea_filt, '1_sea_filt');
stlwrite('2_boat_filt.stl', F_shipOrigin,V_shipFilt);
%% 纵摇数据模型
% kesi_theta = 0.0657;
% w_theta = 0.3;
% dt = 1;
% time = 400;
% A = [0, 1; -(w_theta^2), -2*kesi_theta*w_theta];
% B = [0, w_theta^2]';
% C = [1, 0];
% D = 0;
% G_phi = ss(A, B, C, D);
% tt = 0:dt:time;
% u = he_t_deg;
% [y,t]=lsim(G_phi, u, tt);
% plot(t,y,'-');
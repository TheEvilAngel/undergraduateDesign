%% 1-0du-boat-differentPhi【超算】
% num_sam = [3, 3, 1, 150, 1];
% A = read_Feko_ffe_data(num_sam,1);

load('1-boat-0du-2L-modelMesh-hrrp-lowphi-beside_FarField1.mat');
E_phi = data{6,1};
E_theta = data{7,1};
A = E_theta + E_phi;

fre_Num = 150;
upsample = 50;
band = 9.5e9-9.3e9;
c = 3e8;
Nfft = fre_Num*upsample;
Rmax = c/(2*band)*fre_Num/cos(deg2rad(5));
dataRCS0 = A(:,1)';
dataRCS45 = A(:,2)';
dataRCS90 = A(:,3)';

figure;
plot(abs(dataRCS0(:)),'-o', LineWidth=1.7);
hold on
plot(abs(dataRCS45(:)),'-o', LineWidth=1.7);
plot(abs(dataRCS90(:)),'-o', LineWidth=1.7);
legend('phi = 0°', 'phi = 45°', 'phi = 90°');
title('海上目标雷达回波数据(phi=0°-90°)');


figure;
HRRP = ifftshift(ifft(dataRCS0(:),floor(Nfft)));
% HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-', LineWidth=1.7);
hold on

HRRP = ifftshift(ifft(dataRCS45(:),floor(Nfft)));
% HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-', LineWidth=1.7);

HRRP = ifftshift(ifft(dataRCS90(:),floor(Nfft)));
% HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-', LineWidth=1.7);
legend('phi = 0°', 'phi = 45°', 'phi = 90°');
title('海上目标HRRP(phi=0°-90°)');

%% 1,2-0du-boat-boatsea
load('1-boat-0du-2L-modelMesh-hrrp-lowphi-beside_FarField1.mat');
E_phi_boat = data{6,1};
E_theta_boat = data{7,1};
A_boat = E_theta_boat + E_phi_boat;

load('2-boatSea-0du-2L-modelMesh-hrrp-lowphi-beside_FarField1.mat');
E_phi_boatSea = data{6,1};
E_theta_boatSea = data{7,1};
A_boatSea = E_theta_boatSea + E_phi_boatSea;

fre_Num = 150;
upsample = 50;
band = 9.5e9-9.3e9;
c = 3e8;
Nfft = fre_Num*upsample;

Rmax = c/(2*band)*fre_Num/cos(deg2rad(5));
dataRCS0_boat = A_boat(:,1)';
dataRCS45_boat = A_boat(:,2)';
dataRCS90_boat = A_boat(:,3)';

dataRCS0_boatSea = A_boatSea(:,1)';
dataRCS45_boatSea = A_boatSea(:,2)';
dataRCS90_boatSea = A_boatSea(:,3)';

fre_line = linspace(9.3,9.5,fre_Num);

figure;
plot(fre_line,10*log10(abs(dataRCS0_boat(:))),'-o', LineWidth=1.7);
hold on
plot(fre_line,10*log10(abs(dataRCS0_boatSea(:))),'-.x', LineWidth=1.7);
legend('\fontname{宋体}目标', '\fontname{宋体}目标 + 海面','Location','SouthEast');
xlabel('\fontname{宋体}频率\fontname{Times New Roman}(GHz)');
ylabel('\fontname{Times New Roman}RCS(dBsm)');
set(gca,'FontName','Times New Roman','fontsize',16)
savefig(gcf,'result/4-yawphi0');
exportgraphics(gcf,'result/4-yawphi0.png');

figure;
plot(fre_line,10*log10(abs(dataRCS45_boat(:))),'-o', LineWidth=1.7);
hold on
plot(fre_line,10*log10(abs(dataRCS45_boatSea(:))),'-.x', LineWidth=1.7);
legend('\fontname{宋体}目标', '\fontname{宋体}目标 + 海面','Location','SouthEast');
xlabel('\fontname{宋体}频率\fontname{Times New Roman}(GHz)');
ylabel('\fontname{Times New Roman}RCS(dBsm)');
set(gca,'FontName','Times New Roman','fontsize',16)
savefig(gcf,'result/4-yawphi45');
exportgraphics(gcf,'result/4-yawphi45.png');

figure;
plot(fre_line,10*log10(abs(dataRCS90_boat(:))),'-o', LineWidth=1.7);
hold on
plot(fre_line,10*log10(abs(dataRCS90_boatSea(:))),'-.x', LineWidth=1.7);
legend('\fontname{宋体}目标', '\fontname{宋体}目标 + 海面','Location','SouthEast');
xlabel('\fontname{宋体}频率\fontname{Times New Roman}(GHz)');
ylabel('\fontname{Times New Roman}RCS(dBsm)');
set(gca,'FontName','Times New Roman','fontsize',16)
savefig(gcf,'result/4-yawphi90');
exportgraphics(gcf,'result/4-yawphi90.png');



figure;
HRRP = ifftshift(ifft(dataRCS0_boat(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-', LineWidth=1.7);
hold on
HRRP = ifftshift(ifft(dataRCS0_boatSea(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-.', LineWidth=1.7);
legend('\fontname{宋体}目标', '\fontname{宋体}目标 + 海面');
xlabel('\fontname{宋体}距离\fontname{Times New Roman}(m)');
ylabel('\fontname{宋体}幅度');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/4-yaw-hrrp-phi0');
exportgraphics(gcf,'result/4-yaw-hrrp-phi0.png');

figure;
HRRP = ifftshift(ifft(dataRCS45_boat(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-', LineWidth=1.7);
hold on 
HRRP = ifftshift(ifft(dataRCS45_boatSea(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-.', LineWidth=1.7);
legend('\fontname{宋体}目标', '\fontname{宋体}目标 + 海面');
xlabel('\fontname{宋体}距离\fontname{Times New Roman}(m)');
ylabel('\fontname{宋体}幅度');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/4-yaw-hrrp-phi45');
exportgraphics(gcf,'result/4-yaw-hrrp-phi45.png');

figure;
HRRP = ifftshift(ifft(dataRCS90_boat(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-', LineWidth=1.7);
hold on
HRRP = ifftshift(ifft(dataRCS90_boatSea(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-.', LineWidth=1.7);
legend('\fontname{宋体}目标', '\fontname{宋体}目标 + 海面');
xlabel('\fontname{宋体}距离\fontname{Times New Roman}(m)');
ylabel('\fontname{宋体}幅度');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/4-yaw-hrrp-phi90');
exportgraphics(gcf,'result/4-yaw-hrrp-phi90.png');
%% 3,4-7du-boat-boatsea
load('3-boat-7du-2L-modelMesh-hrrp-lowphi-3D_FarField1.mat');
E_phi_boat = data{6,1};
E_theta_boat = data{7,1};
A_boat = E_theta_boat + E_phi_boat;

load('4-boatSea-7du-2L-modelMesh-hrrp-lowphi-3D_FarField1.mat');
E_phi_boatSea = data{6,1};
E_theta_boatSea = data{7,1};
A_boatSea = E_theta_boatSea + E_phi_boatSea;

fre_Num = 150;
upsample = 50;
band = 9.5e9-9.3e9;
c = 3e8;
Nfft = fre_Num*upsample;

Rmax = c/(2*band)*fre_Num/cos(deg2rad(5));
dataRCS0_boat = A_boat(:,1)';
dataRCS45_boat = A_boat(:,2)';
dataRCS90_boat = A_boat(:,3)';

dataRCS0_boatSea = A_boatSea(:,1)';
dataRCS45_boatSea = A_boatSea(:,2)';
dataRCS90_boatSea = A_boatSea(:,3)';

fre_line = linspace(9.3,9.5,fre_Num);

figure;
plot(fre_line,10*log10(abs(dataRCS0_boat(:))),'-o', LineWidth=1.7);
hold on
plot(fre_line,10*log10(abs(dataRCS0_boatSea(:))),'-.x', LineWidth=1.7);
legend('\fontname{宋体}目标', '\fontname{宋体}目标 + 海面','Location','SouthEast');
xlabel('\fontname{宋体}频率\fontname{Times New Roman}(GHz)');
ylabel('\fontname{Times New Roman}RCS(dBsm)');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/5-yawphi0');
exportgraphics(gcf,'result/5-yawphi0.png');

figure;
plot(fre_line,10*log10(abs(dataRCS45_boat(:))),'-o', LineWidth=1.7);
hold on
plot(fre_line,10*log10(abs(dataRCS45_boatSea(:))),'-.x', LineWidth=1.7);
legend('\fontname{宋体}目标', '\fontname{宋体}目标 + 海面','Location','SouthEast');
xlabel('\fontname{宋体}频率\fontname{Times New Roman}(GHz)');
ylabel('\fontname{Times New Roman}RCS(dBsm)');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/5-yawphi45');
exportgraphics(gcf,'result/5-yawphi45.png');

figure;
plot(fre_line,10*log10(abs(dataRCS90_boat(:))),'-o', LineWidth=1.7);
hold on
plot(fre_line,10*log10(abs(dataRCS90_boatSea(:))),'-.x', LineWidth=1.7);
legend('\fontname{宋体}目标', '\fontname{宋体}目标 + 海面','Location','SouthEast');
xlabel('\fontname{宋体}频率\fontname{Times New Roman}(GHz)');
ylabel('\fontname{Times New Roman}RCS(dBsm)');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/5-yawphi90');
exportgraphics(gcf,'result/5-yawphi90.png');



figure;
HRRP = ifftshift(ifft(dataRCS0_boat(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-', LineWidth=1.7);
hold on
HRRP = ifftshift(ifft(dataRCS0_boatSea(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-.', LineWidth=1.7);
legend('\fontname{宋体}目标', '\fontname{宋体}目标 + 海面');
xlabel('\fontname{宋体}距离\fontname{Times New Roman}(m)');
ylabel('\fontname{宋体}幅度');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/5-yaw-hrrp-phi0');
exportgraphics(gcf,'result/5-yaw-hrrp-phi0.png');

figure;
HRRP = ifftshift(ifft(dataRCS45_boat(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-', LineWidth=1.7);
hold on 
HRRP = ifftshift(ifft(dataRCS45_boatSea(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-.', LineWidth=1.7);
legend('\fontname{宋体}目标', '\fontname{宋体}目标 + 海面');
xlabel('\fontname{宋体}距离\fontname{Times New Roman}(m)');
ylabel('\fontname{宋体}幅度');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/5-yaw-hrrp-phi45');
exportgraphics(gcf,'result/5-yaw-hrrp-phi45.png');

figure;
HRRP = ifftshift(ifft(dataRCS90_boat(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-', LineWidth=1.7);
hold on
HRRP = ifftshift(ifft(dataRCS90_boatSea(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-.', LineWidth=1.7);
legend('\fontname{宋体}目标', '\fontname{宋体}目标 + 海面');
xlabel('\fontname{宋体}距离\fontname{Times New Roman}(m)');
ylabel('\fontname{宋体}幅度');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/5-yaw-hrrp-phi90');
exportgraphics(gcf,'result/5-yaw-hrrp-phi90.png');

%% 1,3-0,7du-boat【论文】
load('1-boat-0du-2L-modelMesh-hrrp-lowphi-beside_FarField1.mat');
E_phi_boat = data{6,1};
E_theta_boat = data{7,1};
A_boat = E_theta_boat + E_phi_boat;

load('3-boat-7du-2L-modelMesh-hrrp-lowphi-3D_FarField1.mat');
E_phi_boatSea = data{6,1};
E_theta_boatSea = data{7,1};
A_boatSea = E_theta_boatSea + E_phi_boatSea;

fre_Num = 150;
upsample = 50;
band = 9.5e9-9.3e9;
c = 3e8;
Nfft = fre_Num*upsample;

Rmax = c/(2*band)*fre_Num/cos(deg2rad(5));
dataRCS0_boat = A_boat(:,1)';
dataRCS45_boat = A_boat(:,2)';
dataRCS90_boat = A_boat(:,3)';

dataRCS0_boatSea = A_boatSea(:,1)';
dataRCS45_boatSea = A_boatSea(:,2)';
dataRCS90_boatSea = A_boatSea(:,3)';

fre_line = linspace(9.3,9.5,fre_Num);

figure;
plot(fre_line,10*log10(abs(dataRCS0_boat(:))),'-o', LineWidth=1.7);
hold on
plot(fre_line,10*log10(abs(dataRCS0_boatSea(:))),'-.x', LineWidth=1.7);
legend('\fontname{宋体}横摇角\fontname{Times New Roman}0°', '\fontname{宋体}横摇角\fontname{Times New Roman}7.3°','Location','SouthEast');
xlabel('\fontname{宋体}频率\fontname{Times New Roman}(GHz)');
ylabel('\fontname{Times New Roman}RCS(dBsm)');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/1-yawphi0');
exportgraphics(gcf,'result/1-yawphi0.png');


figure;
plot(fre_line,10*log10(abs(dataRCS45_boat(:))),'-o', LineWidth=1.7);
hold on
plot(fre_line,10*log10(abs(dataRCS45_boatSea(:))),'-.x', LineWidth=1.7);
legend('\fontname{宋体}横摇角\fontname{Times New Roman}0°', '\fontname{宋体}横摇角\fontname{Times New Roman}7.3°','Location','SouthEast');
xlabel('\fontname{宋体}频率\fontname{Times New Roman}(GHz)');
ylabel('\fontname{Times New Roman}RCS(dBsm)');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/1-yawphi45');
exportgraphics(gcf,'result/1-yawphi45.png');



figure;
plot(fre_line,10*log10(abs(dataRCS90_boat(:))),'-o', LineWidth=1.7);
hold on
plot(fre_line,10*log10(abs(dataRCS90_boatSea(:))),'-.x', LineWidth=1.7);
legend('\fontname{宋体}横摇角\fontname{Times New Roman}0°', '\fontname{宋体}横摇角\fontname{Times New Roman}7.3°','Location','SouthEast');
xlabel('\fontname{宋体}频率\fontname{Times New Roman}(GHz)');
ylabel('\fontname{Times New Roman}RCS(dBsm)');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/1-yawphi90');
exportgraphics(gcf,'result/1-yawphi90.png');
% title('海上目标正面雷达回波数据(phi=90°)');



figure;
HRRP = ifftshift(ifft(dataRCS0_boat(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-', LineWidth=1.7);
hold on
HRRP = ifftshift(ifft(dataRCS0_boatSea(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-.', LineWidth=1.7);
legend('\fontname{宋体}横摇角\fontname{Times New Roman}0°', '\fontname{宋体}横摇角\fontname{Times New Roman}7.3°');
xlabel('\fontname{宋体}距离\fontname{Times New Roman}(m)');
ylabel('\fontname{宋体}幅度');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/2-yaw-hrrp-phi0');
exportgraphics(gcf,'result/2-yaw-hrrp-phi0.png');
% title('海上目标HRRP(phi=0°)');

figure;
HRRP = ifftshift(ifft(dataRCS45_boat(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-', LineWidth=1.7);
hold on 
HRRP = ifftshift(ifft(dataRCS45_boatSea(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-.', LineWidth=1.7);
legend('\fontname{宋体}横摇角\fontname{Times New Roman}0°', '\fontname{宋体}横摇角\fontname{Times New Roman}7.3°');
xlabel('\fontname{宋体}距离\fontname{Times New Roman}(m)');
ylabel('\fontname{宋体}幅度');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/2-yaw-hrrp-phi45');
exportgraphics(gcf,'result/2-yaw-hrrp-phi45.png');
% title('海上目标HRRP(phi=45°)');

figure;
HRRP = ifftshift(ifft(dataRCS90_boat(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-', LineWidth=1.7);
hold on
HRRP = ifftshift(ifft(dataRCS90_boatSea(:),floor(Nfft)));
HRRP = HRRP/max(HRRP);
x = (0:Nfft-1)/Nfft*Rmax;
plot(x,abs(HRRP),'-.', LineWidth=1.7);
legend('\fontname{宋体}横摇角\fontname{Times New Roman}0°', '\fontname{宋体}横摇角\fontname{Times New Roman}7.3°');
xlabel('\fontname{宋体}距离\fontname{Times New Roman}(m)');
ylabel('\fontname{宋体}幅度');
set(gca,'FontName','Times New Roman','fontsize',16);
savefig(gcf,'result/2-yaw-hrrp-phi90');
exportgraphics(gcf,'result/2-yaw-hrrp-phi90.png');
% title('海上目标HRRP(phi=90°)');

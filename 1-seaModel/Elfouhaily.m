function [S,kp]= Elfouhaily(k,U10,age,varargin)
%S = Elfouhaily(k,U10,age)
%S = Elfouhaily(k,U10,age,phi)
%% 根据论文age<5，age=omega_c，逆波龄，与风领域有关，论文2-6下面和2-15，omega==omega_c??
%%
phi = 0;
if nargin == 4 % 输入了4个参数
    phi = varargin{1}; % 输入的第1个参数
end
%constants
g = 9.81; %gravity acceleration
Cd10N = 0.00144; %drag coefficient？？
ustar = sqrt(Cd10N)*U10;%friction velocity at the water surface，论文2-18下面uf，为啥？
km = 370.0; %论文2-4下面
cm = 0.23; %minimum phase speed at wavenumber km，论文2-18下面c(km)
sigma = 0.08*(1+4*age^(-3)); %论文2-14，!!!age=omega_c!!!
alphap = 0.006*age^(0.55); %generalized Phillips-Kitaigorodskii equilibrium range parameter for long waves
%长波的广义philips - kitaigorodskii平衡范围参数,论文2-12上面，原文为0.55次方？？？？
k0 = g/(U10^2);
kp = k0 * age^2; %wavenumber of the spectral peak 无因次峰峰值波数，论文2-7上面
cp = sqrt(g/kp); %phase speed at the spectral peak cp = U10/age，论文2-12上面

% 论文2-18
if (ustar <= cm) %alpham is the generalized Phillips-Kitaigorodskii equilibrium range parameter for short waves
    alpham = 0.01*(1 + log(ustar/cm));
else
    alpham = 0.01*(1 + 3*log(ustar/cm));
end

% 论文2-13
if (age <= 1)
    gamma = 1.7;
else
    gamma = 1.7 + 6*log(age);
end

c = sqrt((g./k).*(1 + (k/km).^2)); %wave phase speed，论文2-12上面
Lpm = exp(-5/4*(kp./k).^2);  %Pierson-Moskowitz shape spectrum，论文2-6下面
Gam = exp(-1/(2*sigma^2)*(sqrt(k/kp) - 1 ).^2 );%论文2-8
Jp = gamma.^Gam; %JONSWAP peak enhancement or "overshoot" factor，论文2-7
Fp = Lpm.*Jp.*exp(-age/sqrt(10)*(sqrt(k/kp) - 1) ); %long-wave side effect function，论文2-12
Fm = Lpm.*Jp.*exp(-0.25*(k/km - 1).^2); %short-wave side effect function，论文2-17
Bl = 0.5*alphap*(cp./c).*Fp; % 论文2-11
Bh = 0.5*alpham*(cm./c).*Fm; % 论文2-16

S = (Bl + Bh)./(k.^3);

% %% compute the spreading function
% %论文2-25
% a0 = log(2)/4;%改动，原本为log2/2
% ap = 4;
% am = 0.13*ustar/cm;
% 
% Delk = tanh(a0 + ap*(c/cp).^(2.5) + am*(cm./c).^(2.5));
% S = 1*S.*1/(2).*(1 + Delk.*cos(2*phi)); % 认为phi=0，实际上为S=S, 此步骤为何意

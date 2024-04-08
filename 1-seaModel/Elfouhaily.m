function [S,kp]= Elfouhaily(k,U10,age,varargin)
%S = Elfouhaily(k,U10,age)
%S = Elfouhaily(k,U10,age,phi)
%% ��������age<5��age=omega_c���沨�䣬��������йأ�����2-6�����2-15��omega==omega_c??
%%
phi = 0;
if nargin == 4 % ������4������
    phi = varargin{1}; % ����ĵ�1������
end
%constants
g = 9.81; %gravity acceleration
Cd10N = 0.00144; %drag coefficient����
ustar = sqrt(Cd10N)*U10;%friction velocity at the water surface������2-18����uf��Ϊɶ��
km = 370.0; %����2-4����
cm = 0.23; %minimum phase speed at wavenumber km������2-18����c(km)
sigma = 0.08*(1+4*age^(-3)); %����2-14��!!!age=omega_c!!!
alphap = 0.006*age^(0.55); %generalized Phillips-Kitaigorodskii equilibrium range parameter for long waves
%�����Ĺ���philips - kitaigorodskiiƽ�ⷶΧ����,����2-12���棬ԭ��Ϊ0.55�η���������
k0 = g/(U10^2);
kp = k0 * age^2; %wavenumber of the spectral peak ����η��ֵ����������2-7����
cp = sqrt(g/kp); %phase speed at the spectral peak cp = U10/age������2-12����

% ����2-18
if (ustar <= cm) %alpham is the generalized Phillips-Kitaigorodskii equilibrium range parameter for short waves
    alpham = 0.01*(1 + log(ustar/cm));
else
    alpham = 0.01*(1 + 3*log(ustar/cm));
end

% ����2-13
if (age <= 1)
    gamma = 1.7;
else
    gamma = 1.7 + 6*log(age);
end

c = sqrt((g./k).*(1 + (k/km).^2)); %wave phase speed������2-12����
Lpm = exp(-5/4*(kp./k).^2);  %Pierson-Moskowitz shape spectrum������2-6����
Gam = exp(-1/(2*sigma^2)*(sqrt(k/kp) - 1 ).^2 );%����2-8
Jp = gamma.^Gam; %JONSWAP peak enhancement or "overshoot" factor������2-7
Fp = Lpm.*Jp.*exp(-age/sqrt(10)*(sqrt(k/kp) - 1) ); %long-wave side effect function������2-12
Fm = Lpm.*Jp.*exp(-0.25*(k/km - 1).^2); %short-wave side effect function������2-17
Bl = 0.5*alphap*(cp./c).*Fp; % ����2-11
Bh = 0.5*alpham*(cm./c).*Fm; % ����2-16

S = (Bl + Bh)./(k.^3);

% %% compute the spreading function
% %����2-25
% a0 = log(2)/4;%�Ķ���ԭ��Ϊlog2/2
% ap = 4;
% am = 0.13*ustar/cm;
% 
% Delk = tanh(a0 + ap*(c/cp).^(2.5) + am*(cm./c).^(2.5));
% S = 1*S.*1/(2).*(1 + Delk.*cos(2*phi)); % ��Ϊphi=0��ʵ����ΪS=S, �˲���Ϊ����

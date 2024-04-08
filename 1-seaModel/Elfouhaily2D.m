function PSI= Elfouhaily2D(k,phi,U10,age,phi_w)
%PSI = Elfouhaily(k,phi, U10,age)
% 
% %constants
g = 9.81; %gravity acceleration
Cd10N = 0.00144; %drag coefficient 阻力系数?
ustar = sqrt(Cd10N)*U10;%friction velocity at the water surface 水面上的摩擦速度?
km = 370.0; % rad/s
cm = 0.23; %minimum phase speed at wavenumber km?
k0 = g/(U10^2);%?
kp = k0 * age^2; %wavenumber of the spectral peak?
cp = sqrt(g/kp); %phase speed at the spectral peak cp = U10/age?
c = sqrt((g./k).*(1 + (k/km).^2)); %wave phase speed 波浪相速度

S = Elfouhaily(k,U10,age);
%% compute the spreading function
%论文2-25
a0 = log(2)/4;%改动，原本为log2/2
ap = 4;
am = 0.13*ustar/cm;

Delk = tanh(a0 + ap*(c/cp).^(2.5) + am*(cm./c).^(2.5));
PSI = 1*S.*1./k.*1/(2*pi).*(1 + Delk.*cos(2*(phi-phi_w)));%论文2-24关系

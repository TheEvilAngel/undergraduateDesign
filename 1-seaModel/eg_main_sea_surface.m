clc
clear all
% 此为海面生成的主程序，为了方便调用sea_surface.m而编写，便于控制主输入。
nx = 1001; xmin =-300; xmax = 300; x = linspace(xmin,xmax,nx); 
ny = 501; ymin =-150; ymax = 150; y = linspace(ymin,ymax,ny);
wind_data.U = 2; wind_data.thetaU = 0; wind_data.X = 1e6;   %wind_data.U = 3,对应二级海情; 5,对应三级海情；2,对应一级海情；7,对应四级海情
[s,Tp,fm,B,Sk,kx,ky] = sea_surface(x,y,wind_data,'PM','none');
figure(1)
%subplot(211)
mesh(x,y,s),ylabel('y(m)'),xlabel('x(m)'); 
axis([-150 150 -50 50 -1 5])% 设置画图范围，可以去掉
%subplot(212),mesh(kx,ky,Sk),ylabel('ky (1/m)'),xlabel('kx (1/m)'); 
%stlwrite('sea5m_601.stl',x,y,s,'mode','ascii')
%surf2stl('sea_surface.stl',x,y,s);

%
%shading interp %
%colormap(gray);
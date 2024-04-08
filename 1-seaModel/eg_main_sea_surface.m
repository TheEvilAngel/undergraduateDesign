clc
clear all
% ��Ϊ�������ɵ�������Ϊ�˷������sea_surface.m����д�����ڿ��������롣
nx = 1001; xmin =-300; xmax = 300; x = linspace(xmin,xmax,nx); 
ny = 501; ymin =-150; ymax = 150; y = linspace(ymin,ymax,ny);
wind_data.U = 2; wind_data.thetaU = 0; wind_data.X = 1e6;   %wind_data.U = 3,��Ӧ��������; 5,��Ӧ�������飻2,��Ӧһ�����飻7,��Ӧ�ļ�����
[s,Tp,fm,B,Sk,kx,ky] = sea_surface(x,y,wind_data,'PM','none');
figure(1)
%subplot(211)
mesh(x,y,s),ylabel('y(m)'),xlabel('x(m)'); 
axis([-150 150 -50 50 -1 5])% ���û�ͼ��Χ������ȥ��
%subplot(212),mesh(kx,ky,Sk),ylabel('ky (1/m)'),xlabel('kx (1/m)'); 
%stlwrite('sea5m_601.stl',x,y,s,'mode','ascii')
%surf2stl('sea_surface.stl',x,y,s);

%
%shading interp %
%colormap(gray);
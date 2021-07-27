clc,clear all,close all

NX=[6]; %an increasing sequence of the number of grid points
E=zeros(size(NX)); H=zeros(size(NX)); %pre-allocation of memory

for kk=1:length(NX)
axy=zeros(NX(kk),NX(kk)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
nx=NX(kk);
ny=nx;             %equal number of grid points in x and y directions
hx=2/(nx-1);       %grid length in x direction
hy=2/(ny-1);       %grid length in y direction
H(kk)=hx;
% for i=1:nx
%     j=1:nx;
%    axy(j,i)=1;%1+sin((i-1)*hx-1)*cos((j-1)*hy-1);
% end

x=(-1:hx:1);
y=(-1:hy:1)';
axy=1+sin(x).*cos(y);
t_final=1;
dt=0.5*hx;
nt = floor(t_final/dt)+1;
dt = t_final/nt;

w=10; kx=6; ky=4;
nxl=nx/2;
hxl=2/(nxl-1); 
%Set up initial data and obtain u0 and u1
u0=u_exact(0,x,y,w,kx,ky);
f2=f2_fun(x,y,w,kx,ky);
L1=laplac(u0(:,1:nxl+1),hx,hy,nx,nxl+1,axy(:,1:nxl+1)); %send an extra column of u0 and axy:send extra 2 if in the middle (0:nxl+1)
L1=L1(:,1:nxl) %but delete the columns for the final L: if end, then delete one row, if middle then delete outside rows
%who is the outside?
%if ix_off = 0 OR ix_off = x_max-nxl ----- what is x_max? it is nx
 L2=laplac(u0(:,nxl:end),hx,hy,nx,nxl+1,axy(:,nxl:end));
 L2=L2(:,2:end)

f=forcing(0,x,y,w,kx,ky);
u1=u0+dt*f2+0.5*(dt^2)*(L+f);

%update BC
u1(1,:)=u_exact(dt,x,-1,w,kx,ky);
u1(ny,:)=u_exact(dt,x,1,w,kx,ky);
u1(:,1)=u_exact(dt,-1,y,w,kx,ky);
u1(:,nx)=u_exact(dt,1,y,w,kx,ky);

for k=0:nt-2
    t=(k+1)*dt;
    %obtin right hand side
    f=forcing(t,x,y,w,kx,ky);
    L=laplac(u1,hx,hy,nx,ny,axy);
    %March in time
    u2=2*u1-u0+(dt^2)*(L+f);
    %update BC
    u2(1,:)=u_exact(t+dt,x,-1,w,kx,ky);
    u2(ny,:)=u_exact(t+dt,x,1,w,kx,ky);
    u2(:,1)=u_exact(t+dt,-1,y,w,kx,ky);
    u2(:,nx)=u_exact(t+dt,1,y,w,kx,ky);
    %switch solution at different time levels
    u0=u1;
    u1=u2;    
end
%compute the error
u_ex_final=u_exact(t_final,x,y,w,kx,ky);
E(kk)=max(max(abs(u2-u_ex_final)));
end

figure
surf(x,y,u2); shading interp

figure
loglog(H,E,'o-',H,H.^2,'r')
legend('\epsilon(h)','h^2')
xlabel('h')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u=u_exact(t,x,y,w,kx,ky)
u=sin(w*t-kx*x).*sin(ky*y);
%u=u';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=f2_fun(x,y,w,kx,ky)
f=w*cos(kx*x).*sin(ky*y);
%f=f';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function force=forcing(t,x,y,w,kx,ky)
auxx=-kx*cos(w*t-kx*x).*cos(x).*cos(y).*sin(ky*y)-(kx^2)*sin(w*t-kx*x).*sin(ky*y).*(1+sin(x).*cos(y));
auyy=ky*sin(w*t-kx*x).*(-sin(x).*sin(y).*cos(ky*y)-ky.*sin(ky*y).*(1+sin(x).*cos(y)));
utt=-w^2*sin(w*t-kx*x).*sin(ky*y);
force=utt-auxx-auyy;
%  force=-(w^2-kx^2-ky^2)*sin(w*t-kx*x).*sin(ky*y);
%force=force';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L=laplac(u,hx,hy,nx,ny,a)
L=zeros(size(u));
for j=2:ny-1
    I=2:nx-1;
    %a(j,I)=1+sin((I-2)*hx-1)'*cos((j-2)*hy-1);
    L(I,j)=(( a(I+1,j)+a(I,j)  )/(2*hx^2)).*u(I+1,j)+(( a(I,j)+a(I-1,j) )/(2*hx^2)).*u(I-1,j)+( ( a(I,j+1)+a(I,j) ) /(2*hy^2)).*u(I,j+1)+( ( a(I,j)+a(I,j-1) ) /(2*hy^2)).*u(I,j-1)-((( a(I+1,j)+2*a(I,j)+a(I-1,j) )/(2*hx^2))+(( a(I,j+1)+2*a(I,j)+a(I,j-1) ) /(2*hy^2))).*u(I,j);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
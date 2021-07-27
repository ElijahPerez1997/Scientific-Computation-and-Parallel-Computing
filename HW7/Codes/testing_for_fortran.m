
format long
NX=[21 41 81 161 321]; %an increasing sequence of the number of grid points
%NX=[6]
E=zeros(size(NX)); H=zeros(size(NX)); %pre-allocation of memory

for kk=1:length(NX)
axy=zeros(1,NX(kk)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
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
t_final=2;
dt=0.5*hx;
nt = floor(t_final/dt)+1;

dt = t_final/nt;

w=10; kx=6; ky=4;

%Set up initial data and obtain u0 and u1
u0=u_exact(0,x,y,w,kx,ky);
f2=f2_fun(x,y,w,kx,ky);
L=laplac(u0,hx,hy,nx,ny,axy);
f=forcing(0,x,y,w,kx,ky);
u1=u0+dt*f2+0.5*(dt^2)*(L+f);

%update BC
u1(1,:)=u_exact(dt,x,-1,w,kx,ky);
u1(ny,:)=u_exact(dt,x,1,w,kx,ky);
u1(:,1)=u_exact(dt,-1,y,w,kx,ky);
u1(:,nx)=u_exact(dt,1,y,w,kx,ky);

for k=0:nt-2 %nt-2 = 7, nt=9
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
% 
% figure
% surf(x,y,u2); shading interp
% 
% figure
% loglog(H,E,'o-',H,H.^2,'r')
% legend('\epsilon(h)','h^2')
% xlabel('h')
% print -dpng error.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=u_exact(t,x,y,w,kx,ky) %NO CHANGES
u=sin(w*t-kx*x).*sin(ky*y);
%u=u';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=f2_fun(x,y,w,kx,ky) %NO CHANGES
f=w*cos(kx*x).*sin(ky*y);
%f=f';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function force=forcing(t,x,y,w,kx,ky)
t1 = sin(ky*y); %n*1
t2 = w*t-kx*x; %1*n
t3 = (1+ (sin(x))'*(cos(y))')'; %n*n

%  x1=sin(x)
%  y1=cos(y)
%   for i=1:5
%     for j=1:5
%       t3new(j,i) = x1(1,i)*y1(j,1);
%     end
%   end
% t3new+1
a1 = (-kx*cos(w*t-kx*x)); %1*n
a2 = cos(x); %1*n
a3 = cos(y); %n*1
a4 = (kx^2)*sin(t2); %1*n

%todo a1.*a2 where both are same dim and output is same dim (1*n)
for i=1:length(a1)
    a1(i)= a1(i)*a2(i);
end

%AUXX SHOULD BE SQUAre from start
% a1size=size(a1) %1*n
% a3size=size(a3) %n*1

% todo: auxx=a1.*a3 where a1=1*n and a3=n*1 and output=n*n
for i=1:length(a1)
    for j=1:length(a1)
    auxx(i,j)= a1(i)*a3(j);
    end
end
auxx=auxx';

%todo auxx.*t1 where auxx=n*n and t1=n*1 and output = n*n
for i=1:length(t1)
    for j=1:length(t1)
    auxx(i,j)= t1(i)*auxx(i,j);
    end
end

%todo auxx2=a4.*t1 where a4=1*n and t1=n*1 and output=n*n
for i=1:length(a4)
    for j=1:length(a4)
    auxx2(i,j)= a4(i)*t1(j);
    end
end

auxx2=auxx2';

%this is n*n * n*n, can use MATMUL
auxx2=auxx2.*t3;

auxx=auxx-auxx2;

t5 = cos(ky*y); %n*1
a2 = -sin(x); %1*n
a3 = sin(y); %n*1
a4 = ky*sin(t2); %1*n

%a2.*a3 = 1*n * n*1
for i=1:length(a2)
    for j=1:length(a2)
    auyy(i,j)= a2(i)*a3(j);
    end
end
auyy=auyy';

%auyy.*t5 = n*n * n*1
for i=1:length(t5)
    for j=1:length(t5)
    auyy(i,j)= t5(i)*auyy(i,j);
    end
end

a3 = ky*t1; %scalar * vector
%a3.*t3 = n*1 * n*n
for i=1:length(t1)
    for j=1:length(t1)
    auyy2(i,j)= a3(i)*t3(i,j);
    end
end
auyy2 = auyy-auyy2;

%auyy=a4.*auyy2 = 1*n * n*n
for i=1:length(t1)
    for j=1:length(t1)
    auyy(i,j)= a4(j)*auyy2(i,j);
    end
end

a2 = -w^2*sin(t2); %1*n
%this is 1*n * n*1, can use MATMUL
utt=times(t1,a2);

force=utt-auxx-auyy;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L=laplac(u,hx,hy,nx,ny,a)
L=zeros(size(u));
%%%%%%%I THINK THIS IS JUST AN ELEMENT WISE LOOP AND WILL WORK WITH REGULAR
%%%%%%%MULTIPLICATION 
for j=2:ny-1
    i=2:nx-1;
    t1 = 2*hx^2; %constant
    t1 = 2*hy^2;
    a1 = ((a(i+1,j)+a(i,j))/t1); %all a's are scalars
    a2 = ((a(i,j)+a(i-1,j))/t1); 
    a3 = ((a(i,j+1)+a(i,j))/t1); 
    a4 = ((a(i,j)+a(i,j-1))/t1);
    a5 = ((a(i+1,j)+2*a(i,j)+a(i-1,j))/t1)+((a(i,j+1)+2*a(i,j)+a(i,j-1))/t1);
    
    %these are (n-2)*1 * (n-2)*1 to output (n-2)*1 
    L(i,j)=a1.*u(i+1,j);
%     t2 = u(i+1,j)
%     %todo a1.*a2 where both are same dim and output is same dim (1*n)
%     for k=1:nx-2
%         t(k)= a1(k)*a2(k);
%     end
%     L(i,j) = t;
    
    L(i,j)=L(i,j)+a2.*u(i-1,j);
    L(i,j)=L(i,j)+a3.*u(i,j+1);
    L(i,j)=L(i,j)+a4.*u(i,j-1);
    L(i,j)=L(i,j)-(a5.*u(i,j));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

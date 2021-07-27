clc,clear,close all
%initalization
%Number of birds
n=50;
%number of repelling birds
L=2;
s=-2;
%Food attraction
gamma1=4;
%Leader attraction
gamma2=5;
%Flock attraction
gamma3=.5;
%bird repulsion
p=0.8;
delta=0.1;


b0=zeros(2*n,1);
%initial condition
for i=1:2*n
    b0(i)=10*rand();
end
b_prime1=@(t,b) gamma1*([4*sin(2*t);1*sin(2*t)]-b);

[xdif,tdif]=rk4(b_prime1,0,10,100,[b0(1);b0(n+1)]);
ct=zeros(2,101);

for i =1:101
   ct(1,i)=2*sin(2*((i-1)*0.1)); 
   ct(2,i)=2*cos(2*((i-1)*0.1));
end

xfinal=zeros(n-1,101);
yfinal=zeros(n-1,101);

for k=2:n
%fl_force_x=@(t,b) gamma3*((sum(b(1:n))/n) - b(k));

%rep_force=p*sum(sum(( (b(k)-b(1:L))/((b(k)-b(1:L)).^2 +delta))))

%We will pick bird s to be the smelly bird
b_primek_x=@(t,b) gamma2*(xdif(1,round(t*10)+1)-b(k))+ gamma3*((sum(xdif(1,round(t*10)+1)+b(2:n))/n) - b(k)) + p*sum(sum(((b(k)-(xdif(1,round(t*10)+1)+b(2:L))/((b(k)-(xdif(1,round(t*10)+1)+b(2:L))).^2 +delta)))))+s*(b(2)-b(k));
%b_primek_x=@(t,b) gamma2*(b(1)-b(k))+ gamma3*((sum(b(1:n))/n) - b(k)) + p*sum(sum(((b(k)-b(1:L))/((b(k)-b(1:L)).^2 +delta))));
%b_primek_y=@(t,b) gamma2*(b(n)-b(k))+ gamma3*((sum(b(n+1:end))/n) - b(k)) + p*sum(sum(((b(k)-b(1+n:L+n))/((b(k)-b(1+n:L+n)).^2 +delta))));
b_primek_y=@(t,b) gamma2*(xdif(2,round(t*10)+1)-b(k))+ gamma3*((sum(xdif(2,round(t*10)+1)+b(2:n))/n) - b(k)) + p*sum(sum(((b(k)-(xdif(2,round(t*10)+1)+b(2:L))/((b(k)-(xdif(2,round(t*10)+1)+b(1:L))).^2 +delta)))))+s*(b(2)-b(k));

[xdif2,tdif2]=rk4(b_primek_x,0,10,100,b0(1:n));
xfinal(k,:)=xdif2(k,:);
[ydif2,tdif3]=rk4(b_primek_y,0,10,100,b0(n+1:n+n));
yfinal(k,:)=ydif2(k,:);

%  plot(xdif(1,:),tdif)
%  hold on
%  plot(xdif2(k,:),tdif2)
%  plot(ct(1,:),tdif2, 'o')

%     plot(xdif22(k,:),ydif22(k,:))
%     hold on
%     plot(xdif(1,:),xdif(2,:))
%     plot(ct(1,:),ct(2,:),'-o')

end

% Initialize video
myVideo = VideoWriter('cSine'); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me the higher the faster
open(myVideo)
% Plot in a loop and grab frames
for i=1:101
    plot(xdif(1,i),xdif(2,i),'r*');
    hold on
    plot(ct(1,i),ct(2,i),'go')
    plot(xfinal(2,i),yfinal(2,i),'d')
    for j=3:n
    plot(xfinal(j,i), yfinal(j,i),'b*');
    end
    ylim([-10,10])
    xlim([-10,10])
    pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    hold off
end
close(myVideo)
max(xfinal,[],'all')
max(yfinal,[],'all')
%xdif
%useful for later
%bx=[b_prime(1); b_primek_x];
%by=[b_prime(2); b_primek_y];
clc,clear,close all
cores=[1,2,5,8,16];
time=[12.3127582073212,6.51575088500977,2.80028295516968,1.83450508117676,1.16910195350647];
speed_up=zeros(1,length(time));
speed_up(1)=0;
for i=1:length(time)-1
    speed_up(i+1)=time(1)-time(i+1);
end

plot(cores,time,'-o')
hold on
plot(cores,speed_up,'-o')
legend('Runtime','Speed up')
title('Runtime and Speed Up for Strong Scaling')
xlabel('Cores')
ylabel('Time (in seconds)')
print -dpng strong_time.png
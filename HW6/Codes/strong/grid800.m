clear all
close all
load time800.txt 
load error800.txt
load s_error800.txt

thread = 1:16;
error800 = round(error800,7)
s_error800 = round(s_error800,7)
speedup = zeros(1,15);
%Compute speedup
for i=1:15
    speedup(i) = time800(1) - time800(i+1);
end

semilogy(thread,time800(1:16),'b-o', 'Markersize',5,'Linewidth',1)
hold on
semilogy(thread(2:16),speedup,'r-*','Markersize',5,'Linewidth',1)
xlabel('Num of Threads')
ylabel('Elapsed Time(s)')
title('Time for Grid Size = 800')
legend('Avg Parallel','Speedup')

print -dpng timeVthread800.png

hold off

semilogy(thread,error800(1:16),'b-o')
hold on
semilogy(thread,s_error800(1:16),'r-*')
xlabel('Num of Threads')
ylabel('Error')
title('Error for Grid Size = 800')
legend('Parallel','Serial')
print -dpng errorVthread800.png

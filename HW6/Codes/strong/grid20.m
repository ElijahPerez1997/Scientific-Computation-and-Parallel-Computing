clear all
close all
load time20.txt 
load error20.txt
load s_error20.txt

thread = 1:16;
error20 = round(error20,7)
s_error20 = round(s_error20,7)
speedup = zeros(1,15);
%Compute speedup
for i=1:15
    speedup(i) = time20(1) - time20(i+1);
end

semilogy(thread,time20(1:16),'b-o', 'Markersize',5,'Linewidth',1)
hold on
semilogy(thread(2:16),speedup,'r-*','Markersize',5,'Linewidth',1)
xlabel('Num of Threads')
ylabel('Elapsed Time(s)')
title('Time for Grid Size = 20')
legend('Avg Parallel','Speedup')

print -dpng timeVthread20.png

hold off

semilogy(thread,error20(1:16),'b-o')
hold on
semilogy(thread,s_error20(1:16),'r-*')
xlabel('Num of Threads')
ylabel('Error')
title('Error for Grid Size = 20')
legend('Parallel','Serial')
print -dpng errorVthread20.png

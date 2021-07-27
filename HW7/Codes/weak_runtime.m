clc,clear,close all
cores=[1,2,5,8,16];
times=[0.289858818054199,0.176476955413818,0.814860105514526,1.07249307632446,1.61291003227234];
plot(cores,times,'-o')
title('Runtime for Weak Scaled Problem')
xlabel('Cores')
ylabel('Runtime (in seconds)')
print -dpng weak_time.png

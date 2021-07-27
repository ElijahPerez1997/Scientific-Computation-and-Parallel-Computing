clear all
close all
load h.txt 
load error.txt 
load errorLaplace.txt


loglog(h,error,'LineWidth',2)
title('H vs Error for Derivatives')
ylabel('error')
xlabel('h')
print -dpng error.png
loglog(h,errorLaplace,'r','LineWidth',2)
title('H vs Error for Laplace')
ylabel('error')
xlabel('h')
print -dpng errorLaplace.png

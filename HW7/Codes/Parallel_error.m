clc,clear,close all
y=[5.7735583462839690E-002,1.5656569620089433E-002,4.0824234039578133E-003,1.0547128832899544E-003, 2.6264596246050598E-004];
h=[0.1,0.05,0.025,0.0125,0.00625];
loglog(h,y,'bo-')
hold on
loglog(h,h.^2,'k')
legend('\epsilon(h)','h^2')
xlabel('h')
ylabel('Error')
title('Error in Parallel Code')
print -dpng error_parallel

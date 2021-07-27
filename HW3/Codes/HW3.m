xt1=(1:length(1));
xt2=(1:length(2));
xg1=(1:length(3));
xg2=(1:length(4));


loglog(xt1,Trap_pi1,'LineWidth',2)
title('Error For Numerical Integration')
ylabel('error')
xlabel('n')
hold on
loglog(xt2,Trap_pi2,'LineWidth',2)
loglog(xg1,Gauss_pi1,'LineWidth',2)
loglog(xg2,Gauss_pi2,'LineWidth',2)
legend('Trap: k=\pi','Trap: k=\pi^2','Gauss: k=\pi','Gauss: k=\pi^2')
saveas(gcf,'HW_plot.png')
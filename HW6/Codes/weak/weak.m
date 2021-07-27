clear all
close all

%thread VS time/gridsize
error = [3.0264977100469928E-005,1.5095998978667389E-005,1.0093342443876056E-005,7.5494846646479528E-006,6.0439564234061329E-006,5.0288523097925204E-006,4.3141177051550715E-006,3.7681058247036378E-006,3.3528573500402834E-006,3.0216981653126210E-006,2.7455530141264900E-006,2.5128433071525104E-006,2.3213460315151302E-006,2.1566915837933680E-006,2.0089540282990944E-006,1.8852892468372988E-006];
time = [6.7313499999999998E-002,0.13682775000000003,0.17810000474855769,0.21716982423095033,0.30293667776277289,0.29169603502668906,0.33606764399155509,0.36300336298882030,0.39572087724809535,0.43016748725494836,0.47698137900442816,0.52926912123803049,0.56982174822769593,0.62376860198855866,0.63565871299942955,0.72528917000454385];
grid = [200,283,346,400,447,490,529,566,600,632,663,693,721,748,775,800];
thread = 1:16;
%scale = time./grid;

%serial results
s_error=[3.0264977100469928E-005,1.5095998978667385E-005,1.0093342443876054E-005,7.5494846646479528E-006,6.0439564234061320E-006, 5.0288523097925204E-006,4.3141177051550715E-006,3.7681058247036378E-006,3.3528573500402825E-006,3.0216981653126206E-006,2.7455530141264904E-006,2.5128433071525104E-006,2.3213460315151302E-006,2.1566915837933676E-006,2.0089540282990948E-006,1.8852892468372986E-006];
s_time=[6.7313499999999998E-002,0.13682775000000003,0.21344675000000002,0.27709375000000014,0.39371650000000002,0.40152175000000012,0.46290024999999968,0.52714524999999934,0.58367049999999931,0.66867950000000231,0.70069375000000100,0.78823724999999634,0.87331525000000099,1.0002342500000037,1.0247650000000021,1.1184422500000011];
%s_scale = s_time./grid;

speedup = s_time-time;

semilogy(grid,time,'b-o', 'Markersize',5,'Linewidth',1)
hold on
semilogy(grid,s_time,'r-*','Markersize',5,'Linewidth',1)
semilogy(grid,speedup,'g-^','Markersize',5,'Linewidth',1)
xlabel('Grid Size')
ylabel('Elapsed Time(s)')
title('Weak Scaling: Time')
legend('Avg Parallel','Serial','Speedup')
print -dpng weaktime.png

hold off

semilogy(grid,error,'b-o')
hold on
semilogy(grid,s_error,'r-*')
xlabel('Grid Size')
ylabel('Error')
title('Weak Scaling: Error ')
legend('Parallel','Serial')
print -dpng weakerror.png

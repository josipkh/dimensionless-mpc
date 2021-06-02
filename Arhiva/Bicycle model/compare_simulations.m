s1 = load('/home/josip/Chalmers/Diplomski/MATLAB/Simulation results/bicycle_full/bicycle_full.mat');
s2 = load('/home/josip/Chalmers/Diplomski/MATLAB/Simulation results/bicycle_full_pi/bicycle_full_pi.mat');
figure
plot(s1.uhist-s2.uhist)
figure
plot(s1.xhist'-s2.xhist')
legend('y','dy','\theta','d\theta')
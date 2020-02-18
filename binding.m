function [t,z] = Binding

close all;
clear all;
    
% Set time of simulation [s] (0 to 0.5 hour)
tspan = [0 1200];                          

% Set initial concentrations of ligands and proteins
z0 = [0 0 0.790751372 0.395375686]; % F, GL, G, L (Local Concentration)

% Simulate
[t,z] = ode45(@DiffCalc1,tspan,z0); %kcat = 1
[x,y] = ode45(@DiffCalc2,tspan,z0); %kcat = .1
[a,b] = ode45(@DiffCalc3,tspan,z0); %kcat = .01
[c,d] = ode45(@DiffCalc4,tspan,z0); %kcat = .001

figure(1);
hold on;
plot(t,z(:,1));
hold on;
plot(x,y(:,1));
hold on;
plot(a,b(:,1));
hold on;
plot(c,d(:,1));

%Output our results
%disp([t,z(:,1)]);

%Write to a file
dlmwrite('Binding_1000nM.csv', [t,z(:,1)])
dlmwrite('Binding_100nM.csv', [x,y(:,1)])
dlmwrite('Binding_10nM.csv', [a,b(:,1)])
dlmwrite('Binding_1nM.csv', [c,d(:,1)])

return;

function d_dt = DiffCalc1(t,z)

% Set parameters
k.sig = 0.033333333; %Rate of mCherry production dF/dt [sec-1]
k.deg = 0.002128988; %Rate of protein degredation [sec-1]

k.on = 0.25; %On rate of GPCR binding dGL/dt [uM-1 sec-1]
k.off = 0.25; %Off rate of GPCR binding dGL/dt [sec-1]

d_dt = zeros(4,1); %Init a matrix with 0

F = z(1,1);
GL = z(2,1);
G = z(3,1);
L = z(4,1);

d_dt(1) = k.sig*GL - k.deg*F;
d_dt(2) = k.on*G*L - k.off*GL;
d_dt(3) = k.off*GL - k.on*G*L;

return;

function d_dt = DiffCalc2(t,z)

% Set parameters
k.sig = 0.033333333; %Rate of mCherry production dF/dt [sec-1]
k.deg = 0.002128988; %Rate of protein degredation [sec-1]

k.on = 0.25; %On rate of GPCR binding dGL/dt [uM-1 sec-1]
k.off = 0.025; %Off rate of GPCR binding dGL/dt [sec-1]

d_dt = zeros(4,1); %Init a matrix with 0

F = z(1,1);
GL = z(2,1);
G = z(3,1);
L = z(4,1);

d_dt(1) = k.sig*GL - k.deg*F;
d_dt(2) = k.on*G*L - k.off*GL;
d_dt(3) = k.off*GL - k.on*G*L;

return;

function d_dt = DiffCalc3(t,z)

% Set parameters
k.sig = 0.033333333; %Rate of mCherry production dF/dt [sec-1]
k.deg = 0.002128988; %Rate of protein degredation [sec-1]

k.on = 0.25; %On rate of GPCR binding dGL/dt [uM-1 sec-1]
k.off = 0.0025; %Off rate of GPCR binding dGL/dt [sec-1]

d_dt = zeros(4,1); %Init a matrix with 0

F = z(1,1);
GL = z(2,1);
G = z(3,1);
L = z(4,1);

d_dt(1) = k.sig*GL - k.deg*F;
d_dt(2) = k.on*G*L - k.off*GL;
d_dt(3) = k.off*GL - k.on*G*L;

return;

function d_dt = DiffCalc4(t,z)

% Set parameters
k.sig = 0.033333333; %Rate of mCherry production dF/dt [sec-1]
k.deg = 0.002128988; %Rate of protein degredation [sec-1]

k.on = 0.25; %On rate of GPCR binding dGL/dt [uM-1 sec-1]
k.off = 0.00025; %Off rate of GPCR binding dGL/dt [sec-1]

d_dt = zeros(4,1); %Init a matrix with 0

F = z(1,1);
GL = z(2,1);
G = z(3,1);
L = z(4,1);

d_dt(1) = k.sig*GL - k.deg*F;
d_dt(2) = k.on*G*L - k.off*GL;
d_dt(3) = k.off*GL - k.on*G*L;

return;

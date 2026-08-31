%% Computational Activity 3
clear;
close all;
clc;

%% Part one: Constant pb and pinf
% ODEs defined at end (MATLAB syntax)

% set constant parameters
pinf = 101.325;
R0 = 1;
rho = 1.2;

%% pb < pinf
pb = 101;

% collect parameters into vector
params = [pb;pinf;R0;rho];

% set start and end times
tspan = [0,10];

% set initial condition as array
delta = 1e-6;
r0 = R0-delta;

[tp,rp] = ode45(@(t,x)bubble_plus(t,x,params),tspan,r0);
[tm,rm] = ode45(@(t,x)bubble_minus(t,x,params),tspan,r0);

figure(1); hold on; box on;
plot(tp,rp')
plot(tm,rm')
xlabel('t')
ylabel('R(t)')
legend('positive branch','negative branch');

%% pb = pinf
pb = 101.325;

% collect parameters into vector
params = [pb;pinf;R0;rho];

% set start and end times
tspan = [0,10];

% set initial condition as array
delta = 1e-6;
r0 = R0-delta;

[tp,rp] = ode45(@(t,x)bubble_plus(t,x,params),tspan,r0);
[tm,rm] = ode45(@(t,x)bubble_minus(t,x,params),tspan,r0);

figure(2); hold on; box on;
plot(tp,rp)
plot(tm,rm)
xlabel('t')
ylabel('R(t)')
legend('positive branch','negative branch');

%% pb > pinf
pb = 101.5;

% collect parameters into vector
params = [pb;pinf;R0;rho];

% set start and end times
tspan = [0,10];

% set initial condition as array
delta = 1e-6;
r0 = R0+delta;

[tp,rp] = ode45(@(t,x)bubble_plus(t,x,params),tspan,r0);
[tm,rm] = ode45(@(t,x)bubble_minus(t,x,params),tspan,r0);

figure(3); hold on; box on;
plot(tp,rp)
plot(tm,rm)
xlabel('t')
ylabel('R(t)')
legend('positive branch','negative branch');

%% Part two: Verifying true solution
% ODE defined at end (MATLAB syntax)

%% pb < pinf
pb = 101;

% collect parameters into vector
params = [pb;pinf;rho];

% set start and end times
tspan = [0,10];

% set initial condition as array
rv0 = [R0;0];

[t,rv] = ode45(@(t,x)bubble_system(t,x,params),tspan,rv0);

figure(4); hold on; box on;
plot(t,rv(:,1)')
xlabel('t')
ylabel('R(t)')

%% We see that the negative branch gives the correct solution for pb<pinf 

%% pb = pinf
pb = pinf;

% collect parameters into vector
params = [pb;pinf;rho];

% set start and end times
tspan = [0,10];

% set initial condition as array
rv0 = [R0;0];

[t,rv] = ode45(@(t,x)bubble_system(t,x,params),tspan,rv0);

figure(5); hold on; box on;
plot(t,rv(:,1)')
xlabel('t')
ylabel('R(t)')

%% We see that both branches gives the correct solution for pb=pinf

%% pb > pinf
pb = 101.5;

% collect parameters into vector
params = [pb;pinf;rho];

% set start and end times
tspan = [0,10];

% set initial condition as array
rv0 = [R0;0];

[t,rv] = ode45(@(t,x)bubble_system(t,x,params),tspan,rv0);

figure(6); hold on; box on;
plot(t,rv(:,1)')
xlabel('t')
ylabel('R(t)')

%% We see that the positive branch gives the correct solution for pb>pinf

%% Part three: Changing environmental pressure
% ODE defined at end (MATLAB sytax)

pb = 101;

cs = 0.2:0.005:0.25;
labels = [];

for ci = cs
    % collect parameters into vector
    params = [pb;pinf;rho;ci];
    
    % set start and end times
    tspan = [0,10];
    
    % set initial condition as array
    rv0 = [R0;0];
    
    [t,rv] = ode45(@(t,x)bubble_system_t(t,x,params),tspan,rv0);
    
    labels = [labels; "c = " + num2str(ci,3)];

    figure(7); hold on;
    plot(t,rv(:,1)')
end

figure(7); box on;
legend(labels,'Location','northwest');
xlabel('t')
ylabel('R(t)')

%% We see the critical rate is between 0.23 and 0.235

%% functions

function drdt = bubble_plus(t,r,params)
    pb = params(1);
    pinf = params(2);
    R0 = params(3);
    rho = params(4);

    drdt = sqrt(2.*(pb-pinf)./(3.*rho).*(1-(R0./r).^3));
    
    if imag(drdt) < 1e-3
        drdt = real(drdt);
    end

%     if ~isreal(drdt)
%         drdt = NaN;
%     end

end

function drdt = bubble_minus(t,r,params)
    pb = params(1);
    pinf = params(2);
    R0 = params(3);
    rho = params(4);

    drdt = -sqrt(2.*(pb-pinf)./(3.*rho).*(1-(R0./r).^3));
    
    if imag(drdt) < 1e-3
        drdt = real(drdt);
    end

%     if ~isreal(drdt)
%         drdt = NaN;
%     end

end

function drvdt = bubble_system(t,rv,params)
    pb = params(1);
    pinf = params(2);
    rho = params(3);
    
    r = rv(1);
    v = rv(2);

    drdt = v;
    dvdt = (1./r).*((pb-pinf)./rho-3./2.*v.^2);

    drvdt = [drdt;dvdt];
end

function drvdt = bubble_system_t(t,rv,params)
    pb = params(1);
    pinf = params(2);
    rho = params(3);
    c = params(4);
    
    r = rv(1);
    v = rv(2);

    drdt = v;
    dvdt = (1./r).*((pb-pinf+c.*t)./rho-3./2.*v.^2);

    drvdt = [drdt;dvdt];
end
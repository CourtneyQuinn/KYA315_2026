%% Computational Activity 4
clear;
close all;
clc;

%% Part one: Poiseuille and Couette Flow
% Visualising analytical solutions

% set parameters
G = 4;
mu = 0.2;
H = 2;
C = 100;

params = [G,mu,H,C];

x = linspace(0,1,1);
z = linspace(0,2,21);
[X,Z] = meshgrid(x,z);

U_P = Poiseuille(X,Z,params);
U_C = Couette(X,Z,params);

figure(1); hold on; box on;
quiver(X,Z,U_P,U_P*0,5)
title('Poiseuille flow')
xlabel('x')
ylabel('z')

figure(2); hold on; box on;
quiver(X,Z,U_C,U_C*0,5)
title(['Couette flow'])
xlabel('x')
ylabel('z')

%% Part two: Finite difference method
% look for solutions of form Mu=h^2f

f =@(z) -G/mu+zeros(size(z));
A = 0;
B = 0;

U_P_FD = finitediff(z,f,H,A,B);
U_C_FD = finitediff(z,f,H,A,C);

figure(3); hold on; box on;
quiver(X,Z,U_P_FD,U_P_FD*0,5)
title('Poiseuille flow (finite-diff)')
xlabel('x')
ylabel('z')

figure(4); hold on; box on;
quiver(X,Z,U_C_FD,U_C_FD*0,5)
title('Couette flow (finite-diff)')
xlabel('x')
ylabel('z')

% Compare approximation with analytical
error_P = abs(U_P(:,1)-U_P_FD(:,1));
error_C = abs(U_C(:,1)-U_C_FD(:,1));

figure(5); hold on; box on;
plot(z,error_P,'bo--')
plot(z,error_C,'ro--')
legend('Poiseuille','Couette')
xlabel('z')
ylabel('error')

%% Part three: Adding time dependence
% Approximate PDE as solution to forced ODE

% Pouiseuille boundary conditions
A = 0;
B = 0;

% set parameters
params = [G,mu,H,A,B];

% set start and end times
tspan = [0,10];
h = 1e-2;

% set initial condition vector
N = size(z,2);
u0 = zeros(N,1);
u0(1) = u0(1) + A;
u0(end) = u0(end) + B;

[u_eval,t_eval] = MyIVP(@(t,u)ODEsys(t,u,params),u0,tspan,h);

[t_mesh,z2_mesh] = meshgrid(t_eval,z);

figure(6); box on;
contourf(t_mesh,z2_mesh,u_eval,20);
xlabel('t')
ylabel('z')

%% Functions

function u = Poiseuille(x,z,params)
    G = params(1);
    mu = params(2);
    H = params(3);

    u = G/(2.*mu).*z.*(H-z);
end

function u = Couette(x,z,params)
    G = params(1);
    mu = params(2);
    H = params(3);
    C = params(4);

    u = C/H.*z+G/(2*mu).*z.*(H-z);
end

function u = finitediff(z,f,H,A,B)
    N = size(z,2);
    z_int = z(:,2:end-1);
    N_int = size(z_int,2);
    h = H/(N-1);
    
    M = diag(ones(N-3,1),-1) + diag(-2*ones(N-2,1),0) + diag(ones(N-3,1),1);

    b = h.^2.*f(z_int);
    b(1) = b(1)-A;
    b(end) = b(end)-B;
    u_int = M\b';

    u = NaN(size(z));
    u(2:end-1) = u_int;
    u(1) = A;
    u(end) = B;

    u = u';
end

function dudt = ODEsys(t,u,params)
    G = params(1);
    nu = params(2);
    H = params(3);
    A = params(4);
    B = params(5);

    N = size(u,1);
    u_int = u(2:end-1);
    N_int = size(u_int,1);
    h = H/(N-1);
    
    M = diag(ones(N_int-1,1),-1) + diag(-2*ones(N_int,1),0) + diag(ones(N_int-1,1),1);

    b = -G/nu + zeros(size(u_int)) + 10.*sin(t.*pi);
    b(1) = b(1)-A.*nu/h^2;
    b(end) = b(end)-B*nu/h^2;

    dudt = NaN(size(u));
    
    dudt_int = M.*nu/(h.^2)*u_int-b;
    dudt(2:end-1) = dudt_int;
    
    dudt(1) = 0;
    dudt(end) = 0;

end
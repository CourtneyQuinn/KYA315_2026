%% Computational Activity 2
clear;
close all;
clc;

%% Part one: Velocity potential function
% Potential function defined at end (MATLAB syntax)

% set discretization for grid
x = linspace(-2,2,101);
y = linspace(-2,2,101);
z = linspace(-2,2,101);

% calculate mesh grid
[X,Y,Z] = meshgrid(x,y,z);

% calculate for different values of c
Phi_05 = potential(X,Y,Z,0.5);
Phi_1 = potential(X,Y,Z,1);
Phi_2 = potential(X,Y,Z,2);

ind_z0 = 51;

f = figure(1);
f.Position(3:4) = [1500 400];

subplot(1,3,1);
contourf(X(:,:,ind_z0),Y(:,:,ind_z0),Phi_05(:,:,ind_z0),10,'Edgecolor','none');
title('c = 0.5')
xlabel('x')
ylabel('y')
colorbar()

subplot(1,3,2);
contourf(X(:,:,ind_z0),Y(:,:,ind_z0),Phi_1(:,:,ind_z0),10,'Edgecolor','none');
title('c = 1')
xlabel('x')
ylabel('y')
colorbar()

subplot(1,3,3);
contourf(X(:,:,ind_z0),Y(:,:,ind_z0),Phi_2(:,:,ind_z0),10,'Edgecolor','none');
title('c = 2')
xlabel('x')
ylabel('y')
colorbar()

%% mask out some of the centre singularities
for ii = 1:length(x)
    for jj = 1:length(y)
        for kk = 1:length(z)
            if ((x(ii)^2+y(jj)^2+z(kk)^2)<0.1)
                Phi_05(jj,ii,kk) = NaN;
                Phi_1(jj,ii,kk) = NaN;
                Phi_2(jj,ii,kk) = NaN;
            end
        end
    end
end

%% plot again with mask

f = figure(1);clf;
f.Position(3:4) = [1500 400];

subplot(1,3,1);
contour(X(:,:,ind_z0),Y(:,:,ind_z0),Phi_05(:,:,ind_z0),10,'Edgecolor','none');
title('c = 0.5')
xlabel('x')
ylabel('y')
colorbar()

subplot(1,3,2);
contour(X(:,:,ind_z0),Y(:,:,ind_z0),Phi_1(:,:,ind_z0),10)%,'Edgecolor','none');
title('c = 1')
xlabel('x')
ylabel('y')
colorbar()

subplot(1,3,3);
contour(X(:,:,ind_z0),Y(:,:,ind_z0),Phi_2(:,:,ind_z0),10,'Edgecolor','none');
title('c = 2')
xlabel('x')
ylabel('y')
colorbar()

%% Part two: Velocity and streamlines
% calculate velocity at each point in mesh
[U_c05,V_c05,W_c05] = velocity(X,Y,Z,0.5);
[U_c1,V_c1,W_c1] = velocity(X,Y,Z,1);
[U_c2,V_c2,W_c2] = velocity(X,Y,Z,2);


% plot streamlines over velocity potential
figure(1); hold on;

subplot(1,3,1);
p = streamslice(X(:,:,ind_z0),Y(:,:,ind_z0),U_c05(:,:,ind_z0),V_c05(:,:,ind_z0));
set(p,'Color','k');

subplot(1,3,2);
p = streamslice(X(:,:,ind_z0),Y(:,:,ind_z0),U_c1(:,:,ind_z0),V_c1(:,:,ind_z0));
set(p,'Color','k');


subplot(1,3,3);
p = streamslice(X(:,:,ind_z0),Y(:,:,ind_z0),U_c2(:,:,ind_z0),V_c2(:,:,ind_z0));
set(p,'Color','k');


%% Part three: Stagnation points
% Find the points at which velocity is 0.
stag1_c05 = fsolve(@(r)velocity_fsolve(r,0.5),[0.1;0.1;0.1]);
stag2_c05 = fsolve(@(r)velocity_fsolve(r,0.5),[-0.1;-0.1;-0.1]);

stag1_c1 = fsolve(@(r)velocity_fsolve(r,1),[0.1;0.1;0.1]);
stag2_c1 = fsolve(@(r)velocity_fsolve(r,1),[-0.1;-0.1;-0.1]);

stag1_c2 = fsolve(@(r)velocity_fsolve(r,2),[0.1;0.1;0.1]);
stag2_c2 = fsolve(@(r)velocity_fsolve(r,2),[-0.1;-0.1;-0.1]);

% plot stagnation points over velocity potential and streamlines
figure(1); hold on;

subplot(1,3,1); hold on;
plot(stag1_c05(1),stag1_c05(2),'ro','LineWidth',2);
plot(stag2_c05(1),stag2_c05(2),'ro','LineWidth',2);

subplot(1,3,2); hold on;
plot(stag1_c1(1),stag1_c1(2),'ro','LineWidth',2);
plot(stag2_c1(1),stag2_c1(2),'ro','LineWidth',2);

subplot(1,3,3); hold on;
plot(stag1_c2(1),stag1_c2(2),'ro','LineWidth',2);
plot(stag2_c2(1),stag2_c2(2),'ro','LineWidth',2);

%% Part four: Pressure 
% Calculate the pressure through Bernoulli case I
P_c05 = pressure(X,Y,Z,0.5);
P_c1 = pressure(X,Y,Z,1);
P_c2 = pressure(X,Y,Z,2);

%% mask out the body of the sphere
for ii = 1:length(x)
    for jj = 1:length(y)
        for kk = 1:length(z)
            if ((x(ii)^2+y(jj)^2+z(kk)^2)<1)
                P_c05(jj,ii,kk) = NaN;
                P_c1(jj,ii,kk) = NaN;
                P_c2(jj,ii,kk) = NaN;
            end
        end
    end
end

%% plot pressure as contour map

f2 = figure(2);clf;
f2.Position(3:4) = [1500 400];

subplot(1,3,1);
contourf(X(:,:,ind_z0),Y(:,:,ind_z0),P_c05(:,:,ind_z0),10,'Edgecolor','none');
title('c = 0.5')
xlabel('x')
ylabel('y')
colorbar()

subplot(1,3,2);
contourf(X(:,:,ind_z0),Y(:,:,ind_z0),P_c1(:,:,ind_z0),10,'Edgecolor','none');
title('c = 1')
xlabel('x')
ylabel('y')
colorbar()

subplot(1,3,3);
contourf(X(:,:,ind_z0),Y(:,:,ind_z0),P_c2(:,:,ind_z0),10,'Edgecolor','none');
title('c = 2')
xlabel('x')
ylabel('y')
colorbar()
%% functions

function Phi = potential(x,y,z,c)
    Phi = c.*x + (c.*x)./(2.*(x.^2+y.^2+z.^2).^(3/2));
end

function [u,v,w] = velocity(x,y,z,c)
    u = c + c./(2.*(x.^2+y.^2+z.^2).^(3/2))-3/2.*(c.*x.^2)./((x.^2+y.^2+z.^2).^(5/2));
    v = -3/2.*(c.*x.*y)./((x.^2+y.^2+z.^2).^(5/2));
    w = -3/2.*(c.*x.*z)./((x.^2+y.^2+z.^2).^(5/2));
end

function vel = velocity_fsolve(r,c)
    x = r(1);
    y = r(2);
    z = r(3);

    u = c + c./(2.*(x.^2+y.^2+z.^2).^(3/2))-3/2.*(c.*x.^2)./((x.^2+y.^2+z.^2).^(5/2));
    v = -3/2.*(c.*x.*y)./((x.^2+y.^2+z.^2).^(5/2));
    w = -3/2.*(c.*x.*z)./((x.^2+y.^2+z.^2).^(5/2));

    vel = [u;v;w];
end

function P = pressure(x,y,z,c)
    pinf = 101.325;
    rho = 1.2;
    
    [u,v,w] = velocity(x,y,z,c);

    P = pinf + 1./2.*rho.*(c.^2-(u.^2+v.^2+w.^2));
end

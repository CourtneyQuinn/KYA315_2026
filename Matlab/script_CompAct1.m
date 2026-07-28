%% Computational Activity 1
clear;
close all;
clc;

%% Part one: Velocity function
% velocity function defined at end (MATLAB syntax)

[u_test,v_test] = velocity(0,1);

fprintf('Velocity of [0,1]^T is \n\n')
disp([u_test;v_test])

%% Part two: Gridded domain
% set discretization for grid
x = linspace(-5,5,20);
y = linspace(-5,5,20);

% calculate paired points
[X,Y] = meshgrid(x,y);

% plot resulting mesh as individual points
figure;
plot(X,Y,'k.');

%% Part three: Streak lines and streamlines
% calculate velocity at each point in mesh
[U,V] = velocity(X,Y);

% plot streak lines
figure;
quiver(X,Y,U,V);

% plot streamlines
figure;
streamslice(X,Y,U,V)

%% Part four: Approximate acceleration
% finite difference function defined at end (MATLAB syntax)
[ax_approx, ay_approx]= accel_approx(X,Y,1e-6);

figure;
quiver(X,Y,ax_approx,ay_approx)

% function for analytical acceleration defined at end (MATLAB syntax)
[ax_exact, ay_exact]= accel_exact(X,Y);

figure;
quiver(X,Y,ax_exact,ay_exact)

% calculate the error
error_x = ax_exact-ax_approx;
error_y = ax_exact-ax_approx;

error = sqrt(error_x.^2+error_y.^2);

figure;
contourf(error)
colorbar()
%% functions

function [u,v] = velocity(x,y)
    u = -y./(x.^2+y.^2);
    v =  x./(x.^2+y.^2);
end

function [ax,ay] = accel_approx(x,y,h)
    xp = x + h;
    xm = x - h;
    
    [u_xp, v_xp] = velocity(xp,y);
    [u_xm, v_xm] = velocity(xm,y);

    dudx = (u_xp - u_xm) ./ (2*h);
    dvdx = (v_xp - v_xm) ./ (2*h);
    
    yp = y + h;
    ym = y - h;
    
    [u_yp, v_yp] = velocity(x,yp);
    [u_ym, v_ym] = velocity(x,ym);
    
    dudy = (u_yp - u_ym) ./ (2*h);
    dvdy = (v_yp - v_ym) ./ (2*h);

    [u,v] = velocity(x,y);

    ax = u.*dudx + v.*dudy;
    ay = u.*dvdx + v.*dvdy;
end

function [ax,ay] = accel_exact(x,y)
    ax = -x./(x.^2+y.^2).^2;
    ay = -y./(x.^2+y.^2).^2;
end
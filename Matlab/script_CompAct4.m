%% Computational Activity 4
clear;
close all;
clc;

%% Part one: The function
% Joukowski transformation function defined at end (MATLAB syntax)

% set body parameters
k = 1;
xi0 = -1;
eta0 = 1;

[x0,y0] = Joukowski(k,xi0,eta0);

figure(1); hold on; box on;
plot(x0,y0,'k')
title('Joukowski airfoil')
xlabel('x')
ylabel('y')


%% Part two: Calculate streamfunction
% function defined at the end (MATLAB syntax)

xi = linspace(-10,10,150);
eta = linspace(-10,10,150);
[Xi,Eta] = meshgrid(xi,eta);

c=1;

[Psi,z] = Streamfunction(k,xi0,eta0,c,Xi,Eta);

% calculate stagnation points
a=((k-xi0)^2+eta0^2)^(0.5);
tau1 = a+xi0+eta0*1i;
tau2 = -a+xi0+eta0*1i;

% convert back to z-plane using Joukowski map
zs1 = tau1 + k^2/tau1;
zs2 = tau2 + k^2/tau2;

figure(2); hold on; box on;
contourf(real(z),imag(z),Psi,10)
plot(real(zs1),imag(zs1),'r*')
plot(real(zs2),imag(zs2),'r*')
plot(x0,y0,'k')
xlabel('x')
ylabel('y')
title('Streamlines');

%% Part three: Adding circulation
% function defined at the end (MATLAB syntax)

Gamma = 4*pi*c*eta0;

[Psi,z] = Streamfunction_circ(k,xi0,eta0,c,Xi,Eta,Gamma);

% calculate stagnation points
a=((k-xi0)^2+eta0^2)^(0.5);
tau1 = -(1i.*Gamma./(2*pi)+(-Gamma^2./(4*pi^2)+4*c^2*a^2)^0.5)./(2*c)+xi0+eta0*1i;
tau2 = -(1i.*Gamma./(2*pi)-(-Gamma^2./(4*pi^2)+4*c^2*a^2)^0.5)./(2*c)+xi0+eta0*1i;

% convet back to z-plane
zs1 = tau1 + k^2/tau1;
zs2 = tau2 + k^2/tau2;

figure(3); hold on; box on;
contourf(real(z),imag(z),Psi,10)
plot(real(zs1),imag(zs1),'r*')
plot(real(zs2),imag(zs2),'r*')
plot(x0,y0,'k')
xlabel('x')
ylabel('y')
title('Streamlines with circulation');

%% Bonus: Modify angle of attack
% function defined at the end (MATLAB syntax)

Gamma = 8*pi*c*eta0;
alpha = pi/4;

[Psi,z] = Streamfunction_angle(k,xi0,eta0,c,Xi,Eta,Gamma,alpha);

% calculate stagnation points
a=((k-xi0)^2+eta0^2)^(0.5);
tau1 = -(1i.*Gamma./(2*pi)+(-Gamma^2./(4*pi^2)+4*c^2*a^2)^0.5)./(2*c.*exp(-1i*alpha))+xi0+eta0*1i;
tau2 = -(1i.*Gamma./(2*pi)-(-Gamma^2./(4*pi^2)+4*c^2*a^2)^0.5)./(2*c.*exp(-1i*alpha))+xi0+eta0*1i;

% convet back to z-plane
zs1 = tau1 + k^2/tau1;
zs2 = tau2 + k^2/tau2;

figure(4); clf; hold on; box on;
contourf(real(z),imag(z),Psi,10)
plot(real(zs1),imag(zs1),'r*')
plot(real(zs2),imag(zs2),'r*')
plot(x0,y0,'k')
xlabel('x')
ylabel('y')
title('Modified angle of attack');

%% functions

function [x,y] = Joukowski(k,xi0,eta0)
    % calculate radius
    r=sqrt((k-xi0)^2+eta0^2);
    
    % create theta array for 2pi rotation
    theta=linspace(0,2*pi,101);
    
    % convert to xi and eta values
    xi=r*cos(theta);
    eta=r*sin(theta);
    
    % shift to centre
    xi=xi+xi0;
    eta=eta+eta0;
    
    % convert back to x and y
    x=xi+(k^2).*xi./(xi.^2+eta.^2);
    y=eta-(k^2).*eta./(xi.^2+eta.^2);

end

function [Psi,z] = Streamfunction(k,xi0,eta0,c,xi,eta)
    tau0 = xi0+eta0.*1i;
    a2=(k-xi0).^2+eta0.^2;

    tau = xi+eta.*1i;
    
    w = c.*(tau+a2./(tau-tau0));

    % add mask for wing shape
    for ii = 1:size(xi,2)
        for jj = 1:size(eta,1)
            if ((xi(jj,ii)-xi0)^2+(eta(jj,ii)-eta0)^2)<a2
                w(jj,ii) = NaN+NaN*1i;
            end
        end
    end

    Psi = imag(w);

    % convert back to x and y
    x=xi+(k.^2).*xi./(xi.^2+eta.^2);
    y=eta-(k.^2)*eta./(xi.^2+eta.^2);
    z = x+y.*1i;

end

function [Psi,z] = Streamfunction_circ(k,xi0,eta0,c,xi,eta,Gamma)
    tau0 = xi0+eta0.*1i;
    a2=(k-xi0).^2+eta0.^2;

    tau = xi+eta.*1i;
    
    w = c.*(tau+a2./(tau-tau0)) + Gamma.*1i./(2.*pi).*log(tau-tau0);

    % add mask for wing shape
    for ii = 1:size(xi,2)
        for jj = 1:size(eta,1)
            if ((xi(jj,ii)-xi0)^2+(eta(jj,ii)-eta0)^2)<a2
                w(jj,ii) = NaN+NaN*1i;
            end
        end
    end

    Psi = imag(w);

    % convert back to x and y
    x=xi+(k.^2).*xi./(xi.^2+eta.^2);
    y=eta-(k.^2)*eta./(xi.^2+eta.^2);
    z = x+y.*1i;

end

function [Psi,z] = Streamfunction_angle(k,xi0,eta0,c,xi,eta,Gamma,alpha)
    tau0 = xi0+eta0.*1i;
    a2=(k-xi0).^2+eta0.^2;

    tau = xi+eta.*1i;
    
    w = c.*(tau.*exp(-1i*alpha)+a2./(exp(-1i*alpha).*(tau-tau0))) + Gamma.*1i./(2.*pi).*log(tau-tau0);

    % add mask for wing shape
    for ii = 1:size(xi,2)
        for jj = 1:size(eta,1)
            if ((xi(jj,ii)-xi0)^2+(eta(jj,ii)-eta0)^2)<a2
                w(jj,ii) = NaN+NaN*1i;
            end
        end
    end

    Psi = imag(w);

    % convert back to x and y
    x=xi+(k.^2).*xi./(xi.^2+eta.^2);
    y=eta-(k.^2)*eta./(xi.^2+eta.^2);
    z = x+y.*1i;

end
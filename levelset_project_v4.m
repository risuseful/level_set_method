% levelset_project_v4.m

% by Akmal Aulia

% the code uses UPWIND DIFFERENCING method

% clear session

clear all

clc

% vector fields

w1=-0.5;%sqrt(0.25/2);

w2=0;%sqrt(0.25/2);

% magnitude of gradient and it's coefficient

mag_grad=1; % magnitude of gradient, using Signed Distance Function (SDF)

a=0.5; % coefficient to mag_grad

% domains

xmin=0.0;

xmax=1.0;

ymin=0.0;

ymax=1.0;

% number of grids

Nx=50;

Ny=50;

% steplength for x and y axis

dx=(xmax-xmin)/(Nx-1);

dy=(ymax-ymin)/(Ny-1);

% creating x and y vector

x=xmin:dx:xmax;

y=ymin:dy:ymax;

% steplength for t

W=2*dx;

dx_prime=0.5; % from Sun and Beckermann's

dt=0.1*dx; % based on Courant Number dt'/dx'=0.1

% initializing time and center-of-circle coordinate

n=1; % time index

xc=1/2*(xmax-xmin);

yc=1/2*(ymax-ymin);

% initializing phi as a cone

for i=1:Nx

for j=1:Ny

% evaluate phi for the 1st time

phi(i,j,n)=sqrt( (x(i)-xc)^2 + (y(j)-yc)^2 );% - R^2; % cone

end

end

% smoothing the sharp center of the cone with filter

% note: see Scardovelli and Zaleski 1999, "Direct Numerical Simulation

% of Free Surface and Interfacial Flow, page

% tempPhi = smoothfilter( phi(:,:,n) );

% phi(:,:,n) = tempPhi(:,:);

%-------up-to-this-line-everything's-okay-!------------------------%

% set size for phi_x and phi_y (first derivative)

phi_dx=zeros(Nx,Ny);

phi_dy=zeros(Nx,Ny);

% make the extended matrix PHI to ease application of upwind diff.

PHI = zeros( Nx+2, Ny+2 );

for n=1:300 % just for testing

% solve for phi_x = d(phi)/dx excluding LEFT and RIGHT wall

% note: apply Central Difference

phi_dx(:,2:(Ny-1)) = (1/(2*dx))*( phi(:,(2+1):(Ny-1+1),n) - phi(:,(2-1):(Ny-1-1),n)

);


% solve for phi_y = d(phi)/dy excluding TOP and BOTTOM wall

% note: apply Central Difference

phi_dy(2:(Nx-1),:) = (1/(2*dx))*( phi((2+1):(Nx-1+1),:,n) - phi((2-1):(Ny-1-1),:,n)

);


% resolving phi_dy at boundary

%phi_dx(:,1) = (2*dx)\( -3*Phi(:,1) + 4*Phi(:,2) - Phi(:,3) ); % 2nd ord. forward

phi_dx(:,1) = ( phi(:,2,n) - phi(:,1,n) )/dx; % forward

phi_dx(:,Ny) = ( phi(:,Ny,n) - phi(:,Ny-1,n) )/dx; % backward


% resolving phi_dy at boundary

phi_dy(1,:) = ( phi(2,:,n) - phi(1,:,n) )/dx; % forward

phi_dy(Nx,:) = ( phi(Nx,:,n) - phi(Nx-1,:,n) )/dx; % backward


% solve for new phi

temp1 = w1*phi_dx + w2*phi_dy;

temp2 = a*mag_grad;

% solving Hamilton-Jacobi eq. with Euler method

phi(:,:,n+1) = phi(:,:,n) - dt*( temp1 + temp2 );



end

Nt=n;

%plotting

mesh(x,y,phi(:,:,Nt))

title('phi')

hold on

level=zeros(Nx,Ny)

mesh(x,y,level)

view(3) 

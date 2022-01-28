% linSearch.m

% by Akmal Aulia

% purpose: function perform steepest descent method

% on a mesh (not function)

% note : makesure to convert phi(i,j,n) into Phi(i,j), without "n"

function minima = linSearch( v0,Phi, x,y,dx,dy,Nx,Ny )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi=Phi; % phi is phi(i,j), not phi(i,j,n) anymore.

% steepest descent method parameters

alpha=0.01; % step length

Bk = eye(2);% define identity matrix Bk


k=1; %counter

while( k<55 ) % since alpha constant, we set 200 iteration for convergence

% find location of v0 in matrix phi

for i=1:Nx

for j=1:Ny

distR(j) = ( x(i)-v0(1) )^2 + ( y(j)-v0(2) )^2 ;

end

[distR_min_i(i), pos_j(i)] = min( distR );

end

[distR_val, pos_i]=min(distR_min_i);

ind_i = pos_i; % i-th position in matrix phi (row)

ind_j = pos_j(pos_i);% j-th position in matrix phi (column)

% record position

position(:,k) = [ind_i ind_j]; % associated position tracking

assoc_phi(k) = phi(ind_i,ind_j); % associate phi tracking

% evaluate gradient at all point

% solve for phi_x = d(phi)/dx excluding LEFT and RIGHT wall

% note: apply Central Difference

phi_dx(:,2:(Ny-1)) = (1/(2*dx))*( phi(:,(2+1):(Ny-1+1)) - phi(:,(2-1):(Ny-1-1)) );

% solve for phi_y = d(phi)/dy excluding TOP and BOTTOM wall

% note: apply Central Difference

phi_dy(2:(Nx-1),:) = (1/(2*dx))*( phi((2+1):(Nx-1+1),:) - phi((2-1):(Ny-1-1),:) );

% resolving phi_dy at boundary

%phi_dx(:,1) = (2*dx)\( -3*Phi(:,1) + 4*Phi(:,2) - Phi(:,3) ); % 2nd ord. forward

phi_dx(:,1) = ( phi(:,2) - phi(:,1) )/dx; % forward

phi_dx(:,Ny) = ( phi(:,Ny) - phi(:,Ny-1) )/dx; % backward

% resolving phi_dy at boundary

phi_dy(1,:) = ( phi(2,:) - phi(1,:) )/dx; % forward

phi_dy(Nx,:) = ( phi(Nx,:) - phi(Nx-1,:) )/dx; % backward

% evaluate gradient at point v0

grad_v0(1) = phi_dx(ind_i, ind_j);

grad_v0(2) = phi_dy(ind_i, ind_j);

% record v0 at k-th

v(:,k) = v0;

% define steepes descent direction

p(:,k) = -inv(Bk)*grad_v0';

% update new coord. "v"

v(:,k+1) = v(:,k) + alpha*p(:,k);

% define v0 as v(:,k+1) for next calculation

v0 = v(:,k+1);

% increment count k

k=k+1;

end

minima=v0

vt=v';

% plotting matrix

tracLin = [vt(1:(k-1),1) vt(1:(k-1),2) assoc_phi']

plot3( tracLin(:,1), tracLin(:,2), tracLin(:,3),'*k' )

% axis([ 0 1 0 1 0 1])

hold on

mesh(x,y,phi(:,:))

%%%%%%%%%%%%

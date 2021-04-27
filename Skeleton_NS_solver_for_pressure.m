clear all
close all
clc

% This syetm that you need to solve will be singular. Matlab gives you a
% warning at each time step. To switch of this warning, remove the comment
% in the next line

warning on

% 00D#MMXIX#

% This file contains the skeleton of the program which will solve the lid
% driven cavity problem on a unit square. The parts that have to be
% supplemented are described in the assignment.
%
% The pieces that need to be filled in are indicated
%

%
% When running the code, determine a suitable time step. A too small time
% step will make the calculation very long, while for a too large time step
% the solution will blow up due to numerical instability.
%

Re = 1000;              % Reynolds number
N = 3;                 % Number of volumes in the x- and y-direction
Delta = 1/N;            % uniform spacing to be used in the mapping to compute tx
tol = 1e-6;             % tol determines when steady state is reached and the program terminates


% wall velocities
U_wall_top = -1;
U_wall_bot = 0;
U_wall_left = 0;
U_wall_right = 0;
V_wall_top = 0;
V_wall_bot = 0;
V_wall_left = 0;
V_wall_right = 0;

%
%   Generation of a non-uniform mesh
%

%
%   tx are the coordinates of the nodal points on the outer-oriented grid
%
tx = zeros(1,N+1);
for i=1:N+1
    xi = (i-1)*Delta;
    tx(i) = 0.5*(1. - cos(pi*xi));
end

% Local mesh size on outer oriented grid
th = zeros(N,1);
th = tx(2:N+1) - tx(1:N);

%
%  x are the coordinates of the nodal points on the inner-oriented grid (including
%  endpoints 0 and 1)
%  h contains the edge lengths on the inner-oriented grid
%
x = 0.5*(tx(1:N) + tx(2:N+1));
x = [0 x 1];

h = zeros(N+1,1);
h = x(2:N+2) - x(1:N+1);

% Determine a suitable time step and stopping criterion, tol
h_min = min(h);
dt = min(h_min, 0.5*Re*h_min^2);             % time step



%
%   Initial condition u=v=0
%
%   Both u and v will be stored in one big vector called 'u'
%
%   The vector u only contains the true unknowns, not the velocities that
%   are prescribed by the boundary conditions
%
%   The vector u contains the *inner-oriented* circulations as unknowns

u = zeros(2*N*(N-1),1);

% Set up the Incidence matrix 'tE21' which connects the fluxes to the
% volumes. Use the orientation described in the assignment.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    This is where you need to insert code
%
tE21 = zeros(N*N , 2*N*(N+1)); % # rows = # of surfaces, # columns = # of edges

h_strip = 0; % will use to add to the indices depending on which horizontal strip we're at
h_count = 1; % keep track of the horizontal strip

for i=1:N*N % loop over surfaces (row of the matrix)
    % do u
    tE21(i,i+h_strip) = -1;
    tE21(i,i+1+h_strip) = 1;

    % do v
    tE21(i,N*(N+1) + i) = -1;
    tE21(i,N*(N+1) + i+N) = 1;

    % horizontal strip shenanigans (see whether or not we've reached the end of a horizontal strip)
    if (h_count==N)
        h_count = 1;
        h_strip = h_strip+1;
    else
        h_count = h_count+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%  Inserting boundary conditions for normal velocity components
%  Multiplication by h components is done in order to obtain integral flux
%  values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Here you need to modify the incidence matrix tE21 to include the
%    boundary conditions and store the boundary terms in u_norm (see
%    assignment)
%
bc_left = 1:N+1:N*(N+1)-N; % indices of the column in tE21 of the left boundary
bc_right = 1+N:N+1:N*(N+1); % indices of the column in tE21 of the right boundary
bc_bottom = N*(N+1)+1 : N*(N+1)+N; % indices of the column in tE21 of the bottom boundary
bc_top = N*(N+1)+1+N*N : N*(N+1)+N*N+N; % indices of the column in tE21 of the top boundary

bc_col = [bc_left bc_right bc_bottom bc_top]; % Store all BC indices into a sigle row vector

tE21(:,bc_col) = []; % Remove the columns in tE21 that correspond to the boundaries

tE21 = sparse(tE21); %% Store in sparse matrix

u_norm = zeros(N*N, 1); %% In the problem the fluxes are zero at all boundaries. Note dimensions is the same as the number of surfaces on the primal mesh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting up simple Hodge matrix which converts the fluxes on the
% outer-oriented grid to circulation on the inner-oriented grid. Assume
% that the fluxes and circulation are constant over each 1-cell. This will
% give a diagonal Hodge matrix. Call this Hodge matrix 'H1t1'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Here you need to construct the incidence matrix H1t1
%
H1t1 = zeros(2*N*(N+1));
Nd = N+1;

count = 1;

% do u
for i=1:N
    for j=1:Nd
        H1t1(count,count) = h(j)/th(i);
        count = count+1;
    end
end

% do v
for i=1:Nd
    for j=1:N
        H1t1(count,count) = h(i)/th(j);
        count = count+1;
    end
end


% Store in Sparse matrix
H1t1 = sparse(H1t1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hu_norm is the vector which will contain the Hodge of the prescribed
% normal fluxes. Calculate the vector 'Hu_norm'

%Vector ubc_norm containing fluxes of prescribed velocities is reused

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Multiply H1t1 with u_norm to get Hu_norm
%
% hu_norm_matrix = zeros(length(bc_col), length(bc_col)); % matrix that is multiplied into the prescribed velocities
% hu_norm_matrix = sparse(hu_norm_matrix);
% % hunm_size = size(hu_norm_matrix);
% %
% % shift = 0;
% % for i=1:length(bc_col)
% %     hu_norm_matrix(i,i) = H1t1(bc_col(i)-shift,bc_col(i)-shift); %% Put into hu_norm_matrix
% %     H1t1(:,bc_col(i)-shift) = []; % remove that column
% %     H1t1(bc_col(i)-shift,:) = []; % remove that row
% %     shift = shift+1;
% % end
%
% Hu_norm = hu_norm_matrix*u_norm; % In the problem we know that u_norm is vector of zeroes anyway
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  Remove corresponding row and columns from the Hodge matrix and also
%  remove the corresp0onding 'rows' from Hu_norm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Remove rows and colums in H1t1 and Hu_norm for prescribed values
%
%
H1t1(:,bc_col) = [];
H1t1(bc_col,:) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the incidence E^{2,1} between 1-cochain circulation and 2-cochain vorticity on
% the inner-oriented (extended) grid
%
% This incidence matrix will be called 'E21' in the program
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Set up the incidence matrix E21 on the dual grid
%
Nd = N+1;
E21 = zeros(Nd*Nd , 2*Nd*(Nd+1)); % # rows = # of surfaces, # columns = # of edges

h_strip = 0;
h_count = 1;

for i=1:Nd*Nd
    % do u
    E21(i, i) = 1;
    E21(i,i+Nd) = -1;
    % do v
    E21(i,Nd*(Nd+1) + i +h_strip) = -1;
    E21(i,Nd*(Nd+1) + i+1 + h_strip) = 1;
    if (h_count==Nd)
        h_strip = h_strip+1;
        h_count = 1;
    else
        h_count = h_count +1;
    end
end

E21 = sparse(E21);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Inserting prescribed tangential bundary conditions

%Vector ubc of size 2(N+2)(N+1) is constructed. It contains all the
%prescribed velocities (normal and tangential) as represented in the inner
%oriented grid.

% Remove columns from the incidence matrix E21 corresponding to both the
% prescribed tangental velocities and normal velocities


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Remove prescribed values from the matrix (just as was done for tE21)
%    and store the known values in a vector u_pres
%

% Indices of each wall (Tangential BC)
bc_left = Nd*(Nd+1)+1 : Nd+1 : 2*Nd*(Nd+1)-Nd  ;
bc_right = Nd*(Nd+1)+1+Nd : Nd+1 : 2*Nd*(Nd+1) ;
bc_top = Nd*(Nd+1)-N : Nd*(Nd+1);
bc_bottom = 1 : Nd;

bc_col = [bc_left bc_right bc_top bc_bottom];
bc_col = sort(bc_col);

% Indices of each wall (Normal BC)
bcn_left = 1+Nd:Nd:Nd*(Nd+1)-N-Nd;
bcn_right = Nd+Nd:Nd:Nd*(Nd+1) - Nd;
bcn_top = 2*Nd*(Nd+1)-Nd+1:2*Nd*(Nd+1)-1;
bcn_bottom = Nd*(Nd+1)+2:Nd*(Nd+1)+Nd;

bcf = [bc_left bc_right bc_top bc_bottom bcn_left bcn_right bcn_top bcn_bottom]; % Normal and tangential parts of the boundary

E21(:,bcf) = []; % take both normal and tangential part

% Get u_pres

count = 1;
u_pres = zeros(Nd^2,1); % number of rows of u_pres corresponds to number of surfaces in dual mesh
for i =Nd^2-N:Nd^2
    u_pres(i) = -U_wall_top*h(count); % known circulation on the surfaces closest to the top wall
    count = count+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set up the Hodge matrix which maps inner-oriented 2-cochains to
% outer-oriented 0-cochains. Assume that the vorticity is constant in the
% inner-oriented 2-cells. This will give a diagonal Hodge matrix. Call this
% Hodhe matrix 'Ht02'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Set up the Hodge matrix which converts integrated values on the dual
%    grid to point values on the primal grid assuming that the physical
%    quantity is constant over the dual surfaces.
%
Ht02 = zeros(Nd);
h_strip = 1;
count = 1;

for i=1:Nd
    for j=1:Nd
        Ht02(count,count) = 1/(h(h_strip)*h(j));
        count = count+1;
    end
    h_strip = h_strip+1;
end
Ht02 = sparse(Ht02);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the Hodge matrix which maps inner-oriented 1-cochains to
% outer-oriented 1-cochains. Call this Hodge matrix 'Ht11'. Assume again
% that everything is constant over the 1-cells, which will then yield a
% diagonal Hdoge matrix.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Set up the Hodge matrix which converts integrated values along dual
%    edge to integral values over primal edges
%
Ht11 = zeros(2*N*(N+1));

count = 1;

% do u
for i=1:N
    for j=1:Nd
        Ht11(count,count) = th(i)/h(j);
        count = count+1;
    end
end

% do v
for i=1:Nd
    for j=1:N
        Ht11(count,count) = th(j)/h(i);
        count = count+1;
    end
end
% Removing columns and rows that correspond to the boundaries (normal comp.)

bc_left = 1:N+1:N*(N+1)-N; % indices of the column in tE21 of the left boundary
bc_right = 1+N:N+1:N*(N+1); % indices of the column in tE21 of the right boundary
bc_bottom = N*(N+1)+1 : N*(N+1)+N; % indices of the column in tE21 of the bottom boundary
bc_top = N*(N+1)+1+N*N : N*(N+1)+N*N+N; % indices of the column in tE21 of the top boundary

bc_col = [bc_left bc_right bc_bottom bc_top]; % Store all BC indices into a sigle row vector
Ht11(:,bc_col) = [];
Ht11(bc_col,:) = [];

Ht11 = sparse(Ht11);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% The prescribed velocties will play a role in the momentum equation
%
u_pres_tvort=Ht02*u_pres; %U_pres to outer oriented 0 form representing contribution of boundary conditions to point wise vorticity
u_pres = H1t1*E21'*Ht02*u_pres; %U_pres to inner oriented 1 forms


% Now all matrices are set up and the time stepping cam start. 'iter' will
% record the number of time steps. This allows you to give output after a
% preselected number of time steps.
%
% 'diff' will be the maximal du/dt or dv/dt. If 'diff' is sufficiently
% small, steady state has been reached. Determine a suitable value for
% 'tol'
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    From here the code enters the time-stepping loop and no new parts in
%    the code need to be inserted. If done correctly, everything should
%    work now.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diff = 1;
iter = 1;


% Set up the matrix for the Poisson equation

A = -tE21*Ht11*tE21';

% Perform an LU-decomposition for the pressure matrix A

[L,U] = lu(A);

% Abbreviation for some matrix products which are constant in the time loop

VLaplace = H1t1*E21'*Ht02*E21;
DIV = tE21*Ht11;

while diff > tol

    %Vector chi is obtained. It corresponds with the point-wise vorticity
    %at each cell
    chi=Ht02*E21*u+u_pres_tvort;

    %Vectors uxchi and uychi correspond with the multiplications of
    %chi with the horizontal and vertical velocity components at each cell.
    %Only the cells required for vector convective are calculated. The
    %ordering of each vector with respect to the ordering of cells in the
    %grid is different (left to right for uxchi and bottom to top for
    %uychi)

    uxchi=zeros((N+1)*(N-1),1);
    uychi=zeros((N+1)*(N-1),1);

    for i=1:N-1
        for j=1:N+1
            k=j+(i-1)*(N+1); %Number of vector component
            if j==1
                uxchi(k)=U_wall_left*chi(j+i*(N+1));
                uychi(k)=V_wall_bot*chi((i+1)+(j-1)*(N+1));
            elseif j==N+1
                uxchi(k)=U_wall_right*chi(j+i*(N+1));
                uychi(k)=V_wall_top*chi((i+1)+(j-1)*(N+1));
            else
                uxchi(k)=(u(j-1+(i-1)*(N-1))+u(j-1+i*(N-1)))/(2*h(j))*chi(j+i*(N+1));
                uychi(k)=(u(N*(N-1)+i+(j-2)*N)+u(N*(N-1)+i+1+(j-2)*N))/(2*h(j))*chi((i+1)+(j-1)*(N+1));
            end
        end
    end

    %Given vectors uxchi and uychi, vector convective can be constructed
    convective=zeros(2*N*(N-1),1);
    for i=1:N
        for j=1:N-1
            convective(j+(i-1)*(N-1))=-h(j+1)/2*(uychi(i+(j-1)*(N+1))+uychi(i+1+(j-1)*(N+1))); %Components along horizontal cell lines
            convective(N*(N-1)+i+(j-1)*N)=h(j+1)/2*(uxchi(i+(j-1)*(N+1))+uxchi(i+1+(j-1)*(N+1))); %Components along vertical cell lines
        end
    end

    % Set up the right hand side for the Poisson equation for the pressure

    rhs_Poisson  =   DIV*(u/dt  - convective - VLaplace*u/Re - u_pres/Re) + u_norm/dt;

    % Solve for the new pressure

    temp = L\rhs_Poisson;
    p = U\temp;

    % Store the velocity from the previous time step in the vector u_old

    uold = u;

    % Update the velocity field

    u = u - dt* (convective - tE21'*p + VLaplace*u/Re + u_pres/Re);

    %
    %  Every other 1000 iterations check whether you approach steady state
    %  and check whether yopu satisfy conservation of mass. The largest
    %  rate at which mass is destroyed or created is denoted by 'maxdiv'.
    %  This number should be very small, in the order of machine precision.

    if mod(iter,1000) == 0

        maxdiv = max(DIV*u + u_norm)

        diff = max(abs(u-uold))/dt

    end
    iter = iter + 1;
end

% Checking the row sum of A
dimA = size(A);
rsA = zeros(dimA(1),1); % Initiate rowsum of A

for i=1:dimA(1)
    rsA(i) = sum(A(i,:));
end

%% Plotting vorticity: chi. According to the description/definition in line 436-438, these are assigned to the points on the primal mesh.

[Xp,Yp] = meshgrid(tx,tx); % coordinates of the points in the primal mesh


CHI = zeros(N+1,N+1); % vector to re-arrange chi so that contour plot can be made

count = 1;
for i=1:N+1
    for j=1:N+1
        CHI(i,j) = chi(count);
        count = count+1;
    end
end

level_chi = [-3.0 -2.0 -1.0 -0.5 0.0 0.5 1.0 2.0 3.0 4.0 5.0]; % Vorticity levels in Botella and Peyret (figure 2, table 8)
figure('Name','Iso-Vorticity')
% contour(Xp,Yp,CHI,level_chi, 'ShowText','on')
contour(Xp,Yp,CHI,level_chi)

%% Plotting the streamlines/streamfunction. these should actually correspond to the fluxes
% (dpsi = u dy, -v dx)
u = Ht11*u; % on primal mesh

% mid point of tx
txm = zeros(1,N);
for i=1:N
    txm(i) = 0.5*(tx(i)+tx(i+1));
end

% Assign values of stream functions to each cell (so the average of the edges)
[X,Y] = meshgrid(txm, txm);
psi = zeros(N,N);

col = N*(N-1); % some dummy to short hand N(N-1)

% % Get the streamfunction values for the surfaces on the bottom strip

count = 1;
i = 1;
for j=1:N
    if (j==1)
        psi(i,j) = 0.25*(u(count) - u(count + col)); % right, top
    elseif (j==N)
        psi(i,j) = 0.25*(u(count - 1) - sum(u(count+col-j+1:count + col))); % left, top
    else
        psi(i,j) = 0.25*(u(count) + u(count - 1) - sum(u(count+col-j+1:count + col))); % right, left, top
    end
    count = count+1;
end

% Next strips besides the top one

for i=2:N-1
    for j=1:N
        if (j==1)
            psi(i,j) = 0.25*(sum(u(j:N-1:count-(i-1))) - sum(u(count+col-j+1:count+col)) - sum(u(count+col-N-j+1:count+col-N))); %right, top, bottom
        elseif (j==N)
            psi(i,j) = 0.25*(sum(u(j-1:N-1:count-(i-1)-1)) - sum(u(count+col-j+1:count+col)) - sum(u(count+col-N-j+1:count+col-N))); % left, top, bottom
        else
            psi(i,j) = 0.25*(sum(u(j:N-1:count-(i-1))) + sum(u(j-1:N-1:count-(i-1)-1)) - sum(u(count+col-j+1:count+col)) - sum(u(count+col-N-j+1:count+col-N))); % right, left, top, bottom
        end
        count = count+1;
    end
end

% Last strip on the very top
i = N;
for j=1:N
    if (j==1)
        psi(i,j) = 0.25*(sum(u(j:N-1:count-(i-1))) - sum(u(count+col-N-j+1:count+col-N))); % right, bottom
    elseif (j==N)
        psi(i,j) = 0.25*(sum(u(j-1:N-1:count-(i-1)-1)) - sum(u(count+col-N-j+1:count+col-N))); % left, bottom
    else
        psi(i,j) = 0.25*(sum(u(j:N-1:count-(i-1))) + sum(u(j-1:N-1:count-(i-1)-1)) - sum(u(count+col-N-j+1:count+col-N)));% right, left, bottom
    end
    count = count+1;
end

% levels from table 7 in Botella and Peyret
level_psi = [-1.5e-3 -1e-3 -5e-4 -2.5e-4 -1e-4 -5e-5 -1e-5 -1e-6 0.0 1e-10 1e-5 1e-4 1e-2 3e-2 5e-2 7e-2 9e-2 0.1 0.11 0.115 0.1175];
figure('Name','Streamfunction')
contour(X,Y,psi, level_psi)


u = H1t1*u; % back to dual mesh

%% Plotting the pressure. This should be assigned to the points in the dual mesh (excluding the boundaries)
level_p = [-0.002 0.0 0.02 0.05 0.07 0.09 0.11 0.12 0.17 0.3]; % pressure contour levels from Botella and Peyret (table 8)


[Xd, Yd] = meshgrid(x(2:length(x)-1), x(2:length(x)-1));
p_plot = zeros(N, N); %% Pressure just to remake p so that contour plot can be made
col = N*(N-1);

count = 1;
for i=1:N
    for j=1:N
        p_plot(i,j) = p(count); % stat from total pressure
        count = count+1;
    end
end

Ek = zeros(N,N) ;  % Kinematic energy corresponding to the point in the pressure 0.5 vel**2

u = Ht11*u; % Go to primal mesh

count = 1;
i = 1;
for j=1:N
    if (j==1)
        Ek(i,j) = 0.5*( ((u(count-(i-1)))/(2*th(i)))^(2) + ((u(count+col))/(2*th(j)))^(2) ); % right, top
    elseif (j==N)
        Ek(i,j) = 0.5*( ((u(count-(i-1)-1))/(2*th(i)))^(2) + ((u(count+col))/(2*th(j)))^(2) ); % left, top
    else
        Ek(i,j) = 0.5*( ((u(count-(i-1))+u(count-(i-1)-1))/(2*th(i)))^(2) + ((u(count+col))/(2*th(j)))^(2) ); % right, left, top
    end
    count = count+1;
end

% Next strips besides the top one

for i=2:N-1
    for j=1:N
        if (j==1)
            Ek(i,j) = 0.5*( ((u(count-(i-1)))/(2*th(i)))^(2) + ((u(count+col)+u(count+col-N))/(2*th(j)))^(2) ); %right, top, bottom
        elseif (j==N)
            Ek(i,j) = 0.5*( ((u(count-(i-1)-1))/(2*th(i)))^(2) + ((u(count+col)+u(count+col-N))/(2*th(j)))^(2) ); % left, top, bottom
        else
            Ek(i,j) = 0.5*( ((u(count-(i-1))+u(count-(i-1)-1))/(2*th(i)))^(2) + ((u(count+col)+u(count+col-N))/(2*th(j)))^(2) ); % right, left, top, bottom
        end
        count = count+1;
    end
end

% Last strip on the very top
i = N;
for j=1:N
    if (j==1)
        Ek(i,j) = 0.5*( ((u(count-(i-1)))/(2*th(i)))^(2) + ((u(count+col-N))/(2*th(j)))^(2) ); % right, bottom
    elseif (j==N)
        Ek(i,j) = 0.5*( ((u(count-(i-1)-1))/(2*th(i)))^(2) + ((u(count+col-N))/(2*th(j)))^(2) ); % left, bottom
    else
        Ek(i,j) = 0.5*( ((u(count-(i-1))+u(count-(i-1)-1))/(2*th(i)))^(2) + ((u(count+col-N))/(2*th(j)))^(2) ); % right, left, bottom
    end
    count = count+1;
end

p_plot = p_plot-Ek; % Static pressure

xr = x(x>0.5);  %% right half of x
c_approx = size(x)-size(xr); %% approximately the index of the centrelines/centre
c_approx = c_approx(2);
pc1 = p_plot(c_approx, c_approx); %% the pressure to the left of the centre
pc2 = p_plot(c_approx+1, c_approx+1); %% the pressure to the right of the centre
pc = pc1 + (0.5 - x(c_approx))*(pc2 - pc1)/(x(c_approx+1) - x(c_approx)); %% linear interpolation for pressure at the centre through a diagonal line

p_plot = p_plot - pc; %% subtract such that the pressure at the centre is zero as in Botella and Peyret


level_p = [-0.002 0.0 0.02 0.05 0.07 0.09 0.11 0.12 0.17 0.3]; % pressure contour levels from Botella and Peyret (table 8)
figure('Name','Isobars')
contour(Xd,Yd,p_plot,level_p);

x_ref = [0 0.0312 0.0391 0.0469 0.0547 0.0937 0.1406 0.1953 0.5 0.7656 0.7734 0.8437 0.9062 0.9219 0.9297 0.9375 1];
y_ref = [0 0.0547 0.0625 0.0703 0.1016 0.1719 0.2813 0.4531 0.5 0.6172 0.7344 0.8516 0.9531 0.9609 0.9688 0.9766 1];

p_ver = [0.110591 0.109689 0.1092 0.108566 0.104187 0.081925 0.040377 0.004434 0 -0.000827 0.012122 0.034910 0.050329 0.050949 0.051514 0.052009 0.052987];
p_hor = [0.090477 0.088445 0.087653 0.086716 0.084386 0.069511 0.04726 0.044848 0 0.034552 0.049029 0.065816 0.077154 0.078148 0.0078685 0.078837 0.077455];
% p_hor should be flipped!

p_vc =  p_plot(:,c_approx) + (0.5 - x(c_approx))*(p_plot(:,c_approx+1)-p_plot(:,c_approx))/(x(c_approx+1) - x(c_approx)); %% pressure along the vertical centreline by linear interpolation
p_hc =  p_plot(c_approx,:) + (0.5 - x(c_approx))*(p_plot(c_approx+1,:)-p_plot(c_approx,:))/(x(c_approx+1) - x(c_approx)); %% pressure along the horizontal centreline by linear interpolation

figure('Name','Horizontal centerline pressure')%Centerline pressure
hold on
plot(x(2:length(x)-1), p_hc, 'r');
plot(x_ref, fliplr(p_hor),'b');
legend(['Numerical values'],['Botella and Peyret'] , 'location', 'best');

figure('Name','Vertical centerline pressure')%Centerline pressure
hold on
plot(x(2:length(x)-1), p_vc, 'r');
plot(y_ref, p_ver,'b');
legend(['Numerical values'],['Botella and Peyret'], 'location', 'best');



u = H1t1*u; % Back to dual mesh

%% Vertical Centreline u

y_ref = [0 0.0547 0.0625 0.0703 0.1016 0.1719 0.2813 0.4531 0.5 0.6172 0.7344 0.8516 0.9531 0.9609 0.9688 0.9766 1];
u_ref = [0 0.1812881 0.2023300 0.2228955 0.3004561 0.3885691 0.2803696 0.1081999 0.0620561 -0.0570178 -0.1886747 -0.3372212 -0.4723329 -0.5169277 -0.5808359 -0.6644227 -1];

u_vc = u(c_approx-1:N-1:Nd*N)/(h(c_approx));
u_vc = u_vc(1:length(u_vc)-2);
figure('Name','Vertical centerline x-velocity')
hold on
plot(x(2:length(x)-1),u_vc, 'r');
plot(y_ref, u_ref, 'b');
legend(['Numerical values'],['Botella and Peyret'], 'location', 'best');


%% Centerline flow field: u, v, p, chi
centre_index = find(x - 0.5 == min(abs(x - 0.5)));
u_reshape = reshape(u,[N*(N-1),2]); %Size of u = (2*N*(N-1),1), reshape it to get x,y components
p_reshape = reshape(p, [N,N]); %Size of p = (N*N,1), reshape it to get x,y components
u_x = u_reshape(:,1); %x-velocity component
u_y = u_reshape(:,2); %y-velocity component
p_x = p_reshape(centre_index,:); %x-pressure component for centre index
p_y = p_reshape(:,centre_index); %y-pressure component for centre index
%CHI to be reshaped


%Table 9,10 of Botella and Peyret
u_ref = [0 0.1812881 0.2023300 0.2228955 0.3004561 0.3885691 0.2803696 0.1081999 0.0620561 -0.0570178 -0.1886747 -0.3372212 -0.4723329 -0.5169277 -0.5808359 -0.6644227 -1];
v_ref = [0 -0.2279225 -0.2936869 -0.3553213 -0.4103754 -0.5264392 -0.4264545 -0.3202137 0.0257995 0.3253592 0.3339924 0.3769189 0.3330442 0.3099097 0.2962703 0.2807056 0];
x_ref = [0 0.0312 0.0391 0.0469 0.0547 0.0937 0.1406 0.1953 0.5 0.7656 0.7734 0.8437 0.9062 0.9219 0.9297 0.9375 1];
y_ref = [0 0.0547 0.0625 0.0703 0.1016 0.1719 0.2813 0.4531 0.5 0.6172 0.7344 0.8516 0.9531 0.9609 0.9688 0.9766 1];
p_ver = [0.110591 0.109689 0.1092 0.108566 0.104187 0.081925 0.040377 0.004434 0 -0.000827 0.012122 0.034910 0.050329 0.050949 0.051514 0.052009 0.052987];
p_hor = [0.090477 0.088445 0.087653 0.086716 0.084386 0.069511 0.04726 0.044848 0 0.034552 0.049029 0.065816 0.077154 0.078148 0.0078685 0.078837 0.077455];
omega_hor = [-5.46217 -8.4435 -8.24616 -7.58524 -6.50867 0.92291 3.43016 2.21171 2.06722 2.06122 2.00174 0.74207 -0.82398 -1.23991 -1.50306 -1.83308 -7.66369];
omega_ver = [-4.16648 -2.4496 -2.31786 -2.20175 -1.63436 1.05467 2.26772 2.06215 2.06722 2.06539 2.09121 1.762 4.85754 6.95968 9.49496 12.067 14.7534];


figure('Name','Horizontal centerline v')%Centerline v
hold on
plot(X(centre_index,:)',u_y(centre_index*N:(centre_index*N)+N-1,1)/h(1,centre_index), 'r');
% figure('Name','Horizontal centerline y-velocity - Botella and Peyret')
plot(x_ref, v_ref,'b');
legend(['Numerical values'],['Botella and Peyret'], 'location', 'best');

u_centre = zeros(N,1);
u_centre(1:N-1) = u(centre_index:N:N*(N-1)-(N-1-centre_index))/h(1,centre_index);
u_centre(N) = -1; % since the last value at the wall was prescribed



figure('Name','Horizontal centerline vorticity')%Horizontal centerline vorticity
hold on %centre_index:N:(N^2)-(N-centre_index)
plot(X(centre_index,:)',CHI(centre_index,1:end-1),'r');
plot(x_ref, omega_hor,'b');
legend(['Numerical values'],['Botella and Peyret'], 'location', 'best');

figure('Name','Vertical centerline vorticity')%Vertical centerline vorticity
hold on
plot(X(centre_index,:)',CHI(1:end-1,centre_index),'r');
plot(y_ref, omega_ver,'b');
legend(['Numerical values'],['Botella and Peyret'], 'location', 'best');

%% Integrated vorticity
vort_integ = 0;

for i=1:Nd
    for j = 1:Nd
        vort_integ = sum(CHI(i,j));%*1/(Ht02(i,j)));
    end
end

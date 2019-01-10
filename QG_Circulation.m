% "Using a simple numerical QG model for wind-driven ocean circulation,
% find the streamfunction when the shape of the basin is varied and several
% islands are added at various locations."

% % % Written by Luke Kachelein as the final project for SIOC 209 (Fall 2017)

% full equation (w/o simplifications) is equation 2.3.2 in Pedlosky.

% I split it up into 2 equations:
% zeta = laplacian of psi
% d zeta/dt = ... with appropriate substitution.
% 
% The problem is now a) inverse solution to Poisson's equation (hard) and a
% time stepping scheme (relatively easy assuming the first step is correct enough).

%% Build domain

clear all

%%% Build in-line (not recommended)
% I = 50; % # of x points
% J = 40; % # of y points
% Domain = ones(J+2,I+2);
%     Domain(1,:) = 0; % define shoreline (will be more important in more
%                        % complicated shorelines later).
%     Domain(end,:) = 0;
%     Domain(:,1) = 0;
%     Domain(:,end) = 0;
    
%%% continental divide
%     Domain(:,30:35) = 0;
%%%%%%%%%%%%%%%
%%% slanted continent
%     for k = 1:min(1+[I J])
%         Domain(k,k:end) = 0;
%     end
%%%%%%%%%%%%%%%
%%% "Hawaii"
%     BigIslandX = floor(I/2);
%     BigIslandY = floor(J/2);
%     Domain(BigIslandY:(BigIslandY+2),BigIslandX:(BigIslandX+2)) = 0;
%     Domain((BigIslandY-5):(BigIslandY-4),(BigIslandX-5):(BigIslandX-4)) = 0;
%     Domain((BigIslandY-8):(BigIslandY-7),(BigIslandX-8):(BigIslandX-7)) = 0;
%%%%%%%%%%%%%%%

%%% USER GENERATED (preferred):
basinDir = '/Users/lukekachelein/Documents/MATLAB/qg_basins';
% ^ Directory in which basin/land configurations are saved as PNG files.
%   To be changed by user.
disp(['Open up Paintbrush or similar program, make a PNG of the right size',...
    ' (black = land, white = water please), place it in ', ...
    basinDir,'/ ',...
    'and load it up below:'])
cd(basinDir)
userBasin = uigetfile('*.png');
Domain = imread(userBasin);
Domain = squeeze(Domain(:,:,1));
Domain = double(Domain./Domain);
%%%%%%%%%%%%%%%

I = size(Domain,2)-2;
J = size(Domain,1)-2;

DomainNoEdge = Domain(2:(end-1),2:(end-1));
DomainNoEdge_vec = reshape(DomainNoEdge,numel(DomainNoEdge),1);

nanDomain = Domain; % Just a land mask for plotting
    nanDomain(Domain==0) = NaN;
    nanDomain(Domain==1) = 0;

imagesc(Domain)

%% Establish parameters and initial conditions

% Grid size
d = 111000; % 111km ~ 1 degree

% Om = (7.2921150 ± 0.0000001) ×10^-5 radians per SI second
Om = (7.2921150)*10^-5;

f0 = 2*Om*sind(45); % average (or typical?) Coriolis parameter
WE = 10^-6; % m/s characteristic Ekman vertical velocity (p.145 in Pedlosky: 10^-4 cm/s)
    Ly = d*(J+2-1); % length of basin in y-direction
    y_min = 0*Ly/2; % I guess m from equator
    y_max = y_min + Ly;
    y = y_min:d:y_max; y = y';
% "wind" example in equation 2.2.12
w = WE*cos(pi*y/Ly); % Not the wind but rather the Ekman pumping/suction
                     % resulting from the wind. This is done to skip the
                     % step of calculating w_EK from a wind field or
                     % wind-stress field. If you need to do that, do it
                     % elsewhere.
    w = repmat(w,1,I+2);
%     w = w.*linspace(0,1,I+2); % make zonally variable, increasing toward the east
%     w = w.*(0.5 + hanning(I+2)'); % make zonally variable, decreasing near the east and west boundaries
    
% latitude = 45; % in degrees please
latitude = y/d; % because one gridbox = 111km = 1degree
beta = 2*Om*cosd(latitude)/(6.3781*10^6); % df(y)/dy, approximated as a constant. Need not be so here.
                                          % denominator = mean radius of Earth in meters
    beta = beta(2:(end-1));

D = 1000; % active layer thickness
U = WE*f0/(abs(beta(round(length(beta)/2))) * D); % characteristic speed
dE = 50; % depth of Ekman layer (m)
r = f0*dE/(2*D); % 2.2.11 "inverse time scale for vorticity decay due to bottom friction."
% AH = 0.04*beta(round(length(beta)/2))*Ly^3; % see Vallis eq 14.44 (wrong somehow)
AH = 10^4; % see Yang Stockie paper
L = 2000 * 10^3; % basin scale (meters)

% relative importance of terms in the equation
eps = U/(beta(round(length(beta)/2)) * L^2);
mu = r/(beta(round(length(beta)/2)) * L);
E = AH/(beta(round(length(beta)/2)) * L^3); % assume negligible since the del^4 is complicated and
                                            % Pedlosky doesn't seem to define AH

%% Laplacian operator matrix
% NOTE: this assumes that dx = dy, which will be the case for this project.
% This can be modified for cases where dx ~= dy, which I have tentatively,
% but possibly incorrectly, labeled.

DEL2 = speye(I*J); % initialize

diag0 = -4*ones(J,1); % diagonal in a main diagonal cell matrix
diag1 =    ones(J,1); % off-diagonal in a main diagonal cell matrix (+1 and -1)
%^ this would be multiplied by (1/dy^2)
Diag = spdiags([diag1 diag0 diag1],[-1 0 1],speye(J));
%     A = spdiags(B,d,A) replaces the diagonals of A specified by d with
%         the columns of B. The output is sparse.

OffDiag = speye(J); % off-diagonal cell matrix
%^ this would be multiplied by (1/dx^2)

for i1=1:I
    for i2=1:I
        if i1==i2
            DEL2( (1 + (i1-1)*J):(i1*J) , (1 + (i1-1)*J):(i1*J) ) = Diag;
        elseif i1==i2+1
            DEL2( (1 + (i2-1)*J):(i2*J) , (1 + (i1-1)*J):(i1*J) ) = OffDiag;
        elseif i1==i2-1
            DEL2( (1 + i1*J):((i1+1)*J) , (1 + (i1-1)*J):(i1*J) ) = OffDiag;
        else
        end
    end
end

DEL2 = (1/d^2)*DEL2; % would not be necessary if dx ~= dy and the above
                     % appropriate steps were taken

% NOTE: This assumes that the streamfunction phi is forced to be zero at
% the boundaries (i.e. just outside the scope of DEL2, which is to say the
% NaN-region at the edge of Domain). For areas like islands, we need to
% modify it slightly:

DEL2(DomainNoEdge_vec==0,:) = 0;
DEL2(:,DomainNoEdge_vec==0) = 0;


%% Build d/dx and d/dy operators analogous to DEL2

DY = 0*speye(I*J); % initialize
Diag = spdiags([diag1 (-1)*diag1],[-1 1],0*speye(J));
for i1=1:I
    for i2=1:I
        if i1==i2
            DY( (1 + (i1-1)*J):(i1*J) , (1 + (i1-1)*J):(i1*J) ) = Diag;
        else
        end
    end
end
DY = DY/(2*d); % Or "DY = DY/(2*dy);" if dx ~= dy

DX = 0*speye(I*J); % initialize
OffDiag = speye(J);
for i1=1:I
    for i2=1:I
        if i1==i2+1
            DX( (1 + (i2-1)*J):(i2*J) , (1 + (i1-1)*J):(i1*J) ) = OffDiag;
        elseif i1==i2-1
            DX( (1 + i1*J):((i1+1)*J) , (1 + (i1-1)*J):(i1*J) ) = -OffDiag;
        else
        end
    end
end
DX = DX/(2*d); % Or "DX = DX/(2*dx);" if dx ~= dy

% Again to account for shorelines:
DX(DomainNoEdge_vec==0,:) = 0;
DX(:,DomainNoEdge_vec==0) = 0;

DY(DomainNoEdge_vec==0,:) = 0;
DY(:,DomainNoEdge_vec==0) = 0;


%% run the simulation (Vallis 2.2.9)

close all
clear MAX MEAN rmsFzeta rmsw rmsbetsv rmszeta Diff Diff2
reshaped_w = reshape(w(2:(end-1),2:(end-1)),numel(w(2:(end-1),2:(end-1))),1);
reshaped_Domain = reshape(Domain(2:(end-1),2:(end-1)),numel(Domain(2:(end-1),2:(end-1))),1);

dt = 10000; % so far 10000 and N=250 work just fine
III = speye(size(DEL2));
psi = zeros(size(Domain)); % initial condition of phi
%     phi = reshape(phi,numel(phi),1); % reshape for matrix stuff
N = 200;%200; % # of time steps
Diff = zeros(N,1);

% Only needs to be done the first time:
    zeta = DEL2*reshape(psi(2:(end-1),2:(end-1)),numel(psi(2:(end-1),2:(end-1))),1);


% % % % % FORWARD SCHEME (dimensional)
    imagesc(psi + nanDomain);colorbar
    title('0')
    pause(0.5)
for n=1:N % time
    
    psi_last = psi;
    
    rehaped_psi = reshape(psi(2:(end-1),2:(end-1)),numel(psi(2:(end-1),2:(end-1))),1);
    psi_x = DX*rehaped_psi;
    F = -psi_x.*DY + (DY*rehaped_psi.*DX) - r*III;
        F = F + AH*DEL2; % problem term
    zeta = zeta + dt*(F*zeta + (f0/D)*reshaped_w - repmat(beta,I,1).*psi_x); % don't forget to reshape w

    
    psi_chopped = DEL2\zeta; % prepare for next psi
    psi_chopped = reshape(psi_chopped,J,I); % reshape
    psi = zeros(size(Domain)); % initialize next psi
    psi(2:(end-1),2:(end-1)) = psi_chopped; % next psi

    Diff(n)  = sum(sum(abs(psi - psi_last)));

    contourf(flip(psi + nanDomain),20);colorbar

    caxis([-8000 8000])
    title(num2str(n))
    pause(0.01)
end % t
% % % % % 

    figure
    contourf([0:d:(d*(I+2-1))],flip(y),psi + nanDomain,20);c=colorbar;c.Label.String = 'psi (area/time of some sort)';
    title(num2str(n))
    %set(gca,'YDir','reverse')
    title(['Stream Function after ',num2str(n),' time steps'],'interpreter','latex','fontsize',20)
    xlabel('x (meters)','interpreter','latex','fontsize',18)
    ylabel('y (meters)','interpreter','latex','fontsize',18)
    
    
% figure
% plot(MAX)
% hold on
% plot(MEAN)
% legend('max','mean')
% xlabel('n')

% figure
% semilogy(rmsFzeta)
% hold on
% semilogy(rmsw)
% hold on
% semilogy(rmsbetsv)
% legend('rmsFzeta','rmsw','rmsbetsv')
% xlabel('n')

figure
semilogy(Diff)
title('Difference between iterations (convergence of solution)','interpreter','latex','fontsize',20)
xlabel('n','interpreter','latex','fontsize',20)
ylabel('$\Sigma\Sigma |\psi_{n+1} - \psi_n|$','interpreter','latex','fontsize',20)


%% visualize wind

% figure
%     contourf([0:d:(d*(I+2-1))],flip(y),w,20);c=colorbar;c.Label.String = 'm/s';
%     title(num2str(n))
%     %set(gca,'YDir','reverse')
%     title(['$w_E$'],'interpreter','latex','fontsize',20)
%     xlabel('x (meters)','interpreter','latex','fontsize',18)
%     ylabel('y (meters)','interpreter','latex','fontsize',18)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % POST-FUNCTION TEST FACILITY:

%% test case for inverse solution to Poisson's equation:

% clear all
% 
% % phi = psi_xx
% % becomes
% % phi(i) = -(dx^2 /2)*phi(i) + 0.5*psi(i+1) + 0.5*psi(i-1)
% % in matrix form:
% % psi = (eye - A)\(-(dx^2 /2)*phi)
% 
% dx = 1;
% 
% phi = [zeros(1,50), ones(1,1), zeros(1,50)]';
% % phi = [linspace(0,1,50),linspace(1-1/50,0,49)]';
% % phi = [zeros(1,25), linspace(0,1,50),linspace(1-1/50,0,49), zeros(1,25)]';
% % phi = zeros(100,1)
% % phi = phi_exp';
% 
% K = length(phi);
% A = speye(K);
% B = (1/4)*[ones(K,1) zeros(K,1) ones(K,1)];
% d = [-1 0 1];
% A = spdiags(B,d,A);
% bc = zeros(size(phi)); bc(1) = 1; bc(end) = 1;
% 
% % psi = (speye(K) - A)\((-dx^2 / 2)*phi + bc);
% psi = (speye(K) - A)\((-dx^2 / 2)*phi);
% % psi = inv(speye(K) - A)*((-dx^2 / 2)*phi);
% 
% 
% figure
% plot(phi)
% hold on
% plot(psi)
% 
% for i = 1:K
%     if i==1
%         phi_exp(i) = (psi(i+1) + 0 - 2*psi(i))/dx^2;
%     elseif i==K
%         phi_exp(i) = (0 + psi(i-1) - 2*psi(i))/dx^2;
%     else
%         phi_exp(i) = (psi(i+1) + psi(i-1) - 2*psi(i))/dx^2;
%     end
% end
% 
% figure
% plot(phi,'k')
% hold on
% plot(phi_exp,'.-r')
% % plot(diff(diff(psi)),'.-g')

%% test to see if DEL2*phi is correct (verdict: it is, and WAY FASTER!*)
% % % % % * not counting the time to build DEL2, which I didn't investigate

% phi = randn(size(Domain));
%     phi(1,:) = 0;
%     phi(end,:) = 0;
%     phi(:,1) = 0;
%     phi(:,end) = 0;
% %     phi(2:50,2:50) = 0;
% 
% tic
% test1 = DEL2*reshape(phi(2:(end-1),2:(end-1)),numel(phi(2:(end-1),2:(end-1))),1);
% test1 = reshape(test1,J,I);
% toc
% 
% tic
% for i=2:(I+1)
%     for j=2:(J+1)
%         test2(j,i) = phi(j+1,i) + phi(j-1,i) + phi(j,i+1) + phi(j,i-1) - 4*phi(j,i);
%     end
% end
% test2 = test2(2:end,2:end);
% toc
% 
% figure;subplot(121);imagesc(test1);subplot(122);imagesc(test2)
% foo = mean(mean(abs(test1-test2)));
% title(['mean(mean(abs(test1-test2))) = ',num2str(foo)])

%% test inverse for solving Poisson's equation

% f = zeros(size(Domain));
% %     f(10,10) = 1;
% for i=1:(I+2)
%     for j=1:(J+2)
%         f(i,j) = sqrt((i-10)^2 + (j-10)^2);
%     end
% end
% 
% 
% g = DEL2\reshape(f(2:(end-1),2:(end-1)),numel(f(2:(end-1),2:(end-1))),1);
% 
% f_calc = DEL2*g;
%     f_calc = reshape(f_calc,J,I);
% 
% g = reshape(g,J,I);
% 
% figure;imagesc(g);colorbar
% 
% figure
% subplot(131)
% imagesc(f(2:(end-1),2:(end-1)))
% % imagesc(log(f))
% colorbar
% subplot(132)
% imagesc(f_calc)
% colorbar
% % imagesc(real(log(f_calc)))
% subplot(133)
% imagesc(f(2:(end-1),2:(end-1))-f_calc)
% colorbar

%% test DX and DY (again, it works and is much faster)

% clear test1 test2
% 
% phi = randn(size(Domain));
%     phi(1,:) = 0;
%     phi(end,:) = 0;
%     phi(:,1) = 0;
%     phi(:,end) = 0;
% %     phi(2:50,2:50) = 0;
%     phi = phi.*Domain;
% 
% tic
% test1 = DX*reshape(phi(2:(end-1),2:(end-1)),numel(phi(2:(end-1),2:(end-1))),1);
% test1 = reshape(test1,J,I);
% toc
% 
% tic
% for i=2:(I+1)
%     for j=2:(J+1)
%         test2(j,i) = ( phi(j,i+1) - phi(j,i-1) )/(2*d); % for DX
% %         test2(j,i) = ( phi(j-1,i) - phi(j+1,i) )/(2*d); % for DY
%     end
% end
% test2 = test2(2:end,2:end);
% toc
% 
% figure;subplot(131);imagesc(phi);subplot(132);imagesc(test1);subplot(133);imagesc(test2)
% foo = mean(mean(abs(test1-test2)));
% title(['mean(mean(abs(test1-test2))) = ',num2str(foo)])


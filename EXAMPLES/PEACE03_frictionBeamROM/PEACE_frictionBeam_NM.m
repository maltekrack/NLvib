%==========================================================================
% DESCRIPTION: 
% Investigation of the dynamics of an Euler-Bernoulli beam
% with elastic dry friction nonlinearity.
% 
% For the model refinement analysis we truncate the model to a subset of M
% free-interface modes. Specifically, we define M_lo=1 and M_hi=5, for
% the low- and the high-fidelity model, respectively. Simulstaneously a
% harmonic refinement from H_lo=1 to H_hi=11 is performed.
%========================================================================
% This file is part of NLvib - PEACE.
% 
% If you use NLvib, please refer to the book:
%   M. Krack, J. Gross: Harmonic Balance for Nonlinear Vibration
%   Problems. Springer, 2019. https://doi.org/10.1007/978-3-030-14023-6.
% 
% COPYRIGHT AND LICENSING: 
% NLvib Version 1.4 Copyright (C) 2024  Malte Krack  
%										(malte.krack@ila.uni-stuttgart.de) 
%                     					Johann Gross 
%										(johann.gross@ila.uni-stuttgart.de)
%                     					University of Stuttgart
% This program comes with ABSOLUTELY NO WARRANTY. 
% NLvib is free software, you can redistribute and/or modify it under the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% For details on license and warranty, see http://www.gnu.org/licenses
% or gpl-3.0.txt.
%==========================================================================
clearvars;
close all;
clc;
addpath('../../SRC');
addpath('../../SRC/MechanicalSystems');
%% Define the beam FE model
len         = .710;             % length
height      = .05;              % height in the bending direction
thickness   = .06;              % thickness in the third dimension
E           = 141.22e9;         % Young's modulus
rho         = 7850;             % density
BCs         = 'clamped-free';   % constraints

% Setup one-dimensional finite element model of an Euler-Bernoulli beam
n_nodes = 30;           % number of equidistant nodes along length
beam = FE_EulerBernoulliBeam(len,height,thickness,E,rho,...
    BCs,n_nodes);

% Apply elastic Coulomb element at 2/5 length from clamping in
% translational direction
inode   = ceil(0.4*n_nodes); % possible are 4:2:14
dir     = 'trans';
kt      = 3e7;      
muN     = .8*3e3;
add_nonlinear_attachment(beam,inode,dir,'elasticDryFriction',...
    'stiffness',kt,'friction_limit_force',muN,'ishysteretic',1);

% Linear modal analysis for sticking contact
w   = beam.nonlinear_elements{1}.force_direction; 
Kt  = beam.nonlinear_elements{1}.stiffness*(w*w'); 
[PHI_stick,OM2] = eig(beam.K+Kt,beam.M);
om_stick_ref = sort(sqrt(diag(OM2)));

% Additional linear damping (mitigates loops that occur only for very low
% damping)
beta  = 2*1e-2/om_stick_ref(1); 
beam.D = beta*beam.K;
%% Set refinement parameters
M_lo = 1;     % number of retained modes for low-fidelity model
H_lo = 1;     % number of retained harmonics for low-fidelity model
M_hi = 5;     % number of retained modes for high-fidelity model
H_hi = 11;    % number of retained harmonics for high-fidelity model
%% Setup low-fidelity model and continue Nonlinear Mode
rom_lo = computeModalROM(beam,M_lo);

% Linear modal analysis for sticking contact
w  = rom_lo.nonlinear_elements{1}.force_direction; 
Kt = rom_lo.nonlinear_elements{1}.stiffness*(w*w'); 
[PHI_stick_lo,OM2_lo] = eig(rom_lo.K+Kt,rom_lo.M);
[om_stick_lo,ind] = sort(sqrt(diag(OM2_lo)));
PHI_stick_lo = PHI_stick_lo(:,ind);

% Nonlinear modal analysis
analysis = 'NMA';
% Analysis parameters
N = 2^9;                % number of time samples per period
log10a_e = 0;           % start vibration level (log10 of modal mass)
log10a_s = -3.5;        % end vibration level (log10 of modal mass)
imod    = 1;            % mode to be analysed                   
inorm   = 1;            % (modal) coordinate used for phase normalization

% Initial guess vector x0 = [Psi;om;del], where del is the modal
% damping ratio, estimate from underlying linear system
n = size(PHI_stick_lo,1);
om = om_stick_lo(imod); phi = PHI_stick_lo(:,imod);
Psi = zeros((2*H_lo+1)*n,1);
Psi(n+(1:n)) = phi;
x0 = [Psi;om;0];

% Solve and continue
ds = .010;
psiscl = 1e0* max(abs(Psi));
fscl = mean(abs(rom_lo.K*phi));
Sopt = struct('Dscale',...
    [1e-1*psiscl*ones(size(Psi));om_stick_ref(1)/2;1e-1;.1]);
[ X_lo , Solinfo_lo ] = solve_and_continue(x0,...
    @(X) HB_residual(X,rom_lo,H_lo,N,analysis,inorm,fscl),...
    log10a_s,log10a_e,ds,Sopt);

% Interpret solver output
om_NM_lo = X_lo(end-2,:);
del_NM_lo = X_lo(end-1,:);
log10a_lo = X_lo(end,:);
a_NM_lo = 10.^log10a_lo;
%% Setup high-fidelity model and carry out refinement
rom_hi = computeModalROM(beam,M_hi);

% Linear modal analysis for sticking contact
w  = rom_hi.nonlinear_elements{1}.force_direction; 
Kt = rom_hi.nonlinear_elements{1}.stiffness*(w*w'); 
[PHI_stick_hi,OM2_hi] = eig(rom_hi.K+Kt,rom_hi.M);
[om_stick_hi,ind] = sort(sqrt(diag(OM2_hi)));
PHI_stick_hi = PHI_stick_hi(:,ind);

% Define auxilary function handle
func = @(X) HB_residual(X,rom_hi,H_hi,N,analysis,inorm);

% In this example, no subset selection is done (all points are refined)
nSP = size(X_lo,2); % number of solution points

% Determine indices of low-fidelity variables
idxQ = 1:(size(X_lo,1)-3);
idxLambda = (size(X_lo,1)-2):size(X_lo,1);

% Initialize refined solution points and their tangents
X_hi = zeros(M_hi*(2*H_hi+1)+length(idxLambda),nSP);
tangent = zeros(M_hi*(2*H_hi+1)+length(idxLambda),nSP);

% Refine points (could be done in parallel)
for ii = 1:nSP    
    [X_hi(:,ii), tangent(:,ii),peaceOut(ii)] = peace(...
        X_lo(idxQ,ii),...
        X_lo(idxLambda,ii),...
        Solinfo_lo.tangents(:,ii),Solinfo_lo.Dscales(:,ii), ...
        M_lo,H_lo,...
        M_hi,H_hi,...
        func);
end

% Interpret solver output
om_NM_hi = X_hi(end-2,:);
del_NM_hi = X_hi(end-1,:);
log10a_hi = X_hi(end,:);
a_NM_hi = 10.^log10a_hi;
%% Reference: Continuation directly applied to high-fidelity model

% Initial guess
n  = size(PHI_stick_hi,1);
om = om_stick_hi(imod); phi = PHI_stick_hi(:,imod);
Psi = zeros((2*H_hi+1)*n,1);
Psi(n+(1:n)) = phi;
x0 = [Psi;om;0];


% Solve and continue w.r.t. Om
ds = .010;
psiscl = 1e0* max(abs(Psi));
fscl = mean(abs(rom_hi.K*phi));
Sopt = struct('Dscale',...
    [1e-1*psiscl*ones(size(Psi));om/2;1e-1;.1]);
[X_ref,Solinfo_HB] = solve_and_continue(x0,...
    @(X) HB_residual(X,rom_hi,H_hi,N,analysis,inorm,fscl),...
    log10a_s,log10a_e,ds,Sopt);

% Interpret solver output
om_NM_ref = X_ref(end-2,:);
del_NM_ref = X_ref(end-1,:);
log10a_ref = X_ref(end,:);
a_NM_ref = 10.^log10a_ref;
%% Visualize results
figure('Color','white')

% low-fidelity model
semilogx(a_NM_lo ,om_NM_lo./om_stick_ref(1),...
    'k.--','LineWidth',1.2 ); hold on;
lgd_FEP = legend('Interpreter','latex');
lgd_FEP.String{end} = ['$H=$',num2str(H_lo) ', $M=$', num2str(M_lo)];

% high-fidelity model
semilogx( ...
     a_NM_hi, om_NM_hi./om_stick_ref(1),...
    'bo','MarkerSize',6,'MarkerFaceColor','b');
lgd_FEP.String{end} = ['$H=$',num2str(H_hi) ', $M=$', num2str(M_hi)];

% failed points
if any([peaceOut.exitflag]<=0)
    semilogx( ...
        a_NM_lo([peaceOut.exitflag]<=0), om_NM_lo([peaceOut.exitflag]<=0)./om_stick_ref(1),...
        'r*','MarkerSize',6,'MarkerFaceColor','b');
    lgd_FEP.String{end} = 'failed';
end

% reference
semilogx(a_NM_ref ,om_NM_ref./om_stick_ref(1),...
    'm-','LineWidth',1.2 ); hold on;
lgd_FEP.String{end} = 'reference';

% additional plot settings
grid on;
xlabel('modal amplitude $a_{\rm NM}$',...
    'Interpreter','latex');
ylabel('$\omega_1/\omega_{1,\rm lin}$','Interpreter','latex');
lgd_FEP.Location = 'best';
lgd_FEP.EdgeColor = 'none';
lgd_FEP.Color = 'none';
%% Illustrate number of iterations during refinement
xSP = 1:nSP;
figure(999);
stem(xSP,[peaceOut.NIT], ...
    'k-','LineWidth',1);hold on
stem(xSP([peaceOut.exitflag]>0), ...
    [peaceOut([peaceOut.exitflag]>0).NIT], ...
    'r*','LineWidth',1,'LineStyle','none');
xlabel('solution point');
ylabel('NIT');
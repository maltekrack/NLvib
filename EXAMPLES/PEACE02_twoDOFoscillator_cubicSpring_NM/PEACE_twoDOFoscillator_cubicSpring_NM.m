%==========================================================================
% DESCRIPTION: 
% Investigation of the dynamics of a two-degree-of-freedom (2DOF)
% oscillator with cubic spring, which is a common benchmark problem, see
% e.g. Kerschen et al. 2009.
% First we compute the low-fidelity model (small harmonic truncation order)
% and subsequently refine to the high-fidelity model (high harmonic
% truncation order). To capture the relevant topology of the high-fidelity
% solution curve we recommend to start with H_lo=3.
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
close all
clc;
addpath('../../SRC');
addpath('../../SRC/MechanicalSystems');

% Parameters of the underlying linear system
mi = [1 1];         % masses
ki = [1 1 1];       % springs
di = 0.0*ki;        % no dampers

% Cubic spring applied to 1st DOF
nonlinear_elements = struct('type','cubicSpring',...
    'stiffness',.5,'force_direction',[1;0]);

% Define system as chain of oscillators
oscillator = ChainOfOscillators(mi,di,ki,...
    nonlinear_elements);

% Number of DOFs
n = oscillator.n;

% Define harmonic order for the low- and high-fidelity model
H_lo    = 3;      % low-fidelity model    
H_hi    = 11;     % high-fidelity model
%% Linear modal analysis
[PHI_lin,OM2] = eig(oscillator.K,oscillator.M);
[om_lin,ind] = sort(sqrt(diag(OM2)));
PHI_lin = PHI_lin(:,ind);
%% Compute low-fidelity solution curve
analysis = 'NMA';

% Analysis parameters
N = 4*H_lo+1;       % number of time samples per period
imod = 1;           % mode to be analyzed
log10a_s = .8;      % start vibration level (log10 of modal mass)
log10a_e = 2.1;     % end vibration level (log10 of modal mass)
inorm = 2;          % coordinate for phase normalization

% Initial guess vector y0 = [Psi;om;del], where del is the modal
% damping ratio, estimate from underlying linear system
om = om_lin(imod); phi = PHI_lin(:,imod);
Psi = zeros((2*H_lo+1)*n,1);
Psi(n+(1:n)) = phi;
x0 = [Psi;om;0];

% Solve and continue w.r.t. Om
% NOTE: Near the loops related to internal resonances, the 
% predictor-corrector continuation can be very sensitive to step length,
% scaling, parametrization and other parameters.
ds      = .3;
Dscale  = [1e-2*ones(size(x0,1)-2,1);1e0;1e-8;1]; % diagonal precond.
Sopt    = struct('Dscale',Dscale,'dsmin',1e-5);
[X_lo,Solinfo_lo] = solve_and_continue(x0,...
    @(X) HB_residual(X,oscillator,H_lo,N,analysis,inorm),...
    log10a_s,log10a_e,ds,Sopt);

% Determine total energy and frequency
[energy_lo,om_lo] = determineEnergyAndFrequency(X_lo,H_lo,oscillator);
%% H-refinement
N = 4*H_hi+1;

% Define auxilary function handle
func_HB = @(X) HB_residual(X,oscillator,H_hi,N,analysis,inorm);

% In this example, no subset selection is done (all points are refined)
nSP = size(X_lo,2); % number of solution points

% Determine indices of low-fidelity variables
idxQ = 1:(size(X_lo,1)-3);
idxLambda = (size(X_lo,1)-2):size(X_lo,1);

% Initialize refined solution points and their tangents
X_hi = zeros(n*(2*H_hi+1)+length(idxLambda),nSP);
tangent_hi = zeros(n*(2*H_hi+1)+length(idxLambda),nSP);

% Refine points (could be done in parallel)
for ii = 1:nSP
    [X_hi(:,ii),tangent_hi(:,ii),peaceOut(ii)] = peace(...
        X_lo(idxQ,ii),...
        X_lo(idxLambda,ii),...
        Solinfo_lo.tangents(:,ii),Solinfo_lo.Dscales(:,ii), ...
        n,H_lo,...
        n,H_hi,...
        func_HB);
end

% Determine total energy and frequency
[energy_hi,om_hi] = determineEnergyAndFrequency(X_hi,H_hi,oscillator);
%% Reference: Continuation directly applied to high-fidelity model

% Initial guess
Psi = zeros((2*H_hi+1)*n,1);
Psi(n+(1:n)) = phi;
x0 = [Psi;om;0];

% Solve and continue
Dscale  = [1e-2*ones(size(x0,1)-2,1);1e0;1e-8;1];
Sopt    = struct('Dscale',Dscale,'dsmin',1e-5);
[X_ref,Solinfo,Sol] = solve_and_continue(x0,...
    @(X) HB_residual(X,oscillator,H_hi,N,analysis,inorm),...
    log10a_s,log10a_e,ds,Sopt);

% Determine total energy and frequency
[energy_ref,om_ref] = determineEnergyAndFrequency(X_ref,H_hi,oscillator);
%% Plot FEP
figure;

% low-fidelity model
semilogx(energy_lo,om_lo./om_lin(1),...
    'k.--','LineWidth',1.2 ); hold on;
lgd_FEP = legend('Interpreter','latex');
lgd_FEP.String{end} = ['$H=$',num2str(H_lo)];

% high-fidelity model
semilogx(energy_hi,om_hi./om_lin(1),...
    'bo','MarkerSize',6,'MarkerFaceColor','b');
lgd_FEP.String{end} = ['$H=$',num2str(H_hi)];

% failed points (if any)
if any([peaceOut.exitflag]<=0)
    semilogx(energy_hi([peaceOut.exitflag]<=0), ...
        om_hi([peaceOut.exitflag]<=0)./om_lin(1),...
        'r*','MarkerSize',6,'MarkerFaceColor','b');
    lgd_FEP.String{end} = 'failed';
end

% reference
semilogx(energy_ref,om_ref./om_lin(1),...
    'm-');
lgd_FEP.String{end} = 'reference';

% additional plot settings
grid on;
xlabel('energy $E$',...
    'Interpreter','latex');
ylabel('$\omega_1/\omega_{1,\rm lin}$',...
    'Interpreter','latex');
lgd_FEP.Location = 'southeast';
%% Illustrate number of iterations during refinement
figure(999);
xSP = 1:nSP;
stem(xSP,[peaceOut.NIT], ...
    'k-','LineWidth',1);hold on
stem(xSP([peaceOut.exitflag]>0), ...
    [peaceOut([peaceOut.exitflag]>0).NIT], ...
    'r*','LineWidth',1,'LineStyle','none');
xlabel('solution point');
ylabel('NIT');
%========================================================================
% DESCRIPTION: 
% This example demonstrates the applicability of the PEACE algorithm to 
% computation of the forced and damped dynamics of a two-degree-of-freedom 
% (2DOF) oscillator with an extremely stiff unilateral spring. The
% reference response (stiffness=200, H=15) is hard to obtain using the
% conventional continuation procedure. Thus, we begin the analysis with a
% simpler system (stiffness=10, H=1), and successively refine.
%
% 
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
%========================================================================
clearvars;
close all;
clc;
addpath('../../SRC');
addpath('../../SRC/MechanicalSystems');

% Parameters of the underlying linear system
mi = [1 1];             % masses
ki = [0 1 1];           % springs
di = .06*ki;            % dampers
Fex1 = [0;.1];          % fundamental excitation force harmonic


% Initially, a relatively soft unilateral spring is applied to the 1st DOF
nonlinear_elements = struct(...
    'type','unilateralSpring',...
    'stiffness',10,...
    'gap',1,'force_direction',[1;0]); 

% Define system as chain of oscillators
system = ChainOfOscillators(mi,di,ki,...
    nonlinear_elements,Fex1);
n       = system.n;   % number of degrees of freedom
Tresp   = [0 1];      % response degree of freedom

% Determine modal frequencies of underlying linear system
[~,omsq] = eig(system.K,system.M);
omLin = sort(sqrt(diag(omsq)));
%% Compute low-fidelity solution curve
analysis    = 'FRF';

% Analysis parameters
H       = 1;                 % harmonic order of the low fidelity model
N_PC    = 2^11;              % number of time samples per period for continuation
N_PEACE = 2^9;               % number of time samples per period for refinement

Om_s    = 0.8*omLin(1);      % start frequency
Om_e    = 1.4*omLin(1);      % end frequency

% Define diagonal preconditioning as function (to adapt itself to the
% problem dimension)
Dscale = @(y0) [1e-1*ones(size(y0));1];

% Initial guess (solution of underlying linear system)
Q1 = (-Om_s^2*system.M + 1i*Om_s*system.D + system.K)\system.Fex1;
y0 = zeros((2*H+1)*length(Q1),1);
y0(length(Q1)+(1:2*length(Q1))) = [real(Q1);-imag(Q1)];

% Solve and continue w.r.t. Om
ds   = .8e-1;   % optional large ds=30e-1;
Sopt = struct( 'Dscale',Dscale(y0));

[X_HB_init,Solinfo_HB_lo] = solve_and_continue(y0,...
    @(X) HB_residual(X,system,H,N_PC,analysis),...
    Om_s,Om_e,ds,Sopt);

% Identify nonlinear solution points, exploiting that unilateral spring 
% induces a non-zero 0-th harmonic
idxNONLINEAR_lo = X_HB_init(1,:)~=0;
nSP_lo = sum(idxNONLINEAR_lo); % number of selected points

% Determine amplitude and frequency
[a_init,Om_init] = determineAmplitudeAndFrequency(X_HB_init,H,n,Tresp);
%% First H-refinement
% Select nonlinear solution points for refinement, acquire tangents
X_HB_init       = X_HB_init(:,idxNONLINEAR_lo);
tangents        = Solinfo_HB_lo.tangents(:,idxNONLINEAR_lo);

% Determine indices of low-fidelity variables
idxQ = 1:(size(X_HB_init,1)-1);
idxLambda = size(X_HB_init,1);

% set target harmonic order for the 1st refinement step
Htarget_p1 = 3;

% Define auxilary function handle
func_HB    = @(X) HB_residual(X,system,Htarget_p1,N_PEACE,analysis);

% Initialize refined solution points and their tangents
X_HB_p1    = zeros(n*(2*Htarget_p1+1)+length(idxLambda),nSP_lo);
tangent_p1 = zeros(n*(2*Htarget_p1+1)+length(idxLambda),nSP_lo);

disp(' ')
disp('first refinement step')
disp('---')

% allocate solver stats
SPtime_p1  = zeros(nSP_lo,1);
NIT_p1 = zeros(nSP_lo,1);

% Refine selected points
for ii = 1:nSP_lo       
    [X_HB_p1(:,ii), tangent_p1(:,ii),Solinfo_PEACE] = peace(...
        X_HB_init(idxQ,ii),...
        X_HB_init(idxLambda,ii),...
        tangents(:,ii),Dscale(y0), ...
        n,H,...
        n,Htarget_p1,...
        func_HB);

    % store solver stats
    SPtime_p1(ii)  = Solinfo_PEACE.sptime;
    NIT_p1(ii) = Solinfo_PEACE.NIT;
end
% Determine amplitude and frequency
[a_peace_1,Om_peace_1] = ...
    determineAmplitudeAndFrequency(X_HB_p1,Htarget_p1,n,Tresp);

%% Perform contact stiffness 'refinement'
system.nonlinear_elements{1}.stiffness = 200;
% harmonic order is not refined here
Hinit_p2   = Htarget_p1; Htarget_p2 = Htarget_p1;

% Determine indices of low-fidelity variables
idxQ = 1:(size(X_HB_p1,1)-1);
idxLambda = size(X_HB_p1,1);

% Define auxilary function handle
func_HB     = @(X) HB_residual(X,system,Htarget_p2,N_PEACE,analysis);

% Initialize refined solution points and their tangents
X_HB_p2     = zeros(n*(2*Htarget_p2+1)+length(idxLambda),nSP_lo);
tangent_p2  = zeros(n*(2*Htarget_p2+1)+length(idxLambda),nSP_lo);

disp(' ')
disp('second refinement step')
disp('---')

% allocate solver stats
SPtime_p2  = zeros(nSP_lo,1);
NIT_p2 = zeros(nSP_lo,1);

% Refine selected points
for ii = 1:nSP_lo    
    [X_HB_p2(:,ii), tangent_p2(:,ii),Solinfo_PEACE] = peace(...
        X_HB_p1(idxQ,ii),...
        X_HB_p1(idxLambda,ii),...
        tangent_p1(:,ii),Dscale(X_HB_p1(1:end-1,1)), ...
        n,Hinit_p2,...
        n,Htarget_p2,...
        func_HB);

    % store solver stats
    SPtime_p2(ii)  = Solinfo_PEACE.sptime;
    NIT_p2(ii) = Solinfo_PEACE.NIT;
end

% Determine amplitude and frequency
[a_peace_2,Om_peace_2] = ...
    determineAmplitudeAndFrequency(X_HB_p2,Htarget_p1,n,Tresp);
%% Perform second H-refinement
Hinit_p3 = 3; Htarget_p3 = 15;

% Determine indices of low-fidelity variables
idxQ = 1:(size(X_HB_p2,1)-1);
idxLambda = size(X_HB_p2,1);

% Define auxilary function handle
func_HB     = @(X) HB_residual(X,system,Htarget_p3,N_PEACE,analysis);

% Initialize refined solution points and their tangents
X_HB_p3     = zeros(n*(2*Htarget_p3+1)+length(idxLambda),nSP_lo);
tangent_p3  = zeros(n*(2*Htarget_p3+1)+length(idxLambda),nSP_lo);


disp(' ')
disp('third refinement step')
disp('---')

% allocate solver stats
SPtime_p3  = zeros(nSP_lo,1);
NIT_p3 = zeros(nSP_lo,1);

% Refine selected points
for ii = 1:nSP_lo    
    [X_HB_p3(:,ii), tangent_p3(:,ii), Solinfo_PEACE] = peace(...
        X_HB_p2(idxQ,ii),...
        X_HB_p2(idxLambda,ii),...
        tangent_p2(:,ii),Dscale(X_HB_p2(1:end-1,1)),...
        n,Hinit_p3,...
        n,Htarget_p3,...
        func_HB);

    % store solver stats
    SPtime_p3(ii) = Solinfo_PEACE.sptime;
    NIT_p3(ii) = Solinfo_PEACE.NIT;
end

% Determine amplitude and frequency
[a_peace_3,Om_peace_3] = ...
    determineAmplitudeAndFrequency(X_HB_p3,Htarget_p1,n,Tresp);
%% Reference: Continuation directly applied to high-fidelity model
H  = 15;

% Initial guess (solution of underlying linear system)
Q1 = (-Om_s^2*system.M + 1i*Om_s*system.D + system.K)\system.Fex1;
y0 = zeros((2*H+1)*length(Q1),1);
y0(length(Q1)+(1:2*length(Q1))) = [real(Q1);-imag(Q1)];

% Solve and continue w.r.t. Om
Sopt = struct( 'Dscale',Dscale(y0));
[X_HB_reference,Solinfo_HB_hi] = solve_and_continue(y0,...
    @(X) HB_residual(X,system,H,N_PC,analysis),...
    Om_s,Om_e,ds,Sopt);
idxNONLINEAR_hi    = X_HB_reference(1,:)~=0;

fprintf('\n\n')
fprintf('**INFO: number of selected nonlinear solution points=%d (lo)\n',nSP_lo)
fprintf('**INFO: number of reference nonlinear points=%d (hi)\n\n',sum(idxNONLINEAR_hi))

fprintf('== CONTINUATION:\n')
fprintf('NIT per SP lo:\t %d/%d/%.2f (min/max/mean)\n',...
    min(Solinfo_HB_lo.NIT(idxNONLINEAR_lo)), ...
    max(Solinfo_HB_lo.NIT(idxNONLINEAR_lo)), ...
    mean(Solinfo_HB_lo.NIT(idxNONLINEAR_lo)))
fprintf('NIT per SP hi:\t %d/%d/%.2f (min/max/mean)\n\n',...
    min(Solinfo_HB_hi.NIT(idxNONLINEAR_hi)), ...
    max(Solinfo_HB_hi.NIT(idxNONLINEAR_hi)), ...
    mean(Solinfo_HB_hi.NIT(idxNONLINEAR_hi)))

fprintf('time per SP lo:\t %.2e/%.2e/%.2e (min/max/mean)\n', ...
    min(Solinfo_HB_lo.sptime(idxNONLINEAR_lo)), ...
    max(Solinfo_HB_lo.sptime(idxNONLINEAR_lo)), ...
    mean(Solinfo_HB_lo.sptime(idxNONLINEAR_lo)))
fprintf('time per SP hi:\t %.2e/%.2e/%.2e (min/max/mean)\n\n', ...
    min(Solinfo_HB_hi.sptime(idxNONLINEAR_hi)), ...
    max(Solinfo_HB_hi.sptime(idxNONLINEAR_hi)), ...
    mean(Solinfo_HB_hi.sptime(idxNONLINEAR_hi)))

fprintf('== PEACE:\n')
fprintf('time and number of iterations needed per solution point refinement:\n')
fprintf('H=1-->H=3:\t\t %.2e/%.2e/%.2e | %d/%d/%.2f (min/max/mean)\n', ...
    min(SPtime_p1),max(SPtime_p1),mean(SPtime_p1), ...
    min(NIT_p1(NIT_p1~=0)),max(NIT_p1),mean(NIT_p1) )
fprintf('gamma lo--> gamma hi:\t %.2e/%.2e/%.2e | %d/%d/%.2f (min/max/mean)\n', ...
    min(SPtime_p2),max(SPtime_p2),mean(SPtime_p2), ...
    min(NIT_p2(NIT_p2~=0)),max(NIT_p2),mean(NIT_p2))
fprintf('H=3-->H=15:\t\t %.2e/%.2e/%.2e | %d/%d/%.2f (min/max/mean)\n', ...
    min(SPtime_p3),max(SPtime_p3),mean(SPtime_p3), ...
    min(NIT_p3(NIT_p3~=0)),max(NIT_p3),mean(NIT_p3))
fprintf('-------------------------------------------\n')


% Determine amplitude and frequency
[a_reference,Om_reference] = ...
    determineAmplitudeAndFrequency(X_HB_reference,H,n,Tresp);

%% Visualize results
figure; hold on;
plot( Om_init./omLin(1), a_init,...
    'k--','LineWidth',1,'Marker','.','DisplayName','initial');
plot( Om_peace_1./omLin(1),a_peace_1,...
    'b-', ...
    'Marker','.', ...
    'MarkerSize',6, ...
    'MarkerFaceColor','b', ...
    'DisplayName','$H=1 \rightarrow H=3$');
plot( Om_peace_2./omLin(1),a_peace_2,...
    'r', ...
    'Marker','.', ...
    'MarkerSize',6, ...
    'MarkerFaceColor','w', ...
    'DisplayName','$\gamma=10 \rightarrow \gamma=200$');
plot( Om_peace_3./omLin(1),a_peace_3,...
    'k', ...
    'Marker','o', ...
    'MarkerSize',6, ...
    'MarkerFaceColor','w', ...
    'DisplayName','$H=3 \rightarrow H=15$');
plot( Om_reference./omLin(1),a_reference,...
    'm', ...
    'Marker','none', ...
    'DisplayName','reference');
xlabel('$\Omega/\omega_{1,\rm lin}$','Interpreter','latex')
ylabel('$\max(q_1)$','Interpreter','latex')
set(gca,'FontSize', 18);
legend('EdgeColor','none','Color','none','Location','best',...
    'Interpreter','latex')

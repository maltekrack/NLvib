%========================================================================
% DESCRIPTION:
% Matlab function for H- and M-refinement based on the PEACE algorithm
% published in "A new paradigm for multi-fidelity continuation using 
% parallel model refinement" (currently under review).
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
% Q_HB_lo       Fourier coefficients of low-fidelity model
% restX_lo      remaining vector of unknowns of low-fidelity model
% tangent_lo    low-fidelity tangent
% M_lo          low-fidelity modal truncation order
% H_lo          low-fidelity harmonic truncation order
% M_hi          high-fidelity modal truncation order
% H_hi          high-fidelity harmonic truncation order
%========================================================================
function [X_HB, tangent, stats] = peace(Q_HB_lo,restX_lo,tangent_lo,...
    Dscale,...
    M_lo,H_lo,...
    M_hi,H_hi,...
    func_HB)
%% Set default options for Newton-type method (fsolve)
opts_fsolve = optimoptions( ...
    'fsolve',...
    'Display', 'off',...    
    'MaxIterations',50,...
    'Jacobian','on');
%% Handle M-refinement

% Set up expansion matrix
if M_hi~=M_lo
    TexpM = eye(M_hi);
    TexpM(:,1+M_lo:M_hi) = [];
else
    TexpM = eye(M_hi);
end

% Expand tangent, diagonal preconditioning and Fourier coefficients
Q_HB_0          = reshape(TexpM*reshape(Q_HB_lo,M_lo,[]),[],1);
tangentTmp_0    = ...
    reshape(TexpM*reshape(tangent_lo(1:length(Q_HB_lo)),M_lo,[]),[],1);
DscaleTmp_0     = ...
    reshape(TexpM*reshape(Dscale(1:length(Q_HB_lo)),M_lo,[]),[],1);
%% Handle H-refinement

% Set up expnansion matrix
if H_hi~=H_lo
    TexpH       = eye(M_hi*(2*H_hi+1));
    idxH_init   = setdiff(1:H_hi,1:H_lo);
    idxH_del    = zeros(M_hi*2*length(idxH_init),1);
    for dd=1:length(idxH_init)
        idxH_del(1+(dd-1)*2*M_hi:(dd)*2*M_hi) = ...
            1+M_hi+(idxH_init(dd)-1)*(2*M_hi) : ...
            (M_hi+(idxH_init(dd))*(2*M_hi));
    end
    TexpH(:,idxH_del) = [];
else
    TexpH = eye(M_hi*(2*H_hi+1));
end

% Expand tangent, diagonal preconditioning and Fourier coefficients
tangent_0   = [TexpH*tangentTmp_0;tangent_lo(1+length(Q_HB_lo):end)];
Dscale_0    = [TexpH*DscaleTmp_0;Dscale(1+length(Q_HB_lo):end)];
Dscale_0(Dscale_0==0) = max(Dscale_0(1:M_hi));
X_HB_0      = [TexpH*Q_HB_0;restX_lo]./Dscale_0;
tangent_0   = tangent_0./Dscale_0;
%% Solve extended residual

% Define extended residual including specific options
Sopt = struct( ...
    'flag',1,...
    'parametrization','orthogonal',...
    'Dscale',Dscale_0,...
    'jac','full');
func_Ext = @(X) extended_residual(X,...
    X_HB_0, tangent_0,...
    @(X) func_HB(X), Sopt);

% Call Newton-type method (fsolve)
tic;
[ X_new,R,exitflag,outInfo, J] = fsolve( ...
            func_Ext, X_HB_0, opts_fsolve );

SPtime = toc;
if exitflag>0
    fprintf(['Refinement successful (transpose(R)*R' ...
        '=%.4e, number of iterations = %d\n'], R'*R, outInfo.iterations);
else
    fprintf(['Refinement failed (transpose(R)*R' ...
        '=%.4e, number of iterations = %d\n'], R'*R, outInfo.iterations);
end

% Obtain original scaling 
X_HB = diag(Dscale_0)*X_new;

% Collect solver statistics
stats.NIT       = outInfo.iterations;
stats.funcCount = outInfo.funcCount;
stats.exitflag  = exitflag;
stats.sptime    = SPtime;

%% Compute new tangent (for the case of further refinement)
J       = J*diag(1./Dscale_0);
tangent = compute_tangent(J,X_HB,tangent_0);
% obtain original scaling
tangent = diag(Dscale_0)*tangent;
end
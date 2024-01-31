%========================================================================
% DESCRIPTION: 
% Matlab function for the computation of the local tangent. This function
% is used in <solve_and_continue.m> and <peace.m>.
%========================================================================
% This file is part of NLvib.
% 
% If you use NLvib, please refer to the book:
%   M. Krack, J. Gross: Harmonic Balance for Nonlinear Vibration
%   Problems. Springer, 2019. https://doi.org/10.1007/978-3-030-14023-6.
% 
% COPYRIGHT AND LICENSING: 
% NLvib Version 1.4 Copyright (C) 2024  Malte Krack  
%                                       (malte.krack@ila.uni-stuttgart.de) 
%                                       Johann Gross 
%                                       (johann.gross@ila.uni-stuttgart.de)
%                                       University of Stuttgart
% This program comes with ABSOLUTELY NO WARRANTY. 
% NLvib is free software, you can redistribute and/or modify it under the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% For details on license and warranty, see http://www.gnu.org/licenses
% or gpl-3.0.txt.
%========================================================================
function [tangent] = compute_tangent(J,X_0,tangent_0)

[~,kk] = sort(abs(tangent_0./max(abs(X_0),1e-4)),...
                    1,'descend');
% Temporarily switch off warning
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
for ik=1:length(kk)
    k = kk(ik);
    c = zeros(length(X_0),1); c(k) = 1;
    % Determine unit tangent to the solution
    % path (Eq. 4.8)
    tangent = [J(1:end-1,:);c']\...
        [zeros(size(J,1)-1,1);1];
    if ~any(isnan(tangent))
        break;
    end
end
% Switch warning on again
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:singularMatrix');
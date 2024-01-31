% ====================================
% compute new tangent
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
%% extended_residual
function [Rext,dRextdX] = extended_residual(X,Xref,zref,...
    fun_residual,Sopt)
%% Evaluation of the residual function and its derivative
switch Sopt.jac
    case {'full','on'}
        [R,dRdX] = feval(fun_residual,diag(Sopt.Dscale)*X);
    case 'x'
        [R,dRdx] = feval(fun_residual,diag(Sopt.Dscale)*X);
        
        % Approximate dfdlam using finite differences
        dlam = max(Sopt.eps*abs(X(end)),Sopt.eps);
        Xtmp = X; Xtmp(end) = X(end)+dlam;
        dlam = Xtmp(end)-X(end);
        Rp = feval(fun_residual,diag(Sopt.Dscale)*Xtmp);
        dRdlam = (Rp-R)/dlam/Sopt.Dscale(end);
        dRdX = [dRdx dRdlam];
    otherwise
        R = feval(fun_residual,diag(Sopt.Dscale)*X);
        Smyopt = struct('epsrel',Sopt.eps,'epsabs',Sopt.eps, ...
            'ikeydx',1,'ikeyfd',1);
        dRdX = finite_difference_jacobian(fun_residual,...
            diag(Sopt.Dscale)*X,Smyopt);
end
%% Evaluation of the parametrization constraint equation and its derivative
if Sopt.flag
    switch Sopt.parametrization
        case {'arc_length'}
            % iteration on a normal plane, perpendicular to tangent
            p = transpose((X-Xref))*(X-Xref)-Sopt.ds^2;
            dpdX = 2*transpose(X-Xref);
        case 'pseudo_arc_length'
            p = Sopt.pseudoxi*...
                transpose((X(1:end-1)-...
                Xref(1:end-1)))*(X(1:end-1)-...
                Xref(1:end-1)) + (1-Sopt.pseudoxi)*...
                (X(end)-Xref(end))^2 - ...
                Sopt.ds^2;
            dpdX = [Sopt.pseudoxi*...
                2*transpose(X(1:end-1)-Xref(1:end-1)) ...
                2*(1-Sopt.pseudoxi)*...
                (X(end)-Xref(end))] ;
        case {'local','orthogonal'}
            % Solution on hyperplane through Xref, normal to zref
            % Eq. (4.22)
            p       = transpose(zref)*(X-Xref);
            dpdX    = transpose(zref);
        case {'normal'}
            [~,kk] = sort(abs(zref./max(abs(X),1e-4)),...
                1,'descend');
            % Temporarily switch off warning
            warning('off','MATLAB:nearlySingularMatrix');
            warning('off','MATLAB:singularMatrix');
            for ik=1:length(kk)
                k = kk(ik);
                c = zeros(length(X),1); c(k) = 1;
                % Determine unit tangent to the solution
                % path (Eq. 4.8)
                ztmp = [dRdX;c']\...
                    [zeros(size(dRdX,1),1);1];
                if ~any(isnan(ztmp))
                    % Erfolgreich!
                    break;
                end
            end
            % Switch warning on again
            warning('on','MATLAB:nearlySingularMatrix');
            warning('on','MATLAB:singularMatrix');
            z       = ztmp/norm(ztmp);
            p       = z'*(X-Xref);
            dpdX    = transpose(z);
        otherwise
            error('Invalid specification of path continuation constraint');
    end
end
%% Assembly of extended residual and derivative
if Sopt.flag
    Rext = [R;p];
    dRextdX = [dRdX*diag(Sopt.Dscale);dpdX];
else
    Rext = R;
    dRextdX = dRdX(:,1:end-1)*diag(Sopt.Dscale(1:end-1));
end
end
%% ===== AUXILLARY FUNCTION =====
% finite_difference_jacobian
function J = finite_difference_jacobian(Hfuncname, x0, Smyopt, varargin)
%========================================================================
% DESCRIPTION:
% Function for the calculation of a Jacobian matrix. The derivatives of the
% Jacobian matrix are approximated by forward finite differences, by back-
% ward finite differences or by central finite differences. Different ways
% are possible for the calculation of the step size needed for the finite
% differences approximation of the derivatives. The algorithm is based on
% the FORTRAN77 subroutine 'fdjac' presented in the book 'Numerical Recipes
% in FORTRAN77' (chapter 9, section 7, pages 376-386), Cambridge University
% Press, 1992.
%========================================================================
%
% Called functions:
%   Hfuncname (function handle)
%
% Input data dictionary:
%   Hfuncname           FunctionHandle  Function handle for the function
%   Smyopt.epsrel       Struct.Rscalar  Maximum relative error
%   Smyopt.epsabs       Struct.Rscalar  Machine precision tolerance
%   Smyopt.ikeydx       Struct.Iscalar  Key for the step size calculation
%   Smyopt.ikeyfd       Struct.Iscalar  Key for the finite differences
%                                       method
%   x0                  Rarray          Input variable vector
%   varargin            Cell Array      Optional parameters for the
%                                       function
%
% Output data dictionary:
%   J                   Rarray          Computed Jacobian matrix
%% Start of the calculation of the Jacobian matrix based on finite
% differences
for ii = 1:length(x0)
    
    %% Calculation of the actual step size for the finite difference method
    switch (Smyopt.ikeydx)
        
        % First method for the calculation of the finite difference step
        % size
        case (1)
            dx = Smyopt.epsrel*abs(x0(ii));
            if (dx == 0)
                dx = Smyopt.epsrel;
            end
            
            % Second method for the calculation of the finite difference
            % step size
        case (2)
            dx = sqrt(Smyopt.epsabs)*(1+abs(x0(ii)));
            
            % Third method for the calculation of the finite difference
            % step size
        case (3)
            dx = Smyopt.epsrel*abs(x0(ii));
            dx = max(dx,Smyopt.epsabs);
            
    end
    %% Jacobian matrix computation based on different finite difference methods
    switch (Smyopt.ikeyfd)
        
        % Calculation based on a forward finite difference approximation
        case (1)
            
            % Initialize the Jacobian matrix and the function value
            if (ii == 1)
                f0 = feval(Hfuncname,x0,varargin{:});
                J = zeros(length(f0),length(x0));
            end
            
            % Udpate of the input variable vector
            xp = x0;
            xp(ii) = x0(ii)+dx;
            
            % Recompute the finite difference step size to reduce the
            % finite precision error
            dx = xp(ii)-x0(ii);
            
            % Calculation of the actual column of the Jacobian matrix
            fp = feval(Hfuncname,xp,varargin{:});
            J(:,ii) = (fp-f0)/dx;
            %% Calculation based on a backward finite difference approximation
        case (2)
            
            % Initialize the Jacobian matrix and the function value
            if (ii == 1)
                f0 = feval(Hfuncname,x0,varargin{:});
                J = zeros(length(f0),length(x0));
            end
            
            % Udpate of the input variable vector
            xm = x0;
            xm(ii) = x0(ii)-dx;
            
            % Recompute the finite difference step size to reduce the
            % finite precision error
            dx = x0(ii)-xm(ii);
            
            % Evaluate the function value at the updated input variable
            % vector
            fm = feval(Hfuncname,xm,varargin{:});
            
            % Calculation of the actual column of the Jacobian matrix
            J(:,ii) = (f0-fm)/dx;
            %% Calculation based on a central finite difference approximation
        case (3)
            
            % Udpate of the input variable vector
            xp = x0;
            xp(ii) = x0(ii)+dx;
            xm = x0;
            xm(ii) = x0(ii)-dx;
            
            % Recompute the finite difference step size to reduce the
            % finite precision error
            dx = (xp(ii)-xm(ii))/2;
            
            % Evaluate the function value at the updated input variable
            % vector
            fp = feval(Hfuncname,xp,varargin{:});
            fm = feval(Hfuncname,xm,varargin{:});
            
            % Initialize the Jacobian matrix
            if (ii == 1)
                J = zeros(length(fp),length(x0));
            end
            
            % Calculation of the actual column of the Jacobian matrix
            J(:,ii) = (fp-fm)/2/dx;
    end
end
end
function [energy,om] = determineEnergyAndFrequency(X_HB,H,oscillator)
% X_HB          solver output
% H             number of harmonics
% oscillator    two-DOF oscillator object with cubic spring attached to 1st
%               DOF
%% Interpret solver output
Psi = X_HB(1:end-3,:);
om = abs( X_HB(end-2,:) );
log10a = X_HB(end,:);
a = 10.^log10a;
Q = Psi.*repmat(a,size(Psi,1),1);

% Determine total energy in the system from the displacement and velocity
% at t=0
energy = zeros(size(a));
for i=1:size(X_HB,2)
    Qi = reshape(Q(:,i),2,2*H+1);
    q0 = Qi(:,1) + sum(Qi(:,2:2:end),2);
    u0 = sum(Qi(:,3:2:end),2)*om(i);
    energy(i) = 1/2*u0'*oscillator.M*u0 + 1/2*q0'*oscillator.K*q0 + ...
        oscillator.nonlinear_elements{1}.stiffness*q0( 1 )^4/4;
end
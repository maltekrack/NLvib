function [qrespMAX,Om_HB] = determineAmplitudeAndFrequency(X_HB,H,n,Tresp)
% X_HB      solver output
% n         number of DOFs;
% H         number of harmonics
% Tresp     row vector selecting response coordinate
%% Extract Fourier coefficients and frequency
Q_HB        = X_HB(1:end-1,:);
Om_HB       = X_HB(end,:);
%% Convert from real to complex Fourier coefficients
I0          = 1:n; ID = n+(1:H*n);
IC          = n+repmat(1:n,1,H)+n*kron(0:2:2*(H-1),ones(1,n)); IS = IC+n;
Q           = zeros(n*(H+1),size(Q_HB,2));
Q(I0,:)     = Q_HB(I0,:);
Q(ID,:)     = Q_HB(IC,:)-1i*Q_HB(IS,:);
%% Evaluate response coordinate in discrete time
tau         = linspace(0,2*pi,2^10)';
H_iDFT      = exp(1i*tau*(0:H));
qresp       = kron(eye(1+H),Tresp)*Q;
%% Determine maximum in discrete time
qrespMAX    = zeros(size(qresp,2),1);
for nn=1:length(qrespMAX)
    qrespMAX(nn) = max(abs(real(H_iDFT*qresp(:,nn))));
end
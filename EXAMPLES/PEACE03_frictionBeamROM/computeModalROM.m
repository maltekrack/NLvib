function rom = computeModalROM(model,nm)
%   model   mechanical model object
%   nm      modal truncation order
%% Compute mass-normalized free interface normal modes
[PHI,OM2] = eigs(model.K,model.M,nm,'sm'); 
[om,ind] = sort(sqrt(diag(OM2)));
PHI = PHI(:,ind);
PHI = PHI./repmat(sqrt(diag(PHI'*model.M*PHI))',size(PHI,1),1);
%% Apply modal truncation and store reduced order model
rom.PHI     = PHI;
rom.M       = eye(size(PHI,2));
rom.K       = diag(om.^2);
rom.D       = PHI'*model.D*PHI;
rom.Fex1    = PHI'*model.Fex1;
rom.n       = nm;
rom.nonlinear_elements = model.nonlinear_elements;
for ii=1:length(rom.nonlinear_elements)
    % NOTE: This should work for local nonlinear elements only
    rom.nonlinear_elements{ii}.force_direction = ...
        PHI'*model.nonlinear_elements{ii}.force_direction;
end

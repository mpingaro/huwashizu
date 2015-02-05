% by Marco Pingaro & Paolo Venini

function [MCORR] = mccook(EL)
%% CORRESPONDENCE MATRIX
% In the correspondence matrix rows are relative to elements whereas 
% columns are relative to the degrees of freedom for each element.
i = EL(:,1)'; j = EL(:,2)'; z = EL(:,3)';
MCORR = [2*i-1 ; 2*i ; 2*j-1 ; 2*j ; 2*z-1 ; 2*z]; MCORR = MCORR';

NNOD = max(max(EL));
NGDLU = 2*size(EL,1)+2*NNOD;
MCORR(:,7) = 2*NNOD+1:2:NGDLU;
MCORR(:,8) = 2*NNOD+2:2:NGDLU;
end
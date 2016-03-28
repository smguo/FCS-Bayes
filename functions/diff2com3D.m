% FCS model function 
% k = wz/wx, structure factor
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function F = diff2com3D(a,t,k)

F = a(5) + abs(a(1))./(1 + t./abs(a(3))).*(1+ t./(abs(a(3)).*k^2)).^(-0.5)...
        + abs(a(2))./(1 + t./abs(a(4))).*(1+ t./(abs(a(4)).*k^2)).^(-0.5); 
    if any(~isreal(F))
error('MODELFUN has returned complex values.');
    end
% Perform weighted non-linear regression
% Reference: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% http://www.mathworks.com/products/statistics/demos.html?file=/products/demos/shipping/stats/wnlsdemo.html
function [bFitw,resid,J,COVBw,mse,sigma_sqi] = nlinfitw (x,y,modelFun,beta0,w)
yw = sqrt(w).*y;
modelFunw = @(b,x) sqrt(w).*modelFun(b,x);

options = statset('MaxIter' , 2000, 'DerivStep',10^-6);
[bFitw,rw,Jw,COVBw,msew] = nlinfit(x,yw,modelFunw,beta0,options);

% % recaculate resid and mse
wn = w ./ sum(w); % normalize weight function
resid = y - modelFun(bFitw,x);
mse = sum(resid.^2.*wn)*( numel(wn)/(numel(wn)-numel(beta0)) ); % calculate mse according to weight
J = Jw;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % recaculate sigma squared at each xi according to fitting residuals
const = mean(resid.^2.*w);
sigma_sqi = 1./w*const*(numel(w))/(numel(w)-numel(beta0));









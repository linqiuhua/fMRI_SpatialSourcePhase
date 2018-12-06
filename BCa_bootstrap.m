function [p,CI] = BCa_bootstrap(data,loo,boot,varargin)

% Computes the bias-corrected and accelerated bootstrap estimate of  Efron,
% B., & Tibshirani, R. J. (1993). An Introduction to the Bootstrap, Chapman
% & Hall/CRC: New York.
% 
% Usage:
% [p,CI] = BCa_bootstrap(data,loo,boot,null,confidence)
% [p,CI] = BCa_bootstrap(data,loo,boot,null,confidence,adjustment)
% 
% p is the two-tailed p-value and CI is the confidence interval.
% 
% data - the statistic of interest calculated on the data
% loo - a vector of length N of the leave-one-out values of the statistic
%       of interest.
% boot - the bootstrapped values of the statistic of interest
% null - (optional) the value of the statistic of interest under the null
%        hypothesis. Defaults to 0.
% confidence - (optional) the % confidence for upper and lower bounds of
%              the confidence interval. Defaults to 95%
% adjustment - (optional) NOT RECOMMENDED. This is included for a
%              specialized application. The integer value entered here will
%              include that many observations equal to the null value into
%              the bootstrap values AFTER calculating the bias and
%              acceleration.
% 
% Author£ºJared Van Snellenberg, PhD 2009

if any(isnan(boot))
    error('Error, NaNs in the bootstrapped data');
end

if isempty(varargin)||isempty(varargin{1})
    null = 0;
else
    null = varargin{1};
end
if isempty(varargin)||length(varargin)<2||isempty(varargin{2})
    c = 0.025;
else
    c = (1-varargin{2})/2;
end
if isempty(varargin)||length(varargin)<3||isempty(varargin{3})
    adjust = 0;
else
    adjust = varargin{3};
end

boot = sort(boot);

a = sum( (mean(loo) - loo).^3 ) / (6 * sum( (mean(loo) - loo).^2 ) ^(3/2));

z = norminv(sum(boot < data) / length(boot));

alpha(1) = normcdf( z + (z + norminv(c)) / (1 - a*(z + norminv(c))) );
alpha(2) = normcdf( z + (z + norminv(1-c)) / (1 - a*(z + norminv(1-c))) );

if adjust
    boot(end+1:end+adjust) = null;
    boot = sort(boot);
end

CI = [boot( max([fix(alpha(1)*length(boot)) 1]) ) boot(fix(alpha(2)*length(boot)))];

p = sum(boot < null) / length(boot);
if ~p
    p = 1 / length(boot);
elseif p == 1
    p = 1 - 1 / length(boot);
end
%p = normcdf( ( 1 - a*z - z / (norminv(p) - z) ) / ( a + 1 / (norminv(p) - z) ) );
p = normcdf( ( norminv(p) * (1 - a*z) + z*(a*z - 2) ) / ( 1 + a*norminv(p) - a*z));
if p < .5
    p = p * 2;
else
    p = (1 - p) * 2;
end

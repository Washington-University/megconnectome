function [fitx, goodness] = hcp_fitspectrum(x, y)

% HCP_FITSPECRUM fits a "one-over-f plus constant" model to a power spectrum
%
% The model that is fitted is
%   y = a/x + b
% with positive a and b.
%
% Example:
%   x = (1:100)';
%   y = 2./x + 3 + 0.1*rand(size(x));
%   [fitx, goodness] = hcp_fitspectrum(x, y)
%
% This replaces the following use of the fit() function from
% the Mathworks curve fitting toolbox
%   g = fittype('abs(a)*1/x + abs(b)')
%   [fitx, goodness] = fit(x, y, g)
%
% Note that this implementation does not guarantee that a
% and b are positive.

% Copyright (C) 2011-2014 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)
%
% This file is part of megconnectome.
%
% megconnectome is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% megconnectome is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with megconnectome.  If not, see <http://www.gnu.org/licenses/>.

x = x(:);
y = y(:);
n = numel(x);

if false
  % convert the model into a linear one, i.e. z = a + x*b
  z = y.*x;
  beta = [ones(n,1) x] \ z;
else
  beta = [1./x ones(n,1)] \ y;
end

a = beta(1);
b = beta(2);

if a<0
  warning('the fitted parameter a is not positive');
end
if b<0
  warning('the fitted parameter b is not positive');
end

% compute the estimation
est = a./x + b;

c = corrcoef(y, est);

fitx.a    = a;
fitx.b    = b;
fitx.est  = est;

goodness.rsquare = c(1,2).^2;
goodness.sse     = sum((y-est).^2);
goodness.dfe     = n-2;
goodness.rmse    = sqrt(goodness.sse/goodness.dfe);

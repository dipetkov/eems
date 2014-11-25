% Copyright (C) 2001 Paul Kienzle <pkienzle@users.sf.net>
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, see <http://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn {Function File} {@var{v} =} nanstd (@var{X})
% @deftypefnx{Function File} {@var{v} =} nanstd (@var{X}, @var{opt})
% @deftypefnx{Function File} {@var{v} =} nanstd (@var{X}, @var{opt}, @var{dim})
% Compute the standard deviation while ignoring NaN values.
%
% @code{nanstd} is identical to the @code{std} function except that NaN values are
% ignored.  If all values are NaN, the standard deviation is returned as NaN.
% If there is only a single non-NaN value, the deviation is returned as 0. 
%
% The argument @var{opt} determines the type of normalization to use. Valid values
% are
%
% @table @asis 
% @item 0:
%   normalizes with @math{N-1}, provides the square root of best unbiased estimator of 
%   the variance [default]
% @item 1:
%   normalizes with @math{N}, this provides the square root of the second moment around 
%   the mean
% @end table
%
% The third argument @var{dim} determines the dimension along which the standard
% deviation is calculated.
%
% @seealso{std, nanmin, nanmax, nansum, nanmedian, nanmean}
% @end deftypefn

function v = nanstd (X, varargin)
if nargin < 1
  print_usage;
else
  if nargin < 2
    dim = min(find(size(X)>1));
    if isempty(dim), dim=1; end;
  else
    dim = varargin{1};
  end
	
  % determine the number of non-missing points in each data set
  n = sum (~isnan(X), varargin{:});
    
  % replace missing data with zero and compute the mean
  X(isnan(X)) = 0;
  meanX = sum (X, varargin{:}) ./ n;
    
  % subtract the mean from the data and compute the sum squared
  sz = ones(1,length(size(X)));
  sz(dim) = size(X,dim);
  v = sumsq (X - repmat(meanX,sz), varargin{:});
    
  % because the missing data was set to zero each missing data
  % point will contribute (-meanX)^2 to sumsq, so remove these
  v = v - (meanX .^ 2) .* (size(X,dim) - n);
  
  % compute the standard deviation from the corrected sumsq using
  % max(n-1,1) in the denominator so that the std for a single point is 0
  v = sqrt ( v ./ max(n - 1, 1) );

end


function Xsq = sumsq (X, varargin)
Xsq = sum (X.*X, varargin{:});

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
% @deftypefn {Function File} {@var{v} =} nanmean (@var{X})
% @deftypefnx{Function File} {@var{v} =} nanmean (@var{X}, @var{dim})
% Compute the mean value while ignoring NaN values.
%
% @code{nanmean} is identical to the @code{mean} function except that NaN values
% are ignored.  If all values are NaN, the mean is returned as NaN. 
%
% @seealso{mean, nanmin, nanmax, nansum, nanmedian}
% @end deftypefn

function v = nanmean (X, varargin) 
if nargin < 1
  print_usage;
else 
  n = sum (~isnan(X), varargin{:});
  n(n == 0) = NaN;
  X(isnan(X)) = 0;
  v = sum (X, varargin{:}) ./ n;
end

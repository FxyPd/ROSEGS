function [xfilt, xi, xmedian, xsigma] = hampel(x, k, nsigma)
%HAMPEL Outlier removal via Hampel identifier.
%   Y = HAMPEL(X) replaces any element in vector X that is more than three
%   standard deviations from the median of itself and up to three
%   neighboring elements with that median value.  The standard deviation is
%   estimated by scaling the local median absolute deviation (MAD) by a
%   constant scale factor.  If X is a matrix, HAMPEL operates over the
%   columns of X.
% 
%   Y = HAMPEL(X,K) specifies the number of adjacent samples, K, on either
%   side of each sample in X over which to compute the Hampel identifier.
%   The default value of K is 3.
% 
%   Y = HAMPEL(X,K,NSIGMA) specifies the number of estimated standard
%   deviations above which it will replace elements in X with the local
%   median.
% 
%   [Y,I] = HAMPEL(...) returns a logical matrix, I, of the same shape as X
%   which is true when the corresponding element in X is identified as an
%   outlier.
% 
%   [Y,I,XMEDIAN,XSIGMA] = HAMPEL(...) additionally returns the local medians
%   and estimated standard deviations (scaled MAD) for each element in X.
% 
%   HAMPEL(...) without output arguments plots the filtered signal as well as
%   the outliers removed from the signal.
% 
%   % Example 1:
%   %   Remove spikes from a sinusoid
%   x = sin(2*pi*(0:99)/100);
%   x(6) = 2;
%   x(20) = -2;
%   hampel(x)
%
%   % Example 2:
%   %   plot statistics returned by hampel identifier
%   x = sin(2*pi*(0:99)/100);
%   x(6) = 2;
%   x(20) = -2;
%   [y,i,xmedian,xsigma] = hampel(x);
%   n = numel(x);
%   plot(1:n,xmedian-3*xsigma,'r', ...
%        1:n,xmedian+3*xsigma,'r', ...
%        1:n,x, ...
%        find(i),x(i),'sk ')
%   legend('lower limit','upper limit','original signal','outliers');
%
%   See also MEDFILT1, MEDIAN, FILTER, SGOLAYFILT.

%   Copyright 2015 The MathWorks, Inc.

narginchk(1,3);

if nargin<2
  k = 3;
else
  validateattributes(k,{'numeric'},{'integer','scalar','positive'});
end

% default is three standard deviations of a gaussian distributed input
if nargin<3
  nsigma=3;
else
  validateattributes(nsigma,{'numeric'},{'real','scalar','nonnegative'});
end

validateattributes(x,{'single','double'},{'real','2d','nonsparse'});

needsTranspose = isrow(x);
if needsTranspose
  x = x(:);
end

% compute the median absolute deviations and the corresponding medians over
% the size of the filter:  ignore samples that contain NaN and truncate
% at the borders of the input
[xmad,xmedian] = movmad(x,2*k+1,1,'omitnan','central');

% scale the MAD by ~1.4826 as an estimate of its standard deviation
scale = -1 /(sqrt(2)*erfcinv(3/2));
xsigma = scale*xmad;

% identify points that are either NaN or beyond the desired threshold
xi = ~(abs(x-xmedian) <= nsigma*xsigma);

% replace identified points with the corresponding median value
xf = x;
xf(xi) = xmedian(xi);

if nargout==0
  plotOutliers(x, xf, xi)
else
  xfilt = xf;
  if needsTranspose
    xfilt = xf.';
    xi = xi.';
    xmedian = xmedian.';
    xsigma = xsigma.';
  end    
end

function plotOutliers(x, xf, xi)
t = (1:size(x,1))';
if size(x,2)==1
  hLine = plot(t,x,'.-', ...
               t,xf, ...
               t(xi),x(xi),'ks');
  if numel(hLine)<3
    legend(getString(message('signal:hampel:OriginalSignal')), ...
           getString(message('signal:hampel:FilteredSignal')));
  else
    legend(getString(message('signal:hampel:OriginalSignal')), ...
           getString(message('signal:hampel:FilteredSignal')), ...
           getString(message('signal:hampel:Outliers')));
  end
else
  colors=get(0,'DefaultAxesColorOrder');
  for i=1:size(x,2)
    if i==1
      hLine = plot(t,xf(:,i),'Color',colors(1+mod(i-1,size(colors,1)),:));
    else
      hLine = line(t,xf(:,i),'Color',colors(1+mod(i-1,size(colors,1)),:));
    end
    hLineOutliers = line(t(xi(:,i)),x(xi(:,i),i),'LineStyle','none','Marker','s','MarkerEdgeColor',hLine.Color);
    if i==1 && ~isempty(hLineOutliers)
      legend(getString(message('signal:hampel:FilteredSignal')), ...
             getString(message('signal:hampel:Outliers')));
    end
  end
end

function x = randsmpl(p, m, n, varargin)
%RANDSMPL  Independent sampling from a discrete distribution.
%          
%   x = randsmpl(p, m, n) returns an m-by-n matrix x of random samples
%   drawn independently from the input (discrete) distribution specified
%   with pmf p. Suppose that the sample space comprises K samples, then p
%   must be a (row- or column-) vector containing K probability masses 
%   summing to 1. The output, x(i,j) = k, k = 1, ..., K, describes that
%   the k-th sample is drawn in the (i,j)-th trial, for i = 1, ..., m and
%   j = 1,...,n. The default output data type of x is 'double'.
%
%   x = randsmpl(p, m, n, classname) returns an m-by-n matrix x whose data  
%   type is specified by classname. The classname must be a valid numeric 
%   class name which includes the following
%     'double' (default) | 'single' | 'int64'  | 'int32' | 'int16' | ...
%      'int8' | 'uint64' | 'uint32' | 'uint16' | 'uint8' | 
%   If classname is not provided, the default 'double' type is used. 
%
%   Remarks:
%       - The main idea is to divide interval [0,1] into K disjoint bins, 
%         each with length proportional to the corresponding probability
%         mass. Then, we draw samples from the unfiorm distribution U(0,1)  
%         and determine the indices of the bins containing those samples.    
%
%       - The histc (introduced in R2006a) and histcounts (introduced in 
%         R2014b) functions determine not only the indices of the bins, but 
%         also the histogram over each bin. In contrast, the discretize 
%         (introduced in R2015a) function is solely for finding the indices
%         of the bins, thus is much faster than histc and histcounts. In
%         addition, discretize is considerably more efficient in terms of 
%         memory consumption compared to histc. Thus, for both performance
%         and memory considerations, we used discretize in this function.  
%
%       - For backward compatibility support, we also provide an interp1 
%         alternative to discretize. The interp1 method is slightly slower 
%         than discretize, but its performance penalty is not significant. 
%         Whether to use the discretize or interp1 approach will be handled
%         automatically in this function based on users' MATLAB version.  
% 
%   See also RAND, RANDI, RANDN, RNG.
%
%   Copyright Peng Liu, Nov. 18, 2015
narginchk(3, 4);
if nargin < 4
    classname = 'double'; % Consistent with MATLAB's default double-precision computation
else
    classname = varargin{:};
    if ~ismember(classname,{'int8','int16','int32','int64','uint8','uint16','uint32','uint64','single','double'})
        error('CLASSNAME input must be a valid numeric class name, for example, ''int32'' or ''double''.');
    end  
end
if ~isvector(p)
    error('Input distribution p must be a vector.')
end
% Error-check of the input distribution
if any(imag(p)) || any(~isfinite(p)) || any(p < 0) || any(p > 1)
    error('The probability elements must be real numbers between 0 and 1.');
elseif abs(sum(p) - 1) > sqrt(eps)
    error('Sum of the probability elements must equal 1.');
end
edges = [0; cumsum(p(:))];
% Deal with floating-point errors due to cumulative sum 
if abs(edges(end) - 1) > sqrt(eps)   
    edges = edges/edges(end);
end
edges(end) = 1 + eps(1);
if verLessThan('matlab','8.5')
    % Code to run in MATLAB R2014b and earlier
    x = cast(interp1(edges,1:length(edges),rand(m,n),'previous'),classname);
else
    % Code to run in MATLAB R2015a and later
    x = cast(discretize(rand(m,n), edges),classname);   
end
end
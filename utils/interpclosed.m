function varargout = interpclosed(x,y,varargin)
% INTERPCLOSED Arc-length interpolation, perimeter and area of 2D closed curves defined by points
%
% xyq = INTERPCLOSED(x,y,tq) Interpolates new data points xyq at given
% positions defined by an arc-length parametrization and the query points
% tq, along the closed curve defined by the points specified by x and y.
% The default method 'spline' is used. tq is a n-by-1  array with its
% elements constrained within [0,1], with 0 being the first point of the
% curve and 1 the last one. 
%
% [len,area] = INTERPCLOSED(x,y) If tq is not specified and two output 
% variables are requested, then only the perimeter and area of the
% interpolated curve are returned. Both outputs are obviously scalars.
%
% [len,area,c] = INTERPCLOSED(x,y) If tq is not specified and three output 
% variables are requested, then the perimeter, area and centroid of the 
% interpolated curve are returned. The centroid is a vector with the
% position as (x,y).
%
% [len,area,c,I] = INTERPCLOSED(x,y) If tq is not specified and four output 
% variables are requested, then the perimeter, area, centroid, and second
% moments of area of the interpolated curve are returned. The second moment
% of area is a vector with three elements: (Ixx, Iyy, Ixy).
%
% pp = INTERPCLOSED(x,y,'pp') If only one output argument is defined and
% the string 'pp' is given as input, the returned variable is the piesewise
% polynomial pp, regardless of the definition of tq.
%
% [xyq,len,area] = INTERPCLOSED(x,y,tq) If tq is provided and there are
% three output variables, the perimeter and area are additionally returned.
%
% [xyq,len,area,c,I] = INTERPCLOSED(x,y,tq) If tq is provided and there are
% five output variables, the perimeter, the area, the centroid, and the
% second moments of area are additionally returned.
%
% [___] = INTERPCLOSED(___,method) By specifying the string method it
% is possible to change how the interpolated points are obtained.
%   'linear': Linear interpolation. The interpolated value at a query point
%             is based on linear interpolation of the values at neighboring
%             points in each respective dimension. This is the fastest
%             method.
%   'spline': Spline interpolation using periodic end conditions. The
%             interpolated value at a query point is based on a cubic
%             interpolation of the values at neighboring points in each
%             respective dimension. This is the default method.
%   'pchip':  Shape-preserving piecewise cubic interpolation. The
%             interpolated value at a query point is based on a shape-
%             preserving piecewise cubic periodic interpolation of the
%             values at neighboring points.
%
% [___] = INTERPCLOSED(___,print) By setting the boolean variable print to
% true, more details about the interpolation can be obtained.
%
% Examples:
%   1) Use the spline method to interpolate at 64 query points within the
%      8 points used to sample the original circle:
%
%       t = linspace(0,2*pi,9);
%       x = sin(t) + 0.2; y = cos(t) - 0.2;
%       [len,area,c,I] = interpclosed(x,y);
%       fprintf(['Perimeter: %4.5f,	Area: %4.5f,\n',...
%           'Centroid: [%4.5f %4.5f],  Iz: %4.5f\n',...
%           'To be compared to\n',...
%           '2*pi:      %4.5f,	pi:   %4.5f,\n',...
%           'Centroid: [%4.5f %4.5f], pi/2: %4.5f\n'],...
%           len,area,c,(I(1)-area*c(2)^2+I(2)-area*c(1)^2),...
%           2*pi,pi,[0.2,-0.2],pi/2)
%
%   2) Get the piecewise polyonomial of a linear interpolation of a given
%      set of points and use the polynomial to make a plot:
%
%       x = [0 .82 .92 0 -.92 -.82]; y = [.66 .9 0 -.83 0 .9];
%       pp = interpclosed(x,y,'pp','linear');
%       tq = min(pp.breaks):0.001:max(pp.breaks);
%       xyq = ppval(pp,tq);
%       figure, plot(xyq(1,:),xyq(2,:))
%
%   Find more examples in the File Exchange website.
%
%   See also CSCVN, PCHIP, MKPP, PPVAL, INTERPARC, ARCLENGTH, INTERP1.
%
% Author:                      Santiago M. Benito
%                              Ruhr-Universität Bochum
% -------------------------------------------------------------------------
% Contact:                     santiago.benito@rub.de
% -------------------------------------------------------------------------
% Current version:             3.0
% -------------------------------------------------------------------------
% Last updated:                17.05.2021
%   Changes:
%       - It is now possible to compute the second moments of area of the 
%       fit with this function.
%% Manage input, output and catch eventual problems.
% Check for errors in the given inputs.
if nargin < 2
    error('INTERPCLOSED:insufficientarguments', ...
        'At least x and y must be supplied.')
end
if ~isvector(x) || ~isvector(y) || (length(x) ~= length(y))
    error('INTERPCLOSED:baddimension', ...
        'x and y must be vectors of the same length.')
end
% Set defaults.
method = 'spline';
print = false;
geomcalc = false; tqgiven = false;
pp = false;
% Initialize output variables.
len = 0;
area = 0;
c = zeros(1,2);
Ixx = 0; Iyy = 0; Ixy = 0;
% Check for other input arguments.
if numel(varargin) > 0
    % At least one other argument was supplied.
    for ii = 1:numel(varargin)
        arg = varargin{ii};
        if ischar(arg)
            % It can be the method or the 'pp'-flag.
            validstrings = {'pp','linear' 'pchip' 'spline'};
            ind = strncmp(arg,validstrings,2);
            if isempty(ind) || (sum(ind) == 0) || (sum(ind) > 1)
                error('INTERPCLOSED:invalidmethod', ...
                    ['Invalid method indicated. Only ''linear'',',...
                    '''pchip'',''spline'' allowed.'])
            end
            if ind(1) == 1
                pp = true;
            else
                method = validstrings{ind>0};
            end
        elseif islogical(arg)
            % It must be the print variable, set the print sampling distance.
            if ~tqgiven, tq = 0:1/32:1; end
            print = arg; 
        else
            % It must be tq, defining the parametric arc-length query
            % points
            if ~isnumeric(arg)
                error('INTERPCLOSED:badtq', ...
                    'tq must be numeric.')
            else
                if max(arg) > 1 || min(arg) < 0
                    error('INTERPCLOSED:badtq', ...
                        'tq elements must be bigger than 0 and smaller than 1.')
                end
                tqgiven = true;
                tq = arg;
            end
        end
    end
end
% If only one output variable is requested and the 'pp' flag was given, no
% need to compute the interpolations, regardless of the definition of tq.
% If three are given, geometry computations will be needed. If two or three
% are given, but tq was not provided, also compute the geometry computations.
if nargout == 1 && pp && ~print
    tqgiven = false;
elseif (nargout == 2 || nargout == 3 || nargout == 4) && ~tqgiven
    geomcalc = true;
elseif (nargout == 3 || nargout == 4 || nargout == 5)
    geomcalc = true;
    if ~tqgiven
        error('INTERPCLOSED:badtq', ...
            'tq was not defined and is needed for interpolation.')
    end
end
% Be sure everything is formatted correclty and group it.
x = x(:)'; y = y(:)';
points = [x;y];
% Round to the 15th decimal position to avoid rounding errors. This is
% necessary for the function to recongnize start and ending points
% properly.
points = round(points,15);
% If the set of points does not describe a closed loop, close it.
if sum(points(:,1) ~= points(:,end)) > 0
    points = [points,points(:,1)];
end
% If less than three distinct points are given, no closed curve can be
% formed.
d = sum((diff(points.').^2).');
if numel(x) - sum(d==0)-1 < 2
    error('INTERPCLOSED:baddimension', ...
        'x and y must be vectors describing at least three distinct points.')
end
%% Actual program start
% Compute the coefficients of the fit-polynomials according to the user's
% choice.
if strcmpi(method,'linear')
    % Remove segments with length equal to zero, the linear interpolation
    % has no continuous derivatives anyway.
    points(:,d==0) = [];
    
    % Compute the linear coefficients of the parametric versions of the
    % lines. First compute the lengths of each segment, then the cumulative
    % length and finally use the slope in each direction to get the coefs.
    seglen = sqrt(sum(diff(points,[],2).^2,1));
    cumarc = [0,cumsum(seglen)];
    coefX = [diff(points(1,:))./diff(cumarc);points(1,1:(end-1))];
    coefY = [diff(points(2,:))./diff(cumarc);points(2,1:(end-1))];
    
    % Create a piecewise polynomial with the given coefficients.
    coefs = zeros(size([coefX,coefY]));
    coefs(:,1:2:end) = coefX;
    coefs(:,2:2:end) = coefY;
    curve = mkpp(cumarc,coefs',2);
    
    % Provide the differentiation array for later use.
    diffarray = [0 0 1;0 0 0];
    
    % Since we already have the lenghts of the individual segments, just
    % sum everything up and save some time.
    len = sum(seglen);
    
elseif strcmpi(method,'spline')
    % MATLAB(R) already has a very useful function that makes all the work
    % for us.
    curve = cscvn(points);
    
    % Provide the differentiation array for later use.
    diffarray = [3 0 0;0 2 0;0 0 1;0 0 0];
    
elseif strcmpi(method,'pchip')
    % Like in the function CSCVN, if the user specified a point where the
    % 2nd derivative is equal to zero, we have to be able to handle the
    % situation.
    d = sum((diff(points.').^2).');
    
    if all(d > 0)
        % The fit is periodic. To have the start and end slopes equal to
        % each other, some tricks must be done. Extra points will be added
        % right before the start and right after the end of the data set. 
        % The fit will be performed with these points, and then the extra
        % pieces will be removed from the general fit.
        %pointsNew = [x(end-2:end-1),x,x(2:3);y(end-2:end-1),y,y(2:3)];
        pointsNew = [points(:,end-2:end-1),points,points(:,2:3)];
        
        % We need the arc length of the modified dataset, therefore we will
        % compute it here.
        seglen = sqrt(sum(diff(pointsNew,[],2).^2,1));
        cumarc = [0,cumsum(seglen)];
        
        % Fit coefficients are obtained from the MATLAB(R) original pchip
        % function.
        temp = pchip(cumarc,pointsNew(1,:)); coefX = temp.coefs;
        temp = pchip(cumarc,pointsNew(2,:)); coefY = temp.coefs;
        
        % Here we remove the unnecesary pieces by removing the extra
        % coefficients.
        coefs = zeros(size([coefX;coefY])-[8,0]);
        coefs(1:2:end,:) = coefX(3:end-2,:);
        coefs(2:2:end,:) = coefY(3:end-2,:);
        
        % Compute the actual arc length
        seglen = sqrt(sum(diff(points,[],2).^2,1));
        cumarc = [0,cumsum(seglen)];
        
    else
        % The 1st derivatives at the end points and at the specified points
        % are not equal, while analysed from both sides. Firstly compute
        % the arclength of the point distribution.
        seglen = sqrt(sum(diff(points,[],2).^2,1));
        cumarc = [0,cumsum(seglen)];    
        
        % Fit coefficients are obtained from the MATLAB(R) original pchip
        % function, according to the desired derivative contiguity.
        dp = find(d>0);
        dpbig = find(diff(dp)>1);
        dpbig = [dpbig,length(dp)];
        idx = dp(1):(dp(dpbig(1))+1);
        temp = pchip(cumarc(idx),points(1,idx)); coefX = temp.coefs;
        temp = pchip(cumarc(idx),points(2,idx)); coefY = temp.coefs;
        for j=2:length(dpbig)
            idx = dp(dpbig(j-1)+1):(dp(dpbig(j))+1);
            temp = pchip(cumarc(idx),points(1,idx));
            coefX = [coefX;temp.coefs];
            temp = pchip(cumarc(idx),points(2,idx));
            coefY = [coefY;temp.coefs];
        end
        
        % Compiling the coefficients in a simple array.
        coefs = zeros(size([coefX;coefY]));
        coefs(1:2:end,:) = coefX(1:end,:);
        coefs(2:2:end,:) = coefY(1:end,:);
        
        % Update the cumulative arclength
        cumarc(:,d==0) = [];
    end
    
    % Finally compute the piecewise polynomial.
    curve = mkpp(cumarc,coefs,2);
    
    % Provide the differentiation array for later use.
    diffarray = [3 0 0;0 2 0;0 0 1;0 0 0];
end
% If tq is given (or a print is required), compute the interpolation using
% the piecewise evaluation function provided in MATLAB(R) and then convert 
% the parametrization into an arc-lenght one.
if tqgiven || print
    step = (max(curve.breaks)-min(curve.breaks))/numel(tq)/30;
    auxtq = min(curve.breaks):step:max(curve.breaks);
    xyq = ppval(curve,auxtq);    
    tqp = pdearcl(auxtq,xyq,tq,0,1);
    xyq = ppval(curve,tqp);
end
% If the geometric parameters (perimeter and area) are required, compute
% them using some calculus.
if geomcalc
    for ii = 1:curve.pieces
        % Get the coefficients of the piecewise polynomial expresions of
        % the parametric form.
        coefX = curve.coefs(2*ii-1,:);
        coefY = curve.coefs(2*ii,:);
        
        % Obtain the derivatives of the polynomials.
        difX = coefX*diffarray;
        difY = coefY*diffarray;
        
        % The length in the linear case is already computed, skip this bit.
        if ~strcmpi(method,'linear')
            % Define the function employed in the arc length and integrate
            % it.
            flen = @(t) sqrt(polyval(difX,t-curve.breaks(ii)).^2 ...
                + polyval(difY,t-curve.breaks(ii)).^2);
            len = len + integral(flen,curve.breaks(ii),curve.breaks(ii+1));
        end
        % The area integral is computed here.
        farea = @(t) polyval(conv(coefY,difX),...
            t-curve.breaks(ii));
        area = area + integral(farea,curve.breaks(ii),...
            curve.breaks(ii+1));
        
        % The centroid is computed here
        fcx = @(t) polyval(conv(coefX,conv(coefY,difX)),...
            t-curve.breaks(ii));
        c(1) = c(1) + integral(fcx,curve.breaks(ii),...
            curve.breaks(ii+1));
        fcy = @(t) polyval(conv(coefY,conv(coefX,difY)),...
            t-curve.breaks(ii));
        c(2) = c(2) - integral(fcy,curve.breaks(ii),...
            curve.breaks(ii+1));
        
        % The area moments of inertia
        fIxx = @(t) polyval(conv(coefY,conv(coefY,conv(coefY,difX))),...
            t-curve.breaks(ii));
        Ixx = Ixx + integral(fIxx,curve.breaks(ii),...
            curve.breaks(ii+1));
        fIyy = @(t) polyval(conv(coefX,conv(coefX,conv(coefX,difY))),...
            t-curve.breaks(ii));
        Iyy = Iyy - integral(fIyy,curve.breaks(ii),...
            curve.breaks(ii+1));
        fIxy = @(t) polyval(conv(coefY,conv(coefY,conv(coefX,difX))),...
            t-curve.breaks(ii));
        Ixy = Ixy + integral(fIxy,curve.breaks(ii),...
            curve.breaks(ii+1));
    end
    c = c / area;
    I = [1/3*Ixx,1/3*Iyy,1/2*Ixy]*sign(area);
    area = abs(area);
    
    
end
%% If required, print some figures to show what the algorithm did.
if print
    figure
    subplot(1,2,1)
    plot(xyq(1,:),xyq(2,:),'*')
    hold on
    plot(points(1,:),points(2,:),'o')
    plot(c(1),c(2),'x')
    xlabel('x'), ylabel('y'), hold off, axis equal
    title('Cartesian representation'), legend('Interpolation','Points',...
        'Centroid')
    subplot(1,2,2)
    plot(tq,xyq), hold on
    for ii = (curve.breaks)/max(curve.breaks)
        line([ii ii],ylim,'LineStyle','--','Color','k')
        line([ii ii],ylim,'LineStyle','--','Color','k')
    end
    hold off, xlim([min(tq) max(tq)]), title('Parametric representation')
    xlabel('t'), ylabel('x(t), y(t)'), legend('x(t)','y(t)')
end
%% Process adequately the variables to be returned.
if nargout == 2 && ~tqgiven
    varargout{1} = len; varargout{2} = area;
elseif nargout == 3 && ~tqgiven
    varargout{1} = len; varargout{2} = area; varargout{3} = c;
elseif nargout == 3
    varargout{1} = xyq;
    varargout{2} = len; varargout{3} = area;
elseif nargout == 4 && tqgiven
    varargout{1} = xyq;
    varargout{2} = len; varargout{3} = area; varargout{4} = c;
elseif nargout == 4 && ~tqgiven
    varargout{1} = len; varargout{2} = area; varargout{3} = c;
    varargout{4} = I;
elseif nargout == 5
    varargout{1} = xyq;
    varargout{2} = len; varargout{3} = area; varargout{4} = c;
    varargout{5} = I;
elseif nargout == 1 && pp
    varargout{1} = curve;
else
    varargout{1} = xyq;
end
end
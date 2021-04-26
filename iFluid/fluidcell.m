classdef fluidcell < matlab.mixin.CustomDisplay

    
properties (Access = public)
    % display properties (not used for storage)
    numel_rapidity1
    numel_rapidity2
    numel_type1
    numel_type2
    numel_space
end
    
    
properties (Access = private)
    
    % Underlying data structure 
    matrix  = [];
   
    % list of index labels (in order)
    labels = {'rapid1', 'space', 'type1', 'rapid2', 'type2'};
end


methods ( Static )
    
    
    function C = eye(rapid_size, type_size)
        % Initialize tensor diagonal in rapidity and type indices
        eye_rapid   = eye(rapid_size);
        eye_type    = eye(type_size);
        
        eye_rapid   = permute( eye_rapid, [1 4 3 2] );
        eye_type    = permute( eye_type, [3 5 1 4 2] );
        
        C           = fluidcell( eye_rapid.*eye_type );
    end
    
    
    function C = zeros(varargin)
        % Initialize tensor of zeros.
        % Size can be specified as a list of comma separated dimensions or
        % an array of dimensions.
        
        if ischar(varargin{1}) % input must be parsed
            
            defaultSize = 1;
            validSize   = @(x) isnumeric(x) && isscalar(x) && (x > 0);

            p = fluidcell.parseIndexLabels(defaultSize, validSize, varargin{:});

            dims        = [ p.Results.rapid1, ...
                            p.Results.space, ...
                            p.Results.type1, ...
                            p.Results.rapid2, ...
                            p.Results.type2 ];       
        
            C = fluidcell( zeros(dims) );
            
        else % input is list of sizes
            C = fluidcell( zeros(varargin{:}) );
        end
    end
    
    
    function C = ones(varargin)
        % Initialize tensor of zeros.
        % Size can be specified as a list of comma separated dimensions or
        % an array of dimensions.
        
        if ischar(varargin{1}) % input must be parsed
            
            defaultSize = 1;
            validSize   = @(x) isnumeric(x) && isscalar(x) && (x > 0);

            p = fluidcell.parseIndexLabels(defaultSize, validSize, varargin{:});

            dims        = [ p.Results.rapid1, ...
                            p.Results.space, ...
                            p.Results.type1, ...
                            p.Results.rapid2, ...
                            p.Results.type2  ];       
        
            C = fluidcell( ones(dims) );
            
        else % input is list of sizes
            C = fluidcell( ones(varargin) );
        end
    end
    
    
    function p = parseIndexLabels(defaultValue,condition,varargin)
        % Parse input based on index labels. If the label is not found in
        % varargin, it will be assigned the value defaultValue. Each passed
        % value can be subjected to a logical condition
        
        % replace aliases
        varargin(strcmpi(varargin,'rapid')) = {'rapid1'};
        varargin(strcmpi(varargin,'rapidity')) = {'rapid1'};
        varargin(strcmpi(varargin,'rapidity1')) = {'rapid1'};
        varargin(strcmpi(varargin,'rapididy2')) = {'rapid2'};
        varargin(strcmpi(varargin,'type')) = {'type1'};
        varargin(strcmpi(varargin,'x')) = {'space'};

        % parse inputs
        p           = inputParser;

        addParameter(p, 'rapid1', defaultValue, condition);
        addParameter(p, 'rapid2', defaultValue, condition);
        addParameter(p, 'type1', defaultValue, condition);
        addParameter(p, 'type2', defaultValue, condition);
        addParameter(p, 'space', defaultValue, condition);

        parse(p, varargin{:});
    end
    
    
    
end % end static methods


methods (Access = protected)
    
    function propgrp = getPropertyGroups(obj)
        % This function controls the appearance of the object properties 
        propList = struct( 'size_rapidity1', size(obj.matrix, 1), ...
                           'size_rapidity2', size(obj.matrix, 4), ...
                           'size_type1', size(obj.matrix, 3), ...
                           'size_type2', size(obj.matrix, 5), ...
                           'size_space', size(obj.matrix, 2)  );
        propgrp = matlab.mixin.util.PropertyGroup(propList);
    end
    
    
    function varargin = convertLabelsToIndexList(obj, varargin)
        % delete all labels from varargin and insert an array of the
        % corresponding indices at the start of varargin
        
        % replace aliases
        varargin(strcmpi(varargin,'rapid')) = {'rapid1'};
        varargin(strcmpi(varargin,'rapidity')) = {'rapid1'};
        varargin(strcmpi(varargin,'rapidity1')) = {'rapid1'};
        varargin(strcmpi(varargin,'rapididy2')) = {'rapid2'};
        varargin(strcmpi(varargin,'type')) = {'type1'};
        varargin(strcmpi(varargin,'x')) = {'space'};
        
        
        idx_list = [];
        for i = 1:length(obj.labels)
            match = strcmpi(varargin, obj.labels{i});
            if any(match)
                idx_list = [idx_list i]; % add to index list
                varargin(match) = ''; % delete from varargin
            end
        end
        
        if ~isempty(idx_list)
            varargin= [{idx_list} varargin];
        end

    end
    
   
end % end protected methods


methods (Access = public)
    
    
    %% Constructor
    function obj = fluidcell( x )
        % Initializes an fluidcell wrapping around the input matrix x.
        if isa(x, 'ifluidcell'); x = double(x); end        
        assert( isnumeric(x), 'Input must be a ND-array' )
        assert( ndims(x) <= 5, 'fluidcell can have max 5 dimensions!' )
        
        obj.matrix = x;
    end
    
    
    %% Accessor methods
    function n = numArgumentsFromSubscript(obj,s,~)
        % This function overload is necessary for obj{i,j} to work.
        % Otherwise MATLAB expects multiple outputs for some reason.
        n = 1;
    end
    
    
    function ind = end(obj,k,n)
        % Overload the end method, which returns the final index of the
        % desired dimension.       
       
        error('Indexing using the end statement is currently not supported')
       
    end
   
    
    function varargout = subsref(obj,s)
        switch s(1).type
        case '.'
            % Use built-in for any other expression
            [varargout{1:nargout}] = builtin('subsref',obj,s);
        case '()'
            if length(s) == 1
                
                if ischar(s.subs{1}) && all(s.subs{1} ~= ':')
                    % Reference via labels (via parsing)
                    defaultSubs = ':';
                    condition = @(x) true; % no condition
                    
                    p = fluidcell.parseIndexLabels(defaultSubs,condition,s.subs{:});
 
                    % create new subsref struct
                    S_new.subs      = {p.Results.rapid1, ...
                                       p.Results.space, ...
                                       p.Results.type1, ...
                                       p.Results.rapid2, ...
                                       p.Results.type2 };
                    S_new.type       = '()';
                    
                    % call built-in subsref on obj.matrix
                    x = builtin('subsref', obj.matrix, S_new);
                    varargout{1} = fluidcell(x);

                else 
                    % direct array referencing (via underlying indices)
                    x = builtin('subsref', obj.matrix, s);
                    varargout{1} = fluidcell(x);
                end
            else
                % Use built-in for any other expression
                [varargout{1:nargout}] = builtin('subsref', obj.matrix, s);
            end
        case '{}'
            % Use built-in for any other expression
            [varargout{1:nargout}] = builtin('subsref',obj,s);
        otherwise
            error('Not a valid indexing expression')
        end
    end
   
    
    function obj = subsasgn(obj,s,varargin)

        switch s(1).type
        case '.'
            % Call built-in subsasgn
            obj = builtin('subsasgn',obj,s,varargin{:});
        case '()'
            if length(s) == 1
                % Implement obj(indices) = varargin{:};
                
                if ischar(s.subs{1}) && all(s.subs{1} ~= ':')
                    assert(length(varargin)==1, 'Number of arguments on right side is too many.')
                    
                    % Assign via labels (via parsing)
                    defaultSubs = ':';
                    condition = @(x) true; % no condition
                    
                    p = fluidcell.parseIndexLabels(defaultSubs,condition,s.subs{:});
 
                    % create new subsasgn struct
                    S_new.subs      = {p.Results.rapid1, ...
                                       p.Results.space, ...
                                       p.Results.type1, ...
                                       p.Results.rapid2, ...
                                       p.Results.type2 };
                    S_new.type       = '()';
                    
                    % call built-in subsref on obj.matrix
                    obj.matrix      = subsasgn( obj.matrix, S_new, double(varargin{1}) );

                else 
                    % direct array assignment (via underlying indices)
                    obj.matrix = builtin('subsasgn',obj.matrix,s,varargin{:});
                end
            else
                % Use built-in for any other expression
                obj = builtin('subsasgn',obj,s,varargin{:});
            end       
        case '{}'       
            % Use built-in for any other expression
            obj = builtin('subsasgn',obj,s,varargin{:});
        otherwise
            error('Not a valid indexing expression')
        end
    end
    
    
    function S = size(obj, varargin)
        % Returns size of underlying matrix
        if nargin == 1
            S = size(obj.matrix);
            S = [S ones(1, 5-length(S))]; % include singletons   
        else
            % dim is either the index or the label of the dimension
            varargin = obj.convertLabelsToIndexList(varargin{:});
            
            S = size(obj.matrix, varargin{:});
        end
    end
    
    
    function L = length(obj)
        % Returns length of underlying matrix
        L = length(obj.matrix);
    end
    
    
    function C = double(obj)
        % Cast fluidcell to double, i.e. return underlying matrix
        C = obj.matrix;
    end  
    
    
    function obj = transpose(obj)
        % Permutes the two rapidity indices and the two type indices
        obj.matrix = permute(obj.matrix, [4 2 5 1 3] );
    end
    
    
    function obj = t(obj)
        % Permutes the two rapidity indices and the two type indices
        obj.matrix = permute(obj.matrix, [4 2 5 1 3] );
    end
    
    
    function C = flatten(obj)
        % Rehapes tensor into 2D- (3D with spacial index) matrix, by
        % combining indices (rapid1, type1) into a single, first index and 
        % indices (rapid2, type2) into a single, second index.
        S   = size(obj);
        C   = reshape(permute(obj.matrix, [1 3 4 5 2]), S(1)*S(3), S(4)*S(5), S(2) );
    end
   
    
    %% Standard matrix operator overloads    
    function C = plus(A,B)       
        % overloads A+B operation
        C = fluidcell(plus(double(A),double(B)));
    end
    
    
    function C = minus(A,B)       
        % overloads A-B operation
        C = fluidcell(minus(double(A),double(B)));
    end
    
    
    function obj = uplus(obj)
        % overloads +A operation
        obj.matrix = uplus(obj.matrix);
    end
    
    
    function obj = uminus(obj)
        % overloads -A operation
        obj.matrix = uminus(obj.matrix);
    end
    
    
    function A = power(A,B)
        % Overloads A.^B operation
        assert( isscalar(B) && isnumeric(B) );
        
        A.matrix = power(A.matrix, B);
    end
    
    
    function C = times(A,B)       
        % overloads A.*B operation
        C = fluidcell(times(double(A),double(B)));
    end
    
    
    function C = rdivide(A,B)
        % Overloads A./B operation
        C = fluidcell(rdivide(double(A),double(B)));
    end
    
    
    function C = mrdivide(A,B)
        % Overloads A/B operation
        C = fluidcell(mrdivide(double(A),double(B)));
    end
    
    
    function C = mtimes(A,B)
        % Overloads A*B
        % Works only for fluidcells and scalars
        
        isa_scalar = @(x) isscalar(x) && isnumeric(x);
        
        if isa_scalar(A) && isa(B,'fluidcell')
            C = fluidcell(mtimes(double(A),double(B)));
            
        elseif isa_scalar(B) && isa(A,'fluidcell')
            C = fluidcell(mtimes(double(A),double(B)));
            
        elseif isa(A,'fluidcell') && isa(B,'fluidcell')
            % Check that sizes match
            SA      = size(A);
            SB      = size(B);
            
            assert( SA(2) == SB(2) || SA(2) == 1 || SB(2) == 1, 'Sizes of space indices do not match')
            
            
            % Compine rapid1 and type1 into first index. 
            % Combine rapid2 and type2 into second index.
            a       = reshape(permute(A.matrix, [1 3 4 5 2]), SA(1)*SA(3), SA(4)*SA(5), SA(2) );
            b       = reshape(permute(B.matrix, [1 3 4 5 2]), SB(1)*SB(3), SB(4)*SB(5), SB(2) );
            
            clear A B % free up some memory
            
            % Make matrices have an equal number of space slices
            Dx      = max(size(a,3),size(b,3));            
            a       = repmat(a, [1, 1, Dx/size(a,3)]);
            b       = repmat(b, [1, 1, Dx/size(b,3)]);
            
            
            % Multiply a and b for each space slice
            if verLessThan('MATLAB','9.9')
                xsize   = [ size(a,1), size(b,2), Dx ];
                x       = zeros(xsize);
                for i = 1:Dx
                    x(:,:,i) = mtimes(a(:,:,i),b(:,:,i));
                end
            else
                x   = pagemtimes(a, b);
            end
            
            % First split the type and rapidity indices by reshaping, then
            % permute to right positions
            SC      = [SA(1), SA(3), SB(4), SB(5), Dx];
            C       = fluidcell( permute(reshape( x, SC ), [1 5 2 3 4]) );
            
        else
            error('Types not supported for this operation')
        end
    end
        
    
    function C = mldivide(A,B)
        % Calculates A\B i.e. solves system of linear equations

        assert( isa( A ,'fluidcell' ) );
        assert( isa( B ,'fluidcell' ) );
        
        % Check that sizes match
        SA      = size(A);
        SB      = size(B);

        assert( SA(2) == SB(2) || SA(2) == 1 || SB(2) == 1, 'Sizes of space indices do not match')

        % Compine rapid1 and type1 into first index. 
        % Combine rapid2 and type2 into second index.
        a       = reshape(permute(A.matrix, [1 3 4 5 2]), SA(1)*SA(3), SA(4)*SA(5), SA(2) );
        b       = reshape(permute(B.matrix, [1 3 4 5 2]), SB(1)*SB(3), SB(4)*SB(5), SB(2) );

        clear A B % free up some memory

        % Make matrices have an equal number of space slices
        Dx      = max(size(a,3),size(b,3));
        a       = repmat(a, [1, 1, Dx/size(a,3)]);
        b       = repmat(b, [1, 1, Dx/size(b,3)]);


        % Solve ax=b for each space slice
        xsize   = [ size(a,2), size(b,2), Dx ];
        x       = zeros(xsize);
        for i = 1:Dx
            x(:,:,i) = mldivide(a(:,:,i),b(:,:,i));
        end

        % First split the type and rapidity indices by reshaping, then
        % permute to right positions
        SC      = [ SA(4), SA(5), SB(4), SB(5), Dx];
        C       = fluidcell( permute(reshape( x, SC ), [1 5 2 3 4]) );
    end
    
    
    %% Overload Matlab logical operators
    function C = lt(A,B)
        % Overloads < operator
        C = lt( double(A), double(B) );
    end
    
    function C = gt(A,B)
        % Overloads > operator
        C = gt( double(A), double(B) );
    end
    
    function C = le(A,B)
        % Overloads <= operator
        C = le( double(A), double(B) );
    end
    
    function C = ge(A,B)
        % Overloads >= operator
        C = ge( double(A), double(B) );
    end
    
    function C = ne(A,B)
        % Overloads ~= operator
        C = ne( double(A), double(B) );
    end
    
    function C = eq(A,B)
        % Overloads == operator
        C = eq( double(A), double(B) );
    end
    
    function C = and(A,B)
        % Overloads & operator
        C = and( double(A), double(B) );
    end
    
    function C = or(A,B)
        % Overloads | operator
        C = or( double(A), double(B) );
    end
    
    function C = not(A,B)
        % Overloads ~ operator
        C = not( double(A), double(B) );
    end
    
    %% Overload Matlab matrix operations 
    function obj = inv(obj)
        % Solve Uinv*U = I
        I = fluidcell.eye(size(obj.matrix, 1), size(obj.matrix, 3));
        obj = obj\I;
    end
    
    
    function obj = exp(obj)
        % Takes elementwise exponential
        obj.matrix = exp(obj.matrix);
    end
    
    
    function obj = log(obj)
        % Takes elementwise logarithm
        obj.matrix = log(obj.matrix);
    end
    
    
    function obj = abs(obj)
        % Takes elementwise absolute value
        obj.matrix = abs(obj.matrix);
    end
    
    
    function obj = sqrt(obj)
        % Takes elementwise absolute value
        obj.matrix = sqrt(obj.matrix);
    end
    
    
    function obj = cos(obj)
        obj.matrix = cos(obj.matrix);
    end
    
    
    function obj = acos(obj)
        obj.matrix = acos(obj.matrix);
    end

    
    function obj = cosh(obj)
        obj.matrix = cosh(obj.matrix);
    end
    
    
    function obj = acosh(obj)
        obj.matrix = acosh(obj.matrix);
    end
    

    function obj = sin(obj)
        obj.matrix = sin(obj.matrix);
    end
    
    
    function obj = asin(obj)
        obj.matrix = asin(obj.matrix);
    end

    
    function obj = sinh(obj)
        obj.matrix = sinh(obj.matrix);
    end
    
    
    function obj = asinh(obj)
        obj.matrix = asinh(obj.matrix);
    end
    
    
    function obj = tan(obj)
        obj.matrix = tan(obj.matrix);
    end
    
    
    function obj = atan(obj)
        obj.matrix = atan(obj.matrix);
    end

    
    function obj = tanh(obj)
        obj.matrix = tanh(obj.matrix);
    end
    
    
    function obj = atanh(obj)
        obj.matrix = atanh(obj.matrix);
    end
    
    
    function obj = csc(obj)
        obj.matrix = csc(obj.matrix);
    end
    
    
    function obj = acsc(obj)
        obj.matrix = acsc(obj.matrix);
    end

    
    function obj = csch(obj)
        obj.matrix = csch(obj.matrix);
    end
    
    
    function obj = acsch(obj)
        obj.matrix = acsch(obj.matrix);
    end
    
    
    function obj = sec(obj)
        obj.matrix = sec(obj.matrix);
    end
    
    
    function obj = asec(obj)
        obj.matrix = asec(obj.matrix);
    end

    
    function obj = sech(obj)
        obj.matrix = sech(obj.matrix);
    end
    
    
    function obj = asech(obj)
        obj.matrix = asech(obj.matrix);
    end
    
    
    function obj = cot(obj)
        obj.matrix = cot(obj.matrix);
    end
    
    
    function obj = acot(obj)
        obj.matrix = acot(obj.matrix);
    end

    
    function obj = coth(obj)
        obj.matrix = coth(obj.matrix);
    end
    
    
    function obj = acoth(obj)
        obj.matrix = acoth(obj.matrix);
    end
    
    
    %% Overload other ND-array functions
    
    function obj = sum(obj, varargin)
        varargin = obj.convertLabelsToIndexList(varargin{:});
        obj.matrix = sum(obj.matrix, varargin{:});
    end
    
    function obj = cumsum(obj, varargin)
        varargin = obj.convertLabelsToIndexList(varargin{:});
        obj.matrix = cumsum(obj.matrix, varargin{:});
    end
        
    function obj = prod(obj, varargin)
        varargin = obj.convertLabelsToIndexList(varargin{:});
        obj.matrix = prod(obj.matrix, varargin{:});
    end
    
    function obj = gradient(obj, varargin)
        % First argument of varargin should be index OR label of the
        % desired direction of differentiation.
        % Second (optional) argument should be the grid OR grid spacing.
        varargin = obj.convertLabelsToIndexList(varargin{:});
        
        % Matlabs gradient function takes the derivative along the x
        % (horizontal) direction. Therefore, we must permute the desired
        % direction to the second index.
        idx             = 1:5;
        idx(varargin{1})= 2;
        idx(2)          = varargin{1};
        
        if length(varargin) > 1
            dx      = num2cell(1:length(size(obj.matrix)));
            dx{1}   = varargin{2};
        else
            dx      = {1};
        end
        
        obj.matrix = permute(gradient(permute(obj.matrix,idx),dx{:}),idx);
    end
    
    
    
    %% Old getter functions
    function C = getX(obj, x_idx, t_str)
        % Returns the iFluidTensor at spatial index x_idx
        mat = obj.matrix(:,x_idx,:,:,:);
        
        % Decide if output is double or iFluidTensor
        if nargin == 3
            if strcmp(t_str, 'd') || strcmp(t_str, 'double')
                C   = mat;
            elseif strcmp(t_str, 't') || strcmp(t_str, 'tensor')
                C   = fluidcell(mat);
            else
                error('Output type not recognized!')
            end
        else
            C   = fluidcell(mat); 
        end
    end
    
    function C = getRapid(obj, rapid_idx, t_str)
        % Returns the iFluidTensor at rapidity index rapid_idx
        mat = obj.matrix(rapid_idx,:,:,:,:);
    
        % Decide if output is double or iFluidTensor
        if nargin == 3
            if strcmp(t_str, 'd') || strcmp(t_str, 'double')
                C   = mat;
            elseif strcmp(t_str, 't') || strcmp(t_str, 'tensor')
                C   = fluidcell(mat);
            else
                error('Output type not recognized!')
            end
        else
            C   = fluidcell(mat); 
        end
    end
    
    function C = getType(obj, type_idx, t_str)
        % Returns the iFluidTensor at type index type_idx
        mat = obj.matrix(:,:,type_idx,:,:);
    
        % Decide if output is double or iFluidTensor
        if nargin == 3
            if strcmp(t_str, 'd') || strcmp(t_str, 'double')
                C   = mat;
            elseif strcmp(t_str, 't') || strcmp(t_str, 'tensor')
                C   = fluidcell(mat);
            else
                error('Output type not recognized!')
            end
        else
            C   = fluidcell(mat); 
        end
    end
    
    function C = getSecRapid(obj, rapid2_idx, t_str)
        % Returns the iFluidTensor at auxillary rapidity index rapid2_idx
        mat = obj.matrix(:,:,:,rapid2_idx,:);
    
        % Decide if output is double or iFluidTensor
        if nargin == 3
            if strcmp(t_str, 'd') || strcmp(t_str, 'double')
                C   = mat;
            elseif strcmp(t_str, 't') || strcmp(t_str, 'tensor')
                C   = fluidcell(mat);
            else
                error('Output type not recognized!')
            end
        else
            C   = fluidcell(mat); 
        end
    end
    
    function C = getSecType(obj, type2_idx, t_str)
        % Returns the iFluidTensor at auxillary type index type2_idx
        mat = obj.matrix(:,:,:,:,type2_idx);
    
        % Decide if output is double or iFluidTensor
        if nargin == 3
            if strcmp(t_str, 'd') || strcmp(t_str, 'double')
                C   = mat;
            elseif strcmp(t_str, 't') || strcmp(t_str, 'tensor')
                C   = fluidcell(mat);
            else
                error('Output type not recognized!')
            end
        else
            C   = fluidcell(mat); 
        end
    end
    
        
end % end public methods



    
end % end classdef
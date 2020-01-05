classdef iFluidTensor < handle
    % Main data structure of iFluid to keep track of the indices. The class
    % wraps around a 5-dimensional matrix and overlaods all the standard
    % matrix operations to generalize them to iFluids index structure.
    %
    % ======================== Index convention ===========================
    %   Index 1: main rapidity index
    %   Index 2: auxiliary rapidity index (used for convolutions)
    %   Index 3: main type index
    %   Index 4: auxiliary type index (used for convolutions)
    %   Index 5: spatial index (leave unused for homogeneous quantities)
    %
    % ========================= Initialization ============================
    %   Via matrix: iFluidTensor( mat )
    %       - initializes a tensor with mat as underlying matrix
    %
    %   Via 0 scalar: iFluidTensor( 0 )
    %       - initializes a tensor with mat as underlying matrix
    %
    %   Via specification of index sizes: iFluidTensor( 4 , 2 , 3 )
    %       - initializes a size 4x2x3 tensor of zeroes
    %
    
properties (Access = private)
    
    % Underlying data structure 
    matrix  = [];
    
    % Size of indices 
    dr1     = []
    dr2     = [] 
    dt1     = [] 
    dt2     = [] 
    dx      = [] 
end


methods (Access = public)
    
    %$ Constructor
    function obj = iFluidTensor(varargin)
        % Determine whether input is a matrix of list of dimensions
        if length(varargin) == 1 && ~isscalar(varargin{1}) % input is matrix
            obj.matrix = varargin{1};
        elseif length(varargin) == 1 && varargin{1} == 0 % input is the scalar 0
            obj.matrix = zeros(1);
        else
            obj.matrix = zeros(varargin{:});
        end
        
        [d1, d2, d3, d4, d5] = size(obj.matrix);
        obj.dr1 = d1;
        obj.dx  = d2;
        obj.dt1 = d3;
        obj.dr2 = d4;
        obj.dt2 = d5;
    end
    
    
    %% Accessor methods
    function C = getX(obj, x_idx)
        % Returns the iFluidTensor at spatial index x_idx
        mat = obj.matrix(:,x_idx,:,:,:);
        C   = iFluidTensor(mat); 
    end
    
    function C = getRapid(obj, rapid_idx)
        % Returns the iFluidTensor at rapidity index rapid_idx
        mat = obj.matrix(rapid_idx,:,:,:,:);
        C   = iFluidTensor(mat); 
    end
    
    function C = getType(obj, type_idx)
        % Returns the iFluidTensor at type index type_idx
        mat = obj.matrix(:,:,type_idx,:,:);
        C   = iFluidTensor(mat); 
    end
    
    function C = getAuxRapid(obj, rapid2_idx)
        % Returns the iFluidTensor at auxillary rapidity index rapid2_idx
        mat = obj.matrix(:,:,:,rapid2_idx,:);
        C   = iFluidTensor(mat); 
    end
    
    function C = getAuxType(obj, type2_idx)
        % Returns the iFluidTensor at auxillary type index type2_idx
        mat = obj.matrix(:,:,:,:,type2_idx);
        C   = iFluidTensor(mat); 
    end
    
    function C = size(obj, dim)
        % Returns size of underlying matrix
        if nargin == 1
            C = size(double(obj));
        else
            C = size(double(obj),dim);
        end
    end
    
    function C = double(obj)
        % Cast iFluidTensor to double, i.e. return underlying matrix
        C = obj.matrix;
    end
    
    
    %% iTensorFluid index manioulations    
    
    function C = flatten(obj)
        % Rehapes tensor into 2D- (3D with spacial index) matrix, by
        % combining indices (rapid1, type1) into a single, first index and 
        % indices (rapid2, type2) into a single, second index.
        C = reshape(permute(obj.matrix, [1 3 4 5 2]), obj.dr1*obj.dt1, obj.dr2*obj.dt2, obj.dx );
    end
    
    
    function C = unflatten(obj, A)
        % Reverse operation of the flatten.
        C = permute(reshape( A, size(obj.matrix) ), [1 3 4 5 2]);
    end        
    
    
    function obj = setIdentity(obj)
        % Sets the tensor to identity in both rapidity and type indices
        if obj.dr1 ~= obj.dr2 && obj.dt1 ~= obj.dt2
            error('Tensor is not square!')
        end
        
        I_rapid = repmat(eye(obj.dr1), 1, 1, 1, 1, 1);
        I_type  = repmat(eye(obj.dt1), 1, 1, 1, 1, 1);
        
        I_rapid = permute( I_rapid, [1 4 3 2] );
        I_type  = permute( I_type, [3 5 1 4 2] );
        x       = I_rapid.*I_type;
         
        x       = repmat(x, 1, obj.dx, 1, 1, 1);
        obj.matrix = x;
    end
    
    
    %% Standard matrix operator overloads    
    function C = plus(A,B)
        a = double(A);
        b = double(B);
        x = plus(a,b);
        
        C = iFluidTensor(x);
    end
    
    
    function C = minus(A,B)
        a = double(A);
        b = double(B);
        x = minus(a,b);
        
        C = iFluidTensor(x);
    end
    
    
    function C = uplus(obj)
        x = uplus(obj.matrix);
        
        C = iFluidTensor(x);
    end
    
    
    function C = uminus(obj)
        x = uminus(obj.matrix);
        
        C = iFluidTensor(x);
    end
    
    
    function C = power(A,B)
        % Overloads A.^B
        a = double(A);
        b = double(B);
        x = power(a,b);
        
        C = iFluidTensor(x);
    end
    
    
    function C = times(A,B)
        a = double(A);
        b = double(B);
        x = times(a,b); 
        
        C = iFluidTensor(x);
    end
    
    
    function C = mtimes(A,B)
        % Overloads A*B
        
        if isa(A,'double') || isa(B,'double')
            % If either A or B is a matrix perform "regular" matrix
            % multiplixation
            a = double(A);
            b = double(B);
            x = mtimes(a,b);
        else
            % If both A and B are iFluidTensor contract over 2nd and 4th
            % index
            assert( isa( A ,'iFluidTensor' ) );
            assert( isa( B ,'iFluidTensor' ) );
        
            a   = flatten(A);
            b   = flatten(B);

            xsize = [ size(a,1), size(b,2), max(size(a,3),size(b,3)) ];
            x   = zeros(xsize);

            if size(a,3) == 1 && size(b,3) == 1 % both homogeneous
                x = mtimes(a,b);
            elseif size(a,3) == 1 % A homogeneous
                for i = 1:size(b,3)
                    x(:,:,i) = mtimes(a,b(:,:,i));
                end
            elseif size(b,3) == 1 % B homogeneous
                for i = 1:size(a,3)
                    x(:,:,i) = mtimes(a(:,:,i),b);
                end
            else % neither homogeneous
                assert( size(a,3) == size(b,3) )
                for i = 1:size(a,3)
                    x(:,:,i) = mtimes(a(:,:,i),b(:,:,i));
                end
            end
            
            % Reverse the flatten()
            size_C      = [ size(A,1), size(B,4), size(A,3), size(B,5), max(size(A,2),size(B,2))];
            size_C([2 3])= size_C([3 2]); % make sure to permute indices, to reverse flatten

            % First split the type and rapidity indices by reshaping, then
            % permute to right positions
            x           = permute(reshape( x, size_C ), [1 5 2 3 4]);
        end
        
        C = iFluidTensor(x);
    end
    
    
    function C = rdivide(A,B)
        % Overloads A./B
        a = double(A);
        b = double(B);
        x = rdivide(a,b);
        
        C = iFluidTensor(x);
    end
    
    
    function C = mrdivide(A,B)
        % Overloads A/B 
        a = double(A);
        b = double(B);
        x = mrdivide(a,b);
        
        C = iFluidTensor(x);
    end
    
    
    function C = mldivide(A,B)
        % Calculates A\B i.e. solves system of linear equations

        assert( isa( A ,'iFluidTensor' ) );
        assert( isa( B ,'iFluidTensor' ) );
        
        % Flattens the tensors to create a system of equations for each
        % type
        a = flatten(A);
        b = flatten(B);
         
        xsize = [ size(a,2), size(b,2), max(size(a,3),size(b,3)) ];
        x   = zeros(xsize);
        
        if size(a,3) == 1 && size(b,3) == 1 % both homogeneous
            x = mldivide(a,b);
        elseif size(a,3) == 1 % A homogeneous
            for i = 1:size(b,3)
                x(:,:,i) = mldivide(a,b(:,:,i));
            end
        elseif size(b,3) == 1 % B homogeneous
            for i = 1:size(a,3)
                x(:,:,i) = mldivide(a(:,:,i),b);
            end
        else % neither homogeneous
            assert( size(a,3) == size(b,3) )
            for i = 1:size(a,3)
                x(:,:,i) = mldivide(a(:,:,i),b(:,:,i));
            end
        end
        
        % Reverse the flatten()
        size_C      = [ size(A,4), size(B,4), size(A,5), size(B,5), max(size(A,2),size(B,2))];
        size_C([2 3])= size_C([3 2]); % make sure to permute indices, to reverse flatten
            
        % First split the type and rapidity indices by reshaping, then
        % permute to right positions
        x           = permute(reshape( x, size_C ), [1 5 2 3 4]);
                    
        C           = iFluidTensor(x);
    end
    
    
    %% Overloads Matlab matrix operations 
    function C = inv(obj)
        % Solve Uinv*U = I
        I = iFluidTensor(obj.dr1, obj.dx, obj.dt1, obj.dr1, obj.dt1).setIdentity();
        
        C = obj\I;
    end
    
    
    function C = exp(obj)
        % Takes elementwise exponential
        x = double(obj);
        x = exp(x);
        C = iFluidTensor(x);
    end
    
    
    function C = log(obj)
        % Takes elementwise logarithm
        x = double(obj);
        x = log(x);
        C = iFluidTensor(x);
    end
    
    
    function C = abs(obj)
        % Takes elementwise absolute value
        x = double(obj);
        C = iFluidTensor(abs(x));
    end
    
    
    function C = transpose(obj)
        % Permutes the two rapidity indices and the two type indices
        x = double(obj);
        x = permute(x, [4 2 5 1 3] );
        C = iFluidTensor(x);
    end
    
    
    function C = t(obj)
        % Permutes the two rapidity indices and the two type indices
        x = double(obj);
        x = permute(x, [4 2 5 1 3] );
        C = iFluidTensor(x);
    end
    
    
end
    
end % end classdef
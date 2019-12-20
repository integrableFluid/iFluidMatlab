classdef iFluidTensor < handle
    
properties (Access = private)
    
    % Underlying data structure 
    matrix = [];
    
    % dimensions
    d1 = [] % main rapidity index
    d2 = [] % auxiliary rapidity index
    d3 = [] % main species index
    d4 = [] % auxiliary species index
    d5 = [] % spatial index
end


methods (Access = public)
    % Constructor
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
        obj.d1 = d1;
        obj.d2 = d2;
        obj.d3 = d3;
        obj.d4 = d4;
        obj.d5 = d5;
    end
    
    
    function C = atSpacePoint(obj, x)
        % Gets tensor at certain spatial coordinate x
        mat = obj.matrix(:,:,:,:,x);
        C   = iFluidTensor(repmat(mat, 1, 1, 1, 1)); % makes sure all indices are preserved
    end
    
    
    function C = flatten(obj)
        % Rehapes tensor into 2D- (3D with spacial index) matrix, by
        % combining indices (rapid1, type1) into a single, first index and 
        % indices (rapid2, type2) into a single, second index.
        C = reshape(permute(obj.matrix, [1 3 2 4 5]), obj.d1*obj.d3, obj.d2*obj.d4, obj.d5 );
    end
    
    
    function C = unflatten(obj, A)
        % Reverse operation of the flatten.
        C = permute(reshape( A, size(obj.matrix) ), [1 3 2 4 5]);
    end
    
    
    function C = plt(obj)
        % Returns a plot-able 3D-matrix
        % 1st index is rapid
        % 2nd index is space
        % 3rd index is type
        
        if obj.d2 ~= 1 || obj.d4 ~= 1
            error( 'plt() only works for tensors with up to 3 indices!' )
        end
        
        C = squeeze(permute(obj.matrix, [1 5 3 2 4]));
    end
    
    
    function C = size(obj, dim)
        if nargin == 1
            C = size(double(obj));
        else
            C = size(double(obj),dim);
        end
    end
    
    
    function C = double(obj)
        % Cast iFluidTensor to doulbe, i.e. return underlying matrix
        C = obj.matrix;
    end
    
    
    function C = inv(obj)
        % Solve Uinv*U = I
        I_rapid = eye(obj.d1);
        I_type  = repmat(eye(obj.d3), 1 ,1, 1, obj.d5);
        I_type  = permute(I_type, [3 4 1 2 5]);
        identity= iFluidTensor(I_rapid.*I_type);
        
        C       = obj\identity;
    end
    
    
    function C = exp(obj)
        x = double(obj);
        x = exp(x);
        C = iFluidTensor(x);
    end
    
    
    function C = log(obj)
        x = double(obj);
        x = log(x);
        C = iFluidTensor(x);
    end
    
    
    function C = abs(obj)
        x = double(obj);
        C = iFluidTensor(abs(x));
    end
    
    
    function C = transpose(obj)
        % Permutes the two rapidity indices and the two type indices
        % seperately
        x = double(obj);
        x = permute(x,[2 1 4 3 5]);
        C = iFluidTensor(x);
    end
    
    
    function obj = setIdentity(obj)
        if obj.d1 ~= obj.d2 && obj.d3 ~= obj.d4
            error('Tensor is not square!')
        end
        
        I_rapid = eye(obj.d1);
        I_type  = repmat(eye(obj.d3), 1 ,1, 1, 1);
        I_type  = permute(I_type, [3 4 1 2]);
        x       = I_rapid.*I_type;
         
        x       = repmat(x , 1, 1, 1, 1, obj.d5);
        obj.matrix = x;
    end
    
    
    %% Matrix operation overloads
    function B = subsref(obj,S)
        B = subsref(obj.matrix, S);
    end
    
    
    function obj = subsasgn(obj,S,B)
        newMat = subsasgn(obj.matrix,S,B);
        obj.matrix = newMat;
    end
    
    
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
            a = double(A);
            b = double(B);
            x = mtimes(a,b);
        else
            % Assumes both objects are GHD tensors
        
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
            size_C      = [ size(A,1), size(B,2), size(A,3), size(B,4), max(size(A,5),size(B,5))];
            size_C([2 3])= size_C([3 2]); % make sure to permute indices, to reverse flatten

            % First split the type and rapidity indices by reshaping, then
            % permute to right positions
            x           = permute(reshape( x, size_C ), [1 3 2 4 5]);
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
        % Calculates A\B
        % This function assumes both objects are GHD tensors
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
        size_C      = [ size(A,2), size(B,2), size(A,4), size(B,4), max(size(A,5),size(B,5))];
        size_C([2 3])= size_C([3 2]); % make sure to permute indices, to reverse flatten
            
        % First split the type and rapidity indices by reshaping, then
        % permute to right positions
        x           = permute(reshape( x, size_C ), [1 3 2 4 5]);
                    
        C           = iFluidTensor(x);
    end
    
    
    
end
    
end % end classdef
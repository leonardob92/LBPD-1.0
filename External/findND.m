function varargout=findND(X,varargin)
%Find non-zero elements in ND-arrays. Replicates all behavior from find.
% The syntax is equivalent to the built-in find, but extended to
% multi-dimensional input.
%
% [...] = findND(X,K) returns at most the first K indices. K must be a
% positive  scalar of any type.
%
% [...] = findND(X,K,side) returns either the first K or the last K
% inidices. The input side  must be a char, either 'first' or 'last'. The
% default behavior is 'first'.
%
% [I1,I2,I3,...,In] = findND(X,...) returns indices along all the
% dimensions of X.
%
% [I1,I2,I3,...,In,V] = findND(X,...) returns indices along all the
% dimensions of X, and additionally returns a vector containg the values.
%
% Note for Matlab 6.5:
% The syntax with more than one input is present in the online doc for R14
% (Matlab 7.0), so this might be the latest release without support for
% this syntax.
%
% Compatibility:
% Matlab: should work on all releases (tested on R2017b, R2012b and 6.5)
% Octave: tested on 4.2.1
% OS:     should work cross-platform
%
% Version: 1.1
% Date:    2017-12-29
% Author:  H.J. Wisselink
% Email=  'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})
%
% Changes from v1.0 to v1.1:
% - Added support for Matlab 6.5 (R13).
% - Fixed a minor bug where the orientation of the output vector did not
%   match the orientation of the built-in find function. It also always
%   returned 1 as value, instead of the true value.
% - Added support upper case third input.

%Parse inputs
%validateattributes(X,{'numeric','logical'},{'nonempty'},1)
if ~(isnumeric(X) || islogical(X)) || numel(X)==0
    error('Expected first input (X) to be a non-empty numeric or logical array.')
end
switch nargin
    case 1
        %[...] = findND(X);
        side='first';
        K=inf;
    case 2
        %[...] = findND(X,K);
        side='first';
        K=varargin{1};
        %validateattributes(K,{'numeric','logical'},{'scalar','positive'},2)
        if ~(isnumeric(K) || islogical(K)) || numel(K)~=1 || any(K<0)
            error('Expected second input (K) to be a positive numeric or logical scalar.')
        end
    case 3
        %[...] = FIND(X,K,'first');
        K=varargin{1};
        %validateattributes(K,{'numeric','logical'},{'scalar','positive'},2)
        if ~(isnumeric(K) || islogical(K)) || numel(K)~=1 || any(K<0)
            error('Expected second input (K) to be a positive numeric or logical scalar.')
        end
        side=varargin{2};
        if ~isa(side,'char') ||...
                ~( strcmpi(side,'first') || strcmpi(side,'last'))
            error('Third input must be either ''first'' or ''last''.')
        end
        side=lower(side);
    otherwise
        error('Incorrent number of inputs.')
end

%parse outputs
nDims=length(size(X));
%allowed outputs: 0, 1, nDims, nDims+1
if nargout>1 && nargout<nDims
    error('Incorrect number of output arguments.')
end

varargout=cell(nargout,1);
v=version;v=str2double(v(1:3));
if v<7
    %The find(X,k,side) syntax was introduced between 6.5 and 7
    if nargout>nDims
        [ind,col_index_equal_to_one,val]=find(X(:));%#ok no tilde pre-R2009b
        %X(:) converts X to a column vector. Treating X(:) as a matrix
        %forces val to be the actual value, instead of the column index.
        if length(ind)>K
            if strcmp(side,'first')
                %select first K outputs
                ind=ind(1:K);
                val=val(1:K);
            else
                %select last K outputs
                ind=ind((end-K):end);
                val=val((end-K):end);
            end
        end
        [varargout{1:(end-1)}] = ind2sub(size(X),ind);
        varargout{end}=val;
    else
        ind=find(X);
        if length(ind)>K
            if strcmp(side,'first')
                %select first K outputs
                ind=ind(1:K);
            else
                %select last K outputs
                ind=ind((end-K):end);
            end
        end
        [varargout{:}] = ind2sub(size(X),ind);
    end
else
    if nargout>nDims
        %Tilde (~) to ignore outputs was introduced in R2009b. It is
        %probably faster to ignore the extra output than to use an if.
        [ind,col_index_equal_to_one,val]=find(X(:),K,side);%#ok
        %X(:) converts X to a column vector. Treating X(:) as a matrix
        %forces val to be the actual value, instead of the column index.
        [varargout{1:(end-1)}] = ind2sub(size(X),ind);
        varargout{end}=val;
    else
        ind=find(X,K,side);
        [varargout{:}] = ind2sub(size(X),ind);
    end
end
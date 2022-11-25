classdef SimFun
    % SimFun Wrapper around similarity function comparing two edge vectors.
    %
    % Matlab does not allow us to infer the number of arguments a function
    % handle receives. This SimFun handle wraps a function handle, simfun,
    % with a type, a SimFunType, that describes the number of arguments it
    % receives.
    %
    % See also: SimFunType
    properties
        fun     % Handle for similarity function
        type    % SimFunType, the number of arguments simfun receives
    end
    methods
        function this = SimFun(fun, type)
            % Defaults
            if nargin == 0
                fun = @(b1, b2)(corr(b1(:), b2(:)));
                type = SimFunType.TwoArgs;
            end
            
            % Construct
            this.fun = fun;
            this.type = type;
        end
        function this = set.fun(this, fun)
            assert(isa(fun, 'function_handle'));
            this.fun = fun;
        end
        function this = set.type(this, type)
            assert(isa(type, 'SimFunType'));
            this.type = type;
        end
        function similarity = feval(this, varargin)
            % Overload feval to call the wrapped similarity function.
            similarity = this.fun(varargin{:});
        end
        function similarity = subsref(this, S)
            % Overload subsref to provide syntactic sugar for calling the
            % wrapped similarity function.
            switch S(1).type
                case '.'
                    % Call builtin subsref so we don't break dot notation.
                    similarity = builtin('subsref', this, S);
                case '()'
                    % Invoke simfun with provided arguments
                    similarity = feval(this, S.subs{:});
                otherwise
                    error('Unsupported subscript reference.');
            end
        end
    end
end
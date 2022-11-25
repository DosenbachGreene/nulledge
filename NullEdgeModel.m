classdef NullEdgeModel
    % NullEdgeModel Express relationship between edges and traits.
    % 
    % This is the main entry point for using nulledge. Provide Y, a
    % participant x edge matrix of observed functional connectivity, and
    % two participant x 1 vectors for traits x1 and x2. To test if the
    % forward projections of x1 and x2, b1 and b2 respectively, are more
    % similar than by chance:
    %
    % model = NullEdgeModel;
    % model.Y = Y;
    % model.x1 = x1;
    % model.x2 = x2;
    % r = model.similarity();      % similarity between b1 and b2
    % rnull = model.permute();     % permuted (null) similarity
    % histogram(abs(rnull)); xline(abs(r));
    % p = sum(abs(rnull) > abs(r)) / length(rnull) % p-value
    %
    % If you have some covariates Z then simply:
    %
    % model.Z = Z;
    %
    % NullEdgeModel computes similarity using correlation by default. If
    % you have a custom similarity function you can set it thusly. Note
    % that b1 and b2 are row vectors.
    %
    % simfun = SimFun;
    % simfun.type = SimFunType.TwoArgs;
    % simfun.simfun = @(b1, b2)(corr(b1(:), b2(:)));
    % model.fun = simfun;
    % 
    % If you have some template b from a different data set but don't know
    % the corresponding x, you can approximate it via backprojection
    %
    % model.x1 = model.backproject(b1);
    %
    % If you wish to inspect the forward projections from the model, you
    % can use:
    %
    % B = model.forwardproject();
    % b1 = B(1,:);
    % b2 = B(2,:);
    %
    % Or to obtain randomized forward projections:
    %
    % [similarity, B_rand] = model.randomize();
    % b1_rand = B_rand(1,:);
    % b2_rand = B_rand(2,:);
    
    properties
        x1      % Participant x 1 vector for trait 1.
        x2      % Participant x 1 vector for trait 2.
        Z       % Partizipant x z matrix of z nuisance covariates.
        Znull   % Null space of Z, Znull = null(Z').
        simfun  % Function to compute similarity, default correlation.
    end
    properties(SetAccess=protected, GetAccess=public)
        U       % Singular value decomposition of Y (components).
        s       % Singular values of edges (as column vector).
        V       % Singular value decomposition of Y (axes).
    end
    properties(Access=protected)
        BU      % Beta vectors for x1 and x2.
    end
    properties(Dependent)
        Y       % Participant x edge matrix of observed connectivity.
    end
    methods
        function this = NullEdgeModel(Y, x1, x2, Z, simfun)
            % Constructor for NullEdgeModel. Creates an empty model when
            % called with no arguments.
            
            % Default to empty values.
            this.U = [];
            this.s = [];
            this.V = [];
            this.x1 = [];
            this.x2 = [];
            this.Z = [];
            this.simfun = SimFun;
            this.BU = [];
            
            % Store whichever arguments we are given.
            if nargin >= 1
                this.Y = Y;
            end
            if nargin >= 2
                this.x1 = x1;
            end
            if nargin >= 3
                this.x2 = x2;
            end
            if nargin >= 4
                this.Z = Z;
            end
            if nargin >= 5
                this.simfun = simfun;
            end
        end
        function Y = get.Y(this)
            % Reconstruct Y from its singular value decomposition.
            Y = this.U .* this.s' * this.V';
        end
        function this = set.Y(this, Y)
            % Y must be a 2-d matrix of floating point values.
            % Rows should be observations and columns should be edges.
            if ~isempty(Y)
                assert(ismatrix(Y));
                assert(isfloat(Y));
                assert(length(size(Y)) == 2);
                if size(Y,2) < size(Y,1)
                    warning('Y should have participants in rows and edges in columns.');
                end
                if ~isempty(this.x1)
                    assert(size(Y,1) == size(this.x1,1));
                elseif ~isempty(this.x2)
                    assert(size(Y,1) == size(this.x2,1));
                elseif ~isempty(this.Z)
                    assert(size(Y,1) == size(this.Z,1));
                end
            end
            
            % Compute and store singular value decomposition of Y.
            if isempty(Y)
                this.U = [];
                this.s = [];
                this.V = [];
            else
                [this.U, S, this.V] = svd(Y, 'econ');
                this.s = diag(S);
            end
            
            % Reset BU.
            this.BU = [];
            this.BU = [];
        end
        function this = set.x1(this, x1)
            % x1 should be a participant x 1 vector.
            if ~isempty(x1)
                assert(isvector(x1));
                assert(isfloat(x1));
                if ~isempty(this.U)
                    assert(size(x1,1) == size(this.U,1));
                elseif ~isempty(this.x2)
                    assert(size(x1,1) == size(this.x2,1));
                elseif ~isempty(this.Z)
                    assert(size(x1,1) == size(this.Z,1));
                end
            end
            this.x1 = x1;
            
            % Reset BU.
            this.BU = [];
        end
        function this = set.x2(this, x2)
            % x2 should be a participant x 1 vector.
            if ~isempty(x2)
                assert(isvector(x2));
                assert(isfloat(x2));
                if ~isempty(this.U)
                    assert(size(x2,1) == size(this.U,1));
                elseif ~isempty(this.x1)
                    assert(size(x2,1) == size(this.x1,1));
                elseif ~isempty(this.Z)
                    assert(size(x2,1) == size(this.Z,1));
                end
            end
            this.x2 = x2;
            
            % Reset BU.
            this.BU = [];
        end
        function this = set.Z(this, Z)
            % Z should be a participant x covariate vector.
            if ~isempty(Z)
                assert(ismatrix(Z));
                assert(isfloat(Z));
                if ~isempty(this.U)
                    assert(size(Z,1) == size(this.U,1));
                elseif ~isempty(this.x1)
                    assert(size(Z,1) == size(this.x1,1));
                elseif ~isempty(this.x2)
                    assert(size(Z,1) == size(this.x2,1));
                end
                assert(size(Z,2) < (size(Z,1) + 3));
            end
            this.Z = Z;
            if isempty(Z)
                this.Znull = [];
            else
                this.Znull = null(Z');
            end
            
            % Reset BU.
            this.BU = [];
        end
        function this = set.simfun(this, simfun)
            assert(isa(simfun, 'SimFun'));
            this.simfun = simfun;
        end
        function x = backproject(this, B)
            % Backproject b to recover x orthogonalized to Z (i.e. with Z
            % regresssed out).
            % See also backproject_null and backproject_nonull.
            x = this.backproject_null(B);
            if ~isempty(this.Znull)
                x = this.Znull * x;
            end
        end
        function xnull = backproject_null(this, B)
            % Backproject b into the null space of Z.
            % See also backproject and backproject_nonull.
            assert(ismatrix(B));
            assert(isfloat(B));
            
            % If b is in the space of Y then project b into the
            % (orthonomal) space of U.
            if size(B,2) == size(this.V,1)
                B = B * this.V .* sinv(this.s)';
            end
            assert(size(B,2) == size(this.U,2));
            
            % Backproject.
            if isempty(this.Znull)
                % Backproject without covariates in Z.
                xnull = this.backproject_nonull(B);
            else
                % Backproject with covariates in Z.
                zub = this.Znull' * this.U * B';
                if size(B,1) == 1
                    % Optimized version for vector B.
                    xnull = zub ./ sum(zub .* zub);
                else
                    % General case.
                    zub_pinv = pinv(zub);
                    xnull = zub * (zub_pinv * zub_pinv');
                end
            end
        end
        function xnull = backproject_nonull(this, B)
            % Backproject b into the space of X without regard to nuisance
            % covariates in Z.
            % See also backproject and backproject_null.
            assert(ismatrix(B));
            assert(isfloat(B));
            
            % If b is in the space of Y then project b into the
            % (orthonomal) space of U.
            if size(B,2) == size(this.V,1)
                B = B * this.V .* sinv(this.s)';
            end
            assert(size(B,2) == size(this.U,2));
            
            % Backproject.            
            if size(B,1) == 1
                % Optimized version for vector B
                xnull = this.U * B' ./ sum(B.*B);
            else
                % General case.
                B_pinv = pinv(B);
                xnull = this.U * B' * (B_pinv' * B_pinv);
            end
        end
        function BU = forwardproject_U(this, X)
            % Obtain an edge vector B for X in the space of U, controlling
            % for Z. X may have multiple columns, in which case B will have
            % multiple rows. Default X = [x1, x2] from this model.
            % See also: forwardproject
            if nargin < 2
                assert(~isempty(this.x1));
                assert(~isempty(this.x2));
                X = [this.x1, this.x2];
            end
            assert(~isempty(this.U));
            assert(ismatrix(X));
            assert(isfloat(X));
            assert(size(X,1) == size(this.U,1) || size(X,1) == size(this.Znull,2));
            assert(size(X,2) < size(X,1));
            if isempty(this.Znull)
                % No covariates Z to account for.
                BU = pinv(X) * this.U;
            elseif size(X,1) == size(this.U,1)
                % Put X in the null space of Z.
                BU = pinv(this.Znull' * X) * this.Znull' * this.U;
            elseif size(X,1) == size(this.Znull,2)
                % X is already in the null space of Z.
                BU = pinv(X) * this.Znull' * this.U;
            else
                error([num2str(size(X,1)), ' x ', num2str(size(X,2)), ' matrix X should have the same number of rows as ', num2str(size(this.U,1)), ' x ', num2str(size(this.V,1)), ' matrix Y.']);
            end
        end
        function B = forwardproject(this, X)
            % Obtain an edge vector B for X in the space of Y, controlling
            % for Z. X may have multiple columns, in which case B will have
            % multiple rows. Default X = [x1, x2] from this model.
            % See also: connectivity_XU
            if nargin < 2
                assert(~isempty(this.x1));
                assert(~isempty(this.x2));
                X = [this.x1, this.x2];
            end
            B = this.forwardproject_U(X) .* this.s' * this.V';
        end
        function [similarity, B] = similarity(this)
            % Compute the similarity of b1 and b2 in this model using
            % simfun. Default similarity function is Pearson linear
            % correlation.
            
            % Forward project into space of Y.
            B = this.forwardproject();
            
            % Compute similarity.
            similarity = this.simfun(B(1,:), B(2,:));
        end
        function [similarity, B_rand] = randomize(this)
            % Randomize edge vectors to obtain a null value for similarity.
            % Optionally returns the two randomized edge vectors b1 and b2.
            
            % Regenerate BU if needed.
            if isempty(this.BU)
                this.BU = this.forwardproject_U();
            end
            
            % Randomize BU by sign flipping columns.
            % Sign flipping preserves the underlying edge structure.
            BU_rand = signflip(this.BU);
            
            % Backproject into the null space.
            Xnull_rand = this.backproject_null(BU_rand);
            
            % Forward project back into space of Y.
            % This step accounts for correlation between the x variables.
            B_rand = this.forwardproject(Xnull_rand);
            
            % Compute similarity.
            similarity = this.simfun(B_rand(1,:), B_rand(2,:));
        end
        function similarities = permute(this, n)
            % Perform n randomizations to obtain a null distribution for
            % similarity. Default n = 1024.
            if nargin == 1
                n = 1024; % default number of permutations
            end
            assert(isscalar(n));
            assert(rem(n,1) == 0); % must be integer
            assert(n > 0); % must be positive
            
            % Preallocate memory.
            similarities = zeros(n,1);
            
            % Perform n randomizations and store resultant similarity
            % values.
            parfor i=1:n
                similarities(i) = this.randomize();
            end
        end
    end
end
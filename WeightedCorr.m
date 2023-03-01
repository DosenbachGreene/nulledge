classdef WeightedCorr
    % Compute the correlation between x and y with observation weights w.
    %
    % Returns the correlation r and covariance rho.
    %
    % Example with frequency weights:
    %     x = [1;1; 2; 5;5;5];
    %     y = [3;3; 2; 4;4;4];
    %     corr(x,y) % 0.8
    %     xx = [1;2;5];
    %     yy = [3;2;4];
    %     WeightedCorr.Frequency(xx, yy, [2;1;3]) % 0.8
    %
    % Example with reliability weights. The last observation seems like an
    % outlier, so we assign it a low reliability.
    %     x = [1;2;3; 10];
    %     y = [1;2;3; -2];
    %     corr(x,y) % -0.8315, outlier makes positive correlation negative
    %     WeightedCorr.Reliability(x, y, [1;1;1; 0.01]) % 0.7043
    %
    % Hint: If you need a function handle for, e.g. reliability weights:
    %     fhandle = @(x,y,w)(feval(WeightedCorr.Reliability, x, y, w));
    enumeration
        Reliability, Frequency
    end
    methods
        function [r, rho] = feval(this, x, y, w)
            arguments
                this WeightedCorr
                x {mustBeNumeric,mustBeVector}
                y {mustBeNumeric,mustBeVector}
                w {mustBeNumeric,mustBeVector,mustBePositive}
            end
            % Overloading feval to compute weighted correlation.
            % See class level documentation for details.
            
            % Additional argument validation.
            assert(length(x) > 1, "Cannot compute correlation between scalars.");
            assert(length(x) == length(y), "x, y, and w must be the same length");
            assert(length(x) == length(w), "x, y, and w must be the same length");

            % Convert all arguments to column vectors to avoid Matlab
            % getting creative with broadcasting.
            x = x(:);
            y = y(:);
            w = w(:);
        
            % Find sums of w.
            sumw = sum(w);
            if this == WeightedCorr.Reliability
                sumw2 = sum(w.*w);
            end
        
            % Center x and y by weighted means.
            meanxw = sum(x.*w)/sumw;
            meanyw = sum(y.*w)/sumw;
            x0w = x - meanxw;
            y0w = y - meanyw;
        
            % Compute weighted sums.
            sumxxw = sum(w.*(x0w.*x0w));
            sumyyw = sum(w.*(y0w.*y0w));
            sumxyw = sum(w.*(x0w.*y0w));
            
            % Compute variances/covariances with appropriate denominator.
            if this == WeightedCorr.Frequency
                % Bessel's correction for sample variance with frequency weights.
                denominator = sumw - 1;
            elseif this == WeightedCorr.Reliability
                % Reliability weights.
                denominator = sumw - sumw2/sumw;
            end
            varxw = sumxxw / denominator;
            varyw = sumyyw / denominator;
            rho = sumxyw / denominator;
        
            % Compute weighted correlation.
            r = rho / sqrt(varxw * varyw);
        end
        function [r, rho] = subsref(this, S)
            % Overload subsref to provide syntactic sugar for feval.
            switch S(1).type
                case '.'
                    % Call builtin subsref so we don't break dot notation.
                    r = builtin('subsref', this, S);
                case '()'
                    % Delegate to feval.
                    [r, rho] = feval(this, S.subs{:});
                otherwise
                    error('Unsupported subscript reference.');
            end
        end
    end
end
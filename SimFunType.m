classdef SimFunType
    % SimFunType Expresses the number of arguments a function receives.
    %
    % A similarity function computes the similarity between two vectors of
    % estimated edges. The simplest type of similarity function receives
    % two arguments, the two edge vectors to compare. This is a
    % SimFunType.TwoArgs.
    %
    % Some similarity functions also account for the standard error of the
    % edges. This type of similarity function receives four arguments. The
    % first two arguments are the edge vectors to compare. The next two
    % arguments are vectors of the standard errors for the first and second
    % edge vector, respectively.  This is a SimFunType.FourArgs.
    %
    % See also: SimFun
    enumeration
        TwoArgs, FourArgs
    end
end
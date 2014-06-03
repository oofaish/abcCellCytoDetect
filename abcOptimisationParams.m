function optimisationParams = abcOptimisationParams()
    optimisationParams = struct();
    optimisationParams.minRadius = 15;
    optimisationParams.maxRadius = 100;
    %optimisationParams.maxAngle = 0;
    %optimisationParams.minAngle = pi;

    optimisationParams.checkUniformIntensity = false;%makes no difference

    optimisationParams.x0         = [ optimisationParams.minRadius * 2, 1, 0 ];
    %optimisationParams.x0         = [ 60, 1.34, 2.0364 / 2 ];
    optimisationParams.lowerBound = [ optimisationParams.minRadius, 1, 0 ];
    optimisationParams.upperBound = [ optimisationParams.maxRadius, 2, pi ];

    optimisationParams.visuallyVerbose = false;
    %optimisationParams.algorithm = 'gridsearch';
    optimisationParams.algorithm = 'simplex';
    %optimisationParams.lossFunction = 'exactNonZeros';
    optimisationParams.lossFunction = 'truncatedSquareLoss';
    optimisationParams.truncation = optimisationParams.maxRadius;
    optimisationParams.epsilon = 1e-6;
end    
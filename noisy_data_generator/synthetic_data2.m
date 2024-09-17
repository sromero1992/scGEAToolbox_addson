function [p,n,X,ylinear,ynlinear] = synthetic_data2(p, n, linear_source, ...
                                        linear_target, nlinear_source, nlinear_target)

    num_linear = size(linear_source, 2);
    num_nlinear = size(nlinear_source, 2);
    B = randn(n,n);
    R = corrcov(B'*B);
    X = mvnrnd(zeros(n,1),R,p);
    % Source features are highly correlated to target features 
    corrStd = 0.1;
    % Linear correlation
    t = 3*X(:, linear_source(1)) - 5*X(:, linear_target(2)) + ...
        7*X(linear_source(3))+ corrStd*randn(p, 1);
    ylinear = rescale(t,0,1);
    % All cells from source multiplied by source predictor
    X(:, linear_target) = X(:, linear_target) + ylinear;

    noiseStd = 0.1;
    % Nonlinear predictor
    t = sin( X(:, nlinear_source(3)).*X(:, nlinear_source(3)) ) + ...
        0.2*exp( X(:, nlinear_source(1))).* ...
        log2( abs(10*X(:, nlinear_source(2)) + 1 ) ) + ...
        noiseStd*randn(p,1);
    ynlinear = rescale(t,0,1);
    % All cells from source multiplied by source predictor
    X(:, nlinear_target) = X(:, nlinear_source).*t;

    X = zscore(X);
end

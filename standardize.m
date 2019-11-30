% function standardize.m standardizes a data matrix such that each column
% has zero mean and unit variance

function X_std = standardize(X_nonstd);

[rows, cols] = size(X_nonstd);
X_std = (X_nonstd-kron(mean(X_nonstd),ones(rows,1)))./(kron(std(X_nonstd),ones(rows,1)));


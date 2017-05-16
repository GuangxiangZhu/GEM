function [ rmse,A2 ] = SVD3D( A,B)

% Conformation alignment using the singular value decomposition (SVD) algorithm
% References:
% Sorkine,O. Least-squares rigid motion using SVD. (2009) Technical notes, 120, 3.
n=size(A,1);
[ret_R, ret_t] = RigidTransform3D(A, B);
A2 = (ret_R*A') + repmat(ret_t, 1, n);
A2 = A2';

% Find the error
err = A2 - B;
err = err .* err;
err = sum(err(:));
rmse = sqrt(err/n);

end


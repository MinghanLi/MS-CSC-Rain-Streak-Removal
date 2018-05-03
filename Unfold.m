%% video to low-rank matrix
function [X] = Unfold( X, dim, i )
X = reshape(shiftdim(X,i-1), dim(i), []);
end
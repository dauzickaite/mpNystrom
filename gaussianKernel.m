function K = gaussianKernel(dataset,sigma)

n = size(dataset,1);
xynorms = zeros(n);

% compute the norms of the upper triangular matrix
% loop through the rows
for rowind = 1:n
    for colind = rowind+1:n
        xynorms(rowind,colind) = ...
            norm(dataset(rowind,:) - dataset(colind,:))^2;
        xynorms(colind,rowind) = xynorms(rowind,colind);
    end
end

K = exp(-xynorms/(2*sigma^2));
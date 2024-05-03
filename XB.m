
function XieBeni = XB(X, U, P, m)
    % X = Data
    % U = Membership matrix
    % H = Centroid matrix
    % m = parameter of fuzziness

    [I, ~] = size(X);
    [~, K] = size(U);
    D = zeros(I, K);
    
    % epsilon = 1e-10;
    for i = 1:I
        for k = 1:K
            D(i, k) = sum((X(i, :) - P(k, :)).^2);
            % D(i, k) = max(sum((X(i, :) - P(k, :)).^2), epsilon);
        end
    end

    distH = Inf;
    for k1 = 1:(K - 1)
        for k2 = (k1 + 1):K
            if sum((P(k1, :) - P(k2, :)).^2) < distH
                distH = sum((P(k1, :) - P(k2, :)).^2);
            end
        end
    end

    XieBeni = sum(sum((U.^m) .* D)) / (I * distH);
end


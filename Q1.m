% INPUT: L is the link matrix that is stochastic and has no cycles
%        n is the size of the link matrix

% OUTPUT: v is the importance vector

function v = project3_Q1(L, n)
    % initialize the first distribution vector
    v = zeros(n, n);
    for i = 1:n
        for j = 1:n
            v(i, j) = 1/n;
        end
    end
    
    x = L*v;
    while x ~= v
        v = x;
        x = L*v;
    end   
end

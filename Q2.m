size = 12;
alpha = 0.1;
u = zeros(size, 1);
e = ones(size, 1);
% for i = 1:size
%     u(i, 1) = 1/size;
% end
u(12) = 1;

A = [0 0 1/5 0 0 0 0 0 0 0 0 0;
     1 0 1/5 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1/2 0 0 0 0 0;
     0 0 1/5 0 1/2 0 0 0 0 0 0 0;
     0 0 1/5 0 0 0 0 0 0 0 0 0;
     0 0 1/5 0 0 0 0 0 0 0 0 0;
     0 0 0 1 1/2 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1/2 0 1 0 0 0;
     0 0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 1 0 0]
 
 B = [0 1/12 1/5 0 0 1/12 0 0 0 0 0 0;
      1 1/12 1/5 0 0 1/12 0 0 0 0 0 0;
      0 1/12 0 0 0 1/12 1/2 0 0 0 0 0;
      0 1/12 1/5 0 1/2 1/12 0 0 0 0 0 0;
      0 1/12 1/5 0 0 1/12 0 0 0 0 0 0;
      0 1/12 1/5 0 0 1/12 0 0 0 0 0 0;
      0 1/12 0 1 1/2 1/12 0 0 0 0 0 0;
      0 1/12 0 0 0 1/12 1/2 0 1 0 0 0;
      0 1/12 0 0 0 1/12 1 0 0 0 0 0;
      0 1/12 0 0 0 1/12 0 0 0 0 1 0;
      0 1/12 0 0 0 1/12 0 0 0 0 0 1;
      0 1/12 0 0 0 1/12 0 0 0 1 0 0]

 C = alpha*B + (1-alpha)*u*transpose(e);
 
 [v, c] = Page_Rank(C, size);
 
 display(v);
%  plot(c);
 
% INPUT: L is the link matrix that is stochastic and has no cycles
%        n is the size of the link matrix

% OUTPUT: v is the importance vector
%         c is a list of rate of convergence

function [v, c] = Page_Rank(L, n)
    % initialize the first distribution vector
    v = zeros(n, 1);
    c = [];
    for i = 1:n
        v(i, 1) = 1/n;
    end
    
    x = L*v;
    while x ~= v
        c = [c, norm(x-v, 2)];
        v = x;
        x = L*v;
    end   
end

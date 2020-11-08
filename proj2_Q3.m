% Q3
% (a)
% load matrix M of 1000 webpage network
load('incidencematrix.mat')
% total number of rows and columns in M
rowNum = size(M,1);
colNum = size(M,2);
% sum of 1's in each column of M as a 1*n matrix
cSum = zeros(1,colNum);
for c = 1:colNum
   column = M(:,c);
   colTotal = sum(column);
   cSum(1,c) = colTotal;
end
% take 1/n for each element of cSum and make Ns
Ns = cSum;
for c = 1:colNum
    if(Ns(1,c) ~= 0)
        Ns(1,c) = 1./Ns(1,c);
    end
end
% use M and Ns to make matrix A
A = zeros(rowNum,colNum);
for i = 1:rowNum
    for j = 1:colNum
        if(M(i,j) == 1)
           A(i,j) = Ns(1,j); 
        end
    end
end
% finished making matrix A
% now modify A to make matrix B
B = zeros(rowNum,colNum);
for j = 1:colNum
   if(sum(A(:,j))==0)
       N = 1/rowNum;
       Ns(1,j) = N;
       for i = 1:rowNum
           B(i,j) = N;
       end
   else
       for i = 1:rowNum
           B(i,j) = A(i,j);
       end
   end
end
% make Cu,a using function created in previous questions
alpha = 0.95;
u = (1/colNum)*ones(colNum,1);
e = ones(colNum,1);
C = alpha*B + (1-alpha)*u*transpose(e);
[v, c] = Page_Rank(C, colNum);
% now v is the importance vector, we need the top10 of them
imp = v;
topList = zeros(10,2);
for i = 1:10
    [m,I] = max(imp);
    topList(i,1)=I;
    topList(i,2)=m;
    imp(I,1)=0;
end
% display top10 to console
% - first column refers to the index
% - second column refers to their importance
display(topList);

% (b)
% because linking out to other websites would make them more important, so
% it's better not to link to any other webpages
% Therefore, the 1001th column should be all zeros, we only need to
% consider buying how many websites to let them link to ours

% in order to calculate price, we first rank all websites in decending
% order for convenience
imp = v;
rankList = zeros(colNum,2);
for i = 1:colNum
    [m,I] = max(imp);
    rankList(i,1)=I;
    rankList(i,2)=m;
    imp(I,1)=0;
end
% according to the price function, the larger the rank, the less price we
% need to pay, so we should start on the least-ranked websites to create
% inward linkage
% first generate a matrix newV with size 1001*1001
newC = zeros(1001);
for i = 1:1000
    for j = 1:1000
        newC(i,j) = C(i,j);
    end
end
% update rowNum and colNum for newV
rowNum = size(newC,1);
colNum = size(newC,2);
% calculate total spending to get our website into top 10%
cost = 1000;
for i = 1000:-1:1
    % ask the least-ranked website to link to ours, update newC and cost
    index = rankList(i,1);
    newC(1001,index) = 1;
    cost = cost + (1000-i+1)^2;
    % calculate new importance vector, importance boundary, and compare
    % whether our website is in the top10%
    [v,c] = Page_Rank(newC,colNum);
    result = Check_Importance(v,0.1,colNum,1001);
    if(result == true)
       break; 
    end
end
% display the total cost to get our website into top 10%
display(cost)






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

% INPUT: v is the link matrix
%        percent is the percentage looking at, in decimal form (e.g. 0.1)
%        N is the row/column number of link matrix
%        i is the index of our own website

% OUTPUT: importance is the importance value at the given percentage in v

function result = Check_Importance(v,percent,N,index)
    result = false;
    imp = v;
    n = ceil(percent * N);
    list = zeros(n,2);
    for i = 1:n
        [m,I] = max(imp);
        if(I == index)
            % if our website is found in the top percent, return true
            result = true;
        end
        list(i,1)=I;
        list(i,2)=m;
        imp(I,1)=0;
    end
end

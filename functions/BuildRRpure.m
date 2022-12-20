function  A  = BuildRRpure( n,k )

% build pure random regular network

stbs= repelem(1:n,k);
stbs = stbs(randperm(n*k));

s = stbs(1:2:end);
t = stbs(2:2:end);

A = sparse([s,t],[t,s],1,n,n);

% remove self and multiple
A(1:n+1:end) = 0;
A = logical(A);

deg = sum(A);

% add link if degree is less than k, and remove if more than k
while any(deg~=k)
    [~,i] = min(deg);
    miss = find(deg<k);    
    opt = setdiff(miss,[i,find(A(i,:))]);
    if isempty(opt)
        opt = 1:n;
    end
    j = opt( randi(length(opt)) );
    
    A(i,j) = 1;
    A(j,i) = 1;
    
    if deg(j) == k
        jn = find(A(j,:));
        opt = setdiff(jn,i);
        tr = opt(randi(length(opt)));
        A(j,tr) = 0;
        A(tr,j) = 0;
    end
    deg = sum(A);       
end

end
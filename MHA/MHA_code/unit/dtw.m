
function [match_idx_P,k_class,D,Dist]=dtw(t,r,d)
%Dynamic Time Warping Algorithm
%Dist is unnormalized distance between t and r
%D is the accumulated distance matrix
%k is the normalizing factor
%w is the optimal path
%t is the vector you are testing against
%r is the vector you are testing
[rows,N]=size(t);
[rows,M]=size(r);



D=zeros(size(d));
D(1,1)=d(1,1);

for n=2:N
    D(n,1)=d(n,1)+D(n-1,1);
end
for m=2:M
    D(1,m)=d(1,m)+D(1,m-1);
end
for n=2:N
    for m=2:M
        D(n,m)=d(n,m)+min([D(n-1,m),D(n-1,m-1),D(n,m-1)]);
    end
end

Dist=D(N,M);
n=N;
m=M;
k_class = 1;

match_idx=[];
match_idx_P=zeros(N,M);
match_idx(1,:)=[N,M];
match_idx_P(N,M)=1;

while ((n+m)~=2)
    if (n-1)==0
        m=m-1;
    elseif (m-1)==0
        n=n-1;
    else 
        [value,number]=min([D(n-1,m),D(n,m-1),D(n-1,m-1)]);
        switch number
         case 1
           n=n-1;
         case 2
           m=m-1;
         case 3
           n=n-1;
           m=m-1;
           k_class=k_class+1;
         end
    end
    match_idx=cat(1,match_idx,[n,m]);
    match_idx_P(n,m) = 1;
end
 



function [ d ] = myEuclidean( t,r,weight )

[rows,N]=size(t);
[rows,M]=size(r);

 
for n=1:N
   for m=1:M
         d(n,m) = weight.*sqrt((t(1,n)-r(1,m)).^2 )+ sqrt((t(2,n)-r(2,m)).^2);
   end
end


end


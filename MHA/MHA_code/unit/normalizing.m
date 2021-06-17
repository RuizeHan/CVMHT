function Y = normalizing(X,a,b)
    
    if length(X) == 1
        Y = (a+b)/2;
    else
        minX=min(min(X));
        maxX=max(max(X));
        Y=(X-minX)*((b-a)/(maxX-minX))+a;
    end
end
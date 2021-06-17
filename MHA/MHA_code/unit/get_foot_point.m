function proj_point = get_foot_point(point,line_p)
    
    x1 = line_p(1);
    y1 = line_p(2);
    x2 = line_p(3);
    y2 = line_p(4);
    x3 = point(1);
    y3 = point(2);

    yk = ((x3-x2)*(x1-x2)*(y1-y2) + y3*(y1-y2)^2 + y2*(x1-x2)^2) / (norm([x1-x2,y1-y2])^2);
    xk = ((x1-x2)*x2*(y1-y2) + (x1-x2)*(x1-x2)*(yk-y2)) / ((x1-x2)*(y1-y2));
    if x1 == x2
        xk = x1;
    end

    if y1 == y2
        xk = x3;
    end
    proj_point = [xk,yk];

end

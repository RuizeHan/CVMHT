
function [scale_y]=get_scale_y(idx_top,idx_ego,vec1,vec2)

% scale_y = y_top:y_ego
 y_top = vec1(2,:);
 y_ego = vec2(2,:);

 
 y_top = y_top(idx_top);
 y_ego = y_ego(idx_ego);
 
 scales = y_top./y_ego;
 scale_y = mean(scales);
end
 



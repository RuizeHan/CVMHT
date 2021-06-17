 function [vec,idx_ord] = get_distribution_vector_land(cet_x,obj_h,img_w)
 
 [cet_x_ord,idx_ord] = sort(cet_x);
 obj_h_ord = obj_h(idx_ord);
 
 idx_num = length(cet_x);
 vec = zeros(2,idx_num);
 
 if idx_num > 0
     for i = 1 : idx_num
        vec(1,i) = cet_x_ord(i) - img_w/2;
        vec(2,i) = obj_h_ord(i);
     end
 else
     vec = [];
 end
 
 
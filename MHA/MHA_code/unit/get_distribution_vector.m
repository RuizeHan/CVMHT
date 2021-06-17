 function [vec,index] = get_distribution_vector(flag,angle_diff,foot,dis)
 
 index = find(flag == 1);
 foot_eff = foot(:,index);
 dis_eff = dis(index);
 angle_eff = angle_diff(index);
 
 [angle_ord,idx_ord] = sort(angle_eff,'descend');
 foot_ord = foot_eff(:,idx_ord);
 dis_ord = dis_eff(idx_ord);
 
 idx_num = length(index);
 vec = zeros(2,idx_num);
 
 if idx_num > 0
     for i = 1 : idx_num
        
        %  Option 1: the disdance between the foot points to the fist one
        %  vec(1,i) = sqrt((foot_ord(1,1)-foot_ord(1,i))^2+(foot_ord(2,1)-foot_ord(2,i))^2);
        
        %  Option 2 : the difference of the left angles
        %  vec(1,i) = -(cot(angle_ord(1,1))-cot(angle_ord(1,i)));
        
        %  Option 3 : the difference of the center angle
        vec(1,i) = cot(angle_ord(1,i));
        
        vec(2,i) = dis_ord(i);
     end
 else
     vec = [];
 end
 
 
 
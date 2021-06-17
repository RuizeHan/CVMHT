function [data] = UAV2vec_one(img,obj,cam,params,data)

    [h,w,~] = size(img);   
    I =img;
    vis = 0;
    try
    cam = [cam(1)+1/2*cam(3),cam(2)+1/2*cam(4)];
    catch err
    end
    obj_x = obj(:,1)+1/2*obj(:,3);
    obj_y = obj(:,2)+1/2*obj(:,4);
    obj_num = length(obj_x);

    obj_dx = obj_x - cam(1);
    obj_dy = obj_y - cam(2);
    angle = zeros(1,obj_num);
    % get the angle of each object
    for i = 1 : obj_num  
        vec_x = [1,0];
        vec_obj = [obj_dx(i),obj_dy(i)];
        angle(i) = acos(dot(vec_x,vec_obj)/(norm(vec_x)*norm(vec_obj)));
        
        if obj_dy(i) > 0
            angle(i) = 2 * pi - angle(i);
        end
            
    end
    angle_deg = angle * 180/pi; %�����ƣ�ת�Ƕ��Ƴ�180/pi
    search_num = params.search_num;
    
    k = data.max_idx;
    search_angle = pi/(search_num/2) * (k);

    vertical_ang = search_angle;
    right_ang = vertical_ang + pi/4;
    left_ang = vertical_ang + pi/4*3;

    % get the distance of each object in specfic view 
    flag = zeros(1,obj_num);
    dis = zeros(1,obj_num);
    foot = zeros(2,obj_num);

    for i = 1 : obj_num  

        if left_ang > 2 * pi  % !todo : when the leftangle greater than 2*pi      
            ang_gre = find(angle < 3/4 * pi & angle > 0);
            angle(ang_gre) = angle(ang_gre) + 2 * pi;        
        end

       P = [obj_dx(i),obj_dy(i)];
       O = [0,0];    
       V = [cos(vertical_ang),-sin(vertical_ang)];  
       OV = [0,0,cos(vertical_ang),-sin(vertical_ang)];  

        if angle(i) > right_ang && angle(i) < left_ang     % !todo : when the leftangle greater than 2*pi
           flag(i) = 1;
           proj_point = get_foot_point(P,OV);
           foot(1,i) = proj_point(1);
           foot(2,i) = proj_point(2);

           %  Option 1 : the distance of the objects P to the vertical line LV
           dis(i) = abs(det([V-O;P-O]))/norm(V-O);
           %  Option 2 : the distance of the objects P to the camera point O
           %  dis(i) = norm(P-O);               
        end    
    end

        % Occlusion
        occ_angle=params.occ_angle;
        if params.occ == 1
        for m = 1 : obj_num - 1
            for n = m + 1 : obj_num
                if  abs(angle_deg(m) - angle_deg(n)) < occ_angle % obj is occluded when the angle is approciate with another     
                    if dis(m) >= dis(n)
                        flag(m) = 0;
                    else
                        flag(n) = 0; 
                    end
                end
            end
        end
        end
     
        angle_diff = angle - vertical_ang;      
        [data.vec_top_1,data.index_top_1] = get_distribution_vector_one(flag,angle_diff,foot,dis);
        data.angle = angle;
        data.visable_flag = flag;
end

function [data] = UAV2vec(params,data,last_angle)
    
    img = data.img_top;
    obj = data.objs_top;
    cam = data.cam_top;    
    size_img_top = data.size_img_top;
    h = size_img_top(1);
    w = size_img_top(2);
    cam = [cam(1)+1/2*cam(3),cam(2)+1/2*cam(4)];
    obj_x = obj(:,1)+1/2*obj(:,3);
    obj_y = obj(:,2)+1/2*obj(:,4);
    I =img;
    obj_num = length(obj_x);
    vis = params.vis;
    if vis == 1
        figure(1)
        set (gcf,'Position',[300,600,600,400])
        imshow(I);
        hold on;
        plot(cam(1),cam(2),'pg','Linewidth',5);
        plot(obj_x(:),obj_y(:),'or','Linewidth',1.5);
    end

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
    
for k = 1 : search_num
         
        search_angle = last_angle - pi/(search_num/2) * (search_num) + pi/(search_num/2) * (k);
        if search_angle < 0     
            search_angle = search_angle + 2 * pi;        
        end
        
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
%                  if vis ==1
%                       %plot(proj_point(1)+cam(1),proj_point(2)+cam(2),'rs');
%                  nd
               
               %  Option 1 : the distance of the objects P to the vertical line LV
               dis(i) = abs(det([V-O;P-O]))/norm(V-O);
               %  Option 2 : the distance of the objects P to the camera point O
%              dis(i) = norm(P-O);               
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
        
        [data.vec_top{k},data.index_top{k}] = get_distribution_vector(flag,angle_diff,foot,dis);

        data.angle = angle;
    if vis ==1
        % get the equation of vertical line in the original axis
        LV  = get_line_equation(cam(1),cam(2),V(1)+cam(1),V(2)+cam(2));

       
            if LV(2) == 0
                plot([-LV(3)/LV(1),-LV(3)/LV(1)],[0,h],'c','LineWidth',1.5);
            else
                plot([0,w],[-LV(3)/LV(2),(-LV(3)-LV(1)*w)/LV(2)],'c','LineWidth',1.5);
            end
        
        corner_x = [w, 0, 0, w];
        corner_y = [0, 0, h, h];
        cor_dx = corner_x - cam(1);
        cor_dy = corner_y - cam(2);

        % get the angle of the line contected the camera with each corner 
        num_cor =  length(corner_x);
        angle_cor = zeros(1,num_cor);
        for i = 1 : num_cor
            vec_x = [1,0];
            vec_cor = [cor_dx(i),cor_dy(i)];
            angle_cor(i) = acos(dot(vec_x,vec_cor)/(norm(vec_x)*norm(vec_cor)));
            if cor_dy(i) > 0
                angle_cor(i) = 2 * pi - angle_cor(i);
            end
        end

      % get the cross point of the side
      for judge_angle  = [left_ang,right_ang]
            % L = [cos(left_ang),-sin(left_ang)];
            % R = [cos(right_ang),-sin(right_ang)]; 
            % judge which edge to cross  
            edge_index = 0;
            if judge_angle >= angle_cor(1) && judge_angle < angle_cor(2)
                edge_index = 1;
                L1 = [0,1,cam(2)];
            elseif judge_angle >= angle_cor(2) && judge_angle < angle_cor(3)
                edge_index = 2;
                L1 = [1,0,cam(1)];
            elseif judge_angle >= angle_cor(3) && judge_angle < angle_cor(4)
                edge_index = 3; 
                L1 = [0,1,cam(2)-h];
            elseif judge_angle <= angle_cor(1) || judge_angle > angle_cor(4)
                edge_index = 4;
                L1 = [1,0,cam(1)-w];
            end

            L2  = get_line_equation(0,0,cos(judge_angle),-sin(judge_angle));
            [X,Y] = get_cross_point(L1,L2);
            if vis ==1
                plot([ceil(X+cam(1)),cam(1)],[ceil(Y+cam(2)),cam(2)],'b','LineWidth',1.5);
            end
      end

        vec_vis = data.vec_top{k};
        if vis ==1 && ~isempty(vec_vis)
          figure(2)     
          plot(vec_vis(1,:),vec_vis(2,:),'ro');
          axis([-inf inf 0 max(vec_vis(2,:))])
        end
    end

end

end

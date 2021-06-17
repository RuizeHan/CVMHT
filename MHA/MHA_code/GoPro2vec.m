function [vec_ego,idx_ord] = GoPro2vec(mha_data,params)

    img = mha_data.img_ego;
    obj = mha_data.objs_ego;
    size_img_hor = mha_data.size_img_hor;
    vis = params.vis;
    img_w = size_img_hor(2);
    obj_x = obj(:,1);
    obj_y = obj(:,2);
    obj_w = obj(:,3);
    obj_h = obj(:,4);

    I = img;
    obj_num = length(obj_x);
    index_ego = 1:obj_num;
    if vis == 1
        figure(3)
        set (gcf,'Position',[300,100,600,400])
        imshow(I);
        hold on;

        for i = 1:obj_num
            rectangle('Position',[obj_x(i),obj_y(i),obj_w(i),obj_h(i)],'LineWidth',2,'EdgeColor','r');
        end
    end
    cet_x = obj_x + 0.5 * obj_w;
    cet_y = obj_y + 0.5 * obj_h;
    
    if vis == 1
        plot(cet_x,cet_y,'ro','Linewidth',1.5);
    end

    [vec_ego,idx_ord] = get_distribution_vector_land(cet_x,obj_h,img_w);

%   figure(4)
%   plot(vec_ego(1,:),vec_ego(2,:),'ro');
%   axis([-inf inf 0 max(vec_ego(2,:))])

end

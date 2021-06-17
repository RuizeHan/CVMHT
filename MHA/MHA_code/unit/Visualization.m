function Visualization(mha_data,vec_top_all,vec_hor_all)

imgs = {mha_data.img_top,mha_data.img_ego};
objs = {mha_data.objs_top, mha_data.objs_ego};
color = {'r','c'};
for view = 1:2
    
    obj_x = objs{view}(:,1);
    obj_y = objs{view}(:,2);
    obj_w = objs{view}(:,3);
    obj_h = objs{view}(:,4);

    I = imgs{view};
    obj_num = length(obj_x);
    
    cet_x = obj_x + 0.5 * obj_w;
    cet_y = obj_y + 0.5 * obj_h;
         
    figure
    imshow(I);
    set (gcf,'Position',[300,100,600,400])
    hold on;

    for i = 1:obj_num   
        rectangle('Position',[obj_x(i),obj_y(i),obj_w(i),obj_h(i)],'LineWidth',2,'EdgeColor','b');        
    end

    for i = 1 : length(cet_x) 
        text(cet_x(i),cet_y(i),num2str(i),'color',color{view},'fontsize',15);
    end

 end   

    figure(3)
    plot(vec_top_all(1,:),vec_top_all(2,:),'ro','Linewidth',1.5); hold on;
    for i = 1:length(vec_top_all)
        text(vec_top_all(1,i)+0.01,vec_top_all(2,i),num2str(i),'color','r','fontsize',15);
    end
    plot(vec_hor_all(1,:),vec_hor_all(2,:),'ch','Linewidth',1.5);
    for i = 1:length(vec_hor_all)
        text(vec_hor_all(1,i)+0.01,vec_hor_all(2,i),num2str(i),'color','c','fontsize',15);
    end
hold off;

end

function Visualization_vec(data,params,match_scores,max_sco,scale_y,flag)

    [~,img_w,~] = size(data.img_ego);
    vec1 = data.vec_top{data.max_idx};  % normalizaing£¿
    vec2 = data.vec_ego; 

     vec1(1,:) = vec1(1,:);
     vec2(1,:) = vec2(1,:)/(img_w/2);
     
    if flag == 1
        vec1(2,:) = normalizing(vec1(2,:),0,1);
        vec2(2,:) = normalizing(-vec2(2,:),0,1); 
    else   
        vec1(2,:) = vec1(2,:)./scale_y;
        vec2(2,:) = 1./vec2(2,:);
    end

    if params.vis == 1    
    figure(2)
    plot(vec1(1,:),vec1(2,:),'ro','Linewidth',1.5);
    axis([-1 1 0 max(vec1(2,:))])
    set (gcf,'Position',[950,600,600,400])

    figure(4)
    plot(vec2(1,:),vec2(2,:),'ro','Linewidth',1.5);
    axis([-1 1 0 max(vec2(2,:))])
    set (gcf,'Position',[950,100,600,400])

    figure(5)
    x = 1:length(match_scores);
%     score_norm = normalizing(-match_scores,0,1);    
    axis([0 length(match_scores) 0 max_sco])
    plot(2*x,match_scores,'r.');

    end
    

end
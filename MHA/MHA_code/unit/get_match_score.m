function [match_score,match_idx] = get_match_score(vec_s,vec_l,factor)

    len_s = length(vec_s(1,:));
    len_l = length(vec_l(1,:));
    
%     factor = 1.5;
    penalty_factor = factor^(len_l-len_s); 

    col = (1:len_l);
    col_sel = combntns(col,len_s);
    sz = size(col_sel);
    match_scores = zeros(1,sz(1));
    
    
    for i = 1 : sz(1)
        vec1 = vec_s;
        vec2 = vec_l(:,col_sel(i,:));
        match_scores(i) = sum((vec1(:) - vec2(:)).^2)/len_s;
    end
    
    [max_score,max_num] = min(match_scores);
    match_idx = col_sel(max_num,:);
    match_score = penalty_factor * max_score;

end

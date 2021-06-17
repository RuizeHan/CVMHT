
function [match_idx_ego,match_idx_top,score]=single_out(d,match_idx_P,thresh)

% vec1-top vec2-ego; n-top;m-ego
[n,m] = size(match_idx_P);
w = d.*match_idx_P;
match_idx_top = (1:n);
match_idx_ego = (1:m);

w_temp = w;
w_temp(w_temp==0) = Inf;

ego_num = sum(match_idx_P~=0,2);
top_num = sum(match_idx_P~=0,1);
idx_temp_top = zeros(1,m);
idx_temp_ego = zeros(1,n);
for i = 1:n
    if ego_num(i)>=2
        [value, loc]=min(w_temp(i,:));
        idx_temp_top(loc)=1;
        w(i,:) = w(i,:).*idx_temp_top;
    end
end

for i = 1:m
    if top_num(i)>=2
        [value, loc]=min(w_temp(:,i));
        idx_temp_ego(loc)=1;
        w(:,i) = w(:,i).*idx_temp_ego';
    end
end

w(w>thresh)=0;

mask_ego = 1 - all(w==0,1);
match_idx_ego = match_idx_ego.* mask_ego;
match_idx_ego(match_idx_ego==0)=[];
mask_top = 1 - all(w==0,2);
match_idx_top = match_idx_top.* mask_top';
match_idx_top(match_idx_top==0)=[];

score = sum(sum(w));

end
 



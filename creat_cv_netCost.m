function cv_net_cost =  creat_cv_netCost(data,param_cv_netCost,param_tracklet)

rho = param_tracklet.rho;
save_tracklet = 0;
use_app = 0;

tv_setments = data.cv_segments{1,1};
hv_segments = data.cv_segments{2,1};
cv_net_cost = cell(length(tv_setments),1);

if use_app ==1
    app_path = [param_tracklet.data_directory, 'AppCost/'];
    tracklets_dis_App = struct2cell(load(fullfile(app_path,['app_net_cost_',param_tracklet.sceneName,'.mat'])));
    tracklets_dis_App = tracklets_dis_App{1};
else
    param_cv_netCost.Sweight = 1.0;
    param_cv_netCost.Aweight = 0;
end

data.last_angle = 0;

for seg_i = 1 : length(tv_setments)
disp(['cv_cost:segment-',num2str(seg_i)]);
tv_segment_i = tv_setments{1,seg_i};
hv_segment_i = hv_segments{1,seg_i};

tv_tracklets = tv_segment_i.tracklet;
hv_tracklets = hv_segment_i.tracklet;

frame_range = (seg_i-1) * param_tracklet.num_frames + 1 : seg_i * param_tracklet.num_frames;
data.frame_range = frame_range;
[tracklets_sim_S,OverLap,data]= tracklet_associate_spatial_vote(data,tv_tracklets,hv_tracklets,param_tracklet,param_cv_netCost);
% [tracklets_dis_S,OverLap]= tracklet_associate_spatial(frame_range,data,tv_tracklets,hv_tracklets);

% Get the tracklets of two views
if save_tracklet == 1
    tracklets = {tv_tracklets,hv_tracklets};
    path = './Data/Tracklets/';
    savepath = strcat(path,'seg_',num2str(seg_i));
    if(~exist(savepath,'dir'))
        mkdir(savepath);
    end
    views = {'tv','hv'};
    for view_i = 1:2
        images = data.cv_images{view_i};
        sv_tracklet = tracklets{view_i};
        for tra_i = 1 : length(sv_tracklet)
        tracklet =  sv_tracklet(tra_i);
        frames = tracklet.frame;
        detection = tracklet.detection;
        sum_det = 0;
        for frm_i = 1 : length(frames)
            frm = frames(frm_i);
            image = images(frm);
            img = imread([image.folder,'/',image.name]);
            reg = detection(frm_i,:);
            det = imcrop(img,[reg(1),reg(2),reg(3)-reg(1),reg(4)-reg(2)]);
            imshow(det);
            det = imresize(det,[100,100]);
            sum_det = double(det) + sum_det;       
        end
        avg_det = uint8(sum_det/frm_i);
        imshow(avg_det);
        imwrite(avg_det,strcat(savepath,'/',views{view_i},'_tra_',num2str(tra_i),'.jpg'));
        end
    end
end

% tracklets_dis_S = normalizing(tracklets_dis_S,0,1);
if use_app ==1
tracklets_dis_A = tracklets_dis_App{seg_i};
end

if use_app ==1
    cv_net_cost{seg_i} =(param_cv_netCost.Sweight * tracklets_sim_S + param_cv_netCost.Aweight * rho.^(OverLap) .* (1 -tracklets_dis_A));
else
    cv_net_cost{seg_i} = param_cv_netCost.Sweight * tracklets_sim_S;
end

end


end
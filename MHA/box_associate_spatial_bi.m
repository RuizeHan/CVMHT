function [Match_top_hor_seg,invis_top,data] = box_associate_spatial_bi(tv_detection,hv_detection,tv_image_seg,hv_image_seg,param_tracklet,data)

ang = 1;
params.interval = ang;
params.search_num = ceil(360/params.interval);
params.vis_vec = 0;
params.vis = 0;
params.search_gopro = param_tracklet.search_gopro; %1: auto, 0:GT  2:Tracking res
% params.detection_method = 2; %1: GT, 2:Raw Detection  3: Detection with selection
params.Rho = param_tracklet.Rho; % 15; % plenty factor
params.x_thresh = 100; 
params.dis_thresh = 1000; 
params.Lambda = param_tracklet.Lambda; % 0.02; % weight for y defined euclidean distance 0.05
params.occ = 1;
params.occ_angle = param_tracklet.occ_angle; % 2
params.fast_search = 0;
params.search_ite = 6;
params.interval = 1;
% params.search_num = ceil(360/params.interval);
params.invis_dis = sqrt(2); % the distance between invisable object in the hor-view and other visable object

params.temp_angle = 1;
params.min_scl_degree = 0.1;   % the min degree search range
params.gma = 0.05; % Motion explotion Average parameter

tv_imgpath = data.cv_im_directory{1};
hv_imgpath = data.cv_im_directory{2};

Match_top_hor_seg = cell(length(tv_detection),1);
invis_top = cell(length(tv_detection),1);
if params.search_gopro == 2
    gopro_trk = load([param_tracklet.dataset_directory,'/Tra_Gopro/',param_tracklet.sceneName,'.mat']);
    gopro_trk = gopro_trk.results.res;
end

mha_data.img_top_1 =  imread([tv_imgpath,'/',tv_image_seg(1).name]);
mha_data.img_ego_1 = imread([hv_imgpath,'/',hv_image_seg(1).name]);
mha_data.size_img_hor = size(mha_data.img_ego_1);
mha_data.size_img_top = size(mha_data.img_top_1);

    for frm_i = 1 : length(tv_detection)
  
        tv_det = tv_detection{frm_i};
        hv_det = hv_detection{frm_i};
        mha_data.objs_top = tv_det;
        mha_data.objs_ego = hv_det;
        
        if params.vis == 1
            mha_data.img_top =  imread([tv_imgpath,'\',tv_image_seg(frm_i).name]);
            mha_data.img_ego = imread([hv_imgpath,'\',hv_image_seg(frm_i).name]);
        else
            mha_data.img_top  = [];
            mha_data.img_ego  = [];
        end
        
        img_hor_w = mha_data.size_img_hor(2);
        switch params.search_gopro 
            case 1
                cam_serch = tv_det;
            case 0 
                cam_serch = tv_det((tv_det(:,5)==0),:);
            case 2
                cam_serch = gopro_trk(data.frame_range(frm_i),:);
                % abs(sum(tv_det(:,1:2) - cam_serch(1:2),2))
                [~,cam_id] = min(sum(abs(tv_det(:,1:2) - cam_serch(1:2)),2));
        end
        
       	if params.temp_angle == 1
            min_scl = params.min_scl_degree;
            N_deg = floor(360*(((1-params.gma).^(data.frame_range(frm_i) -1) + min_scl)/(1 + min_scl)));
            params.search_num = ceil(N_deg/params.interval);
        else
            params.search_num = ceil(360/params.interval);
        end
        last_angle = data.last_angle;

        mha_data.cam_top = cam_serch;
        mha_data.vec_top = cell(params.search_num,1);
        mha_data.index_top = cell(params.search_num,1);  
        max_sco4cam = zeros(size(mha_data.cam_top,1),2);
        for cam_i = 1 : size(mha_data.cam_top,1)
        
            mha_data.cam_top = cam_serch(cam_i,:);
            
            [mha_data.vec_ego,mha_data.index_ego] = GoPro2vec(mha_data,params);
           
            mha_data = UAV2vec(params,mha_data,last_angle);
            
            [match_scores,mha_data,scale_y] = VecMatching_Ransec_DTW(mha_data,params);
           
            [max_sco,max_idx] = min(match_scores);
            scale_y_one = scale_y(max_idx);

            max_sco4cam(cam_i,:) = [max_sco,max_idx];           
        end
        
        [~,best_cam] = min(max_sco4cam(:,1));
        mha_data.cam_idx = best_cam;
        mha_data.cam_top = cam_serch(best_cam,:);  % the camera in top view
        if isempty(mha_data.cam_top)
            disp(strcat('NULL cam label:',num2str(frm_i)));
        end
        mha_data.angle_idx = max_sco4cam(best_cam,2); % the view angle of horizontal viewer
        mha_data.max_idx = mha_data.angle_idx;
        data.last_angle = pi/(params.search_num/2) * (mha_data.max_idx);
        
        [mha_data.vec_ego,mha_data.index_ego] = GoPro2vec(mha_data,params);

        mha_data = UAV2vec_one(mha_data.img_top,mha_data.objs_top,mha_data.cam_top,params,mha_data); 
                 
        vec1 = mha_data.vec_top_1;
        vec2 = mha_data.vec_ego;     
        
        [match_idx_top,match_idx_ego] = VecMatching_Ransec_DTW_one(mha_data,params,vec1,vec2);
        
%         max_len = max(length(tv_det),length(hv_det));
%         idx_top_all = ones(1,max_len)*(-1);
%         idx_hor_all = ones(1,max_len)*(-1);        
%         idx_top_all(1:length(tv_det)) = 1:length(tv_det);        
%         idx_hor_all(1:length(tv_det) = match_idx_ego;
        match_idx_top_all = mha_data.index_top_1(match_idx_top);
        match_idx_ego_all = mha_data.index_ego(match_idx_ego)';
        mha_data.visable_flag(cam_id) = -1;
        
        % invis_top = setdiff(1:length(tv_det),mha_data.index_top_1);
        % visualization for top_idx & ego_idx
        if params.vis_vec == 1
            Visualization(mha_data,vec_top_all,vec_hor_all);
        end
        
        invis_top{frm_i} = mha_data.visable_flag;
        D = [match_idx_top_all;match_idx_ego_all];
        Match_top_hor_seg{frm_i} =  D; 

    end

end
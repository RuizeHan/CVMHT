function Dis_top_hor_seg = box_associate_spatial(tv_detection,hv_detection,tv_image_seg,hv_image_seg)

ang = 1;
params.interval = ang;
params.search_num = ceil(360/params.interval);
params.vis_vec = 0;
params.vis = 0;
params.search_gopro = 0; %1: auto, 0:GT  2:Tracking res
params.detection_method = 2; %1: GT, 2:Raw Detection  3: Detection with selection
params.Rho = 15; % plenty factor
params.x_thresh = 100; 
params.dis_thresh = 1000; 
params.Lambda = 0.02; % weight for y defined euclidean distance 0.05
params.occ = 1;
params.occ_angle = 2;
params.fast_search = 0;
params.search_ite = 6;
params.interval = 1;
params.search_num = ceil(360/params.interval);
params.invis_dis = sqrt(2); % the distance between invisable object in the hor-view and other visable object
tv_imgpath = tv_image_seg(1).folder;
hv_imgpath = hv_image_seg(1).folder;

% if params.search_gopro == 2
%     gopro_trk = load([filePath,'Tra_Gopro\','.mat']);
%     gopro_trk = gopro_trk.results.res;
% end

Dis_top_hor_seg = cell(length(tv_detection),1);

mha_data.img_top_1 =  imread([tv_imgpath,'\',tv_image_seg(1).name]);
mha_data.img_ego_1 = imread([hv_imgpath,'\',hv_image_seg(1).name]);
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
            % case 3
                %cam_serch = gopro_trk(frame_i,:);
        end

        mha_data.cam_top = cam_serch;
        mha_data.vec_top = cell(params.search_num,1);
        mha_data.index_top = cell(params.search_num,1);  
        max_sco4cam = zeros(size(mha_data.cam_top,1),2);
        for cam_i = 1 : size(mha_data.cam_top,1)
        
            mha_data.cam_top = cam_serch(cam_i,:);
            
            [mha_data.vec_ego,mha_data.index_ego] = GoPro2vec(mha_data,params);
           
            mha_data = UAV2vec(mha_data.img_top,mha_data.objs_top,mha_data.cam_top,params,mha_data,mha_data.size_img_top);
            
            [match_scores,mha_data,scale_y] = VecMatching_Ransec_DTW(mha_data,params);
           
            [max_sco,max_idx] = min(match_scores);
            scale_y_one = scale_y(max_idx);

            max_sco4cam(cam_i,:) = [max_sco,max_idx];           
        end
        
        [~,best_cam] = min(max_sco4cam(:,1));
        mha_data.cam_idx = best_cam;
        mha_data.cam_top = cam_serch(best_cam,:);  % the camera in top view
        if isempty(mha_data.cam_top)
             mha_data.cam_top
        end
        mha_data.angle_idx = max_sco4cam(best_cam,2); % the view angle of horizontal viewer
        mha_data.max_idx = mha_data.angle_idx;
        
        [mha_data.vec_ego,mha_data.index_ego] = GoPro2vec(mha_data,params);

        mha_data = UAV2vec_one(mha_data.img_top,mha_data.objs_top,mha_data.cam_top,params,mha_data); 
             
        vec1 = mha_data.vec_top_1;
        vec2 = mha_data.vec_ego;
        if isempty(vec1)
             mha_data.cam_top
        end
        vec1(1,:) = vec1(1,:);
        vec2(1,:) = vec2(1,:)/(img_hor_w/2);
        %  normalizing
        vec1(2,:) = vec1(2,:)./scale_y_one;
        vec2(2,:) = 1./vec2(2,:);

        vec_top_all = zeros(2,length(tv_det));
        vec_hor_all = zeros(2,length(hv_det));
        vec_top_all(:,mha_data.index_top_1) = vec1;
        vec_hor_all(:,mha_data.index_ego) = vec2;

        invis_top = setdiff(1:length(tv_det),mha_data.index_top_1);
        % visualization for top_idx & ego_idx
        if params.vis_vec == 1
            Visualization(mha_data,vec_top_all,vec_hor_all);
        end
        
        D = pdist([vec_top_all';vec_hor_all']);
        Dis = squareform(D);
        Dis_top_hor = Dis(1:size(vec_top_all,2),size(vec_top_all,2)+1:size(vec_top_all,2)+size(vec_hor_all,2));
        Dis_top_hor(invis_top,:) = params.invis_dis;
        Dis_top_hor_seg{frm_i} =  Dis_top_hor; 

    end

end
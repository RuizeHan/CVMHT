function [tracklets_dis_matrix,OverLap,data] = tracklet_associate_spatial_vote(data,tv_tracklets,hv_tracklets,param_tracklet,param_cv_netCost)

punish = param_cv_netCost.spatial_punish;
frame_range = data.frame_range;
tv_detection = data.cv_detections{1,1}(frame_range);
hv_detection = data.cv_detections{2,1}(frame_range);
tv_image_seg = data.cv_images{1,1}(frame_range);
hv_image_seg = data.cv_images{2,1}(frame_range);

[Match_top_hor_seg,invis_top,data] = box_associate_spatial_bi(tv_detection,hv_detection,tv_image_seg,hv_image_seg,param_tracklet,data);

tracklets_dis_matrix = zeros(length(tv_tracklets),length(hv_tracklets));
OverLap = zeros(length(tv_tracklets),length(hv_tracklets));
for i = 1:length(tv_tracklets)

    for j = 1:length(hv_tracklets)
        
        tv_track = tv_tracklets(i);
        hv_track = hv_tracklets(j);
        
        tv_frms = tv_track.frame;
        hv_hrms = hv_track.frame;

        tv_OriDetectionInd = tv_track.origDetectionInd;
        hv_OriDetectionInd = hv_track.origDetectionInd;
        [overlap_frms, tv_overlap_idx, hv_overlap_idx] = intersect(tv_frms,hv_hrms);
       
        tv_over_det = tv_OriDetectionInd(tv_overlap_idx);
        hv_over_det = hv_OriDetectionInd(hv_overlap_idx);

        overlap_num = length(overlap_frms);
        n_match = 0;
        for over_i = 1 : overlap_num
            frm_i = overlap_frms(over_i)-frame_range(1)+1;
            mat_top_hor = Match_top_hor_seg{frm_i};
            tv_det_i = tv_over_det(over_i);
            hv_det_i = hv_over_det(over_i);
            tv_mat_idx = find(mat_top_hor(1,:) == tv_det_i);
            hv_mat_idx = find(mat_top_hor(2,:) == hv_det_i);
            if tv_mat_idx == hv_mat_idx
                n_match = n_match + 1;
            end
        end
      
        invis_cusum = 0;
        cam_cusum = 0;
        for frm_i = 1 : length(tv_frms)
            frm = tv_frms(frm_i) - frame_range(1)+1;
          
            invis_top_i = invis_top{frm};
            try
            if invis_top_i(tv_OriDetectionInd(frm_i)) == 0
            invis_cusum = invis_cusum + 1;
            elseif invis_top_i(tv_OriDetectionInd(frm_i)) == -1
            cam_cusum = cam_cusum + 1;
            end
            catch err
            end
        end

        cam_cusum_ratio = cam_cusum/length(tv_frms);
        
        match_ratio = n_match/length(tv_over_det);
        gamma = 1;
        tracklets_dis_matrix(i,j) = match_ratio ^ gamma;       
        OverLap(i,j) = overlap_num ;
    end
    if cam_cusum > 0

        tracklets_dis_matrix(i,:) = cam_cusum_ratio * (punish);
    end
end

end


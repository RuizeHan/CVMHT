function [tracklets_dis_matrix,OverLap] = tracklet_associate_spatial(frame_range,data,tv_tracklets,hv_tracklets)

tv_detection = data.cv_detections{1,1}(frame_range);
hv_detection = data.cv_detections{2,1}(frame_range);
tv_image_seg = data.cv_images{1,1}(frame_range);
hv_image_seg = data.cv_images{2,1}(frame_range);
Dis_top_hor_seg = box_associate_spatial(tv_detection,hv_detection,tv_image_seg,hv_image_seg);

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
        dis_track = 0;
        overlap_num = length(overlap_frms);
        for over_i = 1 : overlap_num
            frm_i = overlap_frms(over_i)-frame_range(1)+1;
            tv_det_i = tv_over_det(over_i);
            hv_det_i = hv_over_det(over_i);
            dis_top_hor_frm = Dis_top_hor_seg{frm_i};
            dis_track = dis_track + dis_top_hor_frm(tv_det_i,hv_det_i);
        end
               
        dis_track_avg = dis_track/overlap_num;
        tracklets_dis_matrix(i,j) = dis_track_avg;
        OverLap(i,j) = overlap_num;
    end

end

end
function [cv_midLevelTracklets, trackAll] = stitchTracklets(cv_midLevelTracklets)

maxAllowdDist = 5;
% Loop Over All tracks

tv_midLevelTracklets = cv_midLevelTracklets{1};
hv_midLevelTracklets = cv_midLevelTracklets{2};
% maxID = cell(1,2);

for iTrack = 1:length(tv_midLevelTracklets)-1
    
        tv_oldTracklet = tv_midLevelTracklets(iTrack);    
        hv_oldTracklet = hv_midLevelTracklets(iTrack);   
        tv_numOldTracklets = length(tv_oldTracklet.bbox);
        hv_numOldTracklets = length(hv_oldTracklet.bbox);           
        % Associating the tracklet ID between two views
        tv_flag = tv_midLevelTracklets(iTrack).flag; 
        tv_flag_sum = cumsum(tv_flag);
        hv_flag = hv_midLevelTracklets(iTrack).flag;
        hv_flag_sum = cumsum(hv_flag);
        AND = find(tv_flag .* hv_flag == 1);
        tv_flag_num = tv_flag_sum(AND);
        hv_flag_num = hv_flag_sum(AND);
        % Associating the tracklet ID between two views of next iTrack
        tv_flag_next = tv_midLevelTracklets(iTrack+1).flag; 
        tv_flag_sum_next = cumsum(tv_flag_next);
        hv_flag_next = hv_midLevelTracklets(iTrack+1).flag;
        hv_flag_sum_next = cumsum(hv_flag_next);
        AND = find(tv_flag_next .* hv_flag_next == 1);
        tv_flag_num_next = tv_flag_sum_next(AND);
        hv_flag_num_next = hv_flag_sum_next(AND);
               
        
    if iTrack == 1
        
        cv_midLevelTracklets{1}(1,iTrack).IDS = 1:tv_numOldTracklets;
        hv_IDS_idx = 1:hv_numOldTracklets;
        hv_IDS = 1:hv_numOldTracklets;
        hv_IDS(hv_flag_num) = tv_flag_num;
        other_idx = setdiff(hv_IDS_idx,hv_flag_num);
        other_ID = setdiff(hv_IDS_idx,hv_flag_num);
        hv_IDS(other_idx) = other_ID;      
        cv_midLevelTracklets{2}(1,iTrack).IDS = hv_IDS; 
        maxID = max(tv_numOldTracklets,hv_numOldTracklets);        

    end
    
    if  mod(iTrack,10) == 0 
        hv_IDS = cv_midLevelTracklets{2}(1,iTrack).IDS;
        clear cv_midLevelTracklets{2}(1,iTrack).IDS;
        tv_IDS = cv_midLevelTracklets{1}(iTrack).IDS;  
        hv_IDS_idx = 1:hv_numOldTracklets;
        hv_IDS(hv_flag_num) = tv_IDS(tv_flag_num);
        other_idx = setdiff(hv_IDS_idx,hv_flag_num);
        % other_ID = setdiff(hv_IDS,hv_flag_num);
        hv_IDS(other_idx) = hv_IDS(other_idx);    
        cv_midLevelTracklets{2}(1,iTrack).IDS = hv_IDS; 
    end

    for view_i = 1 : 2
        midLevelTracklets = cv_midLevelTracklets{view_i};     
        oldTracklet = midLevelTracklets(iTrack);
        newTracklet = midLevelTracklets(iTrack+1);    
        numOldTracklets = length(oldTracklet.bbox);
        numNewTracklets = length(newTracklet.bbox);
            
        avg_spatial_dist = zeros(numNewTracklets, numOldTracklets);
        IDSNew = ones(1,numNewTracklets)*-1;
        
        for iNewTracklet = 1:numNewTracklets
            for iOldTracklet = 1:numOldTracklets
                frameMaxOld = max(oldTracklet.frames{iOldTracklet});
                frameMinNew = min(newTracklet.frames{iNewTracklet});
                if(frameMinNew>frameMaxOld)
                    avg_spatial_dist(iNewTracklet, iOldTracklet)= Inf;
                    continue;
                end

                % Find frame interesection between two tracklets
                [~, ia, ib] = intersect(oldTracklet.frames{iOldTracklet},newTracklet.frames{iNewTracklet});
                posOld = [(oldTracklet.bbox{iOldTracklet}(ia,1)+oldTracklet.bbox{iOldTracklet}(ia,3))/2, ...
                    (oldTracklet.bbox{iOldTracklet}(ia,2)+oldTracklet.bbox{iOldTracklet}(ia,4))/2]';
                posNew = [(newTracklet.bbox{iNewTracklet}(ib,1)+newTracklet.bbox{iNewTracklet}(ib,3))/2, ...
                    (newTracklet.bbox{iNewTracklet}(ib,2)+newTracklet.bbox{iNewTracklet}(ib,4))/2]';

                avg_spatial_dist(iNewTracklet, iOldTracklet) = norm(posNew-posOld)/length(ia);

            end
        end   
    
    [val, ind] = min(avg_spatial_dist,[],2);
    [sortedVal, sortedInd] = sort(val);
    matchedIDS = [];
    
        for i = sortedInd'
            if(val(i)<maxAllowdDist)
                % Assoaciation
                matchedIDS = [matchedIDS, i];
                try
                IDSNew(1, i)= midLevelTracklets(iTrack).IDS(ind(i));
                catch err
                end
            end
        end

        nonMatched = setdiff(1:numNewTracklets,matchedIDS);
        % Assign the ID of the unmatched tracklets by the other view or a new
        % ID
        if (~isempty(nonMatched))
            if view_i == 1
                for nonmch_i = nonMatched           
    %                 [~,cv_idx,~] = intersect(tv_flag_num,nonmch_i);
    %                 if cv_idx
    %                     IDSNew(1,nonmch_i)= hv_flag_num(cv_idx);
    %                     if IDSNew(1,nonmch_i) >  maxID
    %                         maxID = maxID+1;    
    %                     end
    %                 else
                        IDSNew(1,nonmch_i)=maxID+1;      
                        maxID = maxID +1;
    %                 end
                end      
            else
                if mod(iTrack,10) ~= 0
                    for nonmch_i = nonMatched
                                                         
                        [~,cv_idx,~] = intersect(hv_flag_num_next,nonmch_i);
                        if cv_idx
                            tv_ID = cv_midLevelTracklets{1}(iTrack+1).IDS;                         
                            IDSNew(1,nonmch_i)= tv_ID(tv_flag_num_next(cv_idx));
                            if IDSNew(1,nonmch_i) >  maxID
                                maxID = maxID+1;    
                            end
                        else
                            IDSNew(1,nonmch_i) = maxID + 1;      
                            maxID = maxID + 1;
                        end
                    end
                end
            end
    %     for nonmch_i = nonMatched        
    %         IDSNew(1,nonmch_i)=maxID{view_i}+1;      
    %         maxID{view_i} = maxID{view_i}+1;
    %     end
        end
        
        cv_midLevelTracklets{view_i}(iTrack+1).IDS = IDSNew; 

%         if  mod(iTrack,10) == 0 && view_i == 2           
%            IDSNew = ones(1,numNewTracklets)*-1;
%            for hor_i = 1 : numNewTracklets           
%                 [~,cv_idx,~] = intersect(hv_flag_num,hor_i);
%                 if cv_idx
%                     tv_ID = cv_midLevelTracklets{1}(iTrack).IDS;
%                     IDSNew(1,hor_i)= tv_ID(tv_flag_num(cv_idx));
%                 else
%                     IDSNew(1,hor_i) = maxID + 1;      
%                     maxID = maxID + 1;
%                 end
%            end
%         end

   end
    
end

%% Interpolation
trackAll = cell(1,2);
for view_i = 1:2
    midLevelTracklets = cv_midLevelTracklets{view_i};
    tracks = cell(1,maxID);
    for iTrack = 1:length(midLevelTracklets)
        for iTracklet = 1:length(midLevelTracklets(1,iTrack).bbox)
            currID    = midLevelTracklets(1,iTrack).IDS(iTracklet);
            currBox   = midLevelTracklets(1,iTrack).bbox{iTracklet};
            currFrame = midLevelTracklets(1,iTrack).frames{iTracklet};
            try
            tracks{1, currID} = [tracks{1, currID}; currFrame currBox];  
            catch err
                continue;
            end
        end
    end

    tracksSmooth = cell(1,length(tracks));
    trackAll{view_i} = zeros(100000,6);
    counter = 1;
    for iTrack = 1:length(tracks)
       currTrack = tracks{1, iTrack};
       if ~isempty(currTrack)
       [uniqueFrames, indUniqueFrames] = unique(currTrack(:,1));
       currTrack = currTrack(indUniqueFrames, :);
       % Interpolation for the tracking res
%        minFrame = min(currTrack(:,1));
%        maxFrame = max(currTrack(:,1));
%        x1 = interp1(currTrack(:,1), currTrack(:,2), minFrame:maxFrame);
%        y1 = interp1(currTrack(:,1), currTrack(:,3), minFrame:maxFrame);
%        x2 = interp1(currTrack(:,1), currTrack(:,4), minFrame:maxFrame);
%        y2 = interp1(currTrack(:,1), currTrack(:,5), minFrame:maxFrame);
%        tracksSmooth{1,iTrack}=[(minFrame:maxFrame)', repmat(iTrack, maxFrame-minFrame+1,1), x1', y1', x2', y2'];
%        trackAll{view_i}(counter:(counter+maxFrame-minFrame),:)=tracksSmooth{1,iTrack};
       % -----------------------------------------
       trackAll{view_i}(counter:(counter+length(currTrack(:,1))-1),:)= [(currTrack(:,1)), repmat(iTrack, length(currTrack(:,1)),1),currTrack(:,2:5)];
       
       counter = counter + length(currTrack(:,1));
       end
    end
    trackAll{view_i}(counter:end,:)=[];
end






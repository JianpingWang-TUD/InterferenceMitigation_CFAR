function [ind_ava_seg, ind_0seg] = data_segmentation(det_map)

ind = find(det_map>0.5);  % indices of zeros
if isempty(ind)
    ind_ava_seg = [1, length(det_map)];
    ind_0seg    = [];
else
    ind_diff = ind(2:end)-ind(1:end-1);
    num_0seg = nnz(ind_diff>1)+1;    % number of segments of zeros
    if det_map(1)==det_map(end)
        if det_map(1)<0.5
            num_ava = num_0seg +1;
        else
            num_ava = num_0seg - 1;
        end
    else
        num_ava = num_0seg;
    end
    
    ind_0seg = zeros(num_0seg,2);
    ind_ava_seg = zeros(num_ava,2);
    
    ind_ava_cnt = 0;
    ind_0seg_cnt = 0;
    if det_map(1)<0.5
        ind_ava_seg(1,1)=1;
        ind_ava_cnt = 1;
    else
        ind_0seg(1,1) = 1;
        ind_0seg_cnt = 1;
    end
    
    if det_map(end)<0.5
        ind_ava_seg(end,end)=length(det_map);
    else
        ind_0seg(end,end) = length(det_map);
    end
    
    for kk = 1:length(det_map)-1
        if det_map(kk) < det_map(kk+1)
            ind_0seg_cnt = ind_0seg_cnt + 1;
            ind_ava_seg(ind_ava_cnt,2)=kk;
            ind_0seg(ind_0seg_cnt,1) = kk+1;
        elseif det_map(kk) > det_map(kk+1)
            ind_ava_cnt = ind_ava_cnt + 1;
            ind_ava_seg(ind_ava_cnt,1)=kk+1;
            ind_0seg(ind_0seg_cnt,2)=kk;
        else
            continue;
        end
    end
end
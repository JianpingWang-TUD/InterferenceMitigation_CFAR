function det_map_dlted = maskDilate(det_map, r)
% dilate the mask obtained with CFAR detection.
%

SE =strel('octagon',r);
% SE = strel('disk',r)
det_map_dlted = imdilate(det_map, SE);
import numpy as np
from typing import Any, List, Tuple, Union
import torch
from torch.nn import functional as F
from detectron2.structures.keypoints import Keypoints,_keypoints_to_heatmap,heatmaps_to_keypoints


def heatmaps_to_keypoints2(maps: torch.Tensor, rois: torch.Tensor, ninstances = None) -> torch.Tensor:

    if rois is None: 
      ninstances = ninstances
    else:
      ninstances = rois.shape[0]

    widths = maps.shape[2]
    heights = maps.shape[3]
    widths_ceil = widths
    heights_ceil = heights

    num_rois, num_keypoints = maps.shape[:2]
    xy_preds = maps.new_zeros(ninstances, num_keypoints, 4)

    keypoints_idx = torch.arange(num_keypoints, device=maps.device)

    for i in range(num_rois):
        outsize = (heights_ceil, widths_ceil)
        roi_map = F.interpolate(maps[[i]], size=outsize, mode="bicubic", align_corners=False)

        # Although semantically equivalent, `reshape` is used instead of `squeeze` due
        # to limitation during ONNX export of `squeeze` in scripting mode
        roi_map = roi_map.reshape(roi_map.shape[1:])  # keypoints x H x W

        # softmax over the spatial region
        max_score, _ = roi_map.view(num_keypoints, -1).max(1)
        max_score = max_score.view(num_keypoints, 1, 1)
        tmp_full_resolution = (roi_map - max_score).exp_()
        tmp_pool_resolution = (maps[i] - max_score).exp_()
        # Produce scores over the region H x W, but normalize with POOL_H x POOL_W,
        # so that the scores of objects of different absolute sizes will be more comparable
        roi_map_scores = tmp_full_resolution / tmp_pool_resolution.sum((1, 2), keepdim=True)

        w = roi_map.shape[2]
        pos = roi_map.view(num_keypoints, -1).argmax(1)

        x_int = pos % w
        y_int = (pos - x_int) // w

        assert (
            roi_map_scores[keypoints_idx, y_int, x_int]
            == roi_map_scores.view(num_keypoints, -1).max(1)[0]
        ).all()

        xy_preds[i, :, 0] = x_int
        xy_preds[i, :, 1] = y_int
        xy_preds[i, :, 2] = roi_map[keypoints_idx, y_int, x_int]
        xy_preds[i, :, 3] = roi_map_scores[keypoints_idx, y_int, x_int]

    return xy_preds

def heatmaps_modify_iterative(maps: torch.Tensor, rois: torch.Tensor, ninstances = None) -> torch.Tensor:
    penalty = -30
    nscale = 3
    w = maps.shape[2]
    h = maps.shape[3]
    nkps = maps.shape[1]
    grid = torch.tensor([(x, y) for x in range(w) for y in range(h)], device = maps.device)
    r2 = (max([h,w])/nscale)**2
    
    first_kp_pred = heatmaps_to_keypoints2(maps, rois, ninstances)
    first_best_kp_pos = torch.argmax(first_kp_pred[0][:,2])
    first_best_kp = first_kp_pred[0,first_best_kp_pos][0:2]

    region_mask = grid[
          ((grid[:,0] - int(first_best_kp[0]))**2 + 
          (grid[:,1] - int(first_best_kp[1]))**2) <  r2,
        :]

    flat_inds = [int(region_mask[i, 0] + region_mask[i, 1] * w) for i in range(region_mask.shape[0])]

    for i in range(nkps):
        if(i == first_best_kp_pos):
            continue
        maps.reshape((3,-1,))[i,[flat_inds]] = torch.tensor(
            penalty,
            device = maps.device,
            dtype = maps.dtype
        )

    second_kp_pred = heatmaps_to_keypoints2(maps, rois, ninstances)
    second_best_kp_pos = second_kp_pred[0][:,2].topk(2).indices[1]
    second_best_kp = second_kp_pred[0,second_best_kp_pos][0:2]

    region_mask = grid[
        ((grid[:,0] - int(second_best_kp[0]))**2 + 
        (grid[:,1] - int(second_best_kp[1]))**2) < r2,
      :]

    flat_inds = [int(region_mask[i, 0] + region_mask[i, 1] * w) for i in range(region_mask.shape[0])]

    for i in range(nkps):
        if(i == first_best_kp_pos or i == second_best_kp_pos):
            continue
        maps.reshape((3,-1,))[i,[flat_inds]] = torch.tensor(
            penalty,
            device = maps.device,
            dtype = maps.dtype
        )
    return maps



def heatmaps_to_keypoints_iterative(maps: torch.Tensor, rois: torch.Tensor, ninstances = None) -> torch.Tensor:
    try: maps = heatmaps_modify_iterative(maps, rois, ninstances)
    except: maps = maps
    kypts = heatmaps_to_keypoints(maps, rois)
    return kypts





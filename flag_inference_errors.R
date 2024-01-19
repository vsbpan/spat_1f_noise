source("helper_functions/init.R")


d <- fetch_events(93)

z <- fetch_data_dict(93)

attr(z, "is_rep") <- TRUE

z[1029] %>% 
  plot_mask(frame = 1)
z[1029] %>% 
  plot_keypoint(frame = 1)

z[1030] %>% 
  plot_mask(frame = 1)
z[1030] %>% 
  plot_keypoint(frame = 1)

plot(z, 1030, mask_col = "white")

for (i in 800:1200){
  plot(z, i, mask_col = "white", add = FALSE, main = i)
  cat(sprintf("Rank: %s \r", i))
  Sys.sleep(0.1)
}


f <- function(i){
  if(missing(i)){
    i <- get("i", envir = globalenv())
  }
  plot(z, i, mask_col = "white", add = FALSE, main = i)
  cat(sprintf("Rank: %s \r", i))
  i <<- i + 1
}

f <- forward_plot(z, 800)

f()


# 803 - 805

get_mask(z, 803:805) %>% imappend("z") %>% mask_info()

plot(z, 803, mask_col = "white")

a <- get_mask(z, 800:820) %>% imappend("z") %>% mask_info()

f()


move_seq(a$centroid_x, a$centroid_y)$r > 2

moved_centroid <- move_seq(a$centroid_x, a$centroid_y)$r > 10
moved_head <- move_seq(b$head_x, b$head_y)$r > 10


moved_head * !moved_centroid 


b <- get_keypoints(z[800:820])

f <- forward_plot(z, 800)
f()

iou <- get_mask(z, 800:820) %>% 
  lapply(seq_len(length(.)-1), function(i, l){
    mask_IOU(l[[i]], l[[i+1]], use_C = FALSE)
  }, l = .) %>% 
  do.call("c", .)



moved_head * !(iou < 0.5)



#spatstat.geom:::
?Area.xypolygon


get_polygon(z[100])[[1]] %>%
  as.data.frame() -> g

spatstat.utils::Area.xypolygon(
  list("x" = rev(g$x), "y" = rev(g$y))
)

mask_area(z)



get_centroid()


mask_area(z[100])




i <- 100
plot(z, i, kp_name = "head")
w <- get_centroid2(get_polygon(z[i])[[1]])
points(w[1], w[2], col = "white", pch = 19)

polygon2mask(get_polygon(z[i])[[1]], mini_mask = T, raw_mat = T)

polygon2mask

z[1200]
p <- get_polygon(z[i])[[1]]

mask_area(p)


get_polygon(z[1:3]) %>% length()
mask_info(z[11121])

polygon2mask(get_polygon(z[1])[[1]], raw_mat = TRUE, mini_mask = TRUE)







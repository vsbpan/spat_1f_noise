
img <- image_example()

dim(img)

plot(img)

q2q_trans(
  img, 
  dest_pts = list(
    list(1,1),
    list(1, 200*3),
    list(400*2, 200*2),
    list(400*2, 1)
  )) %>% 
    plot()



q2q_trans(
  img, 
  dest_pts = list(
    list(1,1),
    list(1, 200),
    list(500, 300),
    list(400, 1)
  )
) %>% 
  q2q_trans(
    ., 
    init_pts = list(
      list(1,1),
      list(1, 200),
      list(500, 300),
      list(400, 1)
    ), 
    dest_pts = list(
      list(1,1),
      list(1, 443*2),
      list(810*2, 443*2),
      list(810*2, 1)
    )
  ) %>% 
  plot()


plot(img)





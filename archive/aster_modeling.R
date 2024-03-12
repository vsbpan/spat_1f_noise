source("spat1f/init_analysis.R")
# 
# 
# data("radish")
# 
# levels(radish$varb)
# pred <- c(0,1,2)
# fam <- c(1,3,2)
# 
# rout <- reaster(resp ~ varb + fit : (Site * Region),
#                 list(block = ~ 0 + fit : Block, pop = ~ 0 + fit : Pop),
#                 pred, fam, varb, id, root, data = radish)
# out <- aster(resp ~ varb + fit:(Site*Region) + fit:Pop + fit:Block,
#                 pred, fam, varb, id, root, data = radish)
# 
# summary(rout)
# radish$fit
# 
# summary(out, show.graph = TRUE)
# 
# 
# 
# 
# z <- ref_data[50,]
# 
# 
# z$date_start
# z$pupation_date
# 
# 
# ref_data$date_start


foo <- function(RGR){
  dead <- is.na(RGR)
  RGR <- ifelse(dead, NA, RGR)
  
  data.frame("node1_l5" = as.numeric(!dead), "node2_RGR" = RGR)
}


a<- ref_data %>% 
  rowwise() %>% 
  mutate(
    foo(RGR)
  ) %>% 
  reshape_life_events()

pred <- c(0, 1)
fam <- c(1, 2)
a$u_j <- as.numeric(grepl("node2_RGR", a$node))
fam_list <- list(
  aster::fam.bernoulli(), 
  aster::fam.normal.location(0.005)
)


b <- a %>% 
  filter(
    var_trt != "constant")
  # ) %>% 
  # group_by(
  #   id
  # ) %>%
  # mutate(
  #   has_na = any(is.na(y))
  # ) %>%
  # filter(
  #   !has_na
  # ) %>%
  # ungroup()


aster_fit <- reaster(
  y ~ 
    node + u_j : (scale(log(cat_pre_wt)) * (var_trt + beta)),
  list(session_id = ~ 0 + u_j : session_id),
  pred = pred,
  fam = fam,
  varvar = node,
  idvar = id, 
  root = root, 
  data = b, 
  famlist = fam_list
)


summary(aster_fit, info.tol = 1e-20)





ceiling(range((na.omit(ref_data$pupation_time)) / 10))
ndays <- 4


a <- ref_data %>% 
  rowwise() %>% 
  mutate(
    life_events_flatten(parse_life_events(date_start, 
                       pupation = pupation_date, 
                       death = death_date, 
                       pup_wt = pupal_weight, 
                       t_range = c(1,ndays)))
  ) %>% 
  reshape_life_events()


# Levels
all.equal(levels(a$node), sprintf("node%03d%s%02d",
        seq_len(ndays*3),
        c(
          rep("l", ndays),
          rep("p", ndays),
          rep("m", ndays)
        ), 
        rep(1:ndays, 3)))


# pred
pred <- c(
  c(0:(ndays - 1)), 
  c(1:ndays), 
  c((ndays+1):(2*ndays))
)

# Family list
fam_list <- list(
  aster::fam.bernoulli(), 
  aster::fam.normal.location(0.05)
)

fam <- c(
  rep(1, ndays),
  rep(1, ndays),
  rep(2, ndays)
)



a$u_j <- as.numeric(grepl("node[0-9][0-9][0-9]m", a$node))


b <- a %>% 
  filter(
    var_trt != "constant"
  ) %>% 
  group_by(
    id
  ) %>% 
  mutate(
    has_na = any(is.na(y))
  ) %>% 
  filter(
    !has_na
  ) %>% 
  ungroup()

aster_fit <- aster(
  y ~ 
    node + u_j : (scale(log(cat_pre_wt)) + (var_trt + beta)),
  #list(session_id = ~ 0 + u_j : session_id),
  pred = pred,
  fam = fam,
  varvar = node,
  idvar = id, 
  root = root, 
  data = b, 
  famlist = fam_list
)

summary(aster_fit, info.tol = 1e-20)












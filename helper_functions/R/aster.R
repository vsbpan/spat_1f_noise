# parse_life_events <- function(start, pupation, death, pup_wt, t_range = c(1, 50)){
#   pup_day <- as.numeric(round(as.difftime(pupation - start, units = "d", format = "%Y-%m-%d")))
#   death_day <- as.numeric(round(as.difftime(death - start, units = "d", format = "%Y-%m-%d")))
#   #censor_day <- as.numeric(round(as.difftime(censor - start, units = "d", format = "%Y-%m-%d")))
#   
#   pup_day <- ceiling(pup_day / 10)
#   death_day <- ceiling(death_day / 10)
#   
#   t <- seq(t_range[1], t_range[2], by = 1)
#   nt <- length(t)
#   
#   if(!is.na(death_day) || !is.na(pup_day)){
#     
#     if(is.na(death_day)){
#       l <- rep(0, nt)
#     } else {
#       l <- ifelse(t < death_day, 1, 0)
#     }
#     
#     if(is.na(pup_day)){
#       p <- rep(0, nt)
#       m <- p
#     } else {
#       p <- ifelse(t == pup_day, 1, 0)
#       l <- ifelse(t <= pup_day, 1, 0)
#       m <- ifelse(t == pup_day, pup_wt, 0)
#     }
#     
#     
#     out <- data.frame(
#       "t" = t,
#       "l" = l, # survival as larva
#       "p" = p, # If survived, pupated?
#       "m" = m
#     )
#     return(out)
#   } else{
#     return(
#       data.frame(
#         "t" = t,
#         "l" = NA,
#         "p" = NA,
#         "m" = NA
#       )
#     )
#   }
# }
# 
# life_events_flatten <- function(x){
#   n <- nrow(x)
#   
#   out <- x %>% 
#     gather(key = node, value = val, l:m) %>% 
#     mutate(node_id = case_when(
#       node == "l" ~ 0,
#       node == "p" ~ 1,
#       node == "m" ~ 2
#     ) * n + t) %>% 
#     mutate(node = sprintf("node%03d%s%02d",node_id,node, t)) %>% 
#     dplyr::select(node, val) %>% 
#     spread(
#       value = val, key = node
#     )
#   
#   # node ID 
#   # node %index% %node type% %time index% 
#   return(out)
# }
# 
# 
# reshape_life_events <- function(x){
#   name <- names(x)
#   name <- name[grepl("node",name)]
#   x <- as.data.frame(x)
#   reshape(x, 
#           varying = list(name), 
#           direction = "long", 
#           "timevar" = "node",
#           "times" = as.factor(name), 
#           v.names = "y") %>% 
#     mutate(
#       root = 1
#     )
# }

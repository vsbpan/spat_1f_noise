source("spat1f/init.R")





d <- read_csv("raw_data/T ni conditioning study.csv") %>% 
  filter(Trial!= "58")
d <- d %>% 
  mutate_if(is.character, .funs = tolower) %>% 
  mutate(
    first_diet = ifelse(is.na(first_diet), "n", first_diet),
    second_diet = ifelse(is.na(second_diet), "n", second_diet),
    cat_wt_log = log(cat_wt),
    cat_wt_log_scale = as.numeric(scale(log(cat_wt)))
  )
d_long <- d %>% 
  gather(
    key = type, value = choice, c(first_diet,second_diet,end_of_day_diet,next_morning_choice)
  ) %>% 
  mutate(
    choice_h = ifelse(choice == "h", 1, 0),
    choice_l = ifelse(choice == "l", 1, 0),
    choice_n = ifelse(choice == "n", 1, 0),
    choice_c = ifelse(choice == "c", 1, 0)
  )





d %>% 
  ggplot(aes(x = log(cat_wt), color = condition_trt, y = choice_2h)) + 
  geom_point() + 
  facet_wrap(~center_block_trt) + 
  geom_smooth(method = "lm")




d %>% 
  ggplot(aes(x = bin(log(cat_wt), by = 2), color = condition_trt, y = first_diet)) + 
  geom_count(alpha = 0.5) + 
  facet_wrap(~center_block_trt)








d %>% 
  ggplot(aes(x = bin(log(cat_wt), by = 2), color = condition_trt, y = end_of_day_diet)) + 
  geom_count(alpha = 0.5) + 
  facet_wrap(~center_block_trt)


d %>% 
  ggplot(aes(x = bin(log(cat_wt), by = 2), color = condition_trt, y = next_morning_choice)) + 
  geom_count(alpha = 0.5) + 
  facet_wrap(~center_block_trt)


names(d)

library(brms)


d



d_long$y <- with(d_long, cbind(choice_h, choice_l, choice_n, choice_c))
names(d_long)

get_prior(
  choice ~
    type * (cat_wt_log_scale * 
              center_block_trt * condition_trt) + 
    (1|Trial),
  data = d_long,
  family = categorical(link = "logit")
)

m <- brm(
  choice ~
    cat_wt_log_scale * 
    center_block_trt * condition_trt,
  data = d_long %>% 
    filter(type == "first_diet"),
  family = categorical(link = "logit"),
  prior = c(
    set_prior("normal(0, 3)", class = "b"),
    set_prior("normal(0, 2.5)", class = "Intercept")
  ),
  cores = 6,
  iter = 4000,
  chains = 6,
  control = list(adapt_delta = 0.99)
)
#saveRDS(m, "invisible/fitted_models/first_diet_m.rds")

plot(m)
check_brms(m)


m2 <- brm(
  choice ~
    cat_wt_log_scale * 
    center_block_trt * condition_trt,
  data = d_long %>% 
    filter(type == "end_of_day_diet"),
  family = categorical(link = "logit"),
  prior = c(
    set_prior("normal(0, 3)", class = "b"),
    set_prior("normal(0, 2.5)", class = "Intercept")
  ),
  cores = 6,
  iter = 4000,
  chains = 6,
  control = list(adapt_delta = 0.99)
)
saveRDS(m2, "invisible/fitted_models/end_of_day_diet_m.rds")

m3 <- brm(
  choice ~
    cat_wt_log_scale * 
    center_block_trt * condition_trt,
  data = d_long %>% 
    filter(type == "next_morning_choice"),
  family = categorical(link = "logit"),
  prior = c(
    set_prior("normal(0, 3)", class = "b"),
    set_prior("normal(0, 2.5)", class = "Intercept")
  ),
  cores = 6,
  iter = 4000,
  chains = 6,
  control = list(adapt_delta = 0.99)
)
saveRDS(m3, "invisible/fitted_models/next_morning_choice_m.rds")


m4 <- brm(
  choice ~
    cat_wt_log_scale * 
    center_block_trt * condition_trt,
  data = d_long %>% 
    filter(type == "second_diet"),
  family = categorical(link = "logit"),
  prior = c(
    set_prior("normal(0, 3)", class = "b"),
    set_prior("normal(0, 2.5)", class = "Intercept")
  ),
  cores = 6,
  iter = 4000,
  chains = 6,
  control = list(adapt_delta = 0.99)
)
saveRDS(m4, "invisible/fitted_models/second_choice_m.rds")


marginal_effects(m, terms = c("cat_wt_log_scale", "center_block_trt","condition_trt"))


sjPlot::plot_model(m, type = "pred", 
                   terms = c("cat_wt_log_scale", "condition_trt"))

sjPlot::plot_model(m, type = "pred", 
                   terms = c("cat_wt_log_scale", "center_block_trt"))

sjPlot::plot_model(m, type = "pred", 
                   terms = c("condition_trt", "center_block_trt"))

















m2

m <- readRDS("invisible/fitted_models/first_diet_m.rds")



m1_pred2 <- get_pred(m, 
                         terms = c("cat_wt_log_scale", "center_block_trt","condition_trt"), 
                         along_n = 100, 
                     ndraws = 4000)

g <- m1_pred2 %>% 
  filter(center_block_trt != "y") %>% 
  filter(resp != "n") %>% 
  group_by(
    iter, cat_wt_log_scale, center_block_trt, condition_trt
  ) %>% 
  mutate(
    y_hat = y_hat / sum(y_hat)
  ) %>% 
  filter(
    resp == "l"
  ) %>% 
  mutate(
    condition_trt = ifelse(condition_trt == "h", "Reared on 2 mg/g", "Reared on 0 mg/g")
  ) %>% 
  ggplot(aes(x = exp(cat_wt_log_scale), y = y_hat)) +
  facet_wrap(~ condition_trt) + 
  tidybayes::stat_lineribbon() + 
  scale_fill_brewer() + 
  theme_bw(base_size = 15) + 
  scale_x_log10(label = fancy_scientificb) + 
  geom_hline(color = "violetred", yintercept = 0.2, size = 2, linetype = "dashed") + 
  labs(x = "Pre-weight (g)", y = "First choice is less toxic diet")




ggsave("graphs/first_choice_conditioning.png", g, dpi = 600, width = 8, height = 5)


m1_pred2 <- get_pred(m, 
                     terms = c("cat_wt_log_scale", "center_block_trt","condition_trt"), 
                     along_n = 50, 
                     ndraws = 2000)

m1_pred2 %>% 
  group_by(
    iter, cat_wt_log_scale, center_block_trt, condition_trt
  ) %>% 
  mutate(
    y_hat = 1 - y_hat / sum(y_hat)
  ) %>% 
  filter(
    resp == "n"
  ) %>% 
  mutate(
    condition_trt = ifelse(condition_trt == "h", "Reared on 2 mg/g", "Reared on 0 mg/g"),
    center_block_trt = ifelse(center_block_trt == "y", "Has center block", "No center block")
  ) %>% 
  group_by(iter, cat_wt_log_scale, center_block_trt) %>% 
  summarise(
    y_hat = diff(y_hat)
  )
  slice_sample(n = 5000) %>% 
  ggplot(aes(x = exp(cat_wt_log_scale), y = y_hat)) +
  facet_grid(center_block_trt ~ condition_trt) + 
  tidybayes::stat_lineribbon() + 
  scale_fill_brewer() + 
  theme_bw(base_size = 15) + 
  scale_x_log10(label = fancy_scientificb) + 
  labs(x = "Pre-weight (g)", y = "First choice is less toxic diet")




a <- m1_pred2 %>% 
  group_by(
    iter, cat_wt_log_scale, center_block_trt, condition_trt
  ) %>% 
  mutate(
    y_hat = 1 - y_hat / sum(y_hat)
  ) %>% 
  filter(
    resp == "n"
  ) %>% 
  mutate(
    condition_trt = ifelse(condition_trt == "h", "Reared on 2 mg/g", "Reared on 0 mg/g"),
    center_block_trt = ifelse(center_block_trt == "y", "Has center block", "No center block")
  ) %>% 
  group_by(iter, cat_wt_log_scale, center_block_trt) %>% 
  summarise(
    y_hat = diff(y_hat)
  )



g2 <- a %>% 
  ggplot(aes(x = exp(cat_wt_log_scale), y = y_hat)) +
  facet_wrap(~center_block_trt) + 
  tidybayes::stat_lineribbon() + 
  scale_fill_brewer() + 
  theme_bw(base_size = 15) + 
  scale_x_log10(label = fancy_scientificb) + 
  geom_hline(color = "violetred", yintercept = 0, size = 2, linetype = "dashed") + 
  labs(x = "Pre-weight (g)", y = "Difference in rate of leaving (Prop.)")

ggsave("graphs/first_choice_difference_conditioning.png", g2, dpi = 600, width = 8, height = 5)





m1_pred <- get_mean_pred(m, 
              terms = c("cat_wt_log_scale", "center_block_trt","condition_trt"), 
              along_n = 100)

g <- m1_pred %>% 
  mutate(
    center_block_trt = paste0("center_block: ", toupper(center_block_trt)),
    condition_trt = paste0("condition: ", toupper(condition_trt))
  ) %>% 
  ggplot(aes(x = exp(cat_wt_log_scale), y = y_hat, fill = resp)) + 
  facet_wrap(~center_block_trt+condition_trt) + 
  geom_area() + 
  scale_fill_viridis_d() + 
  theme_minimal(base_size = 15) + 
  scale_x_continuous(trans = "log10") + 
  labs(x = "Cat weight (g)", y = "Proportion", fill = "Choice", subtitle = "First Choice") + 
  theme(legend.position = "right")



m2_pred <- get_mean_pred(m2, 
                         terms = c("cat_wt_log_scale", "center_block_trt","condition_trt"), 
                         along_n = 100)

g <- m2_pred %>% 
  mutate(
    center_block_trt = paste0("center_block: ", toupper(center_block_trt)),
    condition_trt = paste0("condition: ", toupper(condition_trt))
  ) %>% 
  mutate(
    resp = factor(resp, levels = c("h","l","n","c"))
  ) %>% 
  ggplot(aes(x = exp(cat_wt_log_scale), y = y_hat, fill = resp)) + 
  facet_wrap(~center_block_trt+condition_trt) + 
  geom_area() + 
  scale_fill_viridis_d() + 
  theme_minimal(base_size = 15) + 
  scale_x_continuous(trans = "log10") + 
  labs(x = "Cat weight (g)", y = "Proportion", fill = "Choice", subtitle = "End of Day Diet") + 
  theme(legend.position = "right")


ggsave("graphs/end_of_day_choice.png", g, dpi = 400, width = 7, height = 7, bg = "white")





m3_pred <- get_mean_pred(m3, 
                         terms = c("cat_wt_log_scale", "center_block_trt","condition_trt"), 
                         along_n = 100)

g <- m3_pred %>% 
  mutate(
    center_block_trt = paste0("center_block: ", toupper(center_block_trt)),
    condition_trt = paste0("condition: ", toupper(condition_trt))
  ) %>% 
  mutate(
    resp = factor(resp, levels = c("h","l","n","c"))
  ) %>% 
  ggplot(aes(x = exp(cat_wt_log_scale), y = y_hat, fill = resp)) + 
  facet_wrap(~center_block_trt+condition_trt) + 
  geom_area() + 
  scale_fill_viridis_d() + 
  theme_minimal(base_size = 15) + 
  scale_x_continuous(trans = "log10") + 
  labs(x = "Cat weight (g)", y = "Proportion", fill = "Choice", subtitle = "Next Day Diet") + 
  theme(legend.position = "right")



ggsave("graphs/next_day_choice.png", g, dpi = 400, width = 7, height = 7, bg = "white")




m4_pred <- get_mean_pred(m4, 
                         terms = c("cat_wt_log_scale", "center_block_trt","condition_trt"), 
                         along_n = 100)

g <- m4_pred %>% 
  mutate(
    center_block_trt = paste0("center_block: ", toupper(center_block_trt)),
    condition_trt = paste0("condition: ", toupper(condition_trt))
  ) %>% 
  mutate(
    resp = factor(resp, levels = c("h","l","n","c"))
  ) %>% 
  ggplot(aes(x = exp(cat_wt_log_scale), y = y_hat, fill = resp)) + 
  facet_wrap(~center_block_trt+condition_trt) + 
  geom_area() + 
  scale_fill_viridis_d() + 
  theme_minimal(base_size = 15) + 
  scale_x_continuous(trans = "log10") + 
  labs(x = "Cat weight (g)", y = "Proportion", fill = "Choice", subtitle = "Second Choice") + 
  theme(legend.position = "right")



ggsave("graphs/second_choice.png", g, dpi = 400, width = 7, height = 7, bg = "white")




























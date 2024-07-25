source("spat1f/init.R")

# Load data
d <- read_csv("raw_data/T ni conditioning study.csv") %>% 
  filter(Trial!= "58") %>% # Caterpillar dead on arrival 
  mutate_if(is.character, .funs = tolower) %>% 
  mutate(
    first_diet = ifelse(is.na(first_diet), "n", first_diet),
    second_diet = ifelse(is.na(second_diet), "n", second_diet),
    cat_wt_log = log(cat_wt),
    cat_wt_log_scale = as.numeric(scale(log(cat_wt)))
  )


# Fit categorical glm
m <- brm(
  first_diet ~
    cat_wt_log_scale * center_block_trt * condition_trt,
  data = d,
  family = categorical(link = "logit"),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),
    set_prior("normal(0, 2)", class = "Intercept")
  ),
  cores = 6,
  iter = 4000,
  chains = 6,
  control = list(adapt_delta = 0.99)
)

# saveRDS(m, "invisible/fitted_models/first_diet_m.rds")
m <- readRDS("invisible/fitted_models/first_diet_m.rds")


# Check convergence and residuals
plot(m)
check_brms(m)


# Get posterior draws
m1_pred <- pred_draws2(m, 
                       terms = c("cat_wt_log_scale", "center_block_trt","condition_trt"), 
                       along_n = 100, 
                       ndraws = 3000)

# Compute the proportion of first choice that is the less toxic diet for the group without a center block
b <- m1_pred %>% 
  filter(center_block_trt != "y") %>% 
  filter(resp != "n") %>% 
  group_by(
    draw, cat_wt_log_scale, center_block_trt, condition_trt
  ) %>% 
  mutate(
    y_hat = y_hat / sum(y_hat)
  ) %>% 
  filter(
    resp == "l"
  ) %>% 
  mutate(
    condition_trt = ifelse(condition_trt == "h", "Reared on 2 mg/g", "Reared on 0 mg/g")
  ) 

g <- b %>% 
  ggplot(aes(x = exp(unscale(d$cat_wt_log)(cat_wt_log_scale)), y = y_hat)) +
  facet_wrap(~ condition_trt) + 
  tidybayes::stat_lineribbon(alpha = 0.5) + 
  scale_fill_brewer() + 
  theme_bw(base_size = 15) + 
  scale_x_log10(label = fancy_scientific) + 
  geom_hline(color = "violetred", yintercept = 0.2, size = 2, linetype = "dashed") + 
  labs(x = "Caterpillar weight (g)", y = "First choice is less toxic diet")


ggsave("graphs/first_choice_conditioning.png", g, dpi = 600, width = 8, height = 5)


# Compute the marginal effect of having prior toxin experience on whether the caterpillar leaves the center to find a first choice diet block
a <- m1_pred %>% 
  group_by(
    draw, cat_wt_log_scale, center_block_trt, condition_trt
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
  group_by(draw, cat_wt_log_scale, center_block_trt) %>% 
  summarise(
    y_hat = diff(y_hat)
  )

g2 <- a %>% 
  ggplot(aes(x = exp(unscale(d$cat_wt_log)(cat_wt_log_scale)), y = y_hat)) +
  facet_wrap(~center_block_trt) + 
  tidybayes::stat_lineribbon(alpha = 0.5) + 
  scale_fill_brewer() + 
  theme_bw(base_size = 15) + 
  scale_x_log10(label = fancy_scientific) + 
  geom_hline(color = "violetred", yintercept = 0, size = 2, linetype = "dashed") + 
  labs(x = "Caterpillar weight (g)", 
       y = "Marginal effect of prior toxin experience \non rate of leaving (proportion)")

ggsave("graphs/first_choice_difference_conditioning.png", g2, dpi = 600, width = 8, height = 5)




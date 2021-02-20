

# description -------------------------------------------------------------

# Exploring improvements to the general idea of "predict B/Bmsy in new fisheries as a function of catch history and life history"


# setup -------------------------------------------------------------------



library(tidyverse)
library(rstan)
library(rstanarm)
library(patchwork)
library(scales)
library(here)
library(recipes)
library(kernlab)
library(tidymodels)
library(ranger)
library(FishLife)
library(janitor)
library(RcppRoll)

options(mc.cores = parallel::detectCores() / 2)

rstan_options(auto_write = TRUE)

theme_set(theme_classic() + theme(strip.background = element_rect(color = "transparent")))



# load RAM data -----------------------------------------------------------

min_years_catch <- 20

draws <- 3000

crazy_b <- 5

min_draws <- 2000 # minimum number of unique SIR draws

n_cores <- 6 # number of cores for parallel processing

lookup_fmi_names <- FALSE

future::plan(future::multiprocess, workers = n_cores)

# data(Return)

# return is from here, WARNING, changes rapidly, things break check and make sure this isn't why
# https://drive.google.com/drive/u/0/folders/1J46tM6PYDdPwhx5zGrlHMdxUyGRrky7X?ogsrc=32

functions <- list.files(here::here("R"))

functions <- functions[!functions %in% c("zzz.R", "sysdata.rda")]

purrr::walk(functions, ~ source(here::here("R", .x)))

# load data ---------------------------------------------------------------


# load(here("data-raw","Return.Rdata"))

# readr::write_rds(Return,here::here("data-raw","Return.rds"))
#
# rm(Return)

# Return <- readr::read_rds(here::here("data-raw","Return.rds"))


# FishLifeData<- Return[c("ParentChild_gz","beta_gv","Cov_gvv")]
#
# FishLifeData$metadata <- "emailed from Thorson, beta version with newparameters"


if (!file.exists(here("data-raw", "ram.zip"))) {
  download.file(
    "https://www.dropbox.com/s/jpgz0a5s5of3qev/RAM%20v4.491%20Files%20(1-14-20).zip?dl=1",
    destfile = here::here("data-raw", "ram.zip"),
    mode = "wb"
  )
  
  unzip(here::here("data-raw", "ram.zip"), exdir = "data-raw")
  
}



# process RAM data --------------------------------------------------------

ram_dirs <- list.files("data-raw")

ram_dirs <- ram_dirs[str_detect(ram_dirs, "RAM v\\d")]

ram_files <-
  list.files(file.path("data-raw", ram_dirs), recursive = TRUE)

ram_files <- ram_files[str_detect(ram_files, ".RData")]

ram_files <- ram_files[str_detect(ram_files, "Model Fit")]

load(file.path("data-raw", ram_dirs, ram_files[1]))

stock <- stock %>%
  left_join(area, by = "areaid")
# catches
ram_catches <- tcbest.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_tibble() %>%
  gather(stockid, catch, -year)

# B/Bmsy
ram_b_v_bmsy <- divbpref.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_tibble() %>%
  gather(stockid, b_v_bmsy, -year)

# U/Umsy
ram_u_v_umsy <- divupref.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_tibble() %>%
  gather(stockid, u_v_umsy, -year)

# Effort
ram_effort <- effort.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_tibble() %>%
  gather(stockid, effort, -year)

# biomass


ram_total_biomass <- tbbest.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_tibble() %>%
  gather(stockid, total_biomass, -year)

# ssb

ram_ss_biomass <- ssb.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_tibble() %>%
  gather(stockid, ss_biomass, -year)


ram_exp_rate <- ram_catches %>%
  left_join(ram_total_biomass, by = c("stockid", "year")) %>%
  mutate(exploitation_rate = catch / total_biomass) %>%
  select(-catch,-total_biomass)


# put it together

ram_data <- ram_catches %>%
  left_join(bioparams_values_views, by = "stockid") %>%
  left_join(ram_b_v_bmsy, by = c("stockid", "year")) %>%
  left_join(ram_u_v_umsy, by = c("stockid", "year")) %>%
  left_join(ram_exp_rate, by = c("stockid", "year")) %>%
  left_join(ram_effort, by = c("stockid", "year")) %>%
  left_join(ram_total_biomass, by = c("stockid", "year")) %>%
  left_join(ram_ss_biomass, by = c("stockid", "year")) %>%
  left_join(stock, by = "stockid") %>%
  select(stockid, scientificname, commonname, everything())


# create new variables

ram_data <- ram_data %>%
  mutate(tb_v_tb0 = total_biomass / TB0,
         ssb_v_ssb0 = ss_biomass / SSB0)

# filter data

ram_data <- ram_data %>%
  filter(is.na(catch) == FALSE) %>%
  group_by(stockid) %>%
  mutate(
    has_tb0 = !all(is.na(TB0)),
    has_tb = !all(is.na(total_biomass)),
    first_catch_year = year[which(catch > 0)[1]],
    b_v_bmsy = pmin(b_v_bmsy, crazy_b)
  ) %>%
  mutate(
    pchange_effort = lead(u_v_umsy) / (u_v_umsy + 1e-6),
    cs_effort = (u_v_umsy - mean(u_v_umsy)) / sd(u_v_umsy),
    b_rel = dplyr::case_when(
      has_tb0 ~ total_biomass / max(TB0),
      has_tb ~ total_biomass / max(total_biomass),
      TRUE ~ b_v_bmsy / 2
    )
  ) %>%
  ungroup()

# classify stocks by stock history shape  ----------------------------------------------

ram_catches <- ram_data %>%
  ungroup() %>%
  select(stockid, year, catch) %>%
  group_by(stockid) %>%
  mutate(stock_year = 1:length(catch),
         n_years = length(catch)) %>%
  mutate(scaled_catch = scale(catch)) %>%
  ungroup() %>%
  filter(n_years > 25,
         stock_year <= 25)

ram_catches %>%
  ggplot(aes(stock_year, scaled_catch, color = stockid)) +
  geom_line(show.legend = FALSE)

ram_catches <- ram_catches %>%
  select(stockid, stock_year, scaled_catch) %>%
  pivot_wider(names_from = stock_year, values_from = scaled_catch) %>%
  ungroup()

nstocks <- nrow(ram_catches)

map_dbl(ram_catches, ~ sum(is.na(.x)))

a = ram_catches %>% select(-stockid) %>% as.matrix()
set.seed(42)
catch_pca <- specc(a, centers = 5)

# centers(catch_pca)
# size(catch_pca)
# withinss(catch_pca)

cluster <- as.vector(catch_pca)

ram_catches$cluster <- cluster


ram_catches <- ram_catches  %>%
  pivot_longer(c(-stockid,-cluster),
               names_to = "stock_year",
               values_to = "catch",) %>%
  mutate(stock_year = as.integer(stock_year))

ram_catches %>%
  ggplot(aes(stock_year, catch, group = stockid)) +
  geom_line(alpha = 0.5) +
  facet_wrap( ~ cluster)



cluster_data <- ram_catches %>%
  pivot_wider(names_from = stock_year, values_from = catch) %>%
  ungroup() %>%
  mutate(cluster = as.factor(cluster)) %>%
  janitor::clean_names()

cluster_splits <-
  rsample::initial_split(cluster_data, strata = cluster)
#
# class_model <- parsnip::nearest_neighbor(neighbors = 15,
#                                          mode = "classification") %>%
#   set_engine(engine = 'kknn')

cluster_model <-
  rand_forest(mtry = tune(),
              min_n = tune(),
              trees = 1000) %>%
  set_engine("ranger", num.threads = 8) %>%
  set_mode("classification")


cluster_recipe <-
  recipes::recipe(cluster ~ ., data = training(cluster_splits) %>% select(-stockid)) %>%
  themis::step_upsample(cluster)

cluster_workflow <-
  workflows::workflow() %>%
  workflows::add_model(cluster_model) %>%
  workflows::add_recipe(cluster_recipe)

val_set <- training(cluster_splits) %>% select(-stockid) %>%
  rsample::vfold_cv()

set.seed(345)
cluster_tuning <-
  cluster_workflow %>%
  tune_grid(
    val_set,
    grid = 20,
    control = control_grid(save_pred = TRUE),
    metrics = metric_set(roc_auc)
  )

cluster_tuning %>%
  collect_metrics()

best_forest <- cluster_tuning %>%
  select_best("roc_auc")

final_workflow <-
  cluster_workflow %>%
  finalize_workflow(best_forest)

cluster_fit <-
  final_workflow %>%
  fit(data = training(cluster_splits) %>% select(-stockid))

cluster_fit <- workflows::pull_workflow_fit(cluster_fit)

# cluster_fit <- class_model %>%
#   parsnip::fit(cluster ~ . , data = training(cluster_splits) %>% select(-stockid))

cluster_fit$fit %>% summary()

training_data <- training(cluster_splits) %>%
  bind_cols(predict(cluster_fit, new_data = .)) %>%
  mutate(split = "training")


testing_data <- testing(cluster_splits) %>%
  bind_cols(predict(cluster_fit, new_data = .)) %>%
  mutate(split = "testing")

cluster_predictions <- training_data %>%
  bind_rows(testing_data) %>%
  rename(predicted_cluster = .pred_class)

cluster_predictions %>%
  group_by(cluster) %>%
  count()

cluster_model_performance <- cluster_predictions %>%
  group_by(split, cluster) %>%
  summarise(accuracy = mean(cluster == predicted_cluster))

cluster_model_performance %>%
  ggplot(aes(cluster, accuracy, fill = split)) +
  geom_col(position = "dodge")

cluster_predictions %>%
  group_by(split) %>%
  summarise(accuracy = mean(cluster == predicted_cluster)) %>%
  pivot_wider(names_from = "split", values_from = "accuracy") %>%
  mutate(testing_loss = testing / training - 1)


status_model_data <- ram_data %>%
  filter(stockid %in% unique(cluster_predictions$stockid)) %>%
  left_join(cluster_predictions %>% select(stockid, predicted_cluster),
            by = "stockid") %>%
  # left_join(sraplus::fao_taxa$fao_species %>% select(scientific_name, isscaap_group), by = c("scientificname" = "scientific_name")) %>%
  group_by(stockid) %>%
  mutate(
    c_div_maxc = catch / max(catch, na.rm = TRUE),
    c_div_meanc = catch / mean(catch, na.rm = TRUE),
    c_length = log(length(catch)),
    fishery_year = 1:length(catch)
  ) %>%
  mutate(
    c_roll_meanc = RcppRoll::roll_meanr(c_div_meanc, 5),
    c_roll_maxc = catch / cummax(catch)
  ) %>%
  gather(metric, value, b_v_bmsy, u_v_umsy, exploitation_rate) %>%
  select(stockid,
         year,
         contains('c_'),
         metric,
         value,
         predicted_cluster,
         fishery_year) %>%
  mutate(log_value = log(value + 1e-3)) %>%
  unique() %>%
  na.omit() %>%
  ungroup()


status_model_data %>%
  ggplot(aes(year, c_roll_maxc, group = stockid)) +
  geom_line()

# add in life history -----------------------------------------------------



# fit models --------------------------------------------------------------

b_status_data <- status_model_data %>%
  left_join(ram_data %>% select(stockid, primary_country) %>% unique(), by = "stockid") %>%
  filter(metric == "b_v_bmsy")
# actually add life history here
get_thorson <- function(sciname) {
  hopefully_genus_species <-
    stringr::str_split(sciname, pattern = " ", simplify = TRUE)
  
  loc <- FishLife::Match_species(genus_species = sciname)
  
  lh <-
    as_tibble(t(as.matrix(FishLife::FishBase_and_RAM$beta_gv[loc$GroupNum, ])))
  
  return(lh)
  
}

lh_data <-
  tibble(scientificname = unique(ram_data$scientificname)) %>%
  mutate(lh = map(scientificname, safely(get_thorson)))

lh_data <- lh_data %>%
  mutate(lh_worked = map_lgl(map(lh, "error"), is.null))

lh_data <- lh_data %>%
  filter(lh_worked) %>%
  mutate(lh = map(lh, "result"))

lh_data <- lh_data %>%
  select(-lh_worked) %>%
  unnest(cols = lh)

stock_lookup <- ram_data %>%
  select(stockid, scientificname) %>%
  unique()

b_status_data <- b_status_data %>%
  left_join(stock_lookup, by = "stockid")

b_status_data <- b_status_data %>%
  left_join(lh_data, by = "scientificname")

b_status_data <- b_status_data %>%
  mutate(m_v_k = exp(M) / exp(K),
         lmat_v_linf = exp(Lm) / exp(Loo))




b_status_data %>%
  ggplot(aes(value, fill = predicted_cluster)) +
  geom_histogram(position = "dodge") +
  facet_wrap( ~ predicted_cluster)

b_status_data %>%
  ggplot(aes(c_div_meanc, value, color = predicted_cluster)) +
  geom_point(position = "dodge") +
  geom_smooth(method = "lm") +
  facet_wrap( ~ predicted_cluster)


b_status_data %>%
  ggplot(aes(c_div_maxc, value, color = predicted_cluster)) +
  geom_point(position = "dodge", alpha = 0.25) +
  geom_smooth(method = "lm") +
  facet_wrap( ~ predicted_cluster)

# b_status_data %>%
#   ggplot(aes(fishery_length, value, color = predicted_cluster)) +
#   geom_point(position = "dodge") +
#   geom_smooth(method = "lm") +
#   facet_wrap(~predicted_cluster)


training_stocks <-
  sample(unique(b_status_data$stockid),
         round(n_distinct(b_status_data$stockid) * .8),
         replace = FALSE)

b_status_data$training <-
  !b_status_data$stockid %in% training_stocks

b_status_data <- b_status_data %>%
  arrange(training, stockid, year)

training <- b_status_data %>%
  filter(stockid %in% training_stocks)

testing <- b_status_data %>%
  filter(!stockid %in% training_stocks)


catch_split <-
  initial_time_split(b_status_data, prop = last(which(b_status_data$training == FALSE)) / nrow(b_status_data))

training_data <- training(catch_split)

testing_data <- testing(catch_split)


aa_splits <-
  rsample::group_vfold_cv(training_data, group = primary_country)

tune_grid <-
  parameters(
    min_n(range(2, 10)),
    tree_depth(range(4, 15)),
    learn_rate(range = c(-2, 0)),
    mtry(),
    loss_reduction(),
    sample_prop(range = c(0.25, 1)),
    trees(range = c(500, 2000))
  ) %>%
  dials::finalize(mtry(), x = training_data %>% select(-(1:2)))


xgboost_grid <- grid_latin_hypercube(tune_grid, size = 10)

xgboost_model <-
  parsnip::boost_tree(
    mode = "regression",
    mtry = tune(),
    min_n = tune(),
    loss_reduction = tune(),
    sample_size = tune(),
    learn_rate = tune(),
    tree_depth = tune(),
    trees = tune()
  ) %>%
  parsnip::set_engine("xgboost")


xgboost_workflow <- workflows::workflow() %>%
  add_formula(value ~ c_div_maxc + fishery_year + predicted_cluster + c_div_meanc + c_roll_maxc) %>%
  add_model(xgboost_model)


set.seed(234)
doParallel::registerDoParallel(cores = parallel::detectCores() - 4)
message("this takes a long time to run; go do something else and check on it after a nice lunch")
xgboost_tuning <- tune_grid(
  xgboost_workflow,
  resamples = aa_splits,
  grid = xgboost_grid,
  control = control_grid(save_pred = TRUE)
)

# xgboost_tuning

a = collect_metrics(xgboost_tuning) %>%
  select(mean, mtry:sample_size, .metric) %>%
  pivot_longer(mtry:sample_size, names_to = "dial", values_to = "level") %>%
  ggplot(aes(level, mean)) +
  geom_point() +
  facet_grid(.metric ~ dial, scales = "free")

show_best(xgboost_tuning, "rmse")
#

best_rmse <- tune::select_best(xgboost_tuning, metric = "rmse")

final_workflow <- finalize_workflow(xgboost_workflow,
                                    best_rmse)



catch_model <-
  parsnip::boost_tree(
    mode = "regression",
    mtry = best_rmse$mtry,
    min_n = best_rmse$min_n,
    loss_reduction = best_rmse$loss_reduction,
    sample_size = best_rmse$sample_size,
    learn_rate = best_rmse$learn_rate,
    tree_depth = best_rmse$tree_depth,
    trees = best_rmse$trees
  ) %>%
  parsnip::set_engine("xgboost") %>%
  parsnip::fit(value ~ c_div_maxc + fishery_year + predicted_cluster + c_div_meanc + c_roll_maxc,
               data = training_data)


catch_model %>%
  vip::vi() %>%
  vip::vip(geom = "point")


training_data$boost_fit <-
  predict(catch_model, new_data = training_data)$.pred

testing_data$boost_fit <-
  predict(catch_model, new_data = testing_data)$.pred

catch_model_predictions <- training_data %>%
  mutate(set = "training") %>%
  bind_rows(testing_data %>% mutate(set = "testing"))

catch_model_predictions %>%
  ggplot(aes(value, boost_fit)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap( ~ set)

yardstick::rmse(catch_model_predictions %>% group_by(set),
                truth = value,
                estimate = boost_fit)

yardstick::rsq(catch_model_predictions %>% group_by(set),
               truth = value,
               estimate = boost_fit)


# try again with life history ---------------------------------------------



# catch_split <- initial_time_split(b_status_data, prop = last(which(b_status_data$training == FALSE)) / nrow(b_status_data))
#
# training_data <- training(catch_split)
#
# testing_data <- testing(catch_split)



catch_lh_model <-
  parsnip::boost_tree(
    mode = "regression",
    mtry = best_rmse$mtry,
    min_n = best_rmse$min_n,
    loss_reduction = best_rmse$loss_reduction,
    sample_size = best_rmse$sample_size,
    learn_rate = best_rmse$learn_rate,
    tree_depth = best_rmse$tree_depth,
    trees = best_rmse$trees
  ) %>%
  parsnip::set_engine("xgboost") %>%
  parsnip::fit(
    value ~ c_div_maxc + fishery_year + predicted_cluster + c_div_meanc + c_roll_maxc + m_v_k + lmat_v_linf ,
    data = training_data
  )


catch_lh_model %>%
  vip::vi() %>%
  vip::vip(geom = "point")


training_data$boost_fit <-
  predict(catch_lh_model, new_data = training_data)$.pred

testing_data$boost_fit <-
  predict(catch_lh_model, new_data = testing_data)$.pred

catch_model_predictions <- training_data %>%
  mutate(set = "training") %>%
  bind_rows(testing_data %>% mutate(set = "testing"))

catch_model_predictions %>%
  ggplot(aes(boost_fit, value )) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = .2) +
  geom_smooth(method = "lm") +
  facet_wrap( ~ set) + 
  scale_x_continuous(name = "Predicted B/Bmsy") + 
  scale_y_continuous(name = "RAM B/Bmsy")

yardstick::rmse(catch_model_predictions %>% group_by(set),
                truth = value,
                estimate = boost_fit)

yardstick::rsq(catch_model_predictions %>% group_by(set),
               truth = value,
               estimate = boost_fit)

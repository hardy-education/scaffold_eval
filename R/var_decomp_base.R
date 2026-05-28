library(parallel)
library(data.table)
library(glue)
library(tictoc)
library(here)
library(performance)
library(parameters)
library(modelsummary)
library(modelbased)
library(report)
library(effectsize)
library(insight)
library(energy)
library(see)
library(tidyverse)
library(lme4)
library(lmerTest)
library(fixest)
library(broom.mixed)
library(janitor)
library(naniar)
library(broom)
library(ggthemes)
library(ggrepel)
library(brms)
library(bayestestR)
library(furrr)
library(mirt)
library(future)
library(parallel)
library(fastverse)
fastverse_extend(Rfast,dtplyr, 
                 coop,Rfast2,cheapr,
                 stringi,stringdist,
                 vctrs, dqrng,
                 parallelDist
)


# Load data. ---------
df = fread(here("data/hal_parsed_response_matrix7.csv")) |> 
  mutate(model_name = str_split_i(model_name,":",1),
            model_name = case_when(
              model_name == "DeepSeek-R1" ~ "deepseek-r1",
              model_name == "claude-sonnet-4.5" ~ "claude-sonnet-4-5",
              TRUE ~ model_name
            ),
            model_name = str_c(model_name,reasoning_effort))

## dataset with tokens included for multivariate model (optional)
# dft = fread(here("data/hal_parsed_response_matrix_tok.csv")) |> 
#   mutate(model_name = str_split_i(model_name,":",1),
#          model_name = case_when(
#            model_name == "DeepSeek-R1" ~ "deepseek-r1",
#            model_name == "claude-sonnet-4.5" ~ "claude-sonnet-4-5",
#            TRUE ~ model_name
#          ),
#          model_name = str_c(model_name,reasoning_effort)) 




# Helper functions for extracting Bayesian variance components and calculating Ep2 ---------
get_mod_ep2 = function(drawlist,nitems=100,nbench = 1){
  m_ep2 = drawlist[["sd_model_name__Intercept"]]^2/(drawlist[["sd_model_name__Intercept"]]^2+ drawlist[["sd_benchmark:model_name__Intercept"]]^2+drawlist[["sd_benchmark:task_id:model_name__Intercept"]]^2/(nitems*nbench)+drawlist[["sd_model_name:agent_name__Intercept"]]^2+drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2/(nitems*nbench) + drawlist[["sd_benchmark:model_name:agent_name__Intercept" ]]^2/(nbench) + drawlist[["sigma"]]^2/(nitems*nbench))  
  make_draws_stats(m_ep2)
}
Vectorize(get_mod_ep2)


get_bmod_ep2 = function(drawlist,nitems=100,nbench = 1){
  m_ep2 = drawlist[["sd_benchmark:model_name__Intercept"]]^2/(drawlist[["sd_benchmark:model_name__Intercept"]]^2+drawlist[["sd_benchmark:task_id:model_name__Intercept"]]^2/(nitems*nbench)+drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2/(nitems*nbench) + drawlist[["sd_benchmark:model_name:agent_name__Intercept" ]]^2/(nbench) + drawlist[["sigma"]]^2/(nitems*nbench))  
  make_draws_stats(m_ep2)
}
Vectorize(get_bmod_ep2)


get_agnt_ep2 = function(drawlist,nitems=100,nbench = 1){
  a_ep2 = drawlist[["sd_agent_name__Intercept"]]^2/(drawlist[["sd_agent_name__Intercept"]]^2+ drawlist[["sd_benchmark:agent_name__Intercept"]]^2+drawlist[["sd_benchmark:task_id:agent_name__Intercept"]]^2/(nitems*nbench)+drawlist[["sd_model_name:agent_name__Intercept"]]^2+drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2/(nitems*nbench) + drawlist[["sd_benchmark:model_name:agent_name__Intercept" ]]^2/(nbench) + drawlist[["sigma"]]^2/(nitems*nbench))  
  make_draws_stats(a_ep2)
}
Vectorize(get_agnt_ep2)

get_bagnt_ep2 = function(drawlist,nitems=100,nbench = 1){
  a_ep2 = drawlist[["sd_benchmark:agent_name__Intercept"]]^2/(drawlist[["sd_benchmark:agent_name__Intercept"]]^2+drawlist[["sd_benchmark:task_id:agent_name__Intercept"]]^2/(nitems*nbench)+drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2/(nitems*nbench) + drawlist[["sd_benchmark:model_name:agent_name__Intercept" ]]^2/(nbench) + drawlist[["sigma"]]^2/(nitems*nbench))  
  make_draws_stats(a_ep2)
}
Vectorize(get_bagnt_ep2)


get_agmod_ep2 = function(drawlist,nitems=100,nbench = 1){
  (drawlist[["sd_model_name:agent_name__Intercept"]]^2/(drawlist[["sd_model_name:agent_name__Intercept"]]^2+drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2/(nitems*nbench) + drawlist[["sd_benchmark:model_name:agent_name__Intercept" ]]^2/(nbench) + drawlist[["sigma"]]^2/(nitems*nbench))) |> make_draws_stats()
}
Vectorize(get_agmod_ep2)


get_bagmod_ep2 = function(drawlist,nitems=100,nbench = 1){
  (drawlist[["sd_benchmark:model_name:agent_name__Intercept"]]^2/(drawlist[["sd_benchmark:model_name:agent_name__Intercept" ]]^2/(nbench) +drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2/(nitems*nbench)+ drawlist[["sigma"]]^2/(nitems*nbench))) |> make_draws_stats()
}
Vectorize(get_bagmod_ep2)




# Inverse Logit transformed functions for generalized space "g" to observation space ---------

plogis = function(x) exp(x)/(1+exp(x))


get_mod_ep2g = function(drawlist,nitems=100,nbench = 1){
  m_ep2 = plogis(drawlist[["sd_model_name__Intercept"]]^2)/(plogis(drawlist[["sd_model_name__Intercept"]]^2) + 
                                                              plogis(drawlist[["sd_benchmark:model_name__Intercept"]]^2)+plogis(drawlist[["sd_benchmark:task_id:model_name__Intercept"]]^2)/(nitems*nbench)+plogis(drawlist[["sd_model_name:agent_name__Intercept"]]^2)+plogis(drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2)/(nitems*nbench) + plogis(drawlist[["sd_benchmark:model_name:agent_name__Intercept" ]]^2)/(nbench) + plogis(drawlist[["sigma"]]^2)/(nitems*nbench))  
  make_draws_stats(m_ep2)
}
Vectorize(get_mod_ep2g)


get_bmod_ep2g = function(drawlist,nitems=100,nbench = 1){
  m_ep2 = plogis(drawlist[["sd_benchmark:model_name__Intercept"]]^2)/(plogis(drawlist[["sd_benchmark:model_name__Intercept"]]^2)+plogis(drawlist[["sd_benchmark:task_id:model_name__Intercept"]]^2)/(nitems*nbench)+plogis(drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2)/(nitems*nbench) + plogis(drawlist[["sd_benchmark:model_name:agent_name__Intercept" ]]^2)/(nbench) + plogis(drawlist[["sigma"]]^2)/(nitems*nbench))  
  make_draws_stats(m_ep2)
}
Vectorize(get_bmod_ep2g)


get_agnt_ep2g = function(drawlist,nitems=100,nbench = 1){
  a_ep2 = (plogis(drawlist[["sd_agent_name__Intercept"]]^2))/(plogis(drawlist[["sd_agent_name__Intercept"]]^2)+ plogis(drawlist[["sd_benchmark:agent_name__Intercept"]]^2)+plogis(drawlist[["sd_benchmark:task_id:agent_name__Intercept"]]^2)/(nitems*nbench)+plogis(drawlist[["sd_model_name:agent_name__Intercept"]]^2)+plogis(drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2)/(nitems*nbench) + plogis(drawlist[["sd_benchmark:model_name:agent_name__Intercept" ]]^2)/(nbench) + plogis(drawlist[["sigma"]]^2)/(nitems*nbench))  
  make_draws_stats(a_ep2)
}
Vectorize(get_agnt_ep2g)

get_bagnt_ep2g = function(drawlist,nitems=100,nbench = 1){
  a_ep2 = (plogis(drawlist[["sd_benchmark:agent_name__Intercept"]]^2))/(plogis(drawlist[["sd_benchmark:agent_name__Intercept"]]^2)+plogis(drawlist[["sd_benchmark:task_id:agent_name__Intercept"]]^2)/(nitems*nbench)+plogis(drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2)/(nitems*nbench) + plogis(drawlist[["sd_benchmark:model_name:agent_name__Intercept" ]]^2)/(nbench) + plogis(drawlist[["sigma"]]^2)/(nitems*nbench))  
  make_draws_stats(a_ep2)
}
Vectorize(get_bagnt_ep2g)


get_agmod_ep2g = function(drawlist,nitems=100,nbench = 1){
  (((plogis(drawlist[["sd_model_name:agent_name__Intercept"]]^2)))/(plogis(drawlist[["sd_model_name:agent_name__Intercept"]]^2)+
                                                                      plogis(drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2)/(nitems*nbench) + 
                                                                      plogis(drawlist[["sd_benchmark:model_name:agent_name__Intercept" ]]^2)/(nbench) + plogis(drawlist[["sigma"]]^2)/(nitems*nbench))) |> make_draws_stats()
}
Vectorize(get_agmod_ep2g)


get_bagmod_ep2g = function(drawlist,nitems=100,nbench = 1){
  (((plogis(drawlist[["sd_benchmark:model_name:agent_name__Intercept"]]^2)))/(plogis(drawlist[["sd_benchmark:model_name:agent_name__Intercept" ]]^2/(nbench) +
                                                                                         plogis(drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2)/(nitems*nbench)+ 
                                                                                         plogis(drawlist[["sigma"]]^2)/(nitems*nbench)))) |> make_draws_stats()
}
Vectorize(get_bagmod_ep2g)



## Bayes statistics from posterior draws ---------
### BAYES create combination stats -------
make_draws_stats = function(drawstat,roperange = c(0,0.005)){
  res = summary(drawstat) |> 
    bind_cols(bayestestR::p_significance(drawstat,threshold = 0.05) |> 
                select(-Parameter),
              bayestestR::p_map(drawstat
              ) |> 
                rename(pval=p_MAP) |> 
                select(-Parameter)
              ,
              bayestestR::p_map(drawstat,
                                method="KernSmooth"
              ) |> 
                select(-Parameter),
              bayestestR::map_estimate(drawstat) |> 
                select(-Parameter),
              bayestestR::equivalence_test(
                drawstat,   # effects = "random_variance",
                range=roperange
              ) |> 
                select(-Parameter)
              
    )
  res
}






# Models =========
## LME model-------

tic()
m = df |> lmer(score ~ 1
               + (1|benchmark)
               + (1|benchmark:task_id)
               + (1|model_name)
               + (1|agent_name)
               + (1|benchmark:task_id:model_name)
               + (1|benchmark:task_id:agent_name)
               + (1|model_name:agent_name)
               + (1|benchmark:task_id:model_name:agent_name)
               + (1|benchmark:model_name)
               + (1|benchmark:agent_name)
               + (1|benchmark:model_name:agent_name)
               ,
               na.action = na.exclude,
               data = _,
               control = lmerControl(optimizer = "bobyqa", 
                                     # calc.derivs = FALSE,
                                     optCtrl = list(maxfun = 2e5))
)

m |> summary()

mvar = m  |>
  VarCorr() |>
  janitor::clean_names() |>
  as_tibble() |>
  filter(is.na(var2)) |>
  mutate(
    pct = round(vcov / sum(vcov) * 100, 2),
    agent = str_detect(grp, "agent"),
    bench = str_detect(grp, "benchmark"),
    model = str_detect(grp, "model"),
    task = str_detect(grp, "task"),
    run = str_detect(grp, "run"),
    resid = grp == "Residual",
    model_agent = model & agent,
    across(agent:model_agent, \(x) vcov * as.numeric(x))
  )
mvar
toc()



## GLME model-------

tic()
gm = df |> glmer(score ~ 1
               + (1|benchmark)
               + (1|benchmark:task_id)
               + (1|model_name)
               + (1|agent_name)
               + (1|benchmark:task_id:model_name)
               + (1|benchmark:task_id:agent_name)
               + (1|model_name:agent_name)
               + (1|benchmark:task_id:model_name:agent_name)
               + (1|benchmark:model_name)
               + (1|benchmark:agent_name)
               + (1|benchmark:model_name:agent_name)
               ,
               na.action = na.exclude,
               data = _,
               family = "binomial",
               # nAGQ=0, # keep nAGQ at default for more accurate variance component estimation via Laplace approx, but can set to 0 for faster fitting if needed
               control = glmerControl(optimizer = "bobyqa", 
                                     # calc.derivs = FALSE,
                                     optCtrl = list(maxfun = 2e5))
)

gm |> summary()

gmvar = gm  |>
  VarCorr() |>
  janitor::clean_names() |>
  as_tibble() |>
  filter(is.na(var2)) |>
  bind_rows(as_tibble_row(list(grp="Residual",var1=NA,vcov = pi^2/3))) |> 
  mutate(
    pct = round(vcov / sum(vcov) * 100, 2),
    agent = str_detect(grp, "agent"),
    bench = str_detect(grp, "benchmark"),
    model = str_detect(grp, "model"),
    task = str_detect(grp, "task"),
    run = str_detect(grp, "run"),
    resid = grp == "Residual",
    model_agent = model & agent,
    across(agent:model_agent, \(x) vcov * as.numeric(x))
  )
gmvar
toc()



tic()
gm2 = df |> glmer(score ~ 1
                 + (1|benchmark)
                 + (1|benchmark:task_id)
                 + (1|model_name)
                 + (1|agent_name)
                 + (1|benchmark:task_id:model_name)
                 + (1|benchmark:task_id:agent_name)
                 + (1|model_name:agent_name)
                 + (1|benchmark:task_id:model_name:agent_name)
                 + (1|benchmark:model_name)
                 + (1|benchmark:agent_name)
                 + (1|benchmark:model_name:agent_name)
                 ,
                 na.action = na.exclude,
                 data = _,
                 family = "binomial",
                 # nAGQ=0,
                 control = glmerControl(optimizer = "bobyqa", 
                                        calc.derivs = FALSE,
                                        optCtrl = list(maxfun = 2e5))
)

gm2 |> summary()

gmvar2 = gm2  |>
  VarCorr() |>
  janitor::clean_names() |>
  as_tibble() |>
  filter(is.na(var2)) |>
  bind_rows(as_tibble_row(list(grp="Residual",var1=NA,vcov = pi^2/3))) |> 
  mutate(
    pct = round(vcov / sum(vcov) * 100, 2),
    agent = str_detect(grp, "agent"),
    bench = str_detect(grp, "benchmark"),
    model = str_detect(grp, "model"),
    task = str_detect(grp, "task"),
    run = str_detect(grp, "run"),
    resid = grp == "Residual",
    model_agent = model & agent,
    across(agent:model_agent, \(x) vcov * as.numeric(x))
  )
gmvar2
toc()


# Bayes --------
## linear model ------------

tic()
bm = df |> brm(score ~ 1
               + (1|benchmark)
               + (1|benchmark:task_id)
               + (1|model_name)
               + (1|agent_name)
               + (1|benchmark:task_id:model_name)
               + (1|benchmark:task_id:agent_name)
               + (1|model_name:agent_name)
               + (1|benchmark:task_id:model_name:agent_name)
               + (1|benchmark:model_name)
               + (1|benchmark:agent_name)
               + (1|benchmark:model_name:agent_name)
               ,
               chains = 4,
               cores = 4,
               thin = 5,
               iter = 2000,
               backend = "cmdstanr",
               control = list(adapt_delta = 0.95),
               data = _,
               threads = threading(2),
               save_pars= save_pars(all = T),
               stan_model_args = list(stanc_options = list("O1")),
               file = here::here("data","brmfit_.rds"),
               file_refit = "always"
               
               
)

bm = readRDS(here::here("data","brmfit_.rds"))

bm |> summary()

toc()

# get draws from posterior ---------
bmd = as_draws_rvars(bm)

# make pct summary from draws
bms = lapply(bmd[str_starts(names(bmd),"sd|sig")], function(X) make_draws_stats(X^2/reduce(bmd[str_starts(names(bmd),"sd|sig")],\(acc,nxt) acc + (nxt^2)))) |> list_rbind(names_to = "var")
bms




## Generalized Bayes estimation  -----

tic()
gbm = df |> brm(score ~ 1
               + (1|benchmark)
               + (1|benchmark:task_id)
               + (1|model_name)
               + (1|agent_name)
               + (1|benchmark:task_id:model_name)
               + (1|benchmark:task_id:agent_name)
               + (1|model_name:agent_name)
               + (1|benchmark:task_id:model_name:agent_name)
               + (1|benchmark:model_name)
               + (1|benchmark:agent_name)
               + (1|benchmark:model_name:agent_name)
               ,
               chains = 4,
               cores = 4,
               family = bernoulli(link = "logit"),
               thin = 5,
               iter = 2000,
               backend = "cmdstanr",
               control = list(adapt_delta = 0.95),
               data = _,
               threads = threading(2),
               save_pars= save_pars(all=T),
               stan_model_args = list(stanc_options = list("O1")),
               file = here::here("data","gbrmfit_.rds"),
               file_refit = "always"
               
               
)

gbm = readRDS(here::here("data","gbrmfit_.rds"))

gbm |> summary()

toc()

# get draws
gbmd = as_draws_rvars(gbm)
gbmd$sigma = posterior::rvar(rnorm(800, mean = pi^2/3, sd = 0))

# make pct summary
gbms = lapply(gbmd[str_starts(names(gbmd),"sd|sig")], function(X) make_draws_stats(X^2/reduce(gbmd[str_starts(names(gbmd),"sd|sig")],\(acc,nxt) acc + (nxt^2)))) |> list_rbind(names_to = "var")
gbms




### Bayes MV (simple scaling) --------

tic()
bm2 = dft |> 
  filter(total_tokens!=0) |> 
  brm(bf(mvbind(score,log(total_tokens)) ~ 1
               + (1|benchmark)
               + (1|benchmark:task_id)
               + (1|model_name)
               + (1|agent_name)
               + (1|benchmark:task_id:model_name)
               + (1|benchmark:task_id:agent_name)
               + (1|model_name:agent_name)
               + (1|benchmark:task_id:model_name:agent_name)
               + (1|benchmark:model_name)
               + (1|benchmark:agent_name)
               + (1|benchmark:model_name:agent_name))
               ,
               chains = 4,
               cores = 4,
               thin = 5,
               iter = 2000,
               backend = "cmdstanr",
               control = list(adapt_delta = 0.95),
               data = _,
               # family = list(bernoulli(link = "logit"), 
               #               gaussian(link = "identity")),
               threads = threading(2),
               save_pars = save_pars("all"),
               stan_model_args = list(stanc_options = list("O1")),
               file = here::here("data","brmfit_tok.rds"),
               file_refit = "always"


)

bm2 |> summary()

toc()

# get draws
bmd2 = as_draws_rvars(bm2)

# make pct summary
bms2 = lapply(bmd2[str_starts(names(bmd2),"sd|sig")&str_detect(names(bmd2),"score")], function(X) make_draws_stats(X^2/reduce(bmd2[str_starts(names(bmd2),"sd|sig")&str_detect(names(bmd2),"score")],\(acc,nxt) acc + (nxt^2)))) |> list_rbind(names_to = "var")
bms2

bms2tok = lapply(bmd2[str_starts(names(bmd2),"sd|sig")&!str_detect(names(bmd2),"score")], function(X) make_draws_stats(X^2/reduce(bmd2[str_starts(names(bmd2),"sd|sig")&!str_detect(names(bmd2),"score")],\(acc,nxt) acc + (nxt^2)))) |> list_rbind(names_to = "var")
bms2tok







# Non parametric deviance/distribution decomposition --------
  
n_iter = 3
n_i = 100
prop_m =1
alpha = 1

## full set with 7 benchmarks iteration takes 1hr to run
L = list()

bench_list = df |> pull(benchmark) |> unique()
df

b = bench_list[1]
tic("all iterations")
for (i in 1:1){
  
  model_list = df |> dtplyr::lazy_dt() |> 
    select(model_name) |> 
    distinct() |>
    # slice_sample(prop = 0.1) |>  
    pull(model_name)
  
  # for (b in bench_list){
  tic(glue::glue("Iteration {i} - Benchmark: {b}"))
  cat(glue::glue("Iteration {i} - Benchmark: {b}\n"))
  
  set.seed(i*1234)
  
  item_list = df |> dtplyr::lazy_dt() |> 
    # filter(benchmark==b) |> 
    select(task_id) |> 
    distinct() |> 
    # slice_sample(n = 42) |>  
    pull(task_id)
  
  
  
  
  dat = df |> dtplyr::lazy_dt() |> 
    filter(#benchmark==b,
      task_id %in% item_list,
      model_name %in% model_list) |> 
    drop_na(score) |>
    mutate(resp = as.numeric(score),
           item = as.factor(interaction(benchmark,task_id)),
           agent = as.factor(agent_name),
           model = as.factor(model_name),
           benchmark = as.factor(benchmark),
           item_agent = fct_drop(interaction(benchmark,item,agent)),
           item_model = fct_drop(interaction(benchmark,item,model)),
           model_agent = fct_drop(interaction(model,agent)),
           bench_agent = fct_drop(interaction(benchmark,agent)),
           bench_model = fct_drop(interaction(benchmark,model)),
           bench_model_agent = fct_drop(interaction(benchmark,model,agent)),    
           item_model_agent = fct_drop(interaction(benchmark,item,model,agent))) |> select(-any_of(c("run_id","model_name","agent_name","score","task_id","reasoning_effort")))
  
  disco_fit = disco(dat |> select(resp) |> as.data.table() ,
                    factors = dat |> select(-resp) |> 
                      as.data.table(), 
                    distance = F,
                    index = alpha,
                    R=0)
  
  disco_fit_df = disco_fit$stats |> 
    as_tibble() |> 
    select(-any_of("p-value")) |> 
    mutate(facet = disco_fit$factor.names,Withins = Within, 
           Within = disco_fit$within, Total = disco_fit$total, 
           Ep2 = Trt / (Trt + Withins), SNR = Ep2/(1-Ep2),
           Ep2hat = (Trt/df1) / ((Trt/df1) + Withins/df2), 
           SNRhat = Ep2hat/(1-Ep2hat),
           pct = Trt/(Total-Within),
           mean = Trt/df1, var = Withins/df2,
           n_items = n_i,prop_models = prop_m) |> 
    select(facet,Ep2,SNR,everything())
  
  disco_fit_df = disco_fit_df |> mutate(#benchmark = b,
    n_items = n_i,prop_models = prop_m,
    iteration = i,alpha = alpha)
  
  L[[length(L)+1]] = disco_fit_df
  print(disco_fit_df)
  toc()
  # }
}

disco_results = bind_rows(L)

disco_results |>
  write_csv(file = here::here("data", "results","disco_all_model.csv"))

disco_results = read_csv(here::here("data","results","disco_all_model.csv"))

disco_results
toc()



# Combined table --------
combdf = disco_results |> 
  mutate(facet = str_replace_all(str_replace_all(facet,"item","bench_item"),"bench_","benchmark_"),method = "nonp") |> 
  select(method, facet,pct) |> 
  bind_rows(mvar |> mutate(facet = str_replace_all(str_replace_all(str_remove_all(str_replace_all(grp,"task_id","item"),regex("sd_|__Intercept|_name|_id")),":","_"),"Residual","sigma"),method = "lme",pct=pct/100) |> select(facet, method,pct)) |>
  bind_rows(gmvar |> mutate(facet = str_replace_all(str_replace_all(str_remove_all(str_replace_all(grp,"task_id","item"),regex("sd_|__Intercept|_name|_id")),":","_"),"Residual","sigma"),method = "glme",pct=pct/100) |> select(facet, method,pct)) |> 
  bind_rows(gmvar2 |> mutate(facet = str_replace_all(str_replace_all(str_remove_all(str_replace_all(grp,"task_id","item"),regex("sd_|__Intercept|_name|_id")),":","_"),"Residual","sigma"),method = "glme2",pct=pct/100) |> select(facet, method,pct)) |> 
  bind_rows(bms |> mutate(facet = str_replace_all(str_remove_all(str_replace_all(var,"task_id","item"),regex("sd_|__Intercept|_name|_id")),":","_"),pct = mean,method="bayes") |> select(facet,method,pct)) |>
  bind_rows(gbms |> mutate(facet = str_replace_all(str_remove_all(str_replace_all(var,"task_id","item"),regex("sd_|__Intercept|_name|_id")),":","_"),pct = mean,method="gbayes") |> select(facet,method,pct)) |> mutate(facet=if_else(facet == "benchmark_item_model_agent","sigma",facet))


combdf |> ggplot(aes(x= method, y = pct, fill = facet)) + geom_col() + theme_minimal()

combdf |> group_by(method,facet) |> 
  summarize(pct = sum(pct)) |> 
  ungroup()|> 
  pivot_wider(names_from = method,values_from = pct) |> 
  select(facet,lme,glme2,gbayes,nonp) |> 
  kableExtra::kable(digits = 3,format = "latex", booktabs = T) |> 
  kableExtra::kable_styling(latex_options = c("hold_position","scale_down"))



# combined plots --------
plot_ni = 500
tibble(nitems=1:plot_ni, m_ep2 = get_mod_ep2(bmd,nitems = nitems) |> pull(median), a_ep2 = get_agnt_ep2(bmd,nitems = nitems) |> pull(median)) |> pivot_longer(m_ep2:a_ep2,names_to = "facet",values_to = "bayes")|> inner_join(tibble(nitems=1:plot_ni, m_ep2 = get_mod_ep2(gbmd,nitems = nitems) |> pull(median), a_ep2 = get_agnt_ep2(gbmd,nitems = nitems) |> pull(median)) |> pivot_longer(m_ep2:a_ep2,names_to = "facet",values_to = "gbayes") |> distinct())  |> inner_join(tibble(nitems=1:plot_ni, m_ep2 = mvar[mvar$grp=="model_name","vcov"][[1]]/(mvar[mvar$grp=="model_name","vcov"][[1]]+mvar[mvar$grp=="benchmark_model_name","vcov"][[1]]+mvar[mvar$grp=="model_name_agent_name","vcov"][[1]]+mvar[mvar$grp=="benchmark_model_name_agent_name","vcov"][[1]]+mvar[mvar$grp=="benchmark_task_id_model_name","vcov"][[1]]/(nitems)+mvar[mvar$grp=="benchmark_task_id_model_name_agent_name","vcov"][[1]]/(nitems)+mvar[mvar$grp=="Residual","vcov"][[1]]/(nitems)), a_ep2 =mvar[mvar$grp=="agent_name","vcov"][[1]]/(mvar[mvar$grp=="agent_name","vcov"][[1]]+mvar[mvar$grp=="benchmark_agent_name","vcov"][[1]]+mvar[mvar$grp=="model_name_agent_name","vcov"][[1]]+mvar[mvar$grp=="benchmark_model_name_agent_name","vcov"][[1]]+mvar[mvar$grp=="benchmark_task_id_agent_name","vcov"][[1]]/(nitems)+mvar[mvar$grp=="benchmark_task_id_model_name_agent_name","vcov"][[1]]/(nitems)+mvar[mvar$grp=="Residual","vcov"][[1]]/(nitems)) ) |> pivot_longer(m_ep2:a_ep2,names_to = "facet",values_to = "lmer") |> distinct()) |> inner_join(tibble(nitems=1:plot_ni, m_ep2 = gmvar[gmvar$grp=="model_name","vcov"][[1]]/(gmvar[gmvar$grp=="model_name","vcov"][[1]]+gmvar[gmvar$grp=="benchmark_model_name","vcov"][[1]]+gmvar[gmvar$grp=="model_name_agent_name","vcov"][[1]]+gmvar[gmvar$grp=="benchmark_model_name_agent_name","vcov"][[1]]+gmvar[gmvar$grp=="benchmark_task_id_model_name","vcov"][[1]]/(nitems)+gmvar[gmvar$grp=="benchmark_task_id_model_name_agent_name","vcov"][[1]]/(nitems)+gmvar[gmvar$grp=="Residual","vcov"][[1]]/(nitems)), a_ep2 =gmvar[gmvar$grp=="agent_name","vcov"][[1]]/(gmvar[gmvar$grp=="agent_name","vcov"][[1]]+gmvar[gmvar$grp=="benchmark_agent_name","vcov"][[1]]+gmvar[gmvar$grp=="model_name_agent_name","vcov"][[1]]+gmvar[gmvar$grp=="benchmark_model_name_agent_name","vcov"][[1]]+gmvar[gmvar$grp=="benchmark_task_id_agent_name","vcov"][[1]]/(nitems)+gmvar[gmvar$grp=="benchmark_task_id_model_name_agent_name","vcov"][[1]]/(nitems)+gmvar[gmvar$grp=="Residual","vcov"][[1]]/(nitems)) ) |> pivot_longer(m_ep2:a_ep2,names_to = "facet",values_to = "glmer") |> distinct()) |> inner_join(tibble(nitems=1:plot_ni, m_ep2 = 310/(310+3610/nitems+684+779+889+4392/nitems),a_ep2=311/(311+2022/nitems+684+584+889+4392/nitems)) |> pivot_longer(m_ep2:a_ep2,names_to = "facet",values_to = "nonp") |> distinct()) |> pivot_longer(bayes:nonp, names_to = "method",values_to = "rel") |> ggplot(aes(x=nitems,y=rel,color=facet,linetype=method)) + geom_line() + theme_minimal() + scale_x_log10()



# Bayes plots for paper --------

## Bayes CIs --------
### plot across 1 benchmark -------



bayescidf = tibble(nitems=1:plot_ni, m_ep2 = get_mod_ep2(bmd,nitems = nitems) |> pull(median),   
       m_ep2_low = get_mod_ep2(bmd,nitems = nitems) |> pull(HDI_low),
       m_ep2_hi = get_mod_ep2(bmd,nitems = nitems) |> pull(HDI_high),
       a_ep2 = get_agnt_ep2(bmd,nitems = nitems) |> pull(median),
       a_ep2_low = get_agnt_ep2(bmd,nitems = nitems) |> pull(HDI_low),
       a_ep2_hi = get_agnt_ep2(bmd,nitems = nitems) |> pull(HDI_high),
       am_ep2 = get_agmod_ep2(bmd,nitems = nitems) |> pull(median),
       am_ep2_low = get_agmod_ep2(bmd,nitems = nitems) |> pull(HDI_low),
       am_ep2_hi = get_agmod_ep2(bmd,nitems = nitems) |> pull(HDI_high),
       )


plot_ni = 500
bayescidf = tibble(nitems=c(1:plot_ni), m_ep2 = get_mod_ep2(bmd,nitems = nitems) |> pull(mean),   
                   m_ep2_low = get_mod_ep2(bmd,nitems = nitems) |> pull(q5),
                   m_ep2_hi = get_mod_ep2(bmd,nitems = nitems) |> pull(q95),
                   a_ep2 = get_agnt_ep2(bmd,nitems = nitems) |> pull(mean),
                   a_ep2_low = get_agnt_ep2(bmd,nitems = nitems) |> pull(q5),
                   a_ep2_hi = get_agnt_ep2(bmd,nitems = nitems) |> pull(q95),
                   am_ep2 = get_agmod_ep2(bmd,nitems = nitems) |> pull(mean),
                   am_ep2_low = get_agmod_ep2(bmd,nitems = nitems) |> pull(q5),
                   am_ep2_hi = get_agmod_ep2(bmd,nitems = nitems) |> pull(q95),
)


bayescidf |> pivot_longer(m_ep2:am_ep2_hi,
                          names_to = c("facet","stat","val"),
                          names_sep = "_",
                          values_to = "bayes") |> 
  mutate(val = if_else(is.na(val),"est",val)) |>
  pivot_wider(names_from = val, values_from = bayes) |> 
  filter(facet=="am") |> 
  filter(nitems<=500) |> 
  ggplot(aes(x=nitems,y=est,color = facet,fill=facet)) +
  geom_ribbon(aes(x = nitems,y=est,ymin = low,ymax = hi),alpha=0.12,linewidth = 0) +
  geom_line(linewidth=1.5) +
  scale_x_log10()+
  theme_minimal()+
  ylim(c(0, 1)) +
  labs(x = "Number of Agentic Tasks", y = expression("Estimated Reliability: E"*hat(rho)^2)) +
  theme(text=element_text(family="times",size = 15), legend.position = "none")






## Generalized Bayes CIs --------
### plot across 1 benchmark -------


plot_ni = 350
n_bench = 1
gbayescidf1 = tibble(nitems=1:plot_ni, m_ep2 = get_bmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(mean),   
                    m_ep2_low = get_bmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q5),
                    m_ep2_hi = get_bmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q95),
                    a_ep2 = get_bagnt_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(mean),
                    a_ep2_low = get_bagnt_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q5),
                    a_ep2_hi = get_bagnt_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q95),
                    am_ep2 = get_bagmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(mean),
                    am_ep2_low = get_bagmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q5),
                    am_ep2_hi = get_bagmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q95),
)





gbayescidf1 |> 
  pivot_longer(m_ep2:am_ep2_hi,
                           names_to = c("facet","stat","val"),
                           names_sep = "_",
                           values_to = "bayes") |> 
  mutate(val = if_else(is.na(val),"est",val)) |> 
  pivot_wider(names_from = val, values_from = bayes) |> 
  filter(facet%in%c("m","a","am")) |> 
  ggplot(aes(x=nitems,y=est,color = facet,fill=facet)) +
  geom_ribbon(aes(x = nitems,y=est,ymin = low,ymax = hi),alpha=0.2,linewidth = 0) +
  geom_line(linewidth=1.5) +
  scale_x_log10()+
  theme_minimal()+
  ylim(c(0, 1)) +
  labs(x = "Number of Agentic Tasks", y = expression("Estimated Reliability: E"*hat(rho)^2)) +
  theme(text=element_text(family="times",size = 15), legend.position = "bottom")




### plot across 7 benchmarks -------

plot_ni = 50
n_bench = 7
gbayescidf = tibble(nitems=1:plot_ni, m_ep2 = get_bmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(mean),   
                   m_ep2_low = get_bmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q5),
                   m_ep2_hi = get_bmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q95),
                   a_ep2 = get_bagnt_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(mean),
                   a_ep2_low = get_bagnt_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q5),
                   a_ep2_hi = get_bagnt_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q95),
                   am_ep2 = get_bagmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(mean),
                   am_ep2_low = get_bagmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q5),
                   am_ep2_hi = get_bagmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q95),
)


plot_ni = 14
n_bench = 25
gbayescidf25 = tibble(nitems=1:plot_ni, m_ep2 = get_bmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(mean),   
                    m_ep2_low = get_bmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q5),
                    m_ep2_hi = get_bmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q95),
                    a_ep2 = get_bagnt_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(mean),
                    a_ep2_low = get_bagnt_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q5),
                    a_ep2_hi = get_bagnt_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q95),
                    am_ep2 = get_bagmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(mean),
                    am_ep2_low = get_bagmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q5),
                    am_ep2_hi = get_bagmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q95),
)

plot_ni = 175
n_bench = 2
gbayescidf2 = tibble(nitems=1:plot_ni, m_ep2 = get_bmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(mean),   
                      m_ep2_low = get_bmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q5),
                      m_ep2_hi = get_bmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q95),
                      a_ep2 = get_bagnt_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(mean),
                      a_ep2_low = get_bagnt_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q5),
                      a_ep2_hi = get_bagnt_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q95),
                      am_ep2 = get_bagmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(mean),
                      am_ep2_low = get_bagmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q5),
                      am_ep2_hi = get_bagmod_ep2g(gbmd,nitems = nitems,nbench=n_bench) |> pull(q95),
)






gbayescidf |> pivot_longer(m_ep2:am_ep2_hi,
                          names_to = c("facet","stat","val"),
                          names_sep = "_",
                          values_to = "bayes") |> 
  mutate(val = if_else(is.na(val),"est",val),
         nitems = nitems*7) |> 
  pivot_wider(names_from = val, values_from = bayes) |> 
  filter(facet%in%c("m","a","am")) |> 
  ggplot(aes(x=nitems,y=est,color = facet,fill=facet)) +
  geom_ribbon(aes(x = nitems,y=est,ymin = low,ymax = hi),alpha=0.12,linewidth = 0) +
  geom_line(linewidth=1.5) +
  scale_x_log10()+
  theme_minimal()+
  ylim(c(0, 1)) +
  labs(x = "Total Number of Agentic Tasks (divided across 7 benchmarks)", y = expression("Estimated Reliability: E"*hat(rho)^2)) +
  theme(text=element_text(family="times",size = 15), legend.position = "bottom")



gbayescidf |> pivot_longer(m_ep2:am_ep2_hi,
                           names_to = c("facet","stat","val"),
                           names_sep = "_",
                           values_to = "bayes") |> 
  mutate(val = if_else(is.na(val),"est",val),
         nitems = nitems*7) |> 
  pivot_wider(names_from = val, values_from = bayes) |> 
  filter(facet=="m") |> mutate(Benches_Sampled = "7 benches") |> 
  bind_rows(
    gbayescidf1 |> pivot_longer(m_ep2:am_ep2_hi,
                                names_to = c("facet","stat","val"),
                                names_sep = "_",
                                values_to = "bayes") |> 
      mutate(val = if_else(is.na(val),"est",val)) |> 
      pivot_wider(names_from = val, values_from = bayes) |> 
      filter(facet=="m",nitems>6) |> mutate(Benches_Sampled = "1 bench")
  ) |>
  # bind_rows(
  #   gbayescidf25 |> pivot_longer(m_ep2:am_ep2_hi,
  #                               names_to = c("facet","stat","val"),
  #                               names_sep = "_",
  #                               values_to = "bayes") |> 
  #     mutate(val = if_else(is.na(val),"est",val),nitems = nitems*25) |> 
  #     pivot_wider(names_from = val, values_from = bayes) |> 
  #     filter(facet=="m",nitems>6) |> mutate(Benches_Sampled = "25 benches")
  # ) |>
  bind_rows(
    gbayescidf2 |> pivot_longer(m_ep2:am_ep2_hi,
                                 names_to = c("facet","stat","val"),
                                 names_sep = "_",
                                 values_to = "bayes") |> 
      mutate(val = if_else(is.na(val),"est",val),nitems = nitems*2) |> 
      pivot_wider(names_from = val, values_from = bayes) |> 
      filter(facet=="m",nitems>6) |> mutate(Benches_Sampled = "2 benches")
  ) |>
  ggplot(aes(x=nitems,y=est,color = Benches_Sampled,fill=Benches_Sampled)) +
  geom_ribbon(aes(x = nitems,y=est,ymin = low,ymax = hi),alpha=0.2,linewidth = 0) +
  geom_line(linewidth=1.5) +
  scale_x_log10()+
  theme_minimal()+
  ylim(c(0, 1)) +
  labs(x = "Total Number of Agentic Tasks", y = expression("Estimated Reliability: E"*hat(rho)^2)) +
  theme(text=element_text(family="times",size = 15), legend.position = "bottom")


## Bar plot with all the proportions --------


component_cols = c(
  "agent" = "#009E73",
  "item_agent"="#F0E442",
  "model_agent" = "#0072B2",
  "item_model" = "#D55E00",
  "item" = "#E69F00",
  "model" = "#56B4E9",
  "model_name" = "#56B4E9",
  "sigma" = "#000000",
  "benchmark" = "#CC79A7",
  "benchmark_agent" = "#009E73",
  "benchmark_item_agent"="#F0E442",
  "benchmark_model_agent" = "#0072B2",
  "benchmark_item_model" = "#D55E00",
  "benchmark_item" = "#E69F00",
  "benchmark_model" = "#56B4E9",
  "benchmark_item_model_agent" = "#000000"
)

component_labs = c(
  "agent" = "Agent",
  "benchmark_agent" =  "Benchmark-Agent",
  "benchmark" = "Benchmark",
  "benchmark_model_agent" = "Benchmark-Model-Agent",
  "item_agent"="Item-Agent",
  "benchmark_item_agent"= "Item-Agent",
  "item" = "Item",
  "benchmark_item" = "Item",
  "model_agent" = "Model-Agent",
  "item_model" = "Model-Item",
  "benchmark_item_model" = "Model-Item",
  "model" = "Model",
  "benchmark_model" = "Benchmark-Model",
  "benchmark_item_model_agent" = "Item-Model-Agent",
  "sigma" = "Residual"
)

component_alphas = c(
  "agent" = 1,
  "item_agent"=1,
  "model_agent" = 1,
  "item_model" = 1,
  "item" = 1,
  "model" = 1,
  "sigma" = 1,
  "benchmark" = 1,
  "benchmark_agent" = 0.5,
  "benchmark_item_agent"=1,
  "benchmark_model_agent" = 0.5,
  "benchmark_item_model" = 1,
  "benchmark_item" = 1,
  "benchmark_model" = 0.5,
  "benchmark_item_model_agent" =0.5
)



(bp = gbms |>
  mutate(across(mean:MAP_Estimate,\(x) x / sum(x)),
         facet = str_replace_all(str_remove_all(str_replace_all(var,
                                                                "task_id","item"),
                                                regex("sd_|__Intercept|_name|_id")),
                                 ":","_"),
         method = "Bern. (Bayes)",
         facet = if_else(facet == "benchmark_item_model_agent","sigma",facet),
         facet = factor(facet, levels = names(component_labs),ordered = T)
         ) |>
  inner_join(as_tibble(list(facet=names(component_labs),labels = component_labs))) |> 
  ggplot(aes(x = method, y=median,fill=facet,alpha=facet)) + 
  geom_col() + 
    
    geom_label_repel(
      aes(label = labels),
      position = position_stack(vjust = 0.5),
      direction = "y",
      # hjust = 2, 
      # nudge_x = 0.1,nudge_y = 0.1,
      size =3.5,
      # segment.alpha = 0.3
    ) +
  # geom_label(aes(label = labels), position = position_stack(vjust = 0.5),size = 3.25) + 
  theme_minimal() + coord_flip()+
   labs(y = "Proportion of Variance Explained", fill = "Component", alpha = "Component") + 
  # ggthemes::scale_fill_colorblind() +
  scale_fill_manual(values = component_cols,labels = component_labs) + 
  scale_alpha_manual(values = component_alphas,labels = component_labs) +
  scale_y_continuous(labels = scales::percent) + 
  theme(text=element_text(family="times",size = 14), legend.position = "none", 
        # axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank(),
        axis.minor.ticks.x.top= element_blank(), axis.minor.ticks.x.bottom= element_blank(), 
        axis.line.x = element_blank(), axis.title.y = element_blank(),
        axis.line.y = element_blank(),  axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.direction="vertical") 
  )

bp |> ggsave(filename = here::here("figures","bayes_variance_decomp_full.pdf"),width = 15.6,height = 1.66,device = cairo_pdf)



# BLUPs for rankings ------



blupdf = df |> group_by(benchmark, model_name) |> 
  summarize(score = mean(score,na.rm=T)) |>
  ungroup() |> group_by(benchmark) |> mutate(oldrank = percent_rank(score))|>
  inner_join(
    ranef(gbm, groups= "benchmark:model_name")[[1]] |> 
      as_tibble(rownames="var") |> 
      mutate(benchmark = str_match(var,
                                   paste0("(",paste0(benches,collapse = "|"),")"))[,2],
             model_name = str_split_i(var,str_c(benchmark,"_"),2)) |> 
      clean_names()
  )


blupdf |> group_by(benchmark) |> 
  summarise(kcor = cor(score,estimate_intercept,method="kendall"),
            dcor = dcor(score,estimate_intercept,bc=T)$dcor,
            scor = cor(score,estimate_intercept,method="spearman")) |> 
  kableExtra::kable(digits = 3,format = "latex", booktabs = T) |> 
  kableExtra::kable_styling(latex_options = c("hold_position","scale_down"))







blupdf_ranked |>
  select(benchmark, old_rank, new_rank) |>
  pivot_longer(
    cols = c(old_rank, new_rank),
    names_to = "type",
    values_to = "rank"
  ) |>
  mutate(type = recode(type,
                       old_rank = "Original",
                       new_rank = "Adjusted")) |> 
  ggplot(
       aes(x = type, y = rank, group = benchmark)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2) +
  scale_y_reverse() +
  theme_minimal() +
  labs(
    title = "Benchmark Ranking Changes",
    x = NULL,
    y = "Rank (1 = best)"
  )





blupdf_ranked <- blupdf |>
  mutate(
    old_rank = rank(-score, ties.method = "first"),
    new_rank = rank(-estimate_intercept, ties.method = "first"),
    delta = new_rank - old_rank
  )



plot_df <- blupdf_ranked |>
  select(benchmark,model_name,  old_rank, new_rank) |>
  pivot_longer(
    cols = c(old_rank, new_rank),
    names_to = "type",
    values_to = "rank"
  ) |>
  mutate(type = recode(type,
                       old_rank = "Old",
                       new_rank = "New")) |>
  left_join(blupdf_ranked |> select(benchmark,model_name, delta), by =c( "benchmark","model_name"))

(rp = ggplot(plot_df |> filter(benchmark%in%c("taubench_airline","swebench_verified_mini")),
             aes(x = type, y = rank, group = model_name, color = delta)) +
    geom_line(alpha = 1) +
    geom_point(size = 2) +
    geom_text_repel(
      data = subset(plot_df |> filter(benchmark%in%c("taubench_airline","swebench_verified_mini")), type == "New"),
      aes(label = model_name),
      direction = "y",
      hjust = 0,
      nudge_x = 0.1,
      size = 3,
      segment.alpha = 0.3
    ) +
    scale_y_reverse() +
    # scale_x_discrete(expand = add(mult = c(0.5, 1))) +
    xlim("Old", "New") +
    scale_color_gradient2_tableau() +
    # scale_color_gradient2(midpoint = 0) +
    facet_wrap(~benchmark,ncol=2) +
    ggpubr::theme_pubr() +
    theme(text=element_text(family="times",size = 12), 
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",) +
    labs(
      # title = "Benchmark Ranking Changes",
      y = "Rank (1 = best)",
      color = "Rank Change\n(New - Old)"
    ))

rp |> ggsave(filename = here::here("figures","benchmark_ranking_changes_tau_swe.pdf"),width = 10.5,height = 4.57,device = cairo_pdf)



## Plot ofBLUPs vs raw scores -------
(pointp = blupdf |> 
   filter(benchmark %in% c("taubench_airline","swebench_verified_mini")) |>
   mutate(estimate_intercept = plogis(estimate_intercept),
          est_error_intercept = est_error_intercept* (estimate_intercept*(1-estimate_intercept))) |>
   ggplot(aes(x=score,y=estimate_intercept,color = model_name)) + 
  geom_pointrange(aes(y=estimate_intercept,
                      ymin = estimate_intercept-est_error_intercept,
                      ymax = estimate_intercept+est_error_intercept), position = "jitter") + 
  theme_minimal() + 
   theme(text = element_text(family="times",size = 12)) +
  guides(color = guide_legend(ncol=2) ) + 
  facet_wrap(~benchmark,scales="free_x"))

pointp  |> ggsave(filename = here::here("figures","blups_vs_raw_scoresSE.pdf"),width = 11.5,height = 8.4,device = cairo_pdf)



## Plot ofBLUPs vs raw scores -------
(pointp = blupdf |> 
   mutate(estimate_intercept = plogis(estimate_intercept),
          est_error_intercept = est_error_intercept* (estimate_intercept*(1-estimate_intercept))) |>
   filter(benchmark %in% c("taubench_airline","swebench_verified_mini")) |>
   ggplot(aes(x=score,y=estimate_intercept,color = model_name)) + 
   # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
   geom_pointrange(aes(y=estimate_intercept,
                       ymin = estimate_intercept-est_error_intercept,
                       ymax = estimate_intercept+est_error_intercept), position = "jitter") + 
   #theme_minimal() + 
   theme_minimal() +
   labs(y="Estimated Score", x = "Original Mean Score") +
   theme(text = element_text(family="times",size = 13), 
         legend.position = "bottom",legend.title = element_blank(), #legend.text = element_text(size = 10),
         panel.border = element_rect(color = "grey80", fill = NA, size = 1),
         panel.spacing = unit(0.1, "cm"),
         legend.key.size = unit(0.3, "cm"), # Shrinks the symbol boxes
         legend.text = element_text(size = 8),  # Shrinks the labels
         legend.margin = margin(t = 0, r = 0, b = 0, l = 0), # Removes outer padding
         legend.spacing.y = unit(0, "cm")  # Small gap creates a "line" look
         # panel.background = element_rect(fill = "white"),
         # panel.background = element_rect(fill = "gray90"),
         #panel.spacing = unit(0.5, "lines")
   ) +
   # guides(color = guide_legend(ncol=2) ) + 
   facet_wrap(~benchmark,scales="free_x"))

pointp  |> ggsave(filename = here::here("figures","blups_vs_raw_scoresSE_tau_swe_nolegend.pdf"),width = 8.74,height = 5.32,device = cairo_pdf)




# check max/min ranks for statistical distinguishability

blupdf |> group_by(benchmark) |> 
  mutate(maxdifflo = (max(estimate_intercept)-estimate_intercept) - 
           1.28*sqrt(est_error_intercept^2 + est_error_intercept[max(estimate_intercept)==estimate_intercept]^2), 
         maxdiffhi = (max(estimate_intercept)-estimate_intercept) + 
           1.28*sqrt(est_error_intercept^2 + est_error_intercept[max(estimate_intercept)==estimate_intercept]^2), 
         mindifflo = (estimate_intercept-min(estimate_intercept)) - 
           1.28*sqrt(est_error_intercept^2 + est_error_intercept[min(estimate_intercept)==estimate_intercept]^2),
         mindiffhi = (estimate_intercept-min(estimate_intercept)) + 
           1.28*sqrt(est_error_intercept^2 + est_error_intercept[min(estimate_intercept)==estimate_intercept]^2),
         is_sig_max = if_else(sign(maxdifflo)*sign(maxdiffhi)<0 , 1,0),
         is_sig_min = if_else( sign(mindifflo)*sign(mindiffhi)<0, 1, 0)
         ) |> 
  summarise(prop_insig_max = 1-mean(is_sig_max), prop_insig_min = 1-mean(is_sig_min)) |>
  kableExtra::kable(digits = 3,format = "latex", booktabs = T) |>
  kableExtra::kable_styling(latex_options = c("hold_position","scale_down"))





## ranking correlations for BLUP ranks across draws -------
ranef(gbm, groups= "benchmark:model_name",summary = F)[[1]] |> 
  as_tibble(rownames="draw") |> 
  pivot_longer(-draw,names_to = "var", values_to = "est") |> 
  mutate(benchmark = str_match(var,
                               paste0("(",paste0(benches,collapse = "|"),")"))[,2],
         model_name = str_split_i(var,str_c(benchmark,"_"),2)) |> 
  clean_names() |> 
  pivot_wider(names_from = draw,values_from = est,names_prefix = "d") |> 
  select(-var) |> 
  group_by(benchmark) |> 
  group_modify(~ as_tibble(median((as.vector(cor((.x |> select(-var,-model_name)),method = "spearman")))))) |> 
  pivot_wider(names_from = benchmark, values_from = value) |>
  kableExtra::kable(digits = 3,format = "latex", booktabs = T) |>
  kableExtra::kable_styling(latex_options = c("hold_position","scale_down"))



## ranking correlations for BLUP ranks across draws -------
ranef(gbm, groups= "benchmark:model_name",summary = F)[[1]] |> 
  as_tibble(rownames="draw") |> 
  pivot_longer(-draw,names_to = "var", values_to = "est") |> 
  mutate(benchmark = str_match(var,
                               paste0("(",paste0(benches,collapse = "|"),")"))[,2],
         model_name = str_split_i(var,str_c(benchmark,"_"),2)) |> 
  clean_names() |> 
  pivot_wider(names_from = draw,values_from = est,names_prefix = "d") |> 
  select(-var) |> 
  group_by(benchmark) |> 
  mutate(across(-model_name, \(x) percent_rank(x))) |>ungroup() |> 
  mutate(model_name = str_remove_all(model_name,".Intercept")) |>
  pivot_longer(-c(model_name,benchmark),names_to = "draw",values_to = "rank") |> 
  group_by(benchmark,model_name) |>
  summarize(median_rank = median(rank), 
            mad = mad(rank),
            q5 = quantile(rank, 0.05), 
            q95 = quantile(rank, 0.95), 
            q25 = quantile(rank, 0.25), 
            q75 = quantile(rank, 0.75),
            q2_5 = quantile(rank, 0.025),
            q97_5 = quantile(rank, 0.975) ) |> 
  full_join(df |> group_by(benchmark, model_name) |> 
              summarize(score = mean(score,na.rm=T)) |>
              ungroup() |> group_by(benchmark) |> mutate(oldrank = percent_rank(score))) |> 
  filter(benchmark %in% c("taubench_airline","swebench_verified_mini")) |>
  ggplot(aes(x=oldrank,y=median_rank,color=model_name)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  # geom_point() + 
  # geom_pointrange(aes(ymin=q25,ymax=q75),position = "jitter") + # IQR
  geom_pointrange(aes(ymin=median_rank-1.4826*mad,ymax=median_rank+1.4826*mad),position = "jitter") + # empirical mad estimated as SE
  # geom_pointrange(aes(ymin=q5,ymax=q95),position = "jitter") + 
  facet_wrap(~benchmark) + 
  theme_minimal() + 
  guides(color = guide_legend(ncol=6)) +
  labs(x = "Original/Published Rank", y = "Posterior Rank (Median, 1.48xMAD)") +
  theme(text = element_text(family="times",size = 14), 
        panel.border = element_rect(color = "grey80", fill = NA, size = 1),
        legend.position = "bottom",         legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"), # Shrinks the symbol boxes
        legend.text = element_text(size = 8),  # Shrinks the labels
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0), # Removes outer padding
        legend.spacing.y = unit(0, "cm") # Small gap creates a "line" look
        )


  # 
  # group_modify(~ as_tibble(median((as.vector(cor((.x |> select(-model_name)),method = "spearman")))))) |> 
  # pivot_wider(names_from = benchmark, values_from = value) |>
  # kableExtra::kable(digits = 3,format = "latex", booktabs = T) |>
  # kableExtra::kable_styling(latex_options = c("hold_position","scale_down"))
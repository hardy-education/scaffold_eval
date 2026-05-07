## Benchmark level Variance Decompositions
## 
## 



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

df = fread(here("data/hal_parsed_response_matrix7.csv")) |> 
  mutate(model_name = str_split_i(model_name,":",1),
         model_name = case_when(
           model_name == "DeepSeek-R1" ~ "deepseek-r1",
           model_name == "claude-sonnet-4.5" ~ "claude-sonnet-4-5",
           TRUE ~ model_name
         ),
         model_name = str_c(model_name,reasoning_effort)) # |> 

benches =unique( df$benchmark)
item_counts = df |> select(benchmark,task_id) |> distinct() |> group_by(benchmark) |> summarize(n = n_distinct(task_id))

# Functions -------
## Bayes ---------


## BAYES create combination stats -------
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
                                # null = 0.000001,
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

convert_to_percent_of_var_draws = function(draws,term){
  draws[[term]]^2 / sum(draws[[term]]^2)
}







# Models =========
## LME model-------

lmmods = list()
lmsums = list()
for(b in benches){
  tic(b)
  lmmods[[b]] = df |> 
    filter(benchmark==b) |> 
    lmer(score ~ 1
                 + (1|task_id)
                 + (1|model_name)
                 + (1|agent_name)
                 + (1|task_id:model_name)
                 + (1|task_id:agent_name)
                 + (1|model_name:agent_name)
                 # + (1|task_id:model_name:agent_name)
                 ,
                 na.action = na.exclude,
                 data = _,
                 control = lmerControl(optimizer = "bobyqa", 
                                       calc.derivs = FALSE,
                                       optCtrl = list(maxfun = 2e5))
  )
  
  lmmods[[b]] |> summary()
  
  lmsums[[b]] = lmmods[[b]]  |>
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
  lmsums[[b]]
  toc()
}

lmsums = lmsums |> list_rbind(names_to = "benchmark")







## GLME model-------
glmmods = list()
glmsums = list()
for(b in benches){
  tic(b)
  glmmods[[b]] = df |> 
    filter(benchmark==b) |> 
    glmer(score ~ 1
                 + (1|task_id)
                 + (1|model_name)
                 + (1|agent_name)
                 + (1|task_id:model_name)
                 + (1|task_id:agent_name)
                 + (1|model_name:agent_name)
                 # + (1|task_id:model_name:agent_name)
                 ,
                 na.action = na.exclude,
                 data = _,
                 family = "binomial",
                 nAGQ=0,
                 control = glmerControl(optimizer = "bobyqa", 
                                        calc.derivs = FALSE,
                                        optCtrl = list(maxfun = 2e5))
)

  glmmods[[b]] |> summary()

  glmsums[[b]] = glmmods[[b]]  |>
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
  glmsums[[b]]
toc()
}

glmsums = glmsums |> list_rbind(names_to = "benchmark")

# Bayes --------
## linear model

blmmods = list()
blmdraws = list()
blmsums = list()
for(b in benches){
  tic(b)
  blmmods[[b]] = df |> 
    filter(benchmark==b) |>
    brm(score ~ 1
               + (1|task_id)
               + (1|model_name)
               + (1|agent_name)
               + (1|task_id:model_name)
               + (1|task_id:agent_name)
               + (1|model_name:agent_name)
               # + (1|task_id:model_name:agent_name)
               ,
               chains = 4,
               cores = 4,
               thin = 5,
               iter = 2000,
               backend = "cmdstanr",
               control = list(adapt_delta = 0.95),
               data = _,
               threads = threading(2),
               save_all_pars= TRUE,
               stan_model_args = list(stanc_options = list("O1")),
               file = here::here("data",glue::glue("brmfit_{b}.rds")),
               file_refit = "always"
               
               
)

  blmmods[[b]] = readRDS(here::here("data",glue::glue("brmfit_{b}.rds")))

  blmmods[[b]] |> summary()

toc()

# get draws
blmdraws[[b]] = as_draws_rvars(blmmods[[b]])

# make pct summary
blmsums[[b]] = lapply(blmdraws[[b]][str_starts(names(blmdraws[[b]]),"sd|sig")], function(X) make_draws_stats(X^2/reduce(blmdraws[[b]][str_starts(names(blmdraws[[b]]),"sd|sig")],\(acc,nxt) acc + (nxt^2)))) |> list_rbind(names_to = "var")
blmsums[[b]]
}

blmsums = blmsums |> list_rbind(names_to = "benchmark")

nitems = 100

get_bmod_ep2 = function(drawlist,nitems=100){
  m_ep2 = drawlist[["sd_model_name__Intercept"]]^2/(drawlist[["sd_model_name__Intercept"]]^2+ 
                                                      drawlist[["sd_task_id:model_name__Intercept"]]^2/nitems+
                                                      drawlist[["sd_model_name:agent_name__Intercept"]]^2+
                                                      # drawlist[["sd_task_id:model_name:agent_name__Intercept"]]^2/nitems + 
                                                      drawlist[["sigma"]]^2/nitems)  
  make_draws_stats(m_ep2)
}
Vectorize(get_bmod_ep2)


get_bagnt_ep2 = function(drawlist,nitems=100){
  a_ep2 = drawlist[["sd_agent_name__Intercept"]]^2/(drawlist[["sd_agent_name__Intercept"]]^2+
                                                      drawlist[["sd_task_id:agent_name__Intercept"]]^2/nitems+
                                                      drawlist[["sd_model_name:agent_name__Intercept"]]^2+
                                                      # drawlist[["sd_task_id:model_name:agent_name__Intercept"]]^2/nitems + 
                                                      drawlist[["sigma"]]^2/nitems)  
  make_draws_stats(a_ep2)
}
Vectorize(get_bagnt_ep2)

get_bagmod_ep2 = function(drawlist,nitems=100){
  (drawlist[["sd_model_name:agent_name__Intercept"]]^2/(drawlist[["sd_model_name:agent_name__Intercept"]]^2+
                                                          # drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2/nitems + 
                                                          drawlist[["sigma"]]^2/nitems)) |> 
    make_draws_stats()
}
Vectorize(get_bagmod_ep2)




## reverse logit transformed -----------


plogis = function(x) exp(x)/(1+exp(x))


get_rmod_ep2 = function(drawlist,nitems=100){
  m_ep2 = plogis(drawlist[["sd_model_name__Intercept"]]^2)/(plogis(drawlist[["sd_model_name__Intercept"]]^2)+ 
                                                      plogis(drawlist[["sd_task_id:model_name__Intercept"]]^2)/nitems+
                                                      plogis(drawlist[["sd_model_name:agent_name__Intercept"]]^2)+
                                                      # drawlist[["sd_task_id:model_name:agent_name__Intercept"]]^2/nitems + 
                                                      plogis(drawlist[["sigma"]]^2)/nitems)  
  make_draws_stats(m_ep2)
}
Vectorize(get_rmod_ep2)



get_ragnt_ep2 = function(drawlist,nitems=100){
  a_ep2 = plogis(drawlist[["sd_agent_name__Intercept"]]^2)/(plogis(drawlist[["sd_agent_name__Intercept"]]^2)+
                                                      plogis(drawlist[["sd_task_id:agent_name__Intercept"]]^2)/nitems+
                                                      plogis(drawlist[["sd_model_name:agent_name__Intercept"]]^2)+
                                                      # drawlist[["sd_task_id:model_name:agent_name__Intercept"]]^2/nitems + 
                                                      plogis(drawlist[["sigma"]]^2)/nitems)  
  make_draws_stats(a_ep2)
}
Vectorize(get_ragnt_ep2)

get_ragmod_ep2 = function(drawlist,nitems=100){
  (plogis(drawlist[["sd_model_name:agent_name__Intercept"]]^2)/((plogis(drawlist[["sd_model_name:agent_name__Intercept"]]^2)+
                                                          # drawlist[["sd_benchmark:task_id:model_name:agent_name__Intercept"]]^2/nitems + 
                                                          plogis(drawlist[["sigma"]]^2)/nitems)) )|> 
    make_draws_stats()
}
Vectorize(get_ragmod_ep2)





## quick plot tasks per benchmark needed rough 

## Bayes Generalized -----
# 



gblmmods = list()
gblmdraws = list()
gblmsums = list()

for(b in benches){
  tic(b)
  gblmmods[[b]] = df |> 
    filter(benchmark==b) |> brm(score ~ 1
                + (1|task_id)
                + (1|model_name)
                + (1|agent_name)
                + (1|task_id:model_name)
                + (1|task_id:agent_name)
                + (1|model_name:agent_name)
                # + (1|task_id:model_name:agent_name)
                ,
                chains = 4,
                cores = 4,
                thin = 2,
                iter = 1000,
                backend = "cmdstanr",
                control = list(adapt_delta = 0.95),
                data = _,
                family = bernoulli("logit"),
                threads = threading(2),
                save_pars= save_pars(all=T),
                stan_model_args = list(stanc_options = list("O1")),
                file = here::here("data",glue::glue("gbrmfit_{b}.rds")),
                file_refit = "always"
)

  gblmmods[[b]] = readRDS(here::here("data",glue::glue("gbrmfit_{b}.rds")))

  gblmmods[[b]] |> summary()

toc()

# get draws
gblmdraws[[b]] = as_draws_rvars(gblmmods[[b]])
gblmdraws[[b]]$sigma = posterior::rvar(rnorm(1000, mean = pi^2/3, sd = 0))

# make pct summary
gblmsums[[b]] = lapply(gblmdraws[[b]] [str_starts(names(gblmdraws[[b]] ),"sd|sig")], function(X) make_draws_stats(X^2/reduce(gblmdraws[[b]] [str_starts(names(gblmdraws[[b]] ),"sd|sig")],\(acc,nxt) acc + (nxt^2)))) |> list_rbind(names_to = "var")
gblmsums[[b]]
}

gblms = list_rbind(gblmsums,names_to = "benchmark")

# blmsums = blmsums |> list_rbind(names_to = "benchmark")

# 
# 
# ### Bayes MV --------
# 
# 
# tic()
# bm2 = dft |> brm(mvbind(score,total_tokens,(total_tokens)^2) ~ 1
#                  + (1|benchmark)
#                  + (1|benchmark:task_id)
#                  + (1|model_name)
#                  + (1|agent_name)
#                  + (1|benchmark:task_id:model_name)
#                  + (1|benchmark:task_id:agent_name)
#                  + (1|model_name:agent_name)
#                  + (1|benchmark:task_id:model_name:agent_name)
#                  + (1|benchmark:model_name)
#                  + (1|benchmark:agent_name)
#                  + (1|benchmark:model_name:agent_name)
#                  ,
#                  chains = 4,
#                  cores = 4,
#                  thin = 5,
#                  iter = 2000,
#                  backend = "cmdstanr",
#                  control = list(adapt_delta = 0.95),
#                  data = _,
#                  threads = threading(2),
#                  save_pars = save_pars("all"),
#                  stan_model_args = list(stanc_options = list("O1")),
#                  file = here::here("data","brmfit_tok.rds"),
#                  file_refit = "always"
#                  
#                  
# )
# 
# bm2 |> summary()
# 
# toc()
# 
# # get draws
# bmd2 = as_draws_rvars(bm2)
# 
# # make pct summary
# bms2 = lapply(bmd2[str_starts(names(bmd2),"sd|sig")], function(X) make_draws_stats(X^2/reduce(bmd2[str_starts(names(bmd2),"sd|sig")],\(acc,nxt) acc + (nxt^2)))) |> list_rbind(names_to = "var")
# bms2
# 
# 
# 
# 






n_iter = 3
n_i = 100
prop_m =1
alph = 1

## full set with 7 benchmarks iteration takes 1hr to run
discolist = list()

tic("all iterations")
# for (i in 1:1){
for (b in benches){
  tic(glue::glue("Iteration  Benchmark: {b}"))
  cat(glue::glue("Iteration  Benchmark: {b}\n"))
  
  # set.seed(i*1234)
  

  dat = df |> dtplyr::lazy_dt() |> 
    filter(benchmark==b) |> 
    drop_na(score) |>
    mutate(resp = as.numeric(score),
           item = as.factor(interaction(benchmark,task_id)),
           agent = as.factor(agent_name),
           model = as.factor(model_name),
           item_agent = fct_drop(interaction(benchmark,item,agent)),
           item_model = fct_drop(interaction(benchmark,item,model)),
           model_agent = fct_drop(interaction(model,agent)),
           # item_model_agent = fct_drop(interaction(benchmark,item,model,agent))
           ) |> 
    select(-any_of(c("run_id","model_name","agent_name","score","task_id","reasoning_effort","benchmark")))
  
  disco_fit = disco(dat |> select(resp) |> as.data.table() ,
                    factors = dat |> select(-resp) |> 
                      as.data.table(), 
                    distance = F,
                    index = alph,
                    R=0)
  
  discolist[[b]] = disco_fit$stats |> 
    as_tibble() |> 
    select(-any_of("p-value")) |> 
    mutate(facet = disco_fit$factor.names,Withins = Within, 
           Within = disco_fit$within, Total = disco_fit$total, 
           Ep2 = Trt / (Trt + Withins), SNR = Ep2/(1-Ep2),
           Ep2hat = (Trt/df1) / ((Trt/df1) + Withins/df2), 
           SNRhat = Ep2hat/(1-Ep2hat),
           pct = Trt/(Total-Within),
           mean = Trt/df1, var = Withins/df2) |> 
    select(facet,Ep2,SNR,everything()) |> 
    mutate(benchmark = b,
    n_items = n_unique(dat$item),
    alpha = alph)
  
  # L[[length(L)+1]] = disco_fit_df
  print(discolist[[b]])
  toc()
  # }
}

disco_bench_res = bind_rows(discolist)

disco_bench_res |>
  write_csv(file = here::here("data", "results","disco_all_benches.csv"))

disco_bench_res = read_csv(here::here("data","results","disco_all_benches.csv"))

disco_bench_res
toc()


combbenchdf = disco_bench_res |> mutate(facet = str_replace_all(str_replace_all(facet,"item","bench_item"),"bench_","benchmark_"),method = "nonp") |> select(method,benchmark, facet,pct) |> bind_rows(lmsums |> mutate(facet = str_replace_all(str_replace_all(str_remove_all(str_replace_all(grp,"task_id","item"),regex("sd_|__Intercept|_name|_id")),":","_"),"Residual","sigma"),method = "lme",pct=pct/100) |> select(facet, benchmark,method,pct)) |>bind_rows(glmsums |> mutate(facet = str_replace_all(str_replace_all(str_remove_all(str_replace_all(grp,"task_id","item"),regex("sd_|__Intercept|_name|_id")),":","_"),"Residual","sigma"),method = "glme",pct=pct/100) |> select(facet,benchmark, method,pct)) |> bind_rows(gblms |> mutate(facet = str_replace_all(str_remove_all(str_replace_all(var,"task_id","item"),regex("sd_|__Intercept|_name|_id")),":","_"),pct = mean,method="bayes") |> select(facet,benchmark,method,pct)) |> mutate(facet=if_else(facet == "benchmark_item_model_agent","sigma",facet))
combbenchdf |> ggplot(aes(x= method, y = pct, fill = facet)) + geom_col() + theme_minimal()

combbenchdf |> mutate(facet = str_replace_all(facet, "benchmark_item", "item")) |>
  pivot_wider(names_from = method, values_from = pct) |>
  select(benchmark, facet, lme, glme, bayes, nonp) |>
  kableExtra::kable(digits = 3,format = "latex", booktabs = T) |> 
  kableExtra::kable_styling(latex_options = c("hold_position","scale_down"))

plot_ni = 350
plotdfs = list()
for(b in benches){
  plotdfs[[b]] = tibble(nitems=1:plot_ni,
                        m_ep2 = get_bmod_ep2(blmdraws[[b]],nitems = nitems) |> 
                          pull(median), 
                        a_ep2 = get_bagnt_ep2(blmdraws[[b]],nitems = nitems) |> 
                          pull(median)) |> 
    pivot_longer(m_ep2:a_ep2,names_to = "facet",values_to = "bayes")|> 
    inner_join(tibble(nitems=1:plot_ni,
                      m_ep2 = get_bmod_ep2(gblmdraws[[b]],nitems = nitems) |>
                        pull(median),
                      a_ep2 = get_bagnt_ep2(gblmdraws[[b]],nitems = nitems) |>
                        pull(median)) |>
                 pivot_longer(m_ep2:a_ep2,names_to = "facet",values_to = "gbayes") |>
                 distinct())  |>
    inner_join(tibble(nitems=1:plot_ni, 
                      m_ep2 = lmsums[[b]][lmsums[[b]]$grp=="model_name","vcov"][[1]]/(lmsums[[b]][lmsums[[b]]$grp=="model_name","vcov"][[1]]+lmsums[[b]][lmsums[[b]]$grp=="model_name_agent_name","vcov"][[1]]+lmsums[[b]][lmsums[[b]]$grp=="task_id_model_name","vcov"][[1]]/(nitems)+lmsums[[b]][lmsums[[b]]$grp=="Residual","vcov"][[1]]/(nitems)), 
                      a_ep2 =lmsums[[b]][lmsums[[b]]$grp=="agent_name","vcov"][[1]]/(lmsums[[b]][lmsums[[b]]$grp=="agent_name","vcov"][[1]]+lmsums[[b]][lmsums[[b]]$grp=="model_name_agent_name","vcov"][[1]]+lmsums[[b]][lmsums[[b]]$grp=="task_id_agent_name","vcov"][[1]]/(nitems)+lmsums[[b]][lmsums[[b]]$grp=="Residual","vcov"][[1]]/(nitems)) ) |> 
                 pivot_longer(m_ep2:a_ep2,names_to = "facet",values_to = "lmer") |> 
                 distinct()) |> 
    inner_join(tibble(nitems=1:plot_ni, 
                      m_ep2 = glmsums[[b]][glmsums[[b]]$grp=="model_name","vcov"][[1]]/(glmsums[[b]][glmsums[[b]]$grp=="model_name","vcov"][[1]]+glmsums[[b]][glmsums[[b]]$grp=="model_name_agent_name","vcov"][[1]]+glmsums[[b]][glmsums[[b]]$grp=="task_id_model_name","vcov"][[1]]/(nitems)+glmsums[[b]][glmsums[[b]]$grp=="Residual","vcov"][[1]]/(nitems)), 
                      a_ep2 =glmsums[[b]][glmsums[[b]]$grp=="agent_name","vcov"][[1]]/(glmsums[[b]][glmsums[[b]]$grp=="agent_name","vcov"][[1]]+glmsums[[b]][glmsums[[b]]$grp=="model_name_agent_name","vcov"][[1]]+glmsums[[b]][glmsums[[b]]$grp=="task_id_agent_name","vcov"][[1]]/(nitems)+glmsums[[b]][glmsums[[b]]$grp=="Residual","vcov"][[1]]/(nitems)) ) |> 
                 pivot_longer(m_ep2:a_ep2,names_to = "facet",values_to = "glmer") |> distinct()) |> 
    inner_join(tibble(nitems=1:plot_ni, m_ep2 =  (disco_bench_res |> filter(benchmark==b,facet=="model" ) |> pull(Trt))/((disco_bench_res |> filter(benchmark==b,facet=="model" ) |> pull(Trt))+(disco_bench_res |> filter(benchmark==b,facet=="model_agent" ) |> pull(Trt))+(disco_bench_res |> filter(benchmark==b,facet=="item_model" ) |> pull(Trt))/(nitems)+(disco_bench_res |> filter(benchmark==b,facet=="model" ) |> pull(Withins))/(nitems)),
                      a_ep2=(disco_bench_res |> filter(benchmark==b,facet=="agent" ) |> pull(Trt))/((disco_bench_res |> filter(benchmark==b,facet=="agent" ) |> pull(Trt))+(disco_bench_res |> filter(benchmark==b,facet=="model_agent" ) |> pull(Trt))+(disco_bench_res |> filter(benchmark==b,facet=="item_agent" ) |> pull(Trt))/(nitems)+(disco_bench_res |> filter(benchmark==b,facet=="agent" ) |> pull(Withins))/(nitems))) |> pivot_longer(m_ep2:a_ep2,names_to = "facet",values_to = "nonp") |> distinct()) |> pivot_longer(bayes:nonp, names_to = "method",values_to = "rel")
}

plotdf = list_rbind(plotdfs,names_to = "benchmark")

plotdf |> ggplot(aes(x=nitems,y=rel,color=method,linetype=facet)) + geom_line() + theme_minimal() + scale_x_log10() + facet_wrap(~benchmark)

baydf = list()
for(b in benches){
  nis = min(item_counts[item_counts$benchmark==b,"n"][[1]],plot_ni)
  
  baydf[[b]] = tibble(nitems=1:nis, m_ep2 = get_bmod_ep2(blmdraws[[b]],nitems = nitems) |> pull(median),   
                   m_ep2_low = get_bmod_ep2(blmdraws[[b]],nitems = nitems) |> pull(HDI_low),
                   m_ep2_hi = get_bmod_ep2(blmdraws[[b]],nitems = nitems) |> pull(HDI_high),
                   a_ep2 = get_bagnt_ep2(blmdraws[[b]],nitems = nitems) |> pull(median),
                   a_ep2_low = get_bagnt_ep2(blmdraws[[b]],nitems = nitems) |> pull(HDI_low),
                   a_ep2_hi = get_bagnt_ep2(blmdraws[[b]],nitems = nitems) |> pull(HDI_high),
                   am_ep2 = get_bagmod_ep2(blmdraws[[b]],nitems = nitems) |> pull(median),
                   am_ep2_low = get_bagmod_ep2(blmdraws[[b]],nitems = nitems) |> pull(HDI_low),
                   am_ep2_hi = get_bagmod_ep2(blmdraws[[b]],nitems = nitems) |> pull(HDI_high),
)
}

baydf = list_rbind(baydf,names_to = "benchmark")

baydf |> 
  inner_join(item_counts) |> 
  filter(nitems<=n) |> 
  pivot_longer(m_ep2:am_ep2_hi,
                          names_to = c("facet","stat","val"),
                          names_sep = "_",
                          values_to = "bayes") |> 
  mutate(val = if_else(is.na(val),"est",val)) |>
  pivot_wider(names_from = val, values_from = bayes) |> 
  filter(facet=="m") |> 
  ggplot(aes(x=nitems,y=est,color = facet,fill=facet)) +
  geom_ribbon(aes(x = nitems,y=est,ymin = low,ymax = hi),alpha=0.2,linewidth=0) +
  geom_line() +
  theme_minimal()+
  scale_x_log10()+
  facet_wrap(~benchmark, scales="free_x")



gbaydf = list()
for(b in benches){
  nis = min(item_counts[item_counts$benchmark==b,"n"][[1]],plot_ni)
  
  gbaydf[[b]] = tibble(nitems=1:nis, m_ep2 = get_rmod_ep2(gblmdraws[[b]],nitems = nitems) |> pull(median),   
                      m_ep2_low = get_rmod_ep2(gblmdraws[[b]],nitems = nitems) |> pull(HDI_low),
                      m_ep2_hi = get_rmod_ep2(gblmdraws[[b]],nitems = nitems) |> pull(HDI_high),
                      a_ep2 = get_ragnt_ep2(gblmdraws[[b]],nitems = nitems) |> pull(median),
                      a_ep2_low = get_ragnt_ep2(gblmdraws[[b]],nitems = nitems) |> pull(HDI_low),
                      a_ep2_hi = get_ragnt_ep2(gblmdraws[[b]],nitems = nitems) |> pull(HDI_high),
                      am_ep2 = get_ragmod_ep2(gblmdraws[[b]],nitems = nitems) |> pull(median),
                      am_ep2_low = get_ragmod_ep2(gblmdraws[[b]],nitems = nitems) |> pull(HDI_low),
                      am_ep2_hi = get_ragmod_ep2(gblmdraws[[b]],nitems = nitems) |> pull(HDI_high),
  )
}

gbaydf = list_rbind(gbaydf,names_to = "benchmark")


(mvamp = gbaydf |> 
    inner_join(item_counts) |> 
    filter(nitems<=n) |> 
    pivot_longer(m_ep2:am_ep2_hi,
                 names_to = c("facet","stat","val"),
                 names_sep = "_",
                 values_to = "bayes") |> 
    mutate(val = if_else(is.na(val),"est",val)) |>
    pivot_wider(names_from = val, values_from = bayes) |> 
    filter(facet %in% c("m","am")) |> 
    mutate(facet = if_else(facet=="m","model","agnt-mod")) |> 
    ggplot(aes(x=nitems,y=est,color = facet,fill=facet)) +
    geom_ribbon(aes(x = nitems,y=est,ymin = low,ymax = hi),alpha=0.2,linewidth=0) +
    geom_line() +
    theme_minimal()+
    scale_x_log10()+
    facet_wrap(~benchmark, scales="free_x",nrow=1)+
    ylim(c(0, 1)) +
    labs(x = "Number of Agentic Tasks", y = expression("Estimated Reliability: E"*hat(rho)^2)) +
    theme(text=element_text(family="times",size = 15), legend.position = "right"))

mvamp |> ggsave(filename = here::here("figures","bayes_CI_per_bench_mvamp.pdf"),width = 14.3, height = 3.55,device = cairo_pdf)






(p = gbaydf |> 
  inner_join(item_counts) |> 
  filter(nitems==n) |> 
  pivot_longer(m_ep2:m_ep2_hi,
               names_to = c("facet","stat","val"),
               names_sep = "_",
               values_to = "bayes") |> 
  mutate(val = if_else(is.na(val),"est",val)) |>
  pivot_wider(names_from = val, values_from = bayes) |> 
  filter(facet=="m") |> 
  ggplot(aes(x=benchmark,y=est,color = benchmark,fill=benchmark)) +
  geom_pointrange(aes(x = benchmark,y=est,ymin = low,ymax = hi), 
                  linewidth = 1.4, size =1.4,
                  position = position_dodge(width = 0.1)) +
  # geom_line() +
  theme_minimal()+
  # scale_x_log10()+
  # facet_wrap(~benchmark, scales="free_x",nrow=1)+
  ylim(c(0, 1)) +
  # ggthemes::scale_color_colorblind()+
  labs(y = expression("Estimated Reliability: E"*hat(rho)^2)) +
  theme(text=element_text(family="times",size = 15), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none"))

p |> ggsave(filename = here::here("figures","bayes_CI_per_bench.pdf"),width = 7, height = 5,device = cairo_pdf)

# stacked bar charts that sum to 100% ----- 
# for each benchmark, showing the proportion of variance attributed to each facet (model, agent, item, interactions, residual) according to each method (nonparametric, lme, glme, bayes). Facets on x-axis, percentage of variance on y-axis, fill by method. Facet wrap by benchmark.

gblmsums |> list_rbind(names_to = "benchmark") |> 
  # filter(var != "sigma") |> 
  group_by(benchmark) |>
  mutate(across(mean:MAP_Estimate,\(x) x / sum(x)),facet = str_replace_all(str_remove_all(str_replace_all(var,"task_id","item"),regex("sd_|__Intercept|_name|_id")),":","_")) |>
  ungroup() |>
  ggplot(aes(x=benchmark,y=mean,fill=facet)) + 
  geom_col(stat="identity") + 
  theme_minimal() + 
  labs(x = "Benchmark", y = "Proportion of Variance Explained", fill = "Component") + 
  scale_fill_manual(values = component_cols,labels = component_labs) +
  scale_y_continuous(labels = scales::percent) + 
  guides(fill = guide_legend(title = "Variance Component",nrow=1)) +
  theme(text=element_text(family="times",size = 15), legend.position = "bottom")
  
# disco_bench_res |> group_by(benchmark) |> 
#   mutate(across(Trt,\(x) x / sum(x))) |> ungroup() |> 
#   ggplot(aes(x=benchmark,y=Trt,fill=facet)) +
#   geom_col(stat="identity") +
#   theme_minimal() +
#   labs(x = "Benchmark", y = "Proportion of Variance Explained", fill = "Component") +
#   scale_fill_manual(values = component_cols,labels = component_labs) +
#   scale_y_continuous(labels = scales::percent) +
#   guides(fill = guide_legend(title = "Variance Component",nrow=1)) +
#   theme(text=element_text(family="times",size = 15), legend.position =
#          "bottom")



gblmsums |> list_rbind(names_to = "benchmark") |> 
  # filter(var != "sigma") |> 
  group_by(benchmark) |>
  mutate(across(mean:MAP_Estimate,\(x) x / sum(x)),facet = str_replace_all(str_remove_all(str_replace_all(var,"task_id","item"),regex("sd_|__Intercept|_name|_id")),":","_")) |>
  ungroup() |> select(-variable,-var,-starts_with("ROPE"),-c(ess_bulk:p_MAP),-CI) |> 
  select(benchmark,facet,everything()) |>
  kableExtra::kable(digits = 3,format = "latex", booktabs = T) |> 
  kableExtra::kable_styling(latex_options = c("hold_position","scale_down"))






# BLUP ranking extraction ------
blups = list()
for(b in benches){
  blups[[b]] = glmmods[[b]] |> ranef() |> as_tibble() |> 
    filter(grpvar=="model_name") |> select(grp,condval,condsd) 
}

blupdf = df |> group_by(benchmark, model_name) |> summarize(score = mean(score,na.rm=T)) |> 
  inner_join(list_rbind(blups,names_to = "benchmark"))


# 
# |> 
#                group_by(benchmark, grp) |> summarize(condval = mean(condval))) |> 
#   ggplot(aes(x=condval,y=score)) + geom_point() + facet_wrap(~benchmark,scales="free")
# 
# 
# allblups =  ranef(gbm ) |> as_tibble() 
# 




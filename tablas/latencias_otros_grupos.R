# ---- librerias ----

library(brms)
library(emmeans)
library(tidybayes)
library(modelr)
library(gridExtra)
library(bayesplot)
library(hexbin)
library(patchwork)
library(grid)
library(flextable)
library(webshot)
library(gtsummary)
library(sjPlot)
library(lme4)
library(report)
library(rstatix)
library(ggpubr)
library(performance)


# FLX ----

## ---- dataframes y setup ----

source("./scripts/02_dataframe_water_maze.R")

###
latencias_pre <- df_wm_otros |> 
  select(latencia, id, tratamientos, stage, dia, tiempo, prueba) |> 
  dplyr::filter(str_detect(prueba, "Entr"))


latencias_post <- df_wm_otros |> 
  select(latencia, id, tratamientos, stage, dia, tiempo, prueba) |> 
  dplyr::filter(str_detect(prueba, "Reversa"))


## options ----

options(warn=0)
options( mc.cores=parallel::detectCores() )
Sys.time() # to be used to time the analysis.


## censored ----

# Crear variable 'censored' cuando la latencia = 60 segundos


latencias_pre$censurado = ifelse(latencias_pre$latencia== 60, 1, 0)
latencias_pre$Censurado = ifelse(latencias_pre$latencia==60, "SI", "NO")

latencias_post$censurado = ifelse(latencias_post$latencia== 60, 1, 0)
latencias_post$Censurado = ifelse(latencias_post$latencia==60, "SI", "NO")



## ---- scatters ----



lat_scat_pre <- ggplot(latencias_pre, aes(x = tiempo, y = latencia, color = Censurado)) + 
  geom_jitter(aes(shape = Censurado), size = 5) +
  facet_grid(~ tratamientos) +
  ylim(0, 60) +
  scale_shape_manual(values = c(20, 4)) +
  scale_color_manual(values = c("#A5DEE5", "#850E35")) +
  xlab("Prueba") + 
  ylab("Segundos") +
  # labs(
  #   # title = "Latencia de Escape", 
  # tag = "A") +
  ggthemes::theme_base()

lat_scat_pre


lat_scat_post <- ggplot(latencias_post, aes(x = tiempo, y = latencia, color = Censurado)) + 
  geom_jitter(aes(shape = Censurado), size = 5) +
  facet_grid(~ tratamientos) +
  ylim(0, 60) +
  scale_shape_manual(values = c(20, 4)) +
  scale_color_manual(values = c("#A5DEE5", "#850E35")) +
  xlab("Prueba") + 
  ylab("Segundos") +
  # labs(
  #   # title = "Latencia de Escape", 
  # tag = "A") +
  ggthemes::theme_base()

lat_scat_post

## ---- violin ----



colours <- scales::viridis_pal(option = "viridis")(10)



grad_ungroup <- linearGradient(colours, group = FALSE)



lat_violin_pre <-  ggplot(latencias_pre, aes(x=factor(dia), y= latencia)) +
  # geom_violin(scale="width", aes(fill = latencia)) +
  geom_violin(scale="width", fill = grad_ungroup) +
  facet_grid(~ tratamientos) +
  # facet_wrap(~ tratamientos, scales='free') +
  labs(
    #title = "Latencia de Escap",
    x = "Prueba",
    y = "Segundos") +
  # theme_classic()
  # ggthemes::theme_clean()
  ggthemes::theme_base()
# ggthemes::theme_few()
# 
lat_violin_pre

lat_violin_post <-  ggplot(latencias_post, aes(x=factor(dia), y= latencia)) +
  # geom_violin(scale="width", aes(fill = latencia)) + 
  geom_violin(scale="width", fill = grad_ungroup) + 
  facet_grid(~ tratamientos) +
  # facet_wrap(~ tratamientos, scales='free') +
  labs(
    #title = "Latencia de Escap",
    x = "Prueba",
    y = "Segundos") +
  # theme_classic()
  # ggthemes::theme_clean()
  ggthemes::theme_base()
# ggthemes::theme_few()

lat_violin_post


## ---- patchwork ----

lat_pre <-  (lat_scat_pre / lat_violin_pre) + 
  plot_annotation(
    tag_levels = 'A'
  )

ggsave("./figuras/latencias_scatviolin_pre_otros_grupos.png", plot = lat_pre)


lat_post <-  (lat_scat_post / lat_violin_post) + 
  plot_annotation(
    tag_levels = 'A'
  )

ggsave("./figuras/latencias_scatviolin_post_otros_grupos.png", plot = lat_post)


# Bayes ----

contrasts(latencias_pre$tratamientos) = contr.sum(3) # contrast encoding - represent categorical variable numerically
contrasts(latencias_post$tratamientos) = contr.sum(3) # contrast encoding - represent categorical variable numerically


latencias_pre$tiempo_centrado = latencias_pre$tiempo - mean(latencias_pre$tiempo) # centrar variable de día (0 = mean of variable) coefficients easier, convergence, reduce multicol
latencias_post$tiempo_centrado = latencias_post$tiempo - mean(latencias_post$tiempo) # centrar variable de día (0 = mean of variable) coefficients easier, convergence, reduce multicol

priors = c(
  prior("student_t(3, 4.5, .25)", class = "Intercept"),
  prior("student_t(3, 0, .25)", class = "sd"),
  prior(
    "student_t(3, 0, .1)",
    class = "sd",
    coef = "tiempo_centrado",
    group = "id" ## o sera treatment?
  ),
  prior_string(
    "student_t(3, 0, 0.5)",
    class = "b",
    coef = paste("tratamientos", 1:2, sep = "")
  ),
  prior_string("student_t(3,-.2, .1)", class = "b", coef = "tiempo_centrado"),
  prior_string(
    "student_t(3, 0, 0.1)",
    class = "b",
    coef = paste("tiempo_centrado:tratamientos", 1:2, sep = "")
  )
)

# distribucion gamma
b1.c <-
  brm(
    latencia |
      cens(censurado) ~ tiempo_centrado * tratamientos + (tiempo_centrado |
                                                            id),
    family = Gamma(link = "log"),
    data = latencias_post,
    init = "0",
    iter = 4000,
    prior = priors,
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.95)
  )


b1.ignore <-
  brm(
    latencia ~ tiempo_centrado * tratamientos + (tiempo_centrado |
                                                   id),
    family = Gamma(link = "log"),
    data = latencias_post,
    init = "0",
    iter = 4000,
    prior = priors,
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.95)
  )

b1.c.pre <-
  brm(
    latencia |
      cens(censurado) ~ tiempo_centrado * tratamientos + (tiempo_centrado |
                                                            id),
    family = Gamma(link = "log"),
    data = latencias_pre,
    init = "0",
    iter = 4000,
    prior = priors,
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.95)
  )


b1.ignore.pre <-
  brm(
    latencia ~ tiempo_centrado * tratamientos + (tiempo_centrado |
                                                   id),
    family = Gamma(link = "log"),
    data = latencias_pre,
    init = "0",
    iter = 4000,
    prior = priors,
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.95)
  )

priors = c(
  prior("student_t(3, 90, 45)", class = "Intercept"),
  prior_string(
    "student_t(3, 0, 35)",
    class = "b",
    coef = paste("tratamientos", 1:2, sep = "")
  ),
  prior_string("student_t(3,-5, 5)", class = "b", coef = "tiempo_centrado"),
  prior_string(
    "student_t(3, 0, 5)",
    class = "b",
    coef = paste("tiempo_centrado:tratamientos", 1:2, sep = "")
  )
)

bn <- brm(
  latencia ~ tiempo_centrado * tratamientos + (tiempo_centrado | id),
  data = latencias_post,
  init = "0",
  iter = 4000,
  prior = c(set_prior("student_t(3,-25,25)", class = "b")),
  save_pars = save_pars(all = TRUE)
)


bn.pre <- brm(
  latencia ~ tiempo_centrado * tratamientos + (tiempo_centrado | id),
  data = latencias_pre,
  init = "0",
  iter = 4000,
  prior = c(set_prior("student_t(3,-25,25)", class = "b")),
  save_pars = save_pars(all = TRUE)
)

# ---- Anova ----

## Preparar ----
library(gt)

myData.mean.lat.pre <- aggregate(latencias_pre$latencia,
                                 by = list(latencias_pre$prueba, 
                                           latencias_pre$tratamientos,
                                           latencias_pre$id),
                                 FUN = 'mean')

myData.mean.lat.post <- aggregate(latencias_post$latencia,
                                  by = list(latencias_post$prueba, 
                                            latencias_post$tratamientos,
                                            latencias_post$id),
                                  FUN = 'mean')

colnames(myData.mean.lat.pre) <- c("prueba","tratamientos","id","latencia")

colnames(myData.mean.lat.post) <- c("prueba","tratamientos","id","latencia")


## Modelo -----

aov_latencias_pre <- aov(latencia ~ tratamientos * prueba + Error(id), 
                         data = myData.mean.lat.pre)

aov_latencias_post <- aov(latencia ~ tratamientos * prueba + Error(id), 
                          data = myData.mean.lat.post)

summary(aov_latencias_pre) 

summary(aov_latencias_post) 

anova_summary(aov_latencias_pre) |> 
  as_tibble() |> 
  rename(Efecto = Effect) |> 
  flextable() |> 
  theme_apa()

anova_summary(aov_latencias_post) |> 
  as_tibble() |> 
  rename(Efecto = Effect) |>
  flextable() |> 
  theme_apa()

# lo pasé a Markdown con esta página: https://tabletomarkdown.com/generate-markdown-table/

## Plots ----

aov_lat_pre_plot <- ggstatsplot::grouped_ggwithinstats(
  data = myData.mean.lat.pre,
  x = prueba, 
  y = latencia, 
  grouping.var = tratamientos, 
  type = "parametric",
  # bf.message = F,
  results.subtitle = F,
  xlab = "Entrenamiento",
  p.adjust.method = "none",
  pairwise.display = "all",
  violin.args = list(width = 0, linewidth = 0),
  ggplot.component = scale_y_continuous(breaks = seq(0, 60, 10), limits = c(0, 60))
  # boxplot.args = list(width = 0, linewidth = 0)
)

aov_lat_pre_plot

ggsave("./figuras/aov_lat_pre_plot_otros_grupos.png", plot = aov_lat_pre_plot)

aov_lat_post_plot <- ggstatsplot::grouped_ggwithinstats(
  data = myData.mean.lat.post,
  x = prueba, 
  y = latencia, 
  grouping.var = tratamientos, 
  type = "parametric",
  # bf.message = F,
  results.subtitle = F,
  xlab = "Entrenamiento",
  p.adjust.method = "none",
  pairwise.display = "all",
  violin.args = list(width = 0, linewidth = 0),
  ggplot.component = scale_y_continuous(breaks = seq(0, 60, 10), limits = c(0, 60))
  # boxplot.args = list(width = 0, linewidth = 0)
)
aov_lat_post_plot

ggsave("./figuras/aov_lat_post_plot_otros_grupos.png", plot = aov_lat_post_plot)

# report(aov_latencias_pre)

# report(aov_latencias_post)


### emmeans post ----

latencias_pre_emmeans_aov <- emmeans(aov_latencias_pre, ~ prueba | tratamientos, cov.reduce = F)
latencias_pre_emmeans_aov <- as.data.frame(latencias_pre_emmeans_aov)

latencias_post_emmeans_aov <- emmeans(aov_latencias_post, ~ prueba | tratamientos, cov.reduce = F)
latencias_post_emmeans_aov <- as.data.frame(latencias_post_emmeans_aov)


aov_plot.pre <- ggplot(latencias_pre_emmeans_aov,
                       aes(
                         x = prueba,
                         y = emmean,
                         group = tratamientos,
                         color = tratamientos
                       )) +
  geom_line(linewidth = 0.3, linetype = "dashed", position = position_dodge(0.1)) +
  geom_point(size = 4, position = position_dodge(0.1), aes(shape = tratamientos)) +
  # geom_point(size = 7, shape = 21, position = position_dodge(0.1)) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),
                  size = 0.75, position = position_dodge(0.1)) +
  scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#F2AFEF")) +
  labs(
    caption = "Plotted with SEM",
    color = "Tratamientos",
    shape = "Tratamientos",
    x = "Entrenamiento",
    y = "Segundos") +
  ylim(0, 60) +
  # scale_x_discrete(labels=c("entrenamiento_rev_1"="1", "entrenamiento_rev_2"="2")) +
  # theme_classic()
  # ggthemes::theme_clean() +
  ggthemes::theme_clean() + 
  theme(legend.position='top') 

aov_plot.pre

ggsave("./figuras/aov_plot.pre_otros_grupos.png", plot = aov_plot.pre)

aov_plot.post <- ggplot(latencias_post_emmeans_aov,
                        aes(
                          x = prueba,
                          y = emmean,
                          group = tratamientos,
                          color = tratamientos
                        )) +
  geom_line(linewidth = 0.3, linetype = "dashed", position = position_dodge(0.1)) +
  geom_point(size = 4, position = position_dodge(0.1), aes(shape = tratamientos)) +
  # geom_point(size = 7, shape = 21, position = position_dodge(0.1)) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),
                  size = 0.75, position = position_dodge(0.1)) +
  scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#F2AFEF")) +
  labs(
    caption = "Plotted with SEM",
    color = "Tratamientos",
    shape = "Tratamientos",
    x = "Entrenamiento",
    y = "Segundos") +
  ylim(0, 60) +
  # scale_x_discrete(labels=c("entrenamiento_rev_1"="1", "entrenamiento_rev_2"="2")) +
  # theme_classic()
  # ggthemes::theme_clean() +
  ggthemes::theme_clean() + 
  theme(legend.position='top') 

aov_plot.post

ggsave("./figuras/aov_plot_post_otros_grupos.png", plot = aov_plot.post)

#### contrastes ----

aov_em_post <- emmeans(aov_latencias_post, pairwise ~ tratamientos | prueba, adjust = "tukey")$contrasts |> 
  as_tibble() |> 
  select(-df, -t.ratio, -SE) |>
  rename(diferencia_media_estimada = estimate) |> 
  dplyr::filter(p.value <= 0.05) |>
  separate(contrast, into = c("Tratamiento_1", "Tratamiento_2"), sep = " - ") |> 
  flextable() |> 
  add_footer_lines("Comparaciones múltiples ajustadas con Tukey HSD") |> 
  colformat_double(j = c("diferencia_media_estimada"), digits = 1) |>
  colformat_double(j = c("p.value"), digits = 4) |> 
  theme_vanilla() |> 
  # # add_footer_lines("Tiempo = 0 representa el primer dia, T = 1 al segundo dia \n
  # #                  P ajustada con tasa de descubrimiento falso ( FDR )") |>
  # color(~ Tratamiento_1 == "(Flx-CUMS)", ~ Tratamiento_1, color = "#00ADB5")  |>
  # color(~ Tratamiento_2 == "(Flx-CUMS)", ~ Tratamiento_2, color = "#00ADB5")  |>
  # color(~ Tratamiento_1 == "Flx", ~ Tratamiento_1, color = "#222831")  |>
  # color(~ Tratamiento_2 == "Flx", ~ Tratamiento_2, color = "#222831")  |>
  # color(~ Tratamiento_1 == "(Sal-CUMS-F)", ~ Tratamiento_1, color = "#FF2E63")  |>
  # color(~ Tratamiento_2 == "(Sal-CUMS-F)", ~ Tratamiento_2, color = "#FF2E63")  |>
  colformat_double(j = c("diferencia_media_estimada", "p.value"), digits = 4)


aov_em_post


#### Agregar contrastes a plot

# colours <- scales::viridis_pal(option = "viridis")(10)

# grad_ungroup <- linearGradient(colours, group = FALSE)
# 
# aov_plot.post_comps <-  aov_plot.post + annotate(
#   "text", label = "***", x = "Reversa_1", y = 50, size = 6.7, colour =  "#222831"
# ) +
#   annotate(
#     "text", label = "***", x = "Reversa_1", y = 53, size = 6.7, colour =  "#00ADB5"
#   )

# aov_plot.post_comps

#### Guardar ----

# anova_resultados_post <- aov_plot.post_comps /
#   gen_grob(aov_em_post, fit = "width", just = "centre")


ggsave("./figuras/aov_plot.post_comps_otros_grupos.png", plot = aov_plot.post)
aov_em_post |> save_as_image(path = "./tablas/aov_em_post_otros_grupos.png")


### emmeans pre ----

latencias_pre_emmeans_aov <- emmeans(aov_latencias_pre, ~ prueba | tratamientos, cov.reduce = F)
latencias_pre_emmeans_aov_df <- as.data.frame(latencias_pre_emmeans_aov)

# latencias_post_emmeans_aov <- emmeans(aov_latencias_post, ~ prueba | tratamientos, cov.reduce = F)
# latencias_post_emmeans_aov <- as.data.frame(latencias_post_emmeans_aov)


aov_plot.pre <- ggplot(latencias_pre_emmeans_aov_df,
                       aes(
                         x = prueba,
                         y = emmean,
                         group = tratamientos,
                         color = tratamientos
                       )) +
  geom_line(linewidth = 0.3, linetype = "dashed", position = position_dodge(0.1)) +
  geom_point(size = 4, position = position_dodge(0.1), aes(shape = tratamientos)) +
  # geom_point(size = 7, shape = 21, position = position_dodge(0.1)) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),
                  size = 0.75, position = position_dodge(0.1)) +
  scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#F2AFEF")) +
  labs(
    caption = "Plotted with SEM",
    color = "Tratamientos",
    shape = "Tratamientos",
    x = "Entrenamiento",
    y = "Segundos") +
  ylim(0, 60) +
  scale_x_discrete(labels=c("entrenamiento_rev_1"="1", "entrenamiento_rev_2"="2")) +
  # theme_classic()
  # ggthemes::theme_clean() +
  ggthemes::theme_clean() + 
  theme(legend.position='top') 

aov_plot.pre
ggsave("./figuras/aov_plot.pre_otros_grupos.png", plot = aov_plot.pre)

#### contrastes ----

aov_em_pre <- emmeans(latencias_pre_emmeans_aov, pairwise ~ tratamientos | prueba, adjust = "tukey")$contrasts |> 
  as_tibble() |> 
  select(-df, -t.ratio, -SE) |>
  rename(diferencia_media_estimada = estimate) |> 
  dplyr::filter(p.value <= 0.05) |>
  separate(contrast, into = c("Tratamiento_1", "Tratamiento_2"), sep = " - ") |> 
  flextable() |> 
  add_footer_lines("Comparaciones múltiples ajustadas con Tukey HSD") |> 
  colformat_double(j = c("diferencia_media_estimada"), digits = 1) |>
  colformat_double(j = c("p.value"), digits = 4) |> 
  theme_vanilla() |> 
  # add_footer_lines("Tiempo = 0 representa el primer dia, T = 1 al segundo dia \n
  #                  P ajustada con tasa de descubrimiento falso ( FDR )") |>
  color(~ Tratamiento_1 == "(Flx-CUMS)", ~ Tratamiento_1, color = "#00ADB5")  |>
  color(~ Tratamiento_2 == "(Flx-CUMS)", ~ Tratamiento_2, color = "#00ADB5")  |>
  color(~ Tratamiento_1 == "Flx", ~ Tratamiento_1, color = "#222831")  |>
  color(~ Tratamiento_2 == "Flx", ~ Tratamiento_2, color = "#222831")  |>
  color(~ Tratamiento_1 == "(Sal-CUMS-F)", ~ Tratamiento_1, color = "#FF2E63")  |>
  color(~ Tratamiento_2 == "(Sal-CUMS-F)", ~ Tratamiento_2, color = "#FF2E63")  |>
  colformat_double(j = c("diferencia_media_estimada", "p.value"), digits = 4)


aov_em_pre


#### Agregar contrastes a plot

# colours <- scales::viridis_pal(option = "viridis")(10)

# grad_ungroup <- linearGradient(colours, group = FALSE)

# aov_plot.post_comps <-  aov_plot.post + annotate(
  # "text", label = "***", x = "Reversa_1", y = 50, size = 6.7, colour =  "#222831"
# ) +
  # annotate(
    # "text", label = "***", x = "Reversa_1", y = 53, size = 6.7, colour =  "#00ADB5"
  # )

# aov_plot.post_comps

#### Guardar ----

# anova_resultados_post <- aov_plot.post_comps /
#   gen_grob(aov_em_post, fit = "width", just = "centre")


# ggsave("./figuras/aov_plot.pre_comps_otros_grupos.png", plot = aov_plot.pre)
aov_em_pre |> save_as_image(path = "./tablas/aov_em_pre_otros_grupos.png")

### supuestos ----


# ---- lmer ----

lmer_latencias_post <- lmer(
  latencia ~ tiempo * tratamientos + (1 + tiempo | id),
  data = latencias_post
)

lmer_latencias_post.2 <- lmer(
  latencia ~ tiempo * tratamientos + (1 | id),
  data = latencias_post
)

# equatiomatic::extract_eq(lmer_latencias_post.2)

# model_performance(lmer_latencias_post)
# model_performance(lmer_latencias_post.2)


# report(lmer_latencias_post)
# report(lmer_latencias_post.2)

sjPlot::tab_model(lmer_latencias_post.2,
                  # auto.label = TRUE,
                  show.se = T,
                  # show.obs = T,
                  # show.fstat = T,
                  show.reflvl = F, # nivel de referencia para factores
                  # show.intercept = F,
                  # show.df = T,
                  p.style = "numeric_stars",
                  title = "Latencias Reversa Efectos Mixtos",
                  string.pred = "Predictores",
                  string.est = "Estimados",
                  # pred.labels = c("Intercepto",
                  #                 "Tiempo",
                  #                 "Flx-CUMS", "Sal-CUMS-F",
                  #                 "Tiempo:Flx-CUMS", "Tiempo:Sal-CUMS-F"),
                  # dv.labels = c("Efectos Mixtos lineal"),
                  # file = "./tablas/lmer_latencias_post_otros_grupos.html",
                  string.se = "SEM"
)

# webshot("./tablas/lmer_latencias_post.html", "./tablas/lmer_latencias_post.png")


lmer.latencias_post_emmeans_aov <- emmeans(lmer_latencias_post.2, ~ tiempo | tratamientos, cov.reduce = FALSE)

lmer_latencias_post_emmeans_aov_df <- as.data.frame(lmer.latencias_post_emmeans_aov)


lmer_plot <- ggplot(lmer_latencias_post_emmeans_aov_df,
                    aes(
                      x = tiempo,
                      y = emmean,
                      group = tratamientos,
                      color = tratamientos
                    )) +
  geom_line(linewidth = 0.3, linetype = "dashed", position = position_dodge(0.1)) +
  geom_point(size = 4, position = position_dodge(0.1), aes(shape = tratamientos)) +
  # geom_point(size = 7, shape = 21, position = position_dodge(0.1)) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),
                  size = 0.75, position = position_dodge(0.1)) +
  scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
  labs(
    caption = "Plotted with SEM",
    color = "Tratamientos",
    shape = "Tratamientos",
    x = "Tiempo",
    y = "Segundos") +
  ylim(0, 60) +
  # scale_x_discrete(labels=c("entrenamiento_rev_1"="1", "entrenamiento_rev_2"="2")) +
  # theme_classic()
  # ggthemes::theme_clean() +
  ggthemes::theme_clean() + 
  theme(legend.position='top', axis.text = element_text(size = 12),axis.title = element_text(size = 15)) 

lmer_plot

# ggsave("./figuras/lmer_lat_post.png", plot = lmer_plot)


lmer_em <- emmeans(lmer_latencias_post.2, pairwise ~ tratamientos | tiempo, adjust = "tukey",  cov.reduce = FALSE)$contrasts |> 
  as_tibble() |> 
  select(-df, -t.ratio, -SE) |>
  rename(diferencia_media_estimada = estimate) |> 
  dplyr::filter(p.value <= 0.05) |>
  separate(contrast, into = c("Tratamiento_1", "Tratamiento_2"), sep = " - ") |> 
  flextable() |> 
  add_footer_lines("Comparaciones múltiples ajustadas con Tukey HSD") |> 
  colformat_double(j = c("diferencia_media_estimada"), digits = 1) |>
  colformat_double(j = c("p.value"), digits = 4) |> 
  theme_vanilla() |> 
  # add_footer_lines("Tiempo = 0 representa el primer dia, T = 1 al segundo dia \n
  #                  P ajustada con tasa de descubrimiento falso ( FDR )") |>
  # color(~ Tratamiento_1 == "(Flx-CUMS)", ~ Tratamiento_1, color = "#00ADB5")  |>
  # color(~ Tratamiento_2 == "(Flx-CUMS)", ~ Tratamiento_2, color = "#00ADB5")  |>
  # color(~ Tratamiento_1 == "Flx", ~ Tratamiento_1, color = "#222831")  |>
  # color(~ Tratamiento_2 == "Flx", ~ Tratamiento_2, color = "#222831")  |>
  # color(~ Tratamiento_1 == "(Sal-CUMS-F)", ~ Tratamiento_1, color = "#FF2E63")  |>
  # color(~ Tratamiento_2 == "(Sal-CUMS-F)", ~ Tratamiento_2, color = "#FF2E63")  |>
  colformat_double(j = c("diferencia_media_estimada", "p.value"), digits = 4)


lmer_em

# Bayes_2 ----

## tab model ----

sjPlot::tab_model(b1.c,b1.ignore, bn,
                  # auto.label = TRUE,
                  # show.se = T,
                  # show.obs = T,
                  # show.fstat = T,
                  # show.reflvl = T, # nivel de referencia para factores
                  # show.intercept = F,
                  # show.df = T,
                  p.style = "numeric_stars",
                  title = "Latencias Reversa",
                  string.pred = "Predictores",
                  string.est = "Estimados",
                  # pred.labels = c("Intercepto",
                  #                 "Tiempo",
                  #                 "Flx-CUM", "Sal-CUM-F", 
                  #                 "Tiempo:Flx-CUM", "Tiempo:Sal-CUM-F"),
                  # dv.labels = c("Censurado", "GLM", "Efectos Mixtos lineal"),
                  # file = "./tablas/bayes_latencias_post.html"
                  # string.se = "SEM"
)

# webshot("./tablas/bayes_latencias_post.html", "./tablas/bayes_latencias_post.png")

## plot model ----

bayes_bpe_post_latencias <- plot_model(b1.c,
                                       vline.color = "red", # no effect
                                       bpe = "mean", # mean point estimate
                                       bpe.style = "dot",
                                       colors = "Dark2",
                                       show.values = TRUE,
                                       show.legend = T,
                                       # axis.labels =  rev(c("Tiempo",
                                                            # "Flx-CUMS", "Sal-CUMS-F",
                                                            # "Tiempo:Flx-CUMS", "Tiempo:Sal-CUMS-F"
                                       # )
                                       # ),
                                       transform = NULL # exponentiates coefficients, if appropriate (e.g. for models with log or logit link
) +
  ggthemes::theme_clean() +
  ggplot2::ggtitle("")

ggsave("./figuras/bayes_bpe_post_latencias.png", plot = bayes_bpe_post_latencias)

bayes_bpe_post_latencias_lmer <- plot_model(bn,
                                            vline.color = "red", # no effect
                                            bpe = "mean", # mean point estimate
                                            bpe.style = "dot",
                                            colors = "Dark2",
                                            show.values = TRUE,
                                            show.legend = T,
                                            axis.labels =  rev(c("Tiempo",
                                                                 "Flx-CUMS", "Sal-CUMS-F",
                                                                 "Tiempo:Flx-CUMS", "Tiempo:Sal-CUMS-F"
                                            )
                                            ),
                                            transform = NULL # exponentiates coefficients, if appropriate (e.g. for models with log or logit link
) +
  ggthemes::theme_clean() +
  ggplot2::ggtitle("")

# ggsave("./figuras/bayes_bpe_post_latencias_lmer.png", plot = bayes_bpe_post_latencias_lmer)


# table(latencias_post$tiempo_centrado)
# rstantools::posterior_predict(b1.c)

bayes_predict_post_latencias <- plot_model(b1.c, type = "pred", terms = c("tratamientos", "tiempo_centrado[-0.875, 0.875]")) +
  ggthemes::theme_clean() +
  ggplot2::ggtitle("") +
  ggplot2::theme(legend.position="none", axis.text = element_text(size = 12),axis.title = element_text(size = 15))

# ggsave("./figuras/bayes_predict_post_latencias.png", plot = bayes_predict_post_latencias)


bayes_glm_predict_post_latencias <- plot_model(b1.ignore, type = "pred", terms = c("tratamientos", "tiempo_centrado[-0.875, 0.875]")) +
  ggthemes::theme_clean() +
  ggplot2::ggtitle("") +
  ggplot2::theme(legend.position="none", axis.text = element_text(size = 12),axis.title = element_text(size = 15))

ggsave("./figuras/bayes_glm_predict_post_latencias.png", plot = bayes_glm_predict_post_latencias)

bayes_lmer_predict_post_latencias <- plot_model(bn, type = "pred", terms = c("tratamientos", "tiempo_centrado[-0.875, 0.875]")) +
  ggthemes::theme_clean() +
  ggplot2::ggtitle("") +
  ggplot2::theme(legend.position="none", axis.text = element_text(size = 12),axis.title = element_text(size = 15))

ggsave("./figuras/bayes_lmer_predict_post_latencias.png", plot = bayes_lmer_predict_post_latencias)

# plot_model(b1.c, type = "int")

## report y eq ----

# report::report(b1.c) # NOTA - tarda mucho

# equatiomatic::extract_eq(b1.c)

## emmeans otros extras ----

emmeans(b1.c, ~tratamientos, type="response", at=list(c.dia.comp=-0.8723881)) |> 
  as.tibble() |> 
  flextable() 

emmeans(b1.weib, ~tratamientos, type="response", at=list(c.dia.comp=-0.8723881)) |> 
  as.tibble() |> 
  flextable() 

emmeans(b1.ignore, ~tratamientos, type="response", at=list(c.dia.comp=-0.8723881)) |>
  as.tibble() |> 
  flextable() 

emmeans(bn, ~tratamientos, at=list(c.dia.comp=-0.8723881)) |> 
  as.tibble() |> 
  flextable() 




emmeans(b1.c.pre, ~tratamientos, type="response", at=list(c.dia.comp=-0.8723881)) |> 
  as.tibble() |> 
  flextable() 

emmeans(b1.weib.pre, ~tratamientos, type="response", at=list(c.dia.comp=-0.8723881)) |> 
  as.tibble() |> 
  flextable() 

emmeans(b1.ignore.pre, ~tratamientos, type="response", at=list(c.dia.comp=-0.8723881)) |>
  as.tibble() |> 
  flextable() 

emmeans(bn.pre, ~tratamientos, at=list(c.dia.comp=-0.8723881)) |> 
  as.tibble() |> 
  flextable() 

### slopes ----

# slopes
emtrends(b1.c, ~tratamientos, var="c.dia.comp") |> 
  as.tibble() |> 
  flextable() 


emtrends(b1.weib, ~tratamientos, var="c.dia.comp") |> 
  as.tibble() |> 
  flextable() 


emtrends(b1.ignore, ~tratamientos, var="c.dia.comp") |> 
  as.tibble() |> 
  flextable() 


emtrends(bn, ~tratamientos, var="c.dia.comp") |> 
  as.tibble() |> 
  flextable() 

# slopes
emtrends(b1.c.pre, ~tratamientos, var="c.dia.comp") |> 
  as.tibble() |> 
  flextable() 


emtrends(b1.weib.pre, ~tratamientos, var="c.dia.comp") |> 
  as.tibble() |> 
  flextable() 


emtrends(b1.ignore.pre, ~tratamientos, var="c.dia.comp") |> 
  as.tibble() |> 
  flextable() 


emtrends(bn.pre, ~tratamientos, var="c.dia.comp") |> 
  as.tibble() |> 
  flextable() 



# To plot the emmeans above at a particular value (e.g., dia 0). Can add pairwise, if desired

plot(emmeans(b1.c, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept", type="response")

plot(emmeans(b1.weib, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept", type="response")

plot(emmeans(b1.ignore, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept", type="response")

plot(emmeans(bn, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept")


plot(emmeans(b1.c.pre, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept", type="response")

plot(emmeans(b1.weib.pre, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept", type="response")

plot(emmeans(b1.ignore.pre, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept", type="response")

plot(emmeans(bn.pre, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept")



# The following are for average dia, not dia 0
plot(emmeans(b1.c, ~tratamientos), xlab="Intercept", type="response")

plot(emmeans(b1.weib, ~tratamientos), xlab="Intercept", type="response")

plot(emmeans(b1.ignore, ~tratamientos), xlab="Intercept", type="response")

plot(emmeans(bn, ~tratamientos), xlab="Intercept", type="response")


# For slopes, marginality doesn't matter here because I'm breaking them out across tratamientoss.
# Slope doesn't change as a function on any other variable. 
plot(emtrends(b1.c, ~tratamientos, var="c.dia.comp"), xlab="Slope in Log Scale")
plot(emtrends(b1.weib, ~tratamientos, var="c.dia.comp"), xlab="Slope in Log Scale")
plot(emtrends(b1.ignore, ~tratamientos, var="c.dia.comp"), xlab="Slope in Log Scale")
plot(emtrends(bn, ~tratamientos, var="c.dia.comp"), xlab="Slope in Log Scale")




## plots juntos ----

##  Another way to show results all on the same plot - not shown in manuscript

b1.c_effects_plot <- plot(conditional_effects(b1.c, effects = "tiempo_centrado:tratamientos", prob=.66))[[1]]+
  ggplot2::ylim(0,180) +
  theme_classic() +
  scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
  scale_fill_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
  labs(
    # caption = "Plotted with 66% Credible Intervals",
    color = "Tratamientos",
    fill = "Tratamientos",
    shape = "Tratamientos",
    x = "Tiempo",
    y = "Segundos") +
  theme(legend.position='top', axis.text = element_text(size = 12),axis.title = element_text(size = 15)) 

# b1.weib_effects_plot <- plot(conditional_effects(b1.weib, effects = "tiempo:tratamientos", prob=.66))[[1]]+
# ggplot2::ylim(-10,180) +
#   theme_classic() +
#   scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
#   scale_fill_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
#   labs(
#     # caption = "Plotted with 66% Credible Intervals",
#     color = "Tratamientos",
#     fill = "Tratamientos",
#     shape = "Tratamientos",
#     x = "Tiempo",
#     y = "Segundos") +
#   theme(legend.position='top')

b1.ignore_effects_plot <- plot(conditional_effects(b1.ignore, effects = "tiempo_centrado:tratamientos", prob=.66))[[1]]+
  ggplot2::ylim(0,180)+
  theme_classic() +
  scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
  scale_fill_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
  labs(
    # caption = "Plotted with 66% Credible Intervals",
    color = "Tratamientos",
    fill = "Tratamientos",
    shape = "Tratamientos",
    x = "Tiempo",
    y = "Segundos") +
  theme(legend.position='top', axis.text = element_text(size = 12),axis.title = element_text(size = 15)) 

bn_effects_plot <- plot(conditional_effects(bn, effects = "tiempo_centrado:tratamientos", prob=.66))[[1]]+
  ggplot2::ylim(0,180)+
  theme_classic() +
  scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
  scale_fill_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
  labs(
    # caption = "Plotted with 66% Credible Intervals",
    color = "Tratamientos",
    fill = "Tratamientos",
    shape = "Tratamientos",
    x = "Tiempo",
    y = "Segundos") +
  theme(legend.position='top', axis.text = element_text(size = 12),axis.title = element_text(size = 15)) 


bay_effect_plots <- b1.c_effects_plot + b1.ignore_effects_plot + bn_effects_plot + 
  plot_annotation(tag_levels = 'A') + 
  plot_layout(guides = 'collect', axis_titles = 'collect') &
  theme(legend.position = 'top')

bay_effect_plots

# ggsave("./figuras/bay_effect_plots_post.png", plot = bay_effect_plots)

## Bayes otrs grafs posteriors ----

# summary of posterior distribution
library(bayestestR)
library(insight)

posteriors.b1.c <- describe_posterior(b1.c)
print_md(posteriors.b1.c, digits = 2)

# visualize posterior distribution
posteriors.model.b1.c <- get_parameters(b1.c)

ggplot(posteriors.model.b1.c, aes(x = b_Intercept)) +
  geom_density(fill = "orange") +
  ggthemes::theme_base()

# describe 3 elements of posteriors

# point estimate (one-value summary, similar to beta)
# what single value can represent best my post dist?
# mean, median, mode
mean(posteriors.model.b1.c$b_Intercept)
median(posteriors.model.b1.c$b_Intercept)
map_estimate(posteriors.model.b1.c$b_Intercept) # mode - peak - max A posteriori

ggplot(posteriors.model.b1.c, aes(x = b_Intercept)) +
  geom_density(fill = "orange") +
  # The mean in blue
  geom_vline(xintercept = mean(posteriors.model.b1.c$b_Intercept), color = "blue", linewidth = 1) +
  # The median in red
  geom_vline(xintercept = median(posteriors.model.b1.c$b_Intercept), color = "red", linewidth = 1) +
  # The MAP in purple
  geom_vline(xintercept = as.numeric(map_estimate(posteriors.model.b1.c$b_Intercept)), color = "purple", linewidth = 1)

# como están muy cerca, vamos a usar la mediana
# mediana tiene significado probabilistico directo (50% que efect real esté arriba/abajo)

# credible interval (associated uncertainty)
range(posteriors.model.b1.c$b_Intercept)
# en lugar de incluir todos los valores, vamos a calcular intervalos de cred
# similar a intervalos de confianza 
# los vamos a calcular con highest Density Interval
# give us the range containing the 89% most probable effect values
# 89% level gives more stable results (Kruschke, 2014)
hdi(posteriors.model.b1.c$b_Intercept, ci = 0.89)

# indices of significance (relative importance of this effect)
# si toda la distribucion posterior no está en cero, puede ser evidencia
# para determinar si efecto es sign, podemos ver si el int cre contiene al 0
# si no, puede ser evidencia que nuestro efecto es sign

# se puede hacer con rope rope(posteriors$feedsunflower, range = c(-20, 20), ci = 0.89)

## ROPES ----

library(modelbased)
library(see)

vizdata <- estimate_relation(b1.c)

plot_predicted <- ggplot(vizdata, aes(x = tiempo_centrado, y = Predicted)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.5) +
  geom_line() +
  facet_grid(~tratamientos) +
  theme_modern() +
  # ggthemes::theme_base() +
  # scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9")) +
  # scale_fill_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9")) +
  labs(
    # caption = "Plotted with 66% Credible Intervals",
    # color = "Treatments",
    # fill = "Treatments",
    # shape = "Treatments",
    x = "Tiempo",
    y = "Probabilidad en segundos de Latencia")
# theme(legend.position='top') 

plot_predicted

# ggsave("./figuras/plot_predicted.png", plot = plot_predicted)

plot_rope <- plot(rope(b1.c, ci = .95)) +
  ggthemes::theme_fivethirtyeight()  +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 15)) 

plot_rope

ggsave("./figuras/plot_rope_post_bayes.png", plot = plot_rope)

plot_rope_bn <- plot(rope(bn, ci = .95)) +
  ggthemes::theme_fivethirtyeight()  +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 15)) 

plot_rope_bn

ggsave("./figuras/plot_rope_post_lmer.png", plot = plot_rope_bn)


# PRE ENT ORIGINALES ----

# ---- lmer ----

# lmer_latencias_pre <- lmer(
#   latencia ~ tiempo * tratamientos + (1 + tiempo | id),
#   data = latencias_pre
# )

lmer_latencias_pre.2 <- lmer(
  latencia ~ tiempo * tratamientos + (1 | id),
  data = latencias_pre
)

equatiomatic::extract_eq(lmer_latencias_pre.2)

# model_performance(lmer_latencias_pre)
# model_performance(lmer_latencias_pre.2)


# report(lmer_latencias_pre)
report(lmer_latencias_pre.2)

sjPlot::tab_model(lmer_latencias_pre.2,
                  # auto.label = TRUE,
                  show.se = T,
                  # show.obs = T,
                  # show.fstat = T,
                  show.reflvl = F, # nivel de referencia para factores
                  # show.intercept = F,
                  # show.df = T,
                  p.style = "numeric_stars",
                  title = "Latencias Originales Efectos Mixtos",
                  string.pred = "Predictores",
                  string.est = "Estimados",
                  # pred.labels = c("Intercepto",
                                  # "Tiempo",
                                  # "Flx-CUMS", "Sal-CUMS-F",
                                  # "Tiempo:Flx-CUMS", "Tiempo:Sal-CUMS-F"),
                  # dv.labels = c("Efectos Mixtos lineal"),
                  # file = "./tablas/lmer_latencias_pre_otros_grupos.html",
                  string.se = "SEM"
)

# webshot("./tablas/lmer_latencias_pre.html", "./tablas/lmer_latencias_pre.png")


lmer.latencias_pre_emmeans_aov <- emmeans(lmer_latencias_pre.2, ~ tiempo | tratamientos, cov.reduce = FALSE)

lmer_latencias_pre_emmeans_aov_df <- as.data.frame(lmer.latencias_pre_emmeans_aov)


lmer_plot <- ggplot(lmer_latencias_pre_emmeans_aov_df,
                    aes(
                      x = tiempo,
                      y = emmean,
                      group = tratamientos,
                      color = tratamientos
                    )) +
  geom_line(linewidth = 0.3, linetype = "dashed", position = position_dodge(0.1)) +
  geom_point(size = 4, position = position_dodge(0.1), aes(shape = tratamientos)) +
  # geom_point(size = 7, shape = 21, position = position_dodge(0.1)) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),
                  size = 0.75, position = position_dodge(0.1)) +
  scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
  labs(
    caption = "Plotted with SEM",
    color = "Tratamientos",
    shape = "Tratamientos",
    x = "Tiempo",
    y = "Segundos") +
  ylim(0, 60) +
  # scale_x_discrete(labels=c("entrenamiento_rev_1"="1", "entrenamiento_rev_2"="2")) +
  # theme_classic()
  # ggthemes::theme_clean() +
  ggthemes::theme_clean() + 
  theme(legend.position='top', axis.text = element_text(size = 12),axis.title = element_text(size = 15)) 

lmer_plot

# ggsave("./figuras/lmer_lat_pre_otros_grupos.png", plot = lmer_plot)


lmer_em <- emmeans(lmer_latencias_pre.2, pairwise ~ tratamientos | tiempo, adjust = "tukey",  cov.reduce = FALSE)$contrasts |> 
  as_tibble() |> 
  select(-df, -t.ratio, -SE) |>
  rename(diferencia_media_estimada = estimate) |> 
  dplyr::filter(p.value <= 0.05) |>
  separate(contrast, into = c("Tratamiento_1", "Tratamiento_2"), sep = " - ") |> 
  flextable() |> 
  add_footer_lines("Comparaciones múltiples ajustadas con Tukey HSD") |> 
  colformat_double(j = c("diferencia_media_estimada"), digits = 1) |>
  colformat_double(j = c("p.value"), digits = 4) |> 
  theme_vanilla() |> 
  # add_footer_lines("Tiempo = 0 representa el primer dia, T = 1 al segundo dia \n
  #                  P ajustada con tasa de descubrimiento falso ( FDR )") |>
  # color(~ Tratamiento_1 == "(Flx-CUMS)", ~ Tratamiento_1, color = "#00ADB5")  |>
  # color(~ Tratamiento_2 == "(Flx-CUMS)", ~ Tratamiento_2, color = "#00ADB5")  |>
  # color(~ Tratamiento_1 == "Flx", ~ Tratamiento_1, color = "#222831")  |>
  # color(~ Tratamiento_2 == "Flx", ~ Tratamiento_2, color = "#222831")  |>
  # color(~ Tratamiento_1 == "(Sal-CUMS-F)", ~ Tratamiento_1, color = "#FF2E63")  |>
  # color(~ Tratamiento_2 == "(Sal-CUMS-F)", ~ Tratamiento_2, color = "#FF2E63")  |>
  colformat_double(j = c("diferencia_media_estimada", "p.value"), digits = 4)


# lmer_em |> save_as_image("./tablas/lmer_contrastes_pre.png")

# Bayes_2 ----

## tab model ----

# sjPlot::tab_model(b1.c.pre,b1.ignore.pre, bn.pre,
#                   # auto.label = TRUE,
#                   # show.se = T,
#                   # show.obs = T,
#                   # show.fstat = T,
#                   # show.reflvl = T, # nivel de referencia para factores
#                   # show.intercept = F,
#                   # show.df = T,
#                   p.style = "numeric_stars",
#                   title = "Latencias Reversa",
#                   string.pred = "Predictores",
#                   string.est = "Estimados",
#                   pred.labels = c("Intercepto",
#                                   "Tiempo",
#                                   "Flx-CUM", "Sal-CUM-F", 
#                                   "Tiempo:Flx-CUM", "Tiempo:Sal-CUM-F"),
#                   dv.labels = c("Censurado", "GLM", "Efectos Mixtos lineal"),
#                   file = "./tablas/bayes_latencias_pre.html",
#                   # string.se = "SEM"
# )

# webshot("./tablas/bayes_latencias_pre.html", "./tablas/bayes_latencias_pre.png")

## plot model ----

# bayes_bpe_pre_latencias <- plot_model(b1.c.pre,
#                                        vline.color = "red", # no effect
#                                        bpe = "mean", # mean point estimate
#                                        bpe.style = "dot",
#                                        colors = "Dark2",
#                                        show.values = TRUE,
#                                        show.legend = T,
#                                        axis.labels =  rev(c("Tiempo",
#                                                             "Flx-CUMS", "Sal-CUMS-F",
#                                                             "Tiempo:Flx-CUMS", "Tiempo:Sal-CUMS-F"
#                                        )
#                                        ),
#                                        transform = NULL # exponentiates coefficients, if appropriate (e.g. for models with log or logit link
# ) +
#   ggthemes::theme_clean() +
#   ggplot2::ggtitle("")
# 
# ggsave("./figuras/bayes_bpe_pre_latencias.png", plot = bayes_bpe_pre_latencias)
# 
# bayes_bpe_pre_latencias_lmer <- plot_model(bn.pre,
#                                             vline.color = "red", # no effect
#                                             bpe = "mean", # mean point estimate
#                                             bpe.style = "dot",
#                                             colors = "Dark2",
#                                             show.values = TRUE,
#                                             show.legend = T,
#                                             axis.labels =  rev(c("Tiempo",
#                                                                  "Flx-CUMS", "Sal-CUMS-F",
#                                                                  "Tiempo:Flx-CUMS", "Tiempo:Sal-CUMS-F"
#                                             )
#                                             ),
#                                             transform = NULL # exponentiates coefficients, if appropriate (e.g. for models with log or logit link
# ) +
#   ggthemes::theme_clean() +
#   ggplot2::ggtitle("")

# ggsave("./figuras/bayes_bpe_pre_latencias_lmer.png", plot = bayes_bpe_pre_latencias_lmer)


# table(latencias_pre$tiempo_centrado)
# rstantools::preerior_predict(b1.c.pre)

# bayes_predict_pre_latencias <- plot_model(b1.c.pre, type = "pred", terms = c("tratamientos", "tiempo_centrado[-0.875, 0.875]")) +
#   ggthemes::theme_clean() +
#   ggplot2::ggtitle("") +
#   ggplot2::theme(legend.position="none", axis.text = element_text(size = 12),axis.title = element_text(size = 15))
# 
# # ggsave("./figuras/bayes_predict_pre_latencias.png", plot = bayes_predict_pre_latencias)
# 
# 
# bayes_glm_predict_pre_latencias <- plot_model(b1.ignore.pre, type = "pred", terms = c("tratamientos", "tiempo_centrado[-0.875, 0.875]")) +
#   ggthemes::theme_clean() +
#   ggplot2::ggtitle("") +
#   ggplot2::theme(legend.position="none", axis.text = element_text(size = 12),axis.title = element_text(size = 15))
# 
# ggsave("./figuras/bayes_glm_predict_pre_latencias.png", plot = bayes_glm_predict_pre_latencias)
# 
# bayes_lmer_predict_pre_latencias <- plot_model(bn.pre, type = "pred", terms = c("tratamientos", "tiempo_centrado[-0.875, 0.875]")) +
#   ggthemes::theme_clean() +
#   ggplot2::ggtitle("") +
#   ggplot2::theme(legend.position="none", axis.text = element_text(size = 12),axis.title = element_text(size = 15))
# 
# ggsave("./figuras/bayes_lmer_predict_pre_latencias.png", plot = bayes_lmer_predict_pre_latencias)

# plot_model(b1.c.pre, type = "int")

## report y eq ----

# report::report(b1.c.pre) # NOTA - tarda mucho

# equatiomatic::extract_eq(b1.c.pre)

## emmeans otros extras ----

# emmeans(b1.c.pre, ~tratamientos, type="response", at=list(c.dia.comp=-0.8723881)) |> 
#   as.tibble() |> 
#   flextable() 
# 
# emmeans(b1.weib, ~tratamientos, type="response", at=list(c.dia.comp=-0.8723881)) |> 
#   as.tibble() |> 
#   flextable() 
# 
# emmeans(b1.ignore.pre, ~tratamientos, type="response", at=list(c.dia.comp=-0.8723881)) |>
#   as.tibble() |> 
#   flextable() 
# 
# emmeans(bn.pre, ~tratamientos, at=list(c.dia.comp=-0.8723881)) |> 
#   as.tibble() |> 
#   flextable() 




# emmeans(b1.c.pre.pre, ~tratamientos, type="response", at=list(c.dia.comp=-0.8723881)) |> 
#   as.tibble() |> 
#   flextable() 
# 
# emmeans(b1.weib.pre, ~tratamientos, type="response", at=list(c.dia.comp=-0.8723881)) |> 
#   as.tibble() |> 
#   flextable() 
# 
# emmeans(b1.ignore.pre.pre, ~tratamientos, type="response", at=list(c.dia.comp=-0.8723881)) |>
#   as.tibble() |> 
#   flextable() 
# 
# emmeans(bn.pre.pre, ~tratamientos, at=list(c.dia.comp=-0.8723881)) |> 
#   as.tibble() |> 
#   flextable() 

### slopes ----

# slopes
# emtrends(b1.c.pre, ~tratamientos, var="c.dia.comp") |> 
#   as.tibble() |> 
#   flextable() 
# 
# 
# emtrends(b1.weib, ~tratamientos, var="c.dia.comp") |> 
#   as.tibble() |> 
#   flextable() 
# 
# 
# emtrends(b1.ignore.pre, ~tratamientos, var="c.dia.comp") |> 
#   as.tibble() |> 
#   flextable() 
# 
# 
# emtrends(bn.pre, ~tratamientos, var="c.dia.comp") |> 
#   as.tibble() |> 
#   flextable() 
# 
# # slopes
# emtrends(b1.c.pre.pre, ~tratamientos, var="c.dia.comp") |> 
#   as.tibble() |> 
#   flextable() 
# 
# 
# emtrends(b1.weib.pre, ~tratamientos, var="c.dia.comp") |> 
#   as.tibble() |> 
#   flextable() 
# 
# 
# emtrends(b1.ignore.pre.pre, ~tratamientos, var="c.dia.comp") |> 
#   as.tibble() |> 
#   flextable() 
# 
# 
# emtrends(bn.pre.pre, ~tratamientos, var="c.dia.comp") |> 
#   as.tibble() |> 
#   flextable() 



# To plot the emmeans above at a particular value (e.g., dia 0). Can add pairwise, if desired

# plot(emmeans(b1.c.pre, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept", type="response")
# 
# plot(emmeans(b1.weib, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept", type="response")
# 
# plot(emmeans(b1.ignore.pre, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept", type="response")
# 
# plot(emmeans(bn.pre, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept")
# 
# 
# plot(emmeans(b1.c.pre.pre, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept", type="response")
# 
# plot(emmeans(b1.weib.pre, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept", type="response")
# 
# plot(emmeans(b1.ignore.pre.pre, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept", type="response")
# 
# plot(emmeans(bn.pre.pre, ~tratamientos, at=c(c.dia.comp=-0.8723881)), xlab="Intercept")



# The following are for average dia, not dia 0
# plot(emmeans(b1.c.pre, ~tratamientos), xlab="Intercept", type="response")
# 
# plot(emmeans(b1.weib, ~tratamientos), xlab="Intercept", type="response")
# 
# plot(emmeans(b1.ignore.pre, ~tratamientos), xlab="Intercept", type="response")
# 
# plot(emmeans(bn.pre, ~tratamientos), xlab="Intercept", type="response")


# For slopes, marginality doesn't matter here because I'm breaking them out across tratamientoss.
# Slope doesn't change as a function on any other variable. 
# plot(emtrends(b1.c.pre, ~tratamientos, var="tiempo_centrado"), xlab="Slope in Log Scale")
# plot(emtrends(b1.weib, ~tratamientos, var="c.dia.comp"), xlab="Slope in Log Scale")
# plot(emtrends(b1.ignore.pre, ~tratamientos, var="tiempo_centrado"), xlab="Slope in Log Scale")
# plot(emtrends(bn.pre, ~tratamientos, var="tiempo_centrado"), xlab="Slope in Log Scale")
# 



## plots juntos ----

##  Another way to show results all on the same plot - not shown in manuscript

b1.c.pre_effects_plot <- plot(conditional_effects(b1.c.pre, effects = "tiempo_centrado:tratamientos", prob=.66))[[1]]+
  ggplot2::ylim(0,80) +
  theme_classic() +
  scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
  scale_fill_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
  labs(
    # caption = "Plotted with 66% Credible Intervals",
    color = "Tratamientos",
    fill = "Tratamientos",
    shape = "Tratamientos",
    x = "Tiempo",
    y = "Segundos") +
  theme(legend.position='top', axis.text = element_text(size = 12),axis.title = element_text(size = 15)) 

# b1.weib_effects_plot <- plot(conditional_effects(b1.weib, effects = "tiempo:tratamientos", prob=.66))[[1]]+
# ggplot2::ylim(-10,180) +
#   theme_classic() +
#   scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
#   scale_fill_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
#   labs(
#     # caption = "Plotted with 66% Credible Intervals",
#     color = "Tratamientos",
#     fill = "Tratamientos",
#     shape = "Tratamientos",
#     x = "Tiempo",
#     y = "Segundos") +
#   theme(legend.position='top')

b1.ignore.pre_effects_plot <- plot(conditional_effects(b1.ignore.pre, effects = "tiempo_centrado:tratamientos", prob=.66))[[1]]+
  ggplot2::ylim(0,80)+
  theme_classic() +
  scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
  scale_fill_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
  labs(
    # caption = "Plotted with 66% Credible Intervals",
    color = "Tratamientos",
    fill = "Tratamientos",
    shape = "Tratamientos",
    x = "Tiempo",
    y = "Segundos") +
  theme(legend.position='top', axis.text = element_text(size = 12),axis.title = element_text(size = 15)) 

bn.pre_effects_plot <- plot(conditional_effects(bn.pre, effects = "tiempo_centrado:tratamientos", prob=.66))[[1]]+
  ggplot2::ylim(0,80)+
  theme_classic() +
  scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
  scale_fill_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9", "#862B0D")) +
  labs(
    # caption = "Plotted with 66% Credible Intervals",
    color = "Tratamientos",
    fill = "Tratamientos",
    shape = "Tratamientos",
    x = "Tiempo",
    y = "Segundos") +
  theme(legend.position='top', axis.text = element_text(size = 12),axis.title = element_text(size = 15)) 


bay_effect_plots <- b1.c.pre_effects_plot + b1.ignore.pre_effects_plot + bn.pre_effects_plot + 
  plot_annotation(tag_levels = 'A') + 
  plot_layout(guides = 'collect', axis_titles = 'collect') &
  theme(legend.position = 'top')

bay_effect_plots

# ggsave("./figuras/bay_effect_plots_pre_otros_grupos.png", plot = bay_effect_plots)

## Bayes otrs grafs preeriors ----

# summary of preerior distribution
# library(bayestestR)
# library(insight)
# 
# preeriors.b1.c.pre <- describe_posterior(b1.c.pre)
# print_md(preeriors.b1.c.pre, digits = 2)
# 
# # visualize posterior distribution
# posteriors.model.b1.c.pre <- get_parameters(b1.c.pre)
# 
# ggplot(posteriors.model.b1.c.pre, aes(x = b_Intercept)) +
#   geom_density(fill = "orange") +
#   ggthemes::theme_base()

# describe 3 elements of posteriors

# point estimate (one-value summary, similar to beta)
# what single value can represent best my post dist?
# mean, median, mode
# mean(posteriors.model.b1.c.pre$b_Intercept)
# median(posteriors.model.b1.c.pre$b_Intercept)
# map_estimate(posteriors.model.b1.c.pre$b_Intercept) # mode - peak - max A posteriori
# 
# ggplot(posteriors.model.b1.c.pre, aes(x = b_Intercept)) +
#   geom_density(fill = "orange") +
#   # The mean in blue
#   geom_vline(xintercept = mean(posteriors.model.b1.c.pre$b_Intercept), color = "blue", linewidth = 1) +
#   # The median in red
#   geom_vline(xintercept = median(posteriors.model.b1.c.pre$b_Intercept), color = "red", linewidth = 1) +
#   # The MAP in purple
#   geom_vline(xintercept = as.numeric(map_estimate(posteriors.model.b1.c.pre$b_Intercept)), color = "purple", linewidth = 1)

# como están muy cerca, vamos a usar la mediana
# mediana tiene significado probabilistico directo (50% que efect real esté arriba/abajo)

# credible interval (associated uncertainty)
# range(posteriors.model.b1.c.pre$b_Intercept)
# en lugar de incluir todos los valores, vamos a calcular intervalos de cred
# similar a intervalos de confianza 
# los vamos a calcular con highest Density Interval
# give us the range containing the 89% most probable effect values
# 89% level gives more stable results (Kruschke, 2014)
# hdi(posteriors.model.b1.c.pre$b_Intercept, ci = 0.89)

# indices of significance (relative importance of this effect)
# si toda la distribucion posterior no está en cero, puede ser evidencia
# para determinar si efecto es sign, podemos ver si el int cre contiene al 0
# si no, puede ser evidencia que nuestro efecto es sign

# se puede hacer con rope rope(posteriors$feedsunflower, range = c(-20, 20), ci = 0.89)

## ROPES ----

library(modelbased)
library(see)

# vizdata <- estimate_relation(b1.c.pre)
# 
# plot_predicted <- ggplot(vizdata, aes(x = tiempo_centrado, y = Predicted)) +
#   geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.5) +
#   geom_line() +
#   facet_grid(~tratamientos) +
#   theme_modern() +
#   # ggthemes::theme_base() +
#   # scale_color_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9")) +
#   # scale_fill_manual(values = c("#222831", "#00ADB5", "#FF2E63", "#FF9800", "#8675A9")) +
#   labs(
#     # caption = "Plotted with 66% Credible Intervals",
#     # color = "Treatments",
#     # fill = "Treatments",
#     # shape = "Treatments",
#     x = "Tiempo",
#     y = "Probabilidad en segundos de Latencia")
# # theme(legend.position='top') 
# 
# plot_predicted

# ggsave("./figuras/plot_predicted.png", plot = plot_predicted)

plot_rope <- plot(rope(b1.c.pre, ci = .95)) +
  ggthemes::theme_fivethirtyeight()  +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 15)) 

plot_rope

ggsave("./figuras/plot_rope_pre_bayes.png", plot = plot_rope)

plot_rope_bn.pre <- plot(rope(bn.pre, ci = .95)) +
  ggthemes::theme_fivethirtyeight()  +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 15)) 

plot_rope_bn.pre

ggsave("./figuras/plot_rope_pre_lmer.png", plot = plot_rope_bn.pre)
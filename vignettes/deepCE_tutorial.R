## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----warning = FALSE, message=FALSE-------------------------------------------
# Loading packages required for data handling & visualization
library(ggplot2)
library(tidyr)
library(dplyr)
library(igraph)
library(fgsea)

# Loading keras and tensorflow
library(keras)
library(tensorflow)

# Loading DeepiCE package
library(DeepiCE)
# Loading the demo data
data(HelaCE)
dim(HelaCE)

## -----------------------------------------------------------------------------
dataP <- impute_MissingData(HelaCE,bound_left = 2,bound_right = 2)

## ---- eval=FALSE--------------------------------------------------------------
#  dataP <- RemoveSinglePeak(dataP)

## -----------------------------------------------------------------------------
dataP <- filter_ConsecutivePep(dataP,min_stretch_length=2)

## -----------------------------------------------------------------------------
dataP <- scaling_m(dataP)

## -----------------------------------------------------------------------------
conc_dat <- getPPI(dataP)

## -----------------------------------------------------------------------------
# Load reference set 
data("refcpx")
t_data <- build_trainingData(conc_dat, refcpx)

## ----warning = FALSE, message=FALSE,results="hide"----------------------------
# set the seed to ensure reproducible output
set_random_seed(2)
MLP_interactions <- 
  MLP_m(conc_dat, #concatenated co-elution profiles
      t_data$train_d, #training data matrix
      t_data$train_l, #training data label
      cv_fold =3,
      cutoff = 0.5)

## ----warning = FALSE, message=FALSE, eval=FALSE-------------------------------
#  # set the seed to ensure reproducible output
#  set_random_seed(3)
#  cnn_interactions <-
#    oneD_CNN(conc_dat, #concatenated co-elution profiles
#        t_data$train_d, #training data matrix
#        t_data$train_l, #training data label
#        cv_fold =3,
#        cutoff = 0.5)

## ----results="hide", warning = FALSE, fig.show="hide"-------------------------
set_random_seed(4)
tuning_result <-
    MLP_tuning(t_data$train_d, 
                t_data$train_l,
                nlayers = 2,
                powerto1 = c(4,5,6),
                powerto2 = c(4,5), 
                b_size = 64, 
                drate = 0.8,
                metrics = "accuracy",
                epochs = 20, k = 0.3)

## -----------------------------------------------------------------------------
f<- 
  tuning_result %>%
  filter(metric == "loss") 

min_val_error <- 
  f %>%
  filter(data == "validation") %>%
  group_by(run) %>% summarise(min_value.error = 
                                round(min(value, na.rm = TRUE), 3))
head(min_val_error)

## ----results="hide", warning = FALSE------------------------------------------
ggplot(f, aes(epoch, value, color = data)) +
  geom_point(size = 0.1) +
  geom_line() +
  theme_bw() +
  facet_wrap(.~ run)+
  geom_label(data = min_val_error, aes(label=min_value.error), 
             x = Inf, y = -Inf, hjust=1, vjust=0,
             inherit.aes = FALSE)


## -----------------------------------------------------------------------------
pred_cpx <- get_clusters(csize = 2, 
                         d = 0.3, p = 2,mx_overlap = 0.8,
                         tpath =file.path(system.file("extdata", 
                                                      package = "DeepiCE")))
dim(pred_cpx)

## -----------------------------------------------------------------------------
pred_cpx_mcl <- 
  MCL_clustering(MLP_interactions,pred_cpx, inflation = 8, csize =2)
dim(pred_cpx_mcl)

## -----------------------------------------------------------------------------
# first load the reference complex
data("refcpx")
set.seed(103)
Clust_tuning_result <-
  clusterONE_tuning(refcpx, csize = 3, 
                d = c(0.3,0.4),
                p = c(2, 2.5),
                mx_overlap = c(0.6,0.7),
                tpath =
                  file.path(system.file("extdata", package = "DeepiCE")))

## -----------------------------------------------------------------------------
ggplot(Clust_tuning_result, aes(tune_names, compScore)) +
  geom_point(size = 3) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size=15)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(
    color = "black", # label color and size
    size = 12))  +
  theme(axis.ticks = element_line(
    colour = "black",
    size = 0.75, linetype = "solid"
  )) +
  theme(axis.ticks.length=unit(0.3, "cm"))


## -----------------------------------------------------------------------------
pred_cpx_optimized <- get_clusters(csize = 2, 
                         d = 0.3, p = 2.5,mx_overlap = 0.7,
                         tpath =file.path(system.file("extdata", 
                                                      package = "DeepiCE")))
dim(pred_cpx_optimized)

## -----------------------------------------------------------------------------
mcl_tuning_result <- 
  MCL_tuning(MLP_interactions,
             pred_cpx_optimized, 
             refcpx,
             inflation = c(6,8,9))

## -----------------------------------------------------------------------------
finla_clusters <- 
  MCL_clustering(MLP_interactions,
             pred_cpx_optimized, 
             inflation = 9, 
             csize = 2)
dim(finla_clusters)

## ----results="hide", message = FALSE, warning = FALSE-------------------------
data("HelaCE")
data("refcpx")
MLP_prediction_outputs <- 
  predPPI_MLP(HelaCE,
              refcpx,
              cv_fold = 5,
              tpath = tempdir())

## ----message = FALSE, warning = FALSE, eval=FALSE-----------------------------
#  data("HelaCE")
#  data("refcpx")
#  oneDCNN_prediction_outputs <-
#    predPPI_1D.CNN(HelaCE,
#                refcpx,
#                nlayers = 1,
#                cv_fold = 3)

## ---- eval=FALSE--------------------------------------------------------------
#  ig <-
#     graph_from_data_frame(MLP_prediction_outputs$ppi_input_ClusterOne)
#  
#  createNetworkFromIgraph(ig,"myIgraph")

## ---- warning=FALSE, message=FALSE, results='hide', eval = FALSE--------------
#  # extract the predicted complexes
#  predcpx <- MLP_prediction_outputs$predicted_cpx_MCL
#  
#  enrich_result <-
#    enrichfindCPX(predcpx,
#                  threshold = 0.05,
#                  sources = "GO:BP",
#                  p.corrction.method = "bonferroni",
#                  org = "hsapiens")

## ----message=FALSE, results='hide', warning=FALSE, fig.show="hide"------------
# profile 1
data("m1")
# profile 2
data("m2")
# biological term (here is known complexes)
data("refcpx")
diff_output <- diffPPI(m1, 
                       m2, 
                       term = refcpx,
                       minSize = 2,
                       tpath = tempdir())

## -----------------------------------------------------------------------------
plotEnrichment(diff_output$term_list[["F1F0-ATP synthase, mitochondrial"]],
               diff_output$vec_rank) +
  labs(title="F1F0-ATP synthase, mitochondrial") +
   theme_bw() +
  theme(text = element_text(size=15)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(
    color = "black", # label color and size
    size = 12))  +
  theme(axis.ticks = element_line(
    colour = "black",
    size = 0.75, linetype = "solid"
  )) +
  theme(axis.ticks.length=unit(0.3, "cm"))


## ---- eval=FALSE--------------------------------------------------------------
#  scored_Data <- similarity_score(HelaCE, metric = "pearson")

## -----------------------------------------------------------------------------
sessionInfo()




  #' predPPI_MLP
  #' @title Predict Protein-Protein Interactions via MLP Method and
  #' Putative Complexes via a two-stage clustering procedure.
  #' @param data A data matrix with rows including proteins and fractions
  #' along the columns.
  #' @param refcpx A list of known reference complexes.
  #' @param data_imputed If true, detects and imputes missing values. Defaults
  #' to FALSE. See \code{\link{impute_MissingData}}.
  #' @param bound_left Numeric integer, the minimum number of non-zero
  #' values to the left of a missing value to be replaced with \code{NA}.
  #' Defaults to 2. See \code{\link{impute_MissingData}}.
  #' @param bound_right Numeric integer, the minimum number of non-zero
  #' values to the right of a missing value to be replaced with \code{NA}.
  #' Defaults to 2. See \code{\link{impute_MissingData}}.
  #' @param remove_SinglePeak If TRUE, removes Single Peak. Defaults to FALSE.
  #' see \code{\link{RemoveSinglePeak}}.
  #' @param remove_ConsecutivePep If TRUE, filters peptides by consecutive
  #' peptide detection. Defaults to TRUE.
  #' See \code{\link{filter_ConsecutivePep}}.
  #' @param min_stretch_length Numeric integer, the minimal length a stretch of
  #' continuous identifications has to have in order not to be removed.
  #' Defaults to 2. See \code{\link{filter_ConsecutivePep}}.
  #' @param scaling If TRUE, performs row-wise normalization. Defaults to TRUE.
  #' See \code{\link{scaling_m}}.
  #' @param nlayers Number of hidden layers. Defaults to 2.
  #' See \code{\link{MLP_m}}.
  #' @param powerto1 Integer, the number of neurons in the first hidden layer
  #' as defined by two to the power of this value. Defaults to 6.
  #' See \code{\link{MLP_m}}.
  #' @param powerto2 Integer, the number of neurons in the subsequent hidden
  #' layer as defined by two to the power of this value.Defaults to 7.
  #' See \code{\link{MLP_m}}.
  #' @param drate Numeric, the dropout rates.
  #' Defaults to 0.1. See \code{\link{MLP_m}}.
  #' @param optimizer Name of the optimizer.
  #' For most models, this defaults to "rmsprop".See \code{\link{MLP_m}}.
  #' @param b_size Number of samples per gradient update.
  #' Defaults to 128 See \code{\link{MLP_m}}.
  #' @param epochs Number of epochs to train the model.
  #' Defaults to 50. See \code{\link{MLP_m}}.
  #' @param cv_fold Number of partitions for cross-validation; Defaults to 5.
  #' See \code{\link{MLP_m}}.
  #' @param cutoff An integer range between [0,1] for specifying
  #' the cutoff for classifier confidence score to subset the high-confidence
  #' interactions.Defaults to 0.5.
  #' @param csize  Numerical value, the minimum size of the predicted complexes.
  #' Defaults to 2. See \code{\link{get_clusters}} or
  #' \code{\link{MCL_clustering}}.
  #' @param d A number, density of predicted complexes; defaults to 0.3.
  #' See \code{\link{get_clusters}}.
  #' @param p An integer, penalty value for the inclusion of each node;
  #' defaults to 2. See \code{\link{get_clusters}}.
  #' @param mx_overlap A number, specifies the maximum allowed
  #' overlap between two clusters; defaults to 0.8.
  #' See \code{\link{get_clusters}}.
  #' @param inflation MCL inflation parameter. Defaults to 4.
  #' See \code{\link{MCL_clustering}}.
  #' @param plots Logical value, indicating whether to plot the performance of
  #' the learning algorithm using k-fold cross-validation; Defaults to FALSE.
  #' These plots are :
  #' \itemize{ \item{pr_plot} - Precision-recall PLOT
  #' \item{roc_plot} - ROC plot
  #' \item{radar_plot} - Radar plot showing
  #' accuracy, F1-score , positive predictive value (PPV), sensitivity (SE)
  #' and MCC.}
  #' @param tpath A character string indicating the path to the project
  #' directory. If the directory is
  #' missing, it will be stored in the Temp directory.
  #' @return Return three data sets in the project directory including:
  #'  \itemize{
  #'  \item{ppi_input_ClusterOne.txt} - threshold-based network predicted
  #'  from MLP.
  #'  \item{predicted_cpx_clusterOne} - predicted complexes resulted
  #'  from partitioning the threshold-based network predicted from MLP using
  #'  ClusterOne algorithm.
  #'  \item{predicted_cpx_MCL} - predicted complexes resulted
  #'  from partitioning the threshold-based network predicted from MLP using
  #'  MCL algorithm.}
  #' @description This function first begins by executing several
  #' pre-processing steps to improve the quality of the raw data, followed
  #' by concatenating the co-elution profiles of putative interacting
  #' proteins. Concatenated profiles and class labels generated
  #' from reference complexes are then fed into a multi-layer perceptron (MLP)
  #' deep learning model. These models then generate a weighted protein
  #' interaction network in which edge weights between protein nodes represent
  #' the deep learning model’s probability estimate for interaction.
  #' High-confidence PPIs based on the user-defined threshold is then will be
  #' partitioned via a two-stage approach, clustering first using ClusterOne,
  #' and then further clustering with MCL.
  #' @export
  #' @importFrom tidyr pivot_longer
  #' @examples
  #' # Load the input data
  #' data("HelaCE")
  #' HelaCE <- HelaCE[1:500, ]
  #' # Known reference complexes
  #' data("refcpx")
  #' library(tensorflow)
  #' set_random_seed(2)
  #' # Perform prediction
  #' MLP_prediction_outputs <-
  #' predPPI_MLP(HelaCE,
  #' refcpx,
  #' min_stretch_length =3,
  #' cv_fold = 3,
  #' cutoff = 0.5,
  #' epochs = 5,
  #' tpath = tempdir())


  predPPI_MLP <-
    function(data,
             refcpx,
             data_imputed = FALSE,
             bound_left = 2,
             bound_right = 2,
             remove_SinglePeak = FALSE,
             remove_ConsecutivePep = TRUE,
             min_stretch_length =2,
             scaling = TRUE,
             nlayers = 2,
             powerto1 = 6,
             powerto2 = 7,
             drate = 0.1,
             optimizer = "rmsprop",
             epochs = 50,
             b_size = 128,
             cv_fold = 5,
             plots = FALSE,
             cutoff = 0.5,
             csize =2,
             d= 0.3,
             p = 2,
             mx_overlap = 0.8,
             inflation = 4,
             tpath = tempdir()){

      if (!is.matrix(data)) {
        data <- as.matrix(data)
      }

      if(is.character(data) == TRUE){
        stop("matrix must include numerical variables")
      }

      if(is.null(row.names(data))){
        stop("Please specify the row.names")
      }

      if(!is.list(refcpx)){
        stop("Reference complexes must be list")
      }



      rowname <- NULL
      value <- NULL
      PPI <- NULL
      Prob <- NULL


      data_n <- data

      # Data imputation
      if(data_imputed) {
        data_n <-
          impute_MissingData(data_n,bound_left = bound_left,
                             bound_right = bound_right)
        message("Imputation completes ...")
      }

      # Remove single peaks
      if(remove_SinglePeak) {
        data_n <-
          RemoveSinglePeak(data_n)
        message("Removing single peak completes ...")
      }

      if(remove_ConsecutivePep) {

        data_n <-
          filter_ConsecutivePep(data_n, min_stretch_length=min_stretch_length)

        message("Removing consecutive peptides completes ...")
        message("Number of retained proteins:", nrow(data_n))
      }

      if(scaling) {
        data_n <-
          scaling_m(data_n)
        message("Scaling completes ...")
      }

      # construct PPI profile
      message("Generating co-elution profile ...")
      PPI_prof <- getPPI(data_n)


      # Build training set
      result <- build_trainingData(PPI_prof, refcpx)
      train_d <- result$train_d
      train_l <- result$train_l


      output <- list()
      message("Generating PPIs with MLP ...")
      pred_dat <-
        MLP_m(PPI_prof,
              train_d,
              train_l,
              nlayers = nlayers,
              powerto1 = powerto1,
              powerto2 = powerto2,
              drate = drate,
              optimizer = optimizer,
              epochs = epochs,
              b_size = b_size,
              cv_fold = cv_fold,
              plots = plots,
              cutoff = cutoff,
              tpath = tpath)

      fname <- file.path(tpath, "ppi_input_ClusterOne.txt")
      write.table(pred_dat,
                  file = fname,
                  row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)


      # generate clusters via ClusterOne
      predcpx <-
        get_clusters(csize =csize,
                     d= d,
                     p = p,
                     mx_overlap = mx_overlap,
                     tpath = tpath)

      fname <- file.path(tpath, "predicted_complexes_ClusterOne.txt")
      write.table(predcpx, file = fname,
                  row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)


      # generate clusters via MCL
      MCL_predcpx <-
        MCL_clustering(pred_dat,
                       predcpx,
                       inflation = inflation,
                       csize =csize)

      fname <- file.path(tpath, "predicted_complexes_MCL.txt")
      write.table(MCL_predcpx, file = fname,
                  row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)



      output[["ppi_input_ClusterOne"]] <- pred_dat
      output[["predicted_cpx_clusterOne"]] <- predcpx
      output[["predicted_cpx_MCL"]] <- MCL_predcpx




      return(output)

    }



  #' build_trainingData
  #' @title Construct Training Data
  #' @param conc_data Concatenated co-elution matrix generated from
  #' getPPI function. See \code{\link{getPPI}}.
  #' @param refcpx A list of protein complexes to create class labels for
  #' training set.
  #' @return A list containing training data, along with class labels.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom RcppAlgos comboGeneral
  #' @importFrom dplyr sample_n
  #' @description This function creates a training set with class labels,
  #' where "1" corresponds to positive interactions (i.e., interactions
  #' belongs to  the same complex), and "0" corresponds to negative
  #' interactions (i.e., random interactions).The training set then can be
  #' used as input for model tuning and model training.
  #' @export
  #' @examples
  #' # upload raw data
  #' data("HelaCE")
  #' HelaCE <- HelaCE[1:100,]
  #' # concatenate co-elution profile
  #' conc_dat <- getPPI(HelaCE)
  #' # upload reference data
  #' data("refcpx")
  #' # build labeled data
  #' tdata <- build_trainingData(conc_dat, refcpx)



  build_trainingData <-
    function(conc_data, refcpx){


    if (!is.matrix(conc_data)) {
      conc_data <- as.matrix(conc_data)
    }

    if(is.character(conc_data) == TRUE){
      stop("matrix must include numerical variables")
    }

    if(is.null(row.names(conc_data))){
      stop("Please specify the row.names")
    }

    if(!is.list(refcpx)){
      stop("Refcpx Must Be List")
    }

    # remove redundancy in complexes
    custom_bg <-
        unique(unlist(strsplit(row.names(conc_data), "~")))
    refcpx <-
        RemoveCpxRedundance(refcpx, custom_bg)


    # create a reference set
    PPI <- data.frame(row.names(conc_data))
    colnames(PPI) <- c("PPI")
    PPI <- separate(PPI, PPI, c("p1", "p2"), sep = "~")

    class_labels <-
      generate_refInt(PPI,refcpx)


    pair <- unite(class_labels, PPI, c("p1", "p2"), sep = "~")
    d <- conc_data[match(pair$PPI,row.names(conc_data)),]
    train_d <- subset(d, row.names(d) %in% pair$PPI)


    # create a numeric label set
    train_l <- pair$labels


    output_list <- list()
    output_list[["train_d"]] <- train_d
    output_list[["train_l"]] <- train_l

    return(output_list)


  }

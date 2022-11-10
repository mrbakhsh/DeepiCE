  .int = function(mat) {
    present = !is.na(mat) & mat > 0
    intersect = crossprod(present)
    J = intersect
    return(J)
  }

  .m_df <- function(x, sname = "WCC"){
    rowname <- NULL
    name <- NULL
    datPPI <-
      x %>% as.data.frame %>% rownames_to_column() %>%
      pivot_longer(-rowname) %>% na.omit()%>%
      mutate(`PPI` =paste(`rowname`, `name`, sep = "~")) %>%
      select(4,3)
    colnames(datPPI)[2] <- sname
    return(datPPI)

  }



  #' getPPI
  #' @title Generating Protein-Protein Interaction (PPI)
  #' @param data A Co-elution data matrix with proteins in rows and
  #' fractions in columns.
  #' @return A concatenated matrix.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @description This function first removes proteins pairs for which two
  #' proteins never occur in at least one fraction, then concatenates
  #' co-elution profiles for each possible protein pairs to form a master
  #' matrix for input into the learner.
  #' @importFrom dplyr bind_cols
  #' @export
  #' @examples
  #' data("HelaCE")
  #' data <- HelaCE[1:100, 1:100]
  #' conc_dat <- getPPI(data)



  getPPI <-
    function(data){


      if (!is.matrix(data)) {
        data <- as.matrix(data)
      }

      if(is.character(data) == TRUE){
        stop("matrix must include numerical variables")
      }

      if(is.null(row.names(data))){
        stop("Please specify the row.names")
      }

      rowname <- NULL
      value <- NULL
      intersection <- NULL


      #remove those that observed in at least >= 2 fraction
      data <- data[sort(rownames(data)), ]
      m = .int(t(data))
      m[lower.tri(m, diag=TRUE)] <- NA
      int_d <-
        .m_df(m, sname = "intersection") %>%
        filter(intersection >= 1) %>% select(-intersection)


      pair <- int_d$PPI
      PPI <-
       separate(int_d, PPI, c("p1", "p2"), sep = "~")
      p1 <-
        PPI %>% select("p1")
      p2 <-
        PPI %>% select("p2")


      data <- as.data.frame(data)
      data <- rownames_to_column(data, "p1")

      out_p1<-plyr::join(p1, data, by = "p1")
      out_p1 <- out_p1[, -1]
      colnames(out_p1) <- paste0(colnames(out_p1), "_1")
      colnames(data)[1] <- "p2"
      out_p2<-plyr::join(p2, data, by = "p2")
      out_p2 <- out_p2[, -1]
      colnames(out_p2) <- paste0(colnames(out_p2), "_2")


      co_elutdf <- bind_cols(out_p1,out_p2)
      co_elutdf <- as.matrix(co_elutdf)
      row.names(co_elutdf) <- pair


      return(co_elutdf)
    }

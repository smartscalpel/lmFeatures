#' Custom print table function
#'
#' Set up defaults for printing the LaTeX table with xtable. Function do not
#' return anything it just sends strings to the output. Useful for use with
#' knitr or Sweave.
#'
#' @param mat data to be printed matrix, or data convertable to matrix.
#' @param digits cutoff and rounding parameter. See xtable documentation.
#'
#' @return none
#' @export
print_table <-
function(mat,digits=-5){print_tableL(mat,digits)}

#' Print table in the portrate mode.
#'
#' See print_table.
#'
#' @param mat data to be printed matrix, or data convertable to matrix.
#' @param digits cutoff and rounding parameter. See xtable documentation.
#'
#' @return none
#' @export
#' @import xtable
#'
print_tableP<-function(mat,digits=-5){
  addtorow          <- list()
  addtorow$pos      <- list()
  addtorow$pos[[1]] <- c(0)
  addtorow$command  <- c(
    paste(
      "\\hline \n",
      "\\endhead \n",
      "\\hline \n",
      "\\multicolumn{3}{l}{\\footnotesize Continued on next page} \n",
      "\\endfoot \n",
      "\\endlastfoot \n",sep = ""
    )
  )
  cat(
    sprintf(
      "\\begin{center}\n\\captionof{table}{Wide ranges of continious peaks (width>%d)}\n\\scriptsize",50
    )
  )
  print(
    xtable(
      mat,digits=digits)
    ,size = "small",include.colnames = TRUE,
    tabular.environment = "longtable",
    floating = FALSE,include.rownames = TRUE,
    add.to.row = addtorow,
    hline.after =c(-1)
  )
  cat("\\end{center}\n")
}

#' Print table in the landscape mode.
#'
#' See print_table.
#'
#' @param mat data to be printed matrix, or data convertable to matrix.
#' @param digits cutoff and rounding parameter. See xtable documentation.
#'
#' @return none
#' @export
#' @import xtable
#'
print_tableL<-function(mat,digits=-5){
  addtorow          <- list()
  addtorow$pos      <- list()
  addtorow$pos[[1]] <- c(0)
  addtorow$command  <- c(
    paste(
      "\\hline \n",
      "\\endhead \n",
      "\\hline \n",
      "\\multicolumn{3}{l}{\\footnotesize Continued on next page} \n",
      "\\endfoot \n",
      "\\endlastfoot \n",sep = ""
    )
  )
  cat(
    sprintf(
      "\\newpage\n  \\begin{landscape} \n\\begin{center}\n\\captionof{table}{Wide ranges of continious peaks (width>%d))}\n\\scriptsize",50
    )
  )
  print(
    xtable(
      mat,digits=digits)
    ,size = "small",include.colnames = TRUE,
    tabular.environment = "longtable",
    floating = FALSE,include.rownames = TRUE,
    add.to.row = addtorow,
    hline.after =c(-1)
  )
  cat("\\end{center}\n \\end{landscape}")
}

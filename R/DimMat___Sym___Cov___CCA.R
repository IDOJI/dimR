DimMat___Sym___Cov___CCA = function(Cov.list,
                                    explained_var_prop = 0.9,
                                    epsilon = 1e-6,
                                    method = c("Algorithm_1", "Algorithm_2"),
                                    path_Export = NULL,
                                    file.name = NULL){
  #=============================================================================
  # Arguments
  #=============================================================================
  # Cov.list : 공분산 행렬들의 리스트, N개의 p * p 행렬
  # epsilon: 수렴 기준
  # Set the default method to the first element if not provided
  method <- match.arg(method)






  #=============================================================================
  # Algorithm
  #=============================================================================
  if(method == "Algorithm_1"){

    Results = DimMat___Sym___Cov___CCA___Algorithm1(Cov.list, explained_var_prop, epsilon)

  }else if(method == "Algorithm_2"){

    Results = DimMat___Sym___Cov___CCA___Algorithm2(Cov.list, explained_var_prop, epsilon)

  }








  #=============================================================================
  # Export
  #=============================================================================
  if(!is.null(path_Export)){

    if(!is.null(file.name)){
      file.name = "CCA Results"
    }

    saveRDS(Results, paste0(path_Export, "/", file.name, ".rds"))
    cat("\n", crayon::green("Exporting the results is done!") ,"\n")

  }







  #=============================================================================
  # Return
  #=============================================================================
  return(Results)
}


















# Changing the splitting criterion in a classification tree

The purpose of this repository is to implement a classification tree from scratch in order to parameterise the user's desired splitting criterion. If you want to test your own criterion, you can follow the indications given in the R file and apply the function you want. 

**Mode indication are available in the .Rmd file**

You only have to choose your criterion in the next code. As an example, we implement a very simple criterion equal to the maximum frequency within the node.

```{r}
newcrit <- function(vec, type){
  if(type == 'values'){vec <- table(vec)/sum(table(vec))}
  if(type == 'table'){vec <- vec/sum(vec)}
  res <- max(vec)
  return(1 - res)
}
```

Run the folowing code for all these functions, you don't have to change anything. 

```{r}
delta_criter <- function(vec, vecL, vecR, type){
  if(type == 'values'){
    vec <- table(vec)
    vecL <- table(vecL)
    vecR <- table(vecR)
  }
  delta <- newcrit(vec, "table") - sum(vecL)*newcrit(vecL, "table")/sum(vec) - sum(vecR)*newcrit(vecR, "table")/sum(vec)
  return(delta)
}

delta_gini <- function(vec, vecL, vecR, type){
  if(type == 'values'){
    vec <- table(vec)
    vecL <- table(vecL)
    vecR <- table(vecR)
  }
  delta <- gini(vec, "table") - sum(vecL)*gini(vecL, "table")/sum(vec) - sum(vecR)*gini(vecR, "table")/sum(vec)
  return(delta)
}

indice <- function(to_cut, seuil, y, indic){
  whereL <- which(to_cut <= seuil)
  whereR <- which(to_cut > seuil)
  if(indic == "gini"){
    delta_gini(y[c(whereL,whereR)], y[whereL], y[whereR], type = "values")
  } else if (indic == 'newcrit'){
    delta_criter(y[c(whereL,whereR)], y[whereL], y[whereR], type = "values")
  }
}

all_cut <- function(df, y, to_cut, indic){
  vals <- sort(unique(df[,to_cut]))
  vals <- vals[-length(vals)]
  
  data.frame(to_cut = to_cut ,
             cut = vals,
             indice = sapply(vals, FUN = function(x) indice(df[,to_cut], x, df[, y] ,indic)))
}

clean_df <- function(df, tokeep){
  
  tab_unique <- apply(df[, colnames(df) != tokeep], FUN = function(x) length(unique(x)), MARGIN = 2) 
  col_not_unique <- names(tab_unique[tab_unique != 1])
  df1 <- df1[, c(col_not_unique, tokeep)]
  
}

find_opt <- function(df, yname, indic){
  
  v <- apply(df[, colnames(df) != yname], FUN = function(x) length(unique(x)), MARGIN = 2) 
  col_ <- names(v[v != 1])
  
  result <- lapply(col_, FUN = function(x) all_cut(df[,], yname, x, indic)) %>% bind_rows()
  idx <- which.max(result$indice)
  
  return(list(to_cut = result$to_cut[idx],
                            cut = result$cut[idx]))
  
}

accuracy_table <- function(res, n){
  accuracy <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(accuracy) <- c("depth", "accuracy")
  for(i in unique(res$depth)){
    
    locs <- which(res$depth == i)
    
    n_err <- 0
    for(loc in locs){
      
      tab <- res$table[loc][[1]]
      
      n_err <- n_err + (sum(tab) - max(tab))
      
    }
    
    accuracy <- accuracy %>% add_row(depth = i,
                                     accuracy = 1 - n_err/n)
  }
    return(accuracy)
}

rpart_scratch <- function(df, yname, n_depth, indic){

  df[, yname] <- as.factor(df[, yname])
  
  res <- data.frame(depth = 0,
                    node = 1,
                    type = "to do",
                    where = I(list(1:nrow(df))),
                    table = I(list(table(df[, yname]))),
                    pred_ = which.max(table(df[, yname])) %>% names(),
                    to_cut = NA,
                    cut = NA)
  
  for(i in 1:n_depth){
    
    if(sum(res$type == "to do") == 0){
      break
    }
    
    locs <- which(res$type == "to do")
    
    for(loc in locs){
    
    d <- res$depth[loc] 
    id_node <- res$node[loc] 
    
    where <- res$where[loc][[1]]
    
    if(sum(table(df[where, yname]) != 0) == 1){
      
      res[loc, 'type'] <- "leaf"
      
    } else {
    
      result <- find_opt(df[where,], yname, indic)
      
      res[loc, 'to_cut'] <- result$to_cut
      res[loc, 'cut'] <- result$cut
      res[loc, 'type'] <- "node"
      
      whereL <- where[where %in% which(df[, result$to_cut] <= result$cut)]
      whereR <- where[where %in% which(df[, result$to_cut] > result$cut)]
      
      res <- res %>% add_row(depth = c(d + 1, d + 1),
                      node = c(id_node * 2, id_node *2 + 1),
                      type = c("to do", "to do"),
                      where = c(I(list(whereL)), I(list(whereR))),
                      table = c(I(list(table(df[whereL, yname]))), I(list(table(df[whereR, yname])))),
                      pred_ = c(names(which.max(table(df[whereL, yname]))), names(which.max(table(df[whereR, yname])))))
    
    }
    }
  }
  return(res)
}

rpart_scratch_pred <- function(obs, res){
  node <- 1
  pred_ <- NULL
  for(i in 1:50){
  if(!is.null(pred_)){
    break
  }
  loc <- which(res$node == node)
  if(sum(res$node == 2 * node) == 0){
    pred_ <- res$pred_[loc]
  } else {
    if(obs[, res$to_cut[loc]] <= res$cut[loc]){
      node <- node * 2
    } else {
      node <- node * 2 + 1
    }
  }
  }
  return(pred_)
}

rpart_scratch_pred_all <- function(df, res){
  v <- NULL
  for(i in 1:nrow(df)){
    v <- c(v, rpart_scratch_pred(df[i,], res))
  } 
  return(v)
  
}
```

Example of CART application on the IRIS dataset. 

Application with the Gini index: the following table describes the leaves of the tree obtained as well as the nodes and their optimal splitting. 

```{r}
rpart_scratch(iris, "Species", 6, 'gini')
```

The following table, which shows the results for the newcrit index, allows us to observe that the choices do not change, whatever the indicator.

```{r}
rpart_scratch(iris, "Species", 6, 'newcrit')
```


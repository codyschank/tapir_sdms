modelSpecies = function(inRow){
  
  # Utility functions
  logit = function(pp) { log(pp) - log(1-pp) }
  expit = function(eta) {1/(1+exp(-eta))}
  
  MPA_calc = function(psi,increment,testData){
    
    results = data.frame(threshold=numeric(),percent_correct=numeric(),total_area=numeric())
    thresholds = seq(0,1,by=increment)
    
    for(i in seq(1,length(thresholds))){
      
      thresholded = psi > thresholds[i]
      totalArea = cellStats(thresholded,"sum")*resolution^2
      extracted = extract(thresholded,testData)
      percent_correct = sum(extracted,na.rm=TRUE)/sum(!is.na(extracted))
      
      results[i,"threshold"] = thresholds[i]
      results[i,"total_area"] = totalArea
      results[i,"percent_correct"] = percent_correct
      
    }
    mpa = results[max(which(results$percent_correct>0.9)),"total_area"]
    return(mpa)
  }
  
  BI_calc = function(psi,increment,width,testData){
    
    results = data.frame(threshold1=numeric(),threshold2=numeric(),avg_suitability=numeric(),
                         P_i=numeric(),E_i=numeric(),F_i=numeric())
    thresholds = seq(0,1-width,by=increment)
    
    for(i in seq(1,length(thresholds))){
      
      threshold1 = thresholds[i]
      threshold2 = threshold1+width
      
      b = ( (psi > threshold1) & (psi < threshold2) )
      b_e = extract(b,testData)
      
      P_i = sum(b_e,na.rm=TRUE)/sum(!is.na(b_e))
      E_i = cellStats(b,"sum")/cellStats(!is.na(b),"sum")
      
      results[i,"threshold1"] = threshold1
      results[i,"threshold2"] = threshold2
      results[i,"P_i"] = P_i
      results[i,"E_i"] = E_i
      
    }
    
    results[,"avg_suitability"] = results[,"threshold1"] + width/2
    results[,"F_i"] = results[,"P_i"]/results[,"E_i"]
    filter = !is.na(results$F_i)
    bi = stats::cor(x=results[filter,"threshold1"],y=results[filter,"F_i"])
    return(bi)
  }
  
  ObsInfo.po = function(param) {
    
    beta = param[1:dim(X.po)[2]]
    alpha = param[(dim(X.po)[2]+1):(dim(X.po)[2]+dim(W.po)[2])]
    
    lambda = exp(X.back %*% beta)
    mu = lambda * area.back
    p = expit(W.back %*% alpha)
    
    p.po = expit(W.po %*% alpha)
    
    nxcovs = length(beta)
    nwcovs = length(alpha)
    
    nparams = nxcovs + nwcovs
    Hmat = matrix(nrow=nparams, ncol=nparams)
    
    #  beta partials
    for (i in 1:nxcovs) {
      for (j in 1:i) {
        Hmat[i,j] = sum(X.back[,i] * X.back[,j] * mu * p)
        Hmat[j,i] = Hmat[i,j]
      }
    }
    
    # alpha partials
    for (i in 1:nwcovs) {
      for (j in 1:i) {
        Hmat[nxcovs+i, nxcovs+j] = sum(W.back[,i] * W.back[,j] * mu * p * ((1-p)^3) * (1 - exp(2 * W.back %*% alpha)) ) + sum(W.po[,i] * W.po[,j] * p.po * (1-p.po))
        Hmat[nxcovs+j, nxcovs+i] = Hmat[nxcovs+i, nxcovs+j]
      }
    }
    
    # alpha-beta partials
    for (i in 1:nwcovs) {
      for (j in 1:nxcovs) {
        Hmat[nxcovs+i, j] = sum(X.back[,j] * W.back[,i] * mu * p * (1-p))
        Hmat[j, nxcovs+i] = Hmat[nxcovs+i, j]
      }
    }
    
    Hmat
  }
  
  negLL.po = function(param) {
    
    beta = param[1:dim(X.po)[2]]
    alpha = param[(dim(X.po)[2]+1):(dim(X.po)[2]+dim(W.po)[2])]
    
    lambda = exp(X.back %*% beta)
    mu = lambda * area.back
    p = expit(W.back %*% alpha)
    
    logL.po = sum(X.po %*% beta) + sum(W.po %*% alpha) - sum(log(1 + exp(W.po %*% alpha))) - sum(mu*p)
    
    (-1)*sum(logL.po)
  }
  
  negLL.pa = function(param) {
    beta = param[1:dim(X.pa)[2]]
    lambda.pa = exp(X.pa %*% beta)
    
    alpha = param[(dim(X.pa)[2]+1):(dim(X.pa)[2]+dim(W.pa)[2])]
    p.pa = expit(W.pa %*% alpha)
    
    logL.pa = rep(NA,n.pa)
    
    for (i in 1:n.pa) {
      yvec=y.pa[i,]
      navec=is.na(yvec)
      nd=sum(yvec[!navec])
      nj=sum(!navec)
      pvec=p.pa[i,]
      cp= (pvec^yvec)*((1-pvec)^(1-yvec))
      cp[navec]=1
      logL.pa[i]= log(prod(cp)*(1-exp(-lambda.pa[i]*area.pa[i])) + ifelse(nd==0,1,0)*exp(-lambda.pa[i]*area.pa[i]))
    }
    
    (-1)*sum(logL.pa)
  }
  
  negLL.pa.RN = function(param,Ninfty=50) {
    beta = param[1:dim(X.pa)[2]]
    lambda.pa = exp(X.pa %*% beta)
    alpha = param[(dim(X.pa)[2]+1):(dim(X.pa)[2]+dim(W.pa)[3])]
    rvec = matrix(nrow=n.pa, ncol=J.pa)
    for(j in 1:J.pa){rvec[, j] = expit(as.matrix(W.pa[,j,], nrow=n.pa) %*% alpha)}
    lik = rep(NA,n.pa)
    for(i in 1:n.pa){
      gN = dpois(0:Ninfty, lambda.pa[i]*area.pa[i])
      gN = gN/sum(gN)
      dvec = y.pa[i,]
      naflag = is.na(dvec)
      PMAT = 1-outer( (1-rvec[i,]), 0:Ninfty, "^")
      LIK = t((PMAT^dvec)*((1-PMAT)^(1-dvec)))
      LIK[,naflag] = 1
      likvec = apply(LIK,1,prod)
      lik[i] = sum(likvec*gN)
    }
    (-1)*sum(log(lik))
  }
  
  negLL.poANDpa = function(param)  {
    
    param.po = param[1:(dim(X.po)[2]+dim(W.po)[2])]
    param.pa = param[c(1:dim(X.po)[2], (dim(X.po)[2]+dim(W.po)[2]+1):(dim(X.po)[2]+dim(W.po)[2]+dim(W.pa)[2]))]
    negLL.po(param.po) + negLL.pa(param.pa)
  }
  
  negLL.poANDpa.RN = function(param,Ninfty=50)  {
    
    param.po = param[1:(dim(X.po)[2]+dim(W.po)[2])]
    param.pa = param[c(1:dim(X.po)[2], (dim(X.po)[2]+dim(W.po)[2]+1):(dim(X.po)[2]+dim(W.po)[2]+dim(W.pa)[2]))]
    negLL.po(param.po) + negLL.pa.RN(param.pa)
  }
  
  library(dismo)
  
  minrecipCondNum = 1e-6
  
  # double brackets on first one when inside loop
  # 1 is the thinned data, 2 is the left over data, and 3 is the radius
  #presenceDataTrain = inRow[1][[1]][[1]]
  #presenceDataTest = inRow[1][[1]][[2]]
  presenceData = inRow[1][[1]][[1]]
  r.po = inRow[1][[1]][[3]]
  paDataTrain = inRow[2][[1]][[1]]
  #paDataTest = inRow[2][[1]][[2]]
  r.pa = inRow[2][[1]][[3]]
  # remove presence data sets from list
  inRow = inRow[3:length(inRow)]
  
  # split presenceData Train up into 
  k.index = kfold(presenceData,k=4)
  presenceDataTest = presenceData[k.index==1,]
  presenceDataTrain = presenceData[k.index>1,]
  
  po_test_data = presenceDataTest
  coordinates(po_test_data) = c("X","Y")
  
  ### get info and prep output dataframe
  
  resolution = as.numeric(inRow["resolution"][[1]])
  n.samples = as.numeric(inRow["n.samples"][[1]])
  
  candidate_model= inRow["candidate_model"][[1]]
  outputFolder = paste0(outputFolder.base,"/model_",candidate_model,"/")
  
  outputs = data.frame(inRow)
  outputs = outputs[rep(seq_len(1), each=length(models)),]
  outputs$model = models

  model_iteration = inRow["model.iteration"][[1]]
  
  # remove detection histories with missing data?
  missing.data.option = inRow["missing.data.option"][[1]]
  if(missing.data.option==0){
    paDataTrain = paDataTrain[!rowSums(is.na(paDataTrain[,paste("sample",seq(1,n.samples),sep="_")]))>0,]}
  
  file_name_base = paste0(resolution,"km_",n.samples,"samples_",missing.data.option,"NA_",model_iteration)
  
  n.po = dim(presenceDataTrain)[1]
  n.pa = dim(paDataTrain)[1]
  
  outputs['n.po'] = n.po
  outputs['n.pa'] = n.pa
  outputs['r.pa'] = r.pa
  outputs['r.po'] = r.po
  
  if(!is.na(inRow[["x.po.covs"]])){names.x.po.covs = strsplit(inRow[["x.po.covs"]],",")[[1]]} else{names.x.po.covs=NULL}
  if(!is.na(inRow[["w.po.covs"]])){names.w.po.covs = strsplit(inRow[["w.po.covs"]],",")[[1]]} else{names.w.po.covs=NULL}
  if(!is.na(inRow[["x.pa.covs"]])){names.x.pa.covs = strsplit(inRow[["x.pa.covs"]],",")[[1]]} else{names.x.pa.covs=NULL}
  if(!is.na(inRow[["w.pa.covs"]])){names.w.pa.covs = strsplit(inRow[["w.pa.covs"]],",")[[1]]} else{names.w.pa.covs=NULL}
  names_allCovs = unique(c(names.x.po.covs,names.w.po.covs,names.x.pa.covs,names.w.pa.covs))
  
  ### prep data
  
  allStack.files = list.files(path=paste0("/Users/codyschank/Dropbox/Research/Dissertation/dataPrep_output/",resolution,"km/"),full.names=TRUE)
  allStack = stack(allStack.files)
  
  sgrid = stack(subset(allStack,names_allCovs))

  xycov = na.omit(rasterToPoints(sgrid))
  xy = xycov[,1:2]
  outputs['n.quad.points'] = nrow(xy)
  des = as.data.frame(xycov[,3:dim(xycov)[2]])
  names(des)=names_allCovs
  X.back = as.matrix(cbind(rep(1, dim(des)[1]), des[,names.x.po.covs]))
  W.back = as.matrix(cbind(rep(1, dim(des)[1]), des[,names.w.po.covs]))
  po.des = data.frame(na.omit(extract(sgrid,presenceDataTrain[,c("X","Y")])))
  n.po = dim(po.des)[1]
  X.po = as.matrix(cbind(rep(1, dim(po.des)[1]), po.des[,names.x.po.covs]))
  W.po = as.matrix(cbind(rep(1, dim(po.des)[1]), po.des[,names.w.po.covs]))
  area.back = rep(resolution^2, dim(X.back)[1])
  pa.extract = extract(sgrid,paDataTrain[,c("X","Y")])
  filterNA = !is.na(pa.extract)
  filterNA = apply(filterNA, 1, prod)
  filterNA = filterNA==1
  pa.des = as.data.frame(na.omit(pa.extract))
  paDataSelect = paDataTrain[filterNA,]
  X.pa = as.matrix(cbind(rep(1,dim(pa.des)[1]),pa.des[,names.x.pa.covs]))
  W.pa.extract = as.matrix(pa.des[,names.w.pa.covs])
  W.pa = as.matrix(cbind(rep(1, dim(W.pa.extract)[1]),W.pa.extract))
   
  # detection histories
  y.pa = paDataSelect[,paste("sample",seq(1,n.samples),sep="_")]
  y.pa = as.matrix(y.pa)
  outputs['NA.ratio'] = sum(is.na(y.pa))/(nrow(y.pa)*n.samples)
    
  n.pa = dim(y.pa)[1]
  J.pa =  dim(y.pa)[2]
    
  area.pa = rep(resolution^2, n.pa)
  
  betaGuess = rep(0, dim(X.po)[2])
  alphaGuess.po = c(logit(n.po/nrow(xy)),rep(0, (dim(W.po)[2]-1))) # use naive detectability
  #alphaGuess.po = rep(0, dim(W.po)[2]) 
  alphaGuess.pa = c(logit(sum(y.pa,na.rm=T)/sum(y.pa==0,na.rm=T)),rep(0, (dim(W.pa)[2]-1))) 
  #alphaGuess.pa = rep(0, dim(W.pa)[2])
  
  if("po" %in% models){
    
    paramGuess = c(betaGuess, alphaGuess.po)
    
    possibleError.po = tryCatch(
      expr = (fit.po = optim(par=paramGuess, fn=negLL.po, method='BFGS', hessian=FALSE, control=list(maxit=200))),
      error=function(e) e
    )
    
    if(!inherits(possibleError.po, "error")){
      recipCondNum.po = NA
      se.po = rep(NA, length(fit.po$par))
      if (fit.po$convergence==0) {		
        hess = ObsInfo.po(fit.po$par)
        possibleError.hess = tryCatch(expr = (ev = eigen(hess)$values),error=function(e) e)
        if(!inherits(possibleError.hess, "error")){
          recipCondNum.po = ev[length(ev)]/ev[1]
          if (!is.na(recipCondNum.po)){
            if (recipCondNum.po>minrecipCondNum) {
              vcv = chol2inv(chol(hess))
              se.po = sqrt(diag(vcv))
            }
          }
          if (is.na(recipCondNum.po)) {se.po=rep(NA, length(fit.po$par))}
        }
      }else{se.po=rep(NA, length(fit.po$par))}
      
      # Calculate SDM surface
      linearPredictor.po.sdm = X.back %*% fit.po$par[1:dim(X.back)[2]]
      predict.po.sdm = exp(linearPredictor.po.sdm)
      raster.po.sdm = raster(sgrid)
      cells = cellFromXY(raster.po.sdm, xy)
      raster.po.sdm[cells] = predict.po.sdm
      
      # Predicted surface that includes detetectability, for use with test data to calculate residuals
      raster.po.sdm.psi = 1-exp(-1*raster.po.sdm*resolution^2)
      
      po.k = length(se.po)
      # value is negative log likelihood, that's why we use + not -
      poAIC = (2*po.k) + (2*fit.po$value)
      
      mpa.po = MPA_calc(raster.po.sdm.psi,0.01,po_test_data)
      bi.po = BI_calc(raster.po.sdm.psi,0.01,0.1,po_test_data)
      
      outputs[outputs$model=='po',"mpa.po"]=mpa.po
      outputs[outputs$model=='po',"bi.po"]=bi.po
      
      outputs[outputs$model=='po','conv'] = fit.po$convergence
      outputs[outputs$model=='po','recip'] = recipCondNum.po
      outputs[outputs$model=='po','value'] = fit.po$value
      outputs[outputs$model=='po','AIC'] = poAIC
      outputs[outputs$model=='po','total_population'] = sum(predict.po.sdm*resolution^2,na.rm=TRUE)
      outputs[outputs$model=='po','estimated.se'] = !is.na(se.po[1])
      outputs[outputs$model=='po','converged'] = TRUE
      
      po.coefs = data.frame(t(unlist(fit.po$par)))
      if(length(names.x.po.covs)>0){po.coefs.x.names=paste0('x.',names.x.po.covs)}else{po.coefs.x.names=NULL}
      if(length(names.w.po.covs)>0){
        po.coefs.w.names=paste0('w.po.',names.w.po.covs)
      }else{po.coefs.w.names=NULL}
      names(po.coefs) = c('beta0',po.coefs.x.names,'alpha0.po',po.coefs.w.names)
      outputs[outputs$model=='po',names(po.coefs)]=po.coefs
      
      po.se = data.frame(t(unlist(se.po)))
      if(length(names.x.po.covs)>0){po.se.x.names=paste0('x.',names.x.po.covs,'.se')}else{po.se.x.names=NULL}
      if(length(names.w.po.covs)>0){po.se.w.names=paste0('w.po.',names.w.po.covs,'.se')}else{po.se.w.names=NULL}
      names(po.se) = c('beta0.se',po.se.x.names,'alpha0.po.se',po.se.w.names)
      outputs[outputs$model=='po',names(po.se)]=po.se
      
    }else{
      outputs[outputs$model=='po','converged'] = FALSE
    }
  }
  
  ###Dorazio_pa###
  if("pa" %in% models){
    
    W.pa = as.matrix(cbind(rep(1, dim(W.pa.extract)[1]),W.pa.extract))
    
    paramGuess = c(betaGuess, alphaGuess.pa)
    fit.pa = NA
    
    possibleError.pa = tryCatch(
      expr = (fit.pa = optim(par=paramGuess, fn=negLL.pa, method='BFGS', hessian=TRUE,control=list(maxit=200))),
      error=function(e) e
    )
    
    if(!inherits(possibleError.pa, "error")){
      recipCondNum.pa = NA
      se.pa = rep(NA, length(fit.pa$par))
      if (fit.pa$convergence==0) {
        hess = fit.pa$hessian
        ev = eigen(hess)$values
        recipCondNum.pa = ev[length(ev)]/ev[1]
        if (recipCondNum.pa>minrecipCondNum) {
          vcv = chol2inv(chol(hess))
          se.pa = sqrt(diag(vcv))
        }
      }
      
      # Calculate the SDM surface
      linearPredictor.pa.sdm = X.back %*% fit.pa$par[1:dim(X.back)[2]]
      predict.pa.sdm = exp(linearPredictor.pa.sdm)
      raster.pa.sdm = raster(sgrid)
      cells = cellFromXY(raster.pa.sdm, xy)
      raster.pa.sdm[cells] = predict.pa.sdm
      
      raster.pa.sdm.psi = 1-exp(-1*raster.pa.sdm*resolution^2)
      
      pa.k = length(se.pa)
      paAIC = (2*pa.k) + (2*fit.pa$value)
      outputs[outputs$model=='pa','conv'] = fit.pa$convergence
      outputs[outputs$model=='pa','recip'] = recipCondNum.pa
      outputs[outputs$model=='pa','AIC'] = paAIC
      outputs[outputs$model=='pa','value'] = fit.pa$value
      outputs[outputs$model=='pa','total_population'] = sum(predict.pa.sdm*resolution^2, na.rm=TRUE)
      outputs[outputs$model=='pa','estimated.se'] = !is.na(se.pa[1])
      outputs[outputs$model=='pa','converged'] = TRUE
      
      mpa.po = MPA_calc(raster.pa.sdm.psi,0.01,po_test_data)
      bi.po = BI_calc(raster.pa.sdm.psi,0.01,0.1,po_test_data)
      
      outputs[outputs$model=='pa',"mpa.po"]=mpa.po
      outputs[outputs$model=='pa',"bi.po"]=bi.po
      
      pa.coefs = data.frame(t(unlist(fit.pa$par)))
      if(length(names.x.pa.covs)>0){pa.coefs.x.names=paste0('x.',names.x.pa.covs)}else{pa.coefs.x.names=NULL}
      if(length(names.w.pa.covs)>0){
        pa.coefs.w.names=paste0('w.pa.',names.w.pa.covs)
      }else{pa.coefs.w.names=NULL}
      names(pa.coefs) = c('beta0',pa.coefs.x.names,'alpha0.pa',pa.coefs.w.names)
      outputs[outputs$model=='pa',names(pa.coefs)]=pa.coefs
      
      pa.se = data.frame(t(unlist(se.pa)))
      if(length(names.x.pa.covs)>0){pa.se.x.names=paste0('x.',names.x.pa.covs,'.se')}else{pa.se.x.names=NULL}
      if(length(names.w.pa.covs)>0){pa.se.w.names=paste0('w.pa.',names.w.pa.covs,'.se')}else{pa.se.w.names=NULL}
      names(pa.se) = c('beta0.se',pa.se.x.names,'alpha0.pa.se',pa.se.w.names)
      outputs[outputs$model=='pa',names(pa.se)]=pa.se
    }else{
      outputs[outputs$model=='pa','converged'] = FALSE
    }
  }
    
  ###Dorazio_poANDpa###
  if("poANDpa" %in% models){
    
    W.pa = as.matrix(cbind(rep(1, dim(W.pa.extract)[1]),W.pa.extract))
    
    paramGuess = c(betaGuess, alphaGuess.po, alphaGuess.pa)
    fit.poANDpa = NA
    
    possibleError.poANDpa = tryCatch(
      expr = (fit.poANDpa = optim(par=paramGuess, fn=negLL.poANDpa, method='BFGS', hessian=TRUE, control=list(maxit=200))),
      error=function(e) e
    )
    
    if(!inherits(possibleError.poANDpa, "error")){
      recipCondNum.poANDpa = NA
      se.poANDpa = rep(NA, length(fit.poANDpa$par))
      if (fit.poANDpa$convergence==0) {
        hess = fit.poANDpa$hessian
        ev = eigen(hess)$values
        recipCondNum.poANDpa = ev[length(ev)]/ev[1]
        if (recipCondNum.poANDpa>minrecipCondNum) {
          vcv = chol2inv(chol(hess))
          se.poANDpa = sqrt(diag(vcv))
        }
      }	
      
      # Calculate the SDM surface
      linearPredictor.poANDpa.sdm = X.back %*% fit.poANDpa$par[1:dim(X.back)[2]]
      predict.poANDpa.sdm= exp(linearPredictor.poANDpa.sdm)
      raster.poANDpa.sdm = raster(sgrid)
      cells = cellFromXY(raster.poANDpa.sdm, xy)
      raster.poANDpa.sdm[cells] = predict.poANDpa.sdm
      
      raster.poANDpa.sdm.psi = 1-exp(-1*raster.poANDpa.sdm*resolution^2)
      
      poANDpa.k = length(se.poANDpa)
      poANDpaAIC = (2*poANDpa.k) + (2*fit.poANDpa$value)
      
      mpa.po = MPA_calc(raster.poANDpa.sdm.psi,0.01,po_test_data)
      bi.po = BI_calc(raster.poANDpa.sdm.psi,0.01,0.1,po_test_data)
      
      outputs[outputs$model=='poANDpa',"mpa.po"]=mpa.po
      outputs[outputs$model=='poANDpa',"bi.po"]=bi.po
      
      outputs[outputs$model=='poANDpa','conv'] = fit.poANDpa$convergence
      outputs[outputs$model=='poANDpa','recip'] = recipCondNum.poANDpa
      outputs[outputs$model=='poANDpa','AIC'] = poANDpaAIC
      outputs[outputs$model=='poANDpa','value'] = fit.poANDpa$value
      outputs[outputs$model=='poANDpa','total_population'] = sum(predict.poANDpa.sdm*resolution^2,na.rm=TRUE)
      outputs[outputs$model=='poANDpa','estimated.se'] = !is.na(se.poANDpa[1])
      outputs[outputs$model=='poANDpa','converged'] = TRUE
      
      poANDpa.coefs = data.frame(t(unlist(fit.poANDpa$par)))
      if(length(names.x.pa.covs)>0){poANDpa.coefs.x.names=paste0('x.',names.x.pa.covs)}else{poANDpa.coefs.x.names=NULL}
      if(length(names.w.po.covs)>0){poANDpa.coefs.w.po.names=paste0('w.po.',names.w.po.covs)}else{poANDpa.coefs.w.po.names=NULL}
      if(length(names.w.pa.covs)>0){poANDpa.coefs.w.pa.names=paste0('w.pa.',names.w.pa.covs)}else{poANDpa.coefs.w.pa.names=NULL}
      names(poANDpa.coefs) = c('beta0',poANDpa.coefs.x.names,'alpha0.po',poANDpa.coefs.w.po.names,'alpha0.pa',poANDpa.coefs.w.pa.names)
      outputs[outputs$model=='poANDpa',names(poANDpa.coefs)]=poANDpa.coefs
      
      poANDpa.se = data.frame(t(unlist(se.poANDpa)))
      if(length(names.x.pa.covs)>0){poANDpa.se.x.names=paste0('x.',names.x.pa.covs,'.se')}else{poANDpa.se.x.names=NULL}
      if(length(names.w.po.covs)>0){poANDpa.se.w.po.names=paste0('w.po.',names.w.po.covs,'.se')}else{poANDpa.se.w.po.names=NULL}
      if(length(names.w.pa.covs)>0){poANDpa.se.w.pa.names=paste0('w.pa.',names.w.pa.covs,'.se')}else{poANDpa.se.w.pa.names=NULL}
      names(poANDpa.se) = c('beta0.se',poANDpa.se.x.names,'alpha0.po.se',poANDpa.se.w.po.names,'alpha0.pa.se',poANDpa.se.w.pa.names)
      outputs[outputs$model=='poANDpa',names(poANDpa.se)]=poANDpa.se
      
    }else{
      outputs[outputs$model=='poANDpa','converged'] = FALSE
    }
  }
  
  ###Dorazio_pa.RN###
  if("pa.RN" %in% models){
    
    ptm = proc.time()
    
    W.pa = array(dim=c(n.pa, J.pa, (length(names.w.pa.covs)+1)))
    W.pa[,,1] = 1
    if(length(names.w.pa.covs)>0){
      W.pa[,,2:(length(names.w.pa.covs)+1)] = matrix(rep(W.pa.extract,J.pa), ncol=J.pa)
    }
    
    paramGuess = c(betaGuess, alphaGuess.pa)
    
    possibleError.pa = tryCatch(
      expr = (fit.pa = optim(par=paramGuess, fn=negLL.pa.RN, method='BFGS', hessian=TRUE,control=list(maxit=200))),
      error=function(e) e
    )
    
    if(!inherits(possibleError.pa, "error")){
      recipCondNum.pa = NA
      se.pa = rep(NA, length(fit.pa$par))
      if (fit.pa$convergence==0) {
        hess = fit.pa$hessian
        ev = eigen(hess)$values
        recipCondNum.pa = ev[length(ev)]/ev[1]
        if (recipCondNum.pa>minrecipCondNum) {
          vcv = chol2inv(chol(hess))
          se.pa = sqrt(diag(vcv))
        }
      }
      
      # Calculate the SDM surface
      linearPredictor.pa.sdm = X.back %*% fit.pa$par[1:dim(X.back)[2]]
      predict.pa.sdm = exp(linearPredictor.pa.sdm)
      raster.pa.sdm = raster(sgrid)
      cells = cellFromXY(raster.pa.sdm, xy)
      raster.pa.sdm[cells] = predict.pa.sdm
      
      # Predicted surface that includes detetectability, for use with test data to calculate residuals
      raster.pa.sdm.psi = 1-exp(-1*raster.pa.sdm*resolution^2)
      
      pa.k = length(se.pa)
      paAIC = (2*pa.k) + (2*fit.pa$value)
      
      mpa.po = MPA_calc(raster.pa.sdm.psi,0.01,po_test_data)
      bi.po = BI_calc(raster.pa.sdm.psi,0.01,0.1,po_test_data)
      
      outputs[outputs$model=='pa.RN',"mpa.po"]=mpa.po
      outputs[outputs$model=='pa.RN',"bi.po"]=bi.po
      
      outputs[outputs$model=='pa.RN','conv'] = fit.pa$convergence
      outputs[outputs$model=='pa.RN','recip'] = recipCondNum.pa
      outputs[outputs$model=='pa.RN','value'] = fit.pa$value
      outputs[outputs$model=='pa.RN','AIC'] = paAIC
      outputs[outputs$model=='pa.RN','total_population'] = sum(predict.pa.sdm*resolution^2, na.rm=TRUE)
      outputs[outputs$model=='pa.RN','estimated.se'] = !is.na(se.pa[1])
      outputs[outputs$model=='pa.RN','converged'] = TRUE
      
      pa.coefs = data.frame(t(unlist(fit.pa$par)))
      if(length(names.x.pa.covs)>0){pa.coefs.x.names=paste0('x.',names.x.pa.covs)}else{pa.coefs.x.names=NULL}
      if(length(names.w.pa.covs)>0){
        pa.coefs.w.names=paste0('w.pa.',names.w.pa.covs)
      }else{pa.coefs.w.names=NULL}
      names(pa.coefs) = c('beta0',pa.coefs.x.names,'alpha0.pa',pa.coefs.w.names)
      outputs[outputs$model=='pa.RN',names(pa.coefs)]=pa.coefs
      
      pa.se = data.frame(t(unlist(se.pa)))
      if(length(names.x.pa.covs)>0){pa.se.x.names=paste0('x.',names.x.pa.covs,'.se')}else{pa.se.x.names=NULL}
      if(length(names.w.pa.covs)>0){pa.se.w.names=paste0('w.pa.',names.w.pa.covs,'.se')}else{pa.se.w.names=NULL}
      names(pa.se) = c('beta0.se',pa.se.x.names,'alpha0.pa.se',pa.se.w.names)
      outputs[outputs$model=='pa.RN',names(pa.se)]=pa.se
      
    }else{
      outputs[outputs$model=='pa.RN','converged'] = FALSE
    }
  }
  
  ###Dorazio_poANDpa.RN###
  if("poANDpa.RN" %in% models){
    
    W.pa = array(dim=c(n.pa, J.pa, (length(names.w.pa.covs)+1)))
    W.pa[,,1] = 1
    if(length(names.w.pa.covs)>0){
      W.pa[,,2:(length(names.w.pa.covs)+1)] = matrix(rep(W.pa.extract,J.pa), ncol=J.pa)
    }
    
    paramGuess = c(betaGuess, alphaGuess.po, alphaGuess.pa)
    fit.poANDpa = NA
    
    possibleError.poANDpa = tryCatch(
      expr = (fit.poANDpa = optim(par=paramGuess, fn=negLL.poANDpa.RN, method = "BFGS", hessian=TRUE, control=list(maxit=200))),
      error=function(e) e
    )
    
    if(!inherits(possibleError.poANDpa, "error")){
      recipCondNum.poANDpa = NA
      se.poANDpa = rep(NA, length(fit.poANDpa$par))
      if (fit.poANDpa$convergence==0) {
        hess = fit.poANDpa$hessian
        ev = eigen(hess)$values
        recipCondNum.poANDpa = ev[length(ev)]/ev[1]
        if (recipCondNum.poANDpa>minrecipCondNum) {
          vcv = chol2inv(chol(hess))
          se.poANDpa = sqrt(diag(vcv))
        }
      }	
      
      # Calculate the SDM surface
      linearPredictor.poANDpa.sdm = X.back %*% fit.poANDpa$par[1:dim(X.back)[2]]
      predict.poANDpa.sdm= exp(linearPredictor.poANDpa.sdm)
      raster.poANDpa.sdm = raster(sgrid)
      cells = cellFromXY(raster.poANDpa.sdm, xy)
      raster.poANDpa.sdm[cells] = predict.poANDpa.sdm
      
      # Predicted surface that includes detetectability, for use with test data to calculate residuals
      raster.poANDpa.sdm.psi = 1-exp(-1*raster.poANDpa.sdm*resolution^2)
      
      poANDpa.k = length(se.poANDpa)
      poANDpaAIC = (2*poANDpa.k) + (2*fit.poANDpa$value)
      
      mpa.po = MPA_calc(raster.poANDpa.sdm.psi,0.01,po_test_data)
      bi.po = BI_calc(raster.poANDpa.sdm.psi,0.01,0.1,po_test_data)
      
      outputs[outputs$model=='poANDpa.RN',"mpa.po"]=mpa.po
      outputs[outputs$model=='poANDpa.RN',"bi.po"]=bi.po
      
      outputs[outputs$model=='poANDpa.RN','conv'] = fit.poANDpa$convergence
      outputs[outputs$model=='poANDpa.RN','recip'] = recipCondNum.poANDpa
      outputs[outputs$model=='poANDpa.RN','AIC'] = poANDpaAIC
      outputs[outputs$model=='poANDpa.RN','value'] = fit.poANDpa$value
      outputs[outputs$model=='poANDpa.RN','total_population'] = sum(predict.poANDpa.sdm*resolution^2,na.rm=TRUE)
      outputs[outputs$model=='poANDpa.RN','estimated.se'] = !is.na(se.poANDpa[1])
      outputs[outputs$model=='poANDpa.RN','converged'] = TRUE
      
      poANDpa.coefs = data.frame(t(unlist(fit.poANDpa$par)))
      if(length(names.x.pa.covs)>0){poANDpa.coefs.x.names=paste0('x.',names.x.pa.covs)}else{poANDpa.coefs.x.names=NULL}
      if(length(names.w.po.covs)>0){poANDpa.coefs.w.po.names=paste0('w.po.',names.w.po.covs)}else{poANDpa.coefs.w.po.names=NULL}
      if(length(names.w.pa.covs)>0){poANDpa.coefs.w.pa.names=paste0('w.pa.',names.w.pa.covs)}else{poANDpa.coefs.w.pa.names=NULL}
      names(poANDpa.coefs) = c('beta0',poANDpa.coefs.x.names,'alpha0.po',poANDpa.coefs.w.po.names,'alpha0.pa',poANDpa.coefs.w.pa.names)
      outputs[outputs$model=='poANDpa.RN',names(poANDpa.coefs)]=poANDpa.coefs
      
      poANDpa.se = data.frame(t(unlist(se.poANDpa)))
      if(length(names.x.pa.covs)>0){poANDpa.se.x.names=paste0('x.',names.x.pa.covs,'.se')}else{poANDpa.se.x.names=NULL}
      if(length(names.w.po.covs)>0){poANDpa.se.w.po.names=paste0('w.po.',names.w.po.covs,'.se')}else{poANDpa.se.w.po.names=NULL}
      if(length(names.w.pa.covs)>0){poANDpa.se.w.pa.names=paste0('w.pa.',names.w.pa.covs,'.se')}else{poANDpa.se.w.pa.names=NULL}
      names(poANDpa.se) = c('beta0.se',poANDpa.se.x.names,'alpha0.po.se',poANDpa.se.w.po.names,'alpha0.pa.se',poANDpa.se.w.pa.names)
      outputs[outputs$model=='poANDpa.RN',names(poANDpa.se)]=poANDpa.se
      
    }else{
      outputs[outputs$model=='poANDpa.RN','converged'] = FALSE
    }
  }
  
  removeTmpFiles(h=0)
  gc()
  
  write.csv(outputs,paste0(outputFolder,'outputs_df_',file_name_base,'.csv'),row.names=FALSE)
}

#######################################

library(raster)
library(rgdal)

sr_cea = "+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

source("/Users/codyschank/Dropbox/Research/Dissertation/scripts/rarefy_function.R")

inputs = read.csv("/Users/codyschank/Dropbox/Research/Dissertation/inputs_chp3.csv",stringsAsFactors=FALSE)
presenceData = read.csv("/Users/codyschank/Dropbox/Research/Dissertation/Presence_Data/presence_only_proj.csv")
paData = read.csv(paste0("/Users/codyschank/Dropbox/Research/Dissertation/Presence_Data/detection_histories_proj_10samples.csv"))

n.iterations = 10
n.samples = c(3,6,9)
#resolutions also could include 32, 64
resolutions = c(4)
missing.data.option = c(1) # 0 does not allow NAs. 1 does
#models = c("maxent","hSDM.ZIB","occu","occuRN","po","pa","pa.RN","poANDpa","poANDpa.RN")
models = c("po","pa","poANDpa")

r.po = resolutions*1000*2^0.5
r.pa = resolutions*1000*2^0.5
outputFolder.base = paste0("/Users/codyschank/Dropbox/Research/Dissertation/outputs_chp3_",resolutions,"km_final/")
dir.create(file.path(outputFolder.base), showWarnings = FALSE)

# create folders
candidate_models = unique(inputs$candidate_model)
for(candidate_model in candidate_models){
  candidateFolder = paste0(outputFolder.base,"/model_",candidate_model,"/")
  dir.create(file.path(candidateFolder), showWarnings = FALSE)
}

outputFolder.samples = paste0(outputFolder.base,"samples/")
dir.create(file.path(outputFolder.samples), showWarnings = FALSE)

### thinning data ###

# spatially thin the presence data
presenceDataList = list()
presenceDataTrain.rownames = list()
presenceDataTest.rownames = list()

for(i in seq(1,n.iterations)){
  presenceDataList[[i]] = rarefy(presenceData,r.po)
  
  presenceDataTrain = presenceDataList[[i]][[1]]
  presenceDataTest = presenceDataList[[i]][[2]]
  coordinates(presenceDataTrain) = c("X","Y")
  coordinates(presenceDataTest) = c("X","Y")
  proj4string(presenceDataTrain) = sr_cea
  proj4string(presenceDataTest) = sr_cea
  
  presenceDataTrain.rownames = c(presenceDataTrain.rownames,row.names(presenceDataTrain))
  presenceDataTest.rownames = c(presenceDataTest.rownames,row.names(presenceDataTest))
  
  writeOGR(presenceDataTrain, dsn = paste0(outputFolder.samples,"presenceDataTrain_",i,".shp"), layer= paste0("presenceDataTrain_",i), driver="ESRI Shapefile", overwrite_layer=TRUE)
  writeOGR(presenceDataTest, dsn = paste0(outputFolder.samples,"presenceDataTest_",i,".shp"), layer= paste0("presenceDataTest_",i), driver="ESRI Shapefile", overwrite_layer=TRUE)
  
}

#for all data at once
#presenceDataList[[1]]= list(presenceData,presenceData,r.po)

presenceDataTrain.count = data.frame(table(unlist(presenceDataTrain.rownames)))
presenceData.merge = merge(presenceData,presenceDataTrain.count,by.x="row.names",by.y="Var1",all.x=TRUE)
write.csv(presenceData.merge,paste0(outputFolder.samples,"presenceDataTrain.csv"),row.names=FALSE)

paDataList = list()
paDataTrain.rownames = list()
paDataTest.rownames = list()

for(i in seq(1,n.iterations)){
  paDataList[[i]] = rarefy(paData,r.pa)
  
  paDataTrain = paDataList[[i]][[1]]
  paDataTest = paDataList[[i]][[2]]
  coordinates(paDataTrain) = c("X","Y")
  coordinates(paDataTest) = c("X","Y")
  proj4string(paDataTrain) = sr_cea
  proj4string(paDataTest) = sr_cea
  
  paDataTrain.rownames = c(paDataTrain.rownames,row.names(paDataTrain))
  paDataTest.rownames = c(paDataTest.rownames,row.names(paDataTest))
  
  writeOGR(paDataTrain, dsn = paste0(outputFolder.samples,"paDataTrain_",i,".shp"), layer= paste0("paDataTrain_",i), driver="ESRI Shapefile", overwrite_layer=TRUE)
  writeOGR(paDataTest, dsn = paste0(outputFolder.samples,"paDataTest_",i,".shp"), layer= paste0("paDataTest_",i), driver="ESRI Shapefile", overwrite_layer=TRUE)
  
}

#for all data at once
#paDataList[[1]]= list(paData,paData,r.pa)

paDataTrain.count = data.frame(table(unlist(paDataTrain.rownames)))
paData.merge = merge(paData,paDataTrain.count,by.x="row.names",by.y="Var1",all.x=TRUE)
write.csv(paData.merge,paste0(outputFolder.samples,"paDataTrain.csv"),row.names=FALSE)

### inputs data frame prep ###

#replicate rows by n.samples
n.reps = nrow(inputs)
inputs = inputs[rep(seq_len(nrow(inputs)), each=length(n.samples)),]
inputs['n.samples'] = rep(n.samples,n.reps)

#replicate rows by missing.data.option
n.reps = nrow(inputs)
inputs = inputs[rep(seq_len(nrow(inputs)), each=length(missing.data.option)),]
inputs['missing.data.option'] = rep(missing.data.option,n.reps)

#replicate rows by n.resolutions
n.reps = nrow(inputs)
inputs = inputs[rep(seq_len(nrow(inputs)), each=length(resolutions)),]
inputs['resolution'] = rep(resolutions,n.reps)

#replicate rows by n.iterations
n.reps = nrow(inputs)
inputs = inputs[rep(seq_len(nrow(inputs)), each=n.iterations),]
inputs['model.iteration'] = rep(seq(1,n.iterations),n.reps)

# create data frame to store results
all.x.po = unique(na.omit(unique(unlist(strsplit(inputs$x.po.covs,",")))))
inputs[,paste0('x.',all.x.po)] = NA
inputs[,paste0('x.',all.x.po,'.se')] = NA

if('maxent' %in% models){
  inputs[,paste0(all.x.po,'.permutation.importance')] = NA
  inputs[,paste0(all.x.po,'.contribution')] = NA
}

if(length(inputs$w.po.covs[!is.na(inputs$w.po.covs)])){
  all.w.po = unique(na.omit(unique(unlist(strsplit(inputs$w.po.covs,",")))))
  inputs[,paste0('w.po.',all.w.po)] = NA
  inputs[,paste0('w.po.',all.w.po,'.se')] = NA
}

if(length(inputs$w.pa.covs[!is.na(inputs$w.pa.covs)])){
  all.w.pa = unique(na.omit(unique(unlist(strsplit(inputs$w.pa.covs,",")))))
  inputs[,paste0('w.pa.',all.w.pa)] = NA
  inputs[,paste0('w.pa.',all.w.pa,'.se')] = NA
}else{all.w.pa=NA}

all.pa.vars = c('beta0',paste0('x.',all.x.po),'alpha0',paste0('w.',all.w.pa),'deviance')
if('hSDM.ZIB' %in% models){
  inputs[,paste0(all.pa.vars,'.psrf')] = NA
}

inputs[,c("model","converged","estimated.se","r.po","r.pa","n.po","n.pa","conv","recip","AIC","DIC","value","max.psrf","deviance",
             "total_population","NA.ratio","beta0","beta0.se","alpha0.po","alpha0.po.se","alpha0.pa","alpha0.pa.se","n.quad.points",
              "mpa.po","bi.po")] = NA

presenceDataList.full = rep(presenceDataList,(nrow(inputs)/n.iterations))
paDataList.full = rep(paDataList,(nrow(inputs)/n.iterations))

inputData = cbind(presenceDataList.full,paDataList.full)
inputData = cbind(inputData,inputs)

#random sort inputData, to make sure threads take about the same amount of time 
inputData = inputData[sample(nrow(inputData)),]

library(parallel)
# Initiate cluster
n.cores = 10
cl = makeCluster(spec = rep("localhost",n.cores), mc = getOption("cl.cores", n.cores))
clusterExport(cl=cl, varlist=c("outputFolder.base","rarefy","models"))
parApply(inputData,MARGIN=1,FUN=modelSpecies,cl=cl)
stopCluster(cl)
gc()


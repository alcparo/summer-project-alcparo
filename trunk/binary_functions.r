
# Part of software for the ECML'14 submission "Reliability maps: A tool to enhance probability estimates and improve classification accuracy"
# Author: Meelis Kull, meelis.kull@bristol.ac.uk, meelis.kull@gmail.com
# Written February-April 2014.
# All rights reserved.

library(limSolve)
library(ggplot2)
library(e1071)
library(reshape)
library(RWeka)
library(RColorBrewer)
library(mvtnorm)
library(np)

create_normal_distr = function(name,features,mu,sigma) {
  distr = new.env() 
  distr$family = "normal"
  distr$name = name
  distr$features = features
  distr$mu = mu
  distr$sigma = sigma
  distr$generate = function(n) {
    d = data.frame(rmvnorm(n,distr$mu,distr$sigma))
    colnames(d) = distr$features
    return(d)
  }
  distr$density = function(x) {
    return(dmvnorm(x,distr$mu,distr$sigma))
  }
  return(distr)
}

create_uniwang = function() {
  model = new.env()
  model$name = "uniwang"
  model$n_classes = 5
  model$features = c("x1","x2")
  model$distr_list = list(
    c1 = create_normal_distr("c1",model$features,c(0,0),diag(1,2)),
    c2 = create_normal_distr("c2",model$features,c(3,0),diag(4,2)),
    c3 = create_normal_distr("c3",model$features,c(0,5),diag(9,2)),
    c4 = create_normal_distr("c4",model$features,c(7,0),diag(25,2)),
    c5 = create_normal_distr("c5",model$features,c(0,9),diag(64,2))
  )
  model$classes_factor = factor(names(model$distr_list))
  model$label = "class"
  model$prior = rep(1/model$n_classes,model$n_classes)
  model$generate = function(n,n_per_class=NA,stratified=TRUE) {
    if (!is.na(n_per_class)) {
      if (sum(n_per_class)!=n) {
        stop("sum(n_per_class)!=n")
      }
    } else if (stratified) {
      n_per_class = calc_sizes_towards_target(n,0,model$prior)
    } else {
      n_per_class = count(c(1:model$n_classes,sample(1:model$n_classes,n,replace=TRUE,prob=model$prior)))$freq-1
    }
    mat = data.frame()
    for (i in 1:model$n_classes) {
      if (n_per_class[i]>0) {
        d = model$distr_list[[i]]$generate(n_per_class[i])
        d[[model$label]] = model$classes_factor[i]
        mat = rbind(mat,d)
      }
    }
    colnames(mat) = c(model$features,model$label)
    return(mat)
  }
  model$posterior = function(x) {
    posterior = matrix(NA,nrow=nrow(x),ncol=model$n_classes)
    for (i in 1:model$n_classes) {
      posterior[,i] = model$prior*model$distr_list[[i]]$density(x)
    }
    posterior = posterior/apply(posterior,1,sum)
    colnames(posterior) = model$classes_factor
    rownames(posterior) = rownames(x)
    return(posterior)
  }
  return(model)
}

vlen = function(x) sqrt(sum(x^2))

project_to_line = function(point,normal) {
  nn = normal/vlen(normal)
  return(sign(point[1]*normal[2]-point[2]*normal[1])*vlen(point-sum(point*nn)*nn))
}

signed_dist_to_line = function(point,normal) {
  nn = normal/vlen(normal)
  return(sign(sum(point*nn))*vlen(sum(point*nn)*nn))
}

dist_to_line_with_bias = function(point,normal,bias) {
  # ax+by=bias, normal=c(a,b)
  return(abs(signed_dist_to_line(point,normal)-bias/vlen(normal)))
}

cut_mix = function(w,dist_list) {
  m = vector()
  d = vector()
  s = vector()
  p = vector()
  for (i in 1:length(dist_list)) {
    distr = dist_list[[i]]
    m[i] = project_to_line(distr$mu,w[2:3])
    d[i] = dist_to_line_with_bias(distr$mu,w[2:3],-w[1])
    s[i] = sqrt(distr$sigma[1,1])
    p[i] = 1/s[i]*exp(-d[i]^2/2/s[i]^2)
  }
  p = p/sum(p)
  return(list(m=m,d=d,s=s,p=p))
}

create_true_cal_map = function(lr_weights,model,code) {
  return(function(ss) {
    mm = vector()
    for (j in 1:length(ss)) {
      s = ss[j]
      cm = cut_mix(lr_weights+c(log(1/s-1),0,0),model$distr_list)
      mm[j] = sum(cm$p[code==1])
      if (mm[j]>1) {
        mm[j] = 1
      }
    }
    return(mm)
  })
}

create_true_rel_map = function(lr_weights,model,code,niter=1000) {
  return(function(ss) {
    rr = vector()
    for (j in 1:length(ss)) {
      s = ss[j]
      cm = cut_mix(lr_weights+c(log(1/s-1),0,0),model$distr_list)
      pp = sum(cm$p[code==1])
      minval = min(cm$m-5*cm$s)
      maxval = max(cm$m+5*cm$s)
      values = seq(minval,maxval,length.out=niter)
      fp = 0*values
      ft = 0*values
      for (i in 1:length(code)) {
        if (code[i]==1) {
          fp = fp + cm$p[i]*dnorm(values,cm$m[i],cm$s[i])
        }
        ft = ft + cm$p[i]*dnorm(values,cm$m[i],cm$s[i])
      }
      vv = sum(fp^2/ft)/sum(ft) - pp^2
      if (pp*(1-pp)==0) {
        rr[j] = 1
      } else {
        rr[j] = 1-vv/(pp*(1-pp))
      }
      if (rr[j]>1) {
        rr[j] = 1
      }
      if (rr[j]<0) {
        rr[j] = 0
      }
    }
    return(rr)
  })
}

norm_lower = function(m,d,q=0.025) qnorm(q,m,d)
norm_upper = function(m,d,q=0.975) qnorm(q,m,d)

beta_ab_from_mv = function(mv) {
  m = mv[1]
  s2 = mv[2]
  a_plus_b = m*(1-m)/s2-1
  a = m*a_plus_b
  b = (1-m)*a_plus_b
  return(c(a,b))
}

beta_lower = function(m,d,q=0.025) {
  a_plus_b = m*(1-m)/(d^2)-1
  a = m*a_plus_b
  b = (1-m)*a_plus_b
  res = qbeta(q,a,b)
  res[d==0] = m[d==0]
  return(res)
}

beta_upper = function(m,d,q=0.975) {
  a_plus_b = m*(1-m)/(d^2)-1
  a = m*a_plus_b
  b = (1-m)*a_plus_b
  res = qbeta(q,a,b)
  res[d==0] = m[d==0]
  return(res)
}

kl_dist = function(tmu,ts2,mu,s2) {
  return( (log(s2/ts2)+(ts2+(mu-tmu)^2)/s2-1) / 2 )
}

half_beta_energy = function(tmu,ts2,mu,s2,niter=1000) {
  res = vector()
  for (i in 1:length(tmu)) {
    res[i] = 1/2*beta_energy_mv(tmu[i],ts2[i],mu[i],s2[i],niter=niter)
  }
  return(res)
}

half_norm_energy = function(tmu,ts2,mu,s2,niter=1000) {
  res = vector()
  for (i in 1:length(tmu)) {
    res[i] = 1/2*norm_energy_mv(tmu[i],ts2[i],mu[i],s2[i],niter=niter)
  }
  return(res)
}

norm_energy_mv = function(tmu,ts2,mu,s2,niter=1000) {
  s = (1:niter)/(niter+1)
  ta = pnorm(0,tmu,ts2)
  tb = pnorm(1,tmu,ts2)
  a = pnorm(0,mu,s2)
  b = pnorm(1,mu,s2)
  if (ts2==0) {
    tp = 1*(s>=tmu)
  } else {
    tp = (pnorm(s,tmu,ts2)-ta)/(tb-ta)
  }
  if (s2==0) {
    p = 1*(s>=mu)
  } else {
    p = (pnorm(s,mu,s2)-a)/(b-a)
  }
  d2 = 2*sum((tp-p)^2)/niter
  return(d2)
}

beta_energy_mv = function(tmu,ts2,mu,s2,niter=1000) {
  s = (1:niter)/(niter+1)
  tab = beta_ab_from_mv(c(tmu,ts2))
  if (any(is.na(tab))||(sum(tab)<=0)||(sum(tab)==Inf)) {
    tp = 1*(s>=tmu)
  } else {
    tp = pbeta(s,tab[1],tab[2])
  }
  if (any(is.na(tp))) {
    tp = 1*(s>=tmu)
  }
  ab = beta_ab_from_mv(c(mu,s2))
  if (any(is.na(ab))||(sum(ab)<=0)||(sum(ab)==Inf)) {
    p = 1*(s>=mu)
  } else {
    p = pbeta(s,ab[1],ab[2])
  }
  if (any(is.na(p))) {
    p = 1*(s>=mu)
  }
  d2 = 2*sum((tp-p)^2)/niter
  return(d2)
}

beta_energy_ab = function(a1,b1,a2,b2,niter=1000) {
  s = (1:niter)/(niter+1)
  d2 = 2*sum((pbeta(s,a1,b1)-pbeta(s,a2,b2))^2)/niter
  return(d2)
}

create_datamatrix = function(name) {
  dm = new.env()
  dm$name = name
  dm$mat = data.frame()
  dm$m = 0
  dm$n = 0
  dm$instance_lengths = vector()
  dm$features = list()
  dm$add_instances = function(value_matrix) {
    dm$instance_lengths = c(dm$instance_lengths,rep(ncol(value_matrix),nrow(value_matrix)))
    if (length(dm$features)<ncol(value_matrix)) {
      dm$features = c(dm$features,rep(NA,ncol(value_matrix)-length(dm$features)))
    }
    res = value_matrix
    while (ncol(res)<ncol(dm$mat)) {
      feat = dm$features[[ncol(res)+1]]
      res = cbind(res,feat$calc(res[,feat$inputs]))
    }
    if (ncol(dm$mat)==0) {
      dm$mat = data.frame(res)
    } else {
      dm$mat[nrow(dm$mat)+(1:nrow(res)),1:ncol(res)] = res
    }
    dm$m = ncol(dm$mat)
    dm$n = nrow(dm$mat)
    rownames(dm$mat) = 1:(dm$n)
  }
  dm$add_feature = function(inputs,calc,n_outputs,names) {
    f = new.env(parent=dm)
    if (!is.character(inputs)) {
      inputs = colnames(dm$mat[0,inputs,drop=FALSE]) # get names of input columns as columns can be deleted
    }
    f$inputs = inputs
    f$calc = calc
    dm$features = c(dm$features,f)
    if (n_outputs>=2) {
      dm$features = c(dm$features,rep(NA,n_outputs-1))
    }
    dm$mat = cbind(dm$mat,calc(dm$mat[,inputs]))
    colnames(dm$mat)[(ncol(dm$mat)-n_outputs+1):ncol(dm$mat)] = names
    dm$m = ncol(dm$mat)
    dm$n = nrow(dm$mat)
  }
  dm$delete_features = function(cols) {
    col_names = colnames(dm$mat[0,cols,drop=FALSE])
    col_numbers = which(colnames(dm$mat) %in% col_names)
    rem_col_numbers = setdiff(1:ncol(dm$mat),col_numbers)
    dm$features = dm$features[rem_col_numbers]
    dm$mat = dm$mat[,rem_col_numbers]
    dm$m = ncol(dm$mat)
  }
  return(dm)
}

split_vec_given_sizes = function(vec,sizes) {
  shuffled = vec[sample(1:length(vec))]
  cl = rep(1:length(sizes),sizes)
  folds = list()
  for (i in 1:length(sizes)) {
    folds[[i]] = sort(shuffled[cl==i])
  }
  return(folds)
}

calc_sizes_towards_target = function(n,current,target) {
  ptarget = target/sum(target)
  sizes = ceiling(ptarget*(sum(current)+n))-current
  if (any(sizes<0)) {
    stop("any(sizes<0)")
  }
  decr_sizes = sample(1:length(sizes),sum(sizes)-n)
  sizes[decr_sizes] = sizes[decr_sizes]-1
  if (any(sizes<0)) {
    print('WARNING! Some sizes in calc_sizes_towards_target less than 0')
    print(n)
    print(current)
    print(target)
    print('END-OF-WARNING! Increasing these to 0')
    sizes[sizes<0] = 0
  }
  return(sizes)
}

create_fold = function(datamatrix,rows=1:datamatrix$n,name_suffix="",name=paste(datamatrix$name,name_suffix,sep="")) {
  f = new.env(parent=datamatrix)
  f$datamatrix = datamatrix
  f$rows = rows
  f$name = name
  f$n = length(f$rows)
  f$extract = function(task=NA,code=NA,code_col=NA,label_name=colnames(code$mat)[code_col],posneg=c(1,0),discard_na=FALSE) {
    if (!is.environment(task) && !is.environment(code)) {
      return(f$datamatrix$mat[f$rows,])
    } else if (!is.environment(code)) {
      return(f$datamatrix$mat[f$rows,c(task$features,task$labels)])
    } else {
      task = code$task
      d = f$datamatrix$mat[f$rows,task$features]
      d[[label_name]] = code$calc_labels(f$datamatrix$mat[f$rows,task$labels],code_cols=code_col,posneg=posneg)
      if (discard_na) {
        d = d[!is.na(d[[label_name]]),]
      }
      return(d)
    }
  }
  f$split = function(n_folds,
                     proportions=rep(1,n_folds),
		     sizes=calc_sizes_towards_target(f$n,0,proportions),
		     fold_name_suffices=paste(".fold.",1:n_folds,".of.",n_folds,sep=""),
		     fold_names=paste(f$name,fold_name_suffices,sep=""),
		     stratified=FALSE,task=NA) {
    if (sum(sizes)>f$n) {
      stop("'sizes' must not sum up to more than the total number of instances in the input fold")
    }
    if (stratified) {
      if (length(task$labels)!=1) {
        stop("stratification works currently only on one label")
      }
      classes = unique(f$datamatrix$mat[f$rows,task$labels])
      n_classes = length(classes)
      folds = as.list(sizes)
      for (i in 1:length(folds)) {
        folds[[i]] = vector()
      }
      for (cl in 1:n_classes) {
        vec = f$rows[f$datamatrix$mat[f$rows,task$labels]==classes[cl]]
        newfolds = split_vec_given_sizes(vec,calc_sizes_towards_target(length(vec),sapply(folds,length),sizes))
        for (i in 1:length(newfolds)) {
          folds[[i]] = sort(c(as.numeric(folds[[i]]),newfolds[[i]]))
        }
      }
    } else {
      folds = split_vec_given_sizes(f$rows,sizes)
    }
    names(folds) = fold_names
    for (i in 1:length(folds)) {
      folds[[i]] = create_fold(f$datamatrix,folds[[i]],name=fold_names[i])
    }
    return(folds)
  }
  return(f)
}

create_task = function(datamatrix,features,labels) {
  t = new.env(parent=datamatrix)
  t$datamatrix = datamatrix
  t$features = features
  t$labels = labels
  t$classes_factor = factor(sort(unique(datamatrix$mat[,labels])))
  t$n_classes = length(t$classes_factor)
  return(t)
}

binary_gaussian_optim = function(M,mu,s2) {
  s2 = rep(s2,length.out=length(mu))
  sel = mu<1e-4
  mu[sel] = 1e-4
  s2[sel] = (1e-4)/3
  s2[s2<1e-20] = 1e-20
  s = sqrt(s2)
  M = (M+1)/2 
  A = t(M)/matrix(s,nrow=ncol(M),ncol=nrow(M))
  B = matrix(mu/s)
  E = matrix(rep(1,nrow(M)),nrow=1)
  F = 1
  G = diag(nrow(M))
  H = matrix(rep(0,nrow(M)),ncol=1)
  ep = lsei(A=A,B=B,E=E,F=F,G=G,H=H)
  return(ep$X)
}


# 03/07 Andre: adding type 'binary'
create_ecoc = function(task,type=c("1vsRest","1vs1","custom","binary")[1],custom_mat=NA) {
  if (length(task$labels)!=1) {
    stop("error-correcting codes work currently only on one label")
  }
  code = new.env()
  classes_factor = factor(sort(unique(task$datamatrix$mat[,task$labels])))
  n_classes = length(classes_factor)
  mat = c()
  if (type=="1vsRest") {
    for (i in 1:n_classes) {
      col = 2*((1:n_classes)==i) - rep(1,n_classes)
      mat = cbind(mat,col)
      colnames(mat)[length(colnames(mat))] = paste(i,"R",sep="vs")
    }
  } else if (type=="1vs1") {
    for (i in 1:(n_classes-1)) {
      for (j in (i+1):n_classes) {
        col = 1*((1:n_classes)==i) - 1*((1:n_classes)==j)
        mat = cbind(mat,col)
        colnames(mat)[length(colnames(mat))] = paste(i,j,sep="vs")
      }
    }
  } else if (type=="custom") {
    if (nrow(custom_mat)!=n_classes) {
      stop("number of rows in the code matrix and number of classes must be equal")
    }
    if (!all(unique(as.vector(custom_mat)) %in% c(-1,0,+1))) {
      stop("code matrix can only have values -1,0,+1")
    }
    if ((any(apply(custom_mat,2,max)!=+1)) || (any(apply(custom_mat,2,min)!=-1))) {
      stop("each column in the code matrix must contain both -1 and +1")
    }
    mat = custom_mat
  } else if (type=="binary"){ ### NEW
    n_models=3; ### temp!
    for(i in 1:n_models){
      col = c(-1,+1);
      mat=cbind(mat,col);
      colnames(mat)[length(colnames(mat))] = paste("model",i,sep="")
    }
  } else {
    stop("'1vsRest','1vs1','custom' and 'binary' are the only allowed values for type")
  }
  rownames(mat) = classes_factor
  code$task = task
  code$type = type
  code$classes_factor = classes_factor
  code$n_classes = n_classes
  code$n_cols = ncol(mat)
  code$mat = mat
  code$calc_labels = function(x,code_cols=1:code$n_cols,posneg=c(1,0)) {
    res = matrix(NA,nrow=length(x),ncol=length(code_cols))
    colnames(res) = colnames(code$mat)[code_cols]
    rownames(res) = rownames(x)
    for (col in code_cols) {
      res[code$mat[x,col]==+1,colnames(code$mat)[col]] = posneg[1]
      res[code$mat[x,col]==-1,colnames(code$mat)[col]] = posneg[2]
    }
    return(res)
  }
  code$ls_ecoc = function(cal) {
    res = matrix(NA,nrow=nrow(cal),ncol=code$n_classes)
    for (i in 1:nrow(cal)) {
      res[i,] = binary_gaussian_optim(code$mat,as.numeric(cal[i,]),1)
    }
    return(res)
  }
  code$ls_ecoc_v = function(cal,s2) {
    res = matrix(NA,nrow=nrow(cal),ncol=code$n_classes)
    for (i in 1:nrow(cal)) {
      res[i,] = binary_gaussian_optim(code$mat,as.numeric(cal[i,]),as.numeric(s2[i,]))
    }
    return(res)
  }
  code$ls_ecoc_r = function(cal,rel) {
    s2 = (1-rel)*cal*(1-cal)
    res = matrix(NA,nrow=nrow(cal),ncol=code$n_classes)
    for (i in 1:nrow(cal)) {
      res[i,] = binary_gaussian_optim(code$mat,as.numeric(cal[i,]),as.numeric(s2[i,]))
    }
    return(res)
  }
  return(code)
}

train_lr = function(trainfold,task,code,code_cols=1:code$n_cols) {
  features = list()
  models = list()
  for (i in 1:length(code_cols)) {
    train_data = trainfold$extract(code=code,code_col=code_cols[i],label_name="Y",posneg=c(1,0))
    n_features = ncol(train_data)-1
    feature_names = colnames(train_data)[1:n_features] # paste("X",1:n_features,sep="")
    formula = as.formula(paste("Y",paste(feature_names,collapse="+"),sep="~"))
    model = glm(formula,family=binomial(logit),data=train_data)
    tmp_create_feature = function(model) {
      f = new.env()
      f$model = model
      f$calc = function(x) {
        colnames(x) = feature_names
        predict(f$model,x,type="response")
      }
      f$n_inputs = n_features
      f$n_outputs = 1
      return(f)
    }
    features[[i]] = tmp_create_feature(model)
  }
  return(features)
}

logistic = function(x) 1/(1+exp(-x))

train_svm = function(trainfold,task,code,code_cols=1:code$n_cols,kernel="polynomial",degree=3,type="C-classification") {
  features = list()
  models = list()
  for (i in 1:length(code_cols)) {
    train_data = trainfold$extract(code=code,code_col=code_cols[i],label_name="Y",posneg=factor(c("+","-")))
    n_features = ncol(train_data)-1
    feature_names = colnames(train_data)[1:n_features] # paste("X",1:n_features,sep="")
    model = svm(train_data[,1:n_features],train_data[,ncol(train_data)],kernel=kernel,degree=degree,type=type)
    tmp_create_feature = function(model) {
      f = new.env()
      f$model = model
      f$calc = function(x) {
        colnames(x) = feature_names
        return(logistic((f$model$labels[1]-f$model$labels[2])*attr(predict(f$model,x,decision.values=TRUE),"decision.values")))
      }
      f$n_inputs = n_features
      f$n_outputs = 1
      return(f)
    }
    features[[i]] = tmp_create_feature(model)
  }
  return(features)
}

create_interpolator = function(x_vec,fx_vec) {
  res_func = function(z_vec) {
    lower_neighbour = vector()
    upper_neighbour = vector()
    for (i in 1:length(z_vec)) {
      upper_neighbour[i] = match(TRUE,z_vec[i]<=x_vec)
      if (is.na(upper_neighbour[i])) {
        upper_neighbour[i] = length(x_vec)
      }
      if ((upper_neighbour[i]==1) || (z_vec[i]>=x_vec[upper_neighbour[i]])) {
        lower_neighbour[i] = upper_neighbour[i]
      } else {
        lower_neighbour[i] = upper_neighbour[i]-1
      }
    }
    lower_dist = z_vec-x_vec[lower_neighbour]
    upper_dist = x_vec[upper_neighbour]-z_vec
    c_lower = upper_dist/(lower_dist+upper_dist)
    c_lower[is.nan(c_lower)] = 0.5
    c_lower[abs(c_lower)==Inf] = 0.5
    c_upper = lower_dist/(lower_dist+upper_dist)
    c_upper[is.nan(c_upper)] = 0.5
    c_upper[abs(c_upper)==Inf] = 0.5
    res = c_lower*fx_vec[lower_neighbour]+c_upper*fx_vec[upper_neighbour]
    return(res)
  }
  return(res_func)
}

moving_average = function(x,binsize=3) {
  a = 0*x
  for (i in 1:length(x)) {
    start = i-floor(binsize/2)
    end = start+binsize-1
    if (start<1) {
      start = 1
    }
    if (end>length(x)) {
      end = length(x)
    }
    a[i] = mean(x[start:end])
  }
  return(a)
}

create_my_regression_learner = function(binsize=3) {
  return(function(x,y) create_interpolator(x,moving_average(y,binsize=binsize)))
}

create_np_regression_learner = function(minvalue=0.001,maxvalue=0.999,bwgiven=NA,bwmethod='cv.ls',bwtype='generalized_nn',itmax=10000,ckertype='epanechnikov',regtype="ll") {
  return(function(x,y) {
    if (is.na(bwgiven)) {
      bw = npregbw(formula=y~x,data=data.frame(x=x,y=y),bwmethod=bwmethod,bwtype=bwtype,ckertype=ckertype,regtype=regtype,bandwidth.compute=TRUE,itmax=itmax)
      print(paste('cv bandwidth=',bw$bandwidth))
    } else {
      bw = npregbw(formula=y~x,data=data.frame(x=x,y=y),bwmethod=bwmethod,bwtype=bwtype,ckertype=ckertype,regtype=regtype,bandwidth.compute=FALSE,bws=bwgiven)
      print(paste('given bandwidth=',bw$bandwidth))
    }
    model = npreg(formula=y~x,data=data.frame(x=x,y=y),bws=bw,ckertype=ckertype,regtype=regtype)
    return(function(z) pmax(minvalue,pmin(maxvalue,predict(model,data=data.frame(x=x,y=y),newdata=data.frame(x=z)))))
  })
}

create_ksmooth_regression_learner = function(bandwidth,minvalue,maxvalue) {
  return(function(x,y) function(z) pmax(minvalue,pmin(maxvalue,ksmooth(x,y,kernel="normal",bandwidth=bandwidth,x.points=z)$y))[rank(z)])
}

calibration_map_learner = function(cal_fold,feature_col,target_col,regression_learner) {
  d = cal_fold$extract()
  o = order(d[,feature_col])
  x = d[o,feature_col]
  y = d[o,target_col]
  return(regression_learner(x,y))
}

reliability_from_posterior_learner = function(cal_fold,feature_col,calib_col,post_col,regression_learner) {
  d = cal_fold$extract()
  o = order(d[,feature_col])
  cal = new.env()
  g = d[o,feature_col]
  m = d[o,calib_col]
  q = d[o,post_col]
  r = 1-((q-m)^2)/(m*(1-m))
  return(regression_learner(g,r))
}

# DIANA clustering algorithm from:
#[BOOK] Finding groups in data: an introduction to cluster analysis
#L Kaufman, PJ Rousseeuw - 2009 - books.google.com
# Adapted to give fixed-size clusters

diana_step = function(d,cluster_size,to_clust,labels) {
  if (length(labels)*cluster_size!=length(to_clust)) {
    print("Warning! Cannot produce equal-sized clusters of the required size")
  }
  if (length(to_clust)<=cluster_size) {
    return(rep(labels[1],length(to_clust)))
  }
  n_parts = floor(length(labels)/2)
  half_size = cluster_size*n_parts
  if (length(half_size)>length(to_clust)) {
    stop("This should never happen! length(half_size)>length(to_clust)")
  }
  dist_to_others = apply(d[to_clust,to_clust],1,sum)
  cluster = vector()
  remainder = to_clust
  dist_rem_to_clust = rep(0,length(remainder))
  dist_rem_to_rem = apply(d[remainder,remainder],1,sum)
  while (length(cluster)<half_size) {
    best_in_rem = which.max(dist_rem_to_rem/length(remainder)-dist_rem_to_clust/max(1,length(cluster)))
    best = remainder[best_in_rem]
    dist_diff = d[remainder,best]
    dist_rem_to_clust = dist_rem_to_clust+dist_diff
    dist_rem_to_rem = dist_rem_to_rem-dist_diff
    cluster = c(cluster,best)
    without_best = setdiff(1:length(remainder),best_in_rem)
    remainder = remainder[without_best]
    dist_rem_to_clust = dist_rem_to_clust[without_best]
    dist_rem_to_rem = dist_rem_to_rem[without_best]
  }
  cluster = sort(cluster)
  remainder = sort(remainder)
  res = vector()
  res[cluster] = diana_step(d,cluster_size,cluster,labels[1:n_parts])
  res[remainder] = diana_step(d,cluster_size,remainder,labels[(n_parts+1):length(labels)])
  return(res[to_clust])
}

eucl_distance_matrix_calculator = function(x) {
  return(as.matrix(dist(x)))
}

create_segmentwise_diana_learner = function(approxbinsize=200,clustsize=10,dist_calc=eucl_distance_matrix_calculator) {
  return(function(fold,orig_feature_cols,cal_feature_col) {
    d = fold$extract()
    o = order(d[,cal_feature_col])
    n = nrow(d)
    nbins = round(n/approxbinsize)
    if (nbins==0) {
      nbins = 1
    }
    nclustperbin = floor(n/nbins/clustsize)
    if (nclustperbin==0) {
      clustsize = n/nbins
      nclustperbin = 1
    }
    binsize = nclustperbin*clustsize
    nbins = floor(n/binsize)
    m = clustsize
    clust = vector()
    for (i in 1:nbins) {
      start = (i-1)*binsize+1
      end = start+binsize-1
      distance_matrix = dist_calc(d[o[start:end],orig_feature_cols])
      clust[o[start:end]] = diana_step(distance_matrix,clustsize,1:binsize,(1:nclustperbin)+(i-1)*nclustperbin)
      if (FALSE) {
        df = d[o[start:end],orig_feature_cols]
        df$Y1 = factor((clust-1)%/%7)
        df$Y2 = factor(1+clust%%7)
        p = ggplot(df)
        p = p + geom_point(aes(x=x1,y=x2,colour=Y2,shape=Y1),size=3)
        p = p + scale_colour_brewer(type="qual",palette="Set1")
        print(p)
      }
    }
    remain = n-nbins*binsize
    nclust = floor(remain/clustsize)
    if (nclust>=1) {
      start = nbins*binsize+1
      end = start+nclust*clustsize-1
      distance_matrix = dist_calc(d[o[start:end],orig_feature_cols])
      clust[o[start:end]] = diana_step(distance_matrix,clustsize,1:(nclust*clustsize),(1:nclust)+nclustperbin*nbins)
    }
    remain = n-nbins*binsize-nclust*clustsize
    if (remain>=1) {
      start = nbins*binsize+nclust*clustsize+1
      end = n
      clust[o[start:end]] = nclustperbin*nbins+nclust+1
    }
    return(clust)
  })
}

reliability_from_labels_learner = function(cal_fold,orig_feature_cols,feature_col,labelcal_col,label_col,regression_learner,clustering_learner) {
  d = cal_fold$extract()
  o = order(d[,feature_col])
  clust = clustering_learner(cal_fold,orig_feature_cols,feature_col)
  n_clust = max(clust)
  mu = vector()
  mf = vector()
  z = vector()
  m = vector()
  r = vector()
  for (i in 1:n_clust) {
    m[i] = sum(clust==i)
  }
  for (i in which(m==m[1])) {
    mf[i] = mean(d[clust==i,feature_col])
    mu[i] = mean(d[clust==i,labelcal_col])
    z[i] = sum(d[clust==i,label_col])
    r[i] = 1 - ((z[i]-m[i]*mu[i])^2-m[i]*mu[i]*(1-mu[i]))/(m[i]*(m[i]-1)*mu[i]*(1-mu[i]))
  }
  return(regression_learner(mf,r))
}

mae_func = function(x) apply(abs(x[,1:(ncol(x)/2),drop=FALSE]-x[,(ncol(x)/2+1):ncol(x),drop=FALSE]),1,mean)

mse_func = function(x) apply((x[,1:(ncol(x)/2),drop=FALSE]-x[,(ncol(x)/2+1):ncol(x),drop=FALSE])^2,1,mean)

pea_func = function(x) apply(x,1,function(a) cor(a[1:(length(a)/2)],a[(length(a)/2+1):length(a)]))

neq_func = function(x) 1*(x[,1]!=x[,2])

calc_map_prediction = function(probs) {
  return(colnames(probs)[apply(probs,1,which.max)])
}

rmse = function(x,y) {
  return(sqrt(mean((x-y)^2)))
}

mae = function(x,y) {
  return(mean(abs((x-y)^1)))
}

err = function(x,y) {
  return(sum(x!=y)/length(x))
}

which.unique = function(x) match(unique(x),x)

plot_calrel = function(dm,prefix,ecoc,cols=1:ecoc$n_cols,calrel=".hrel") {
  p = ggplot(dm$mat)
  p = p + xlim(c(0,1))
  p = p + ylim(c(-0.02,1))
  for (i in cols) {
    code = paste(prefix,".model",i,sep="")
    rel = paste(code,calrel,sep="")
    p = p + geom_line(aes_string(x=code,y=rel))
    ymin = min(dm$mat[,rel])-0.02
    xmin = dm$mat[which.min(dm$mat[,rel]),code]
    p = p + annotate("text",x=xmin,y=ymin,label=paste(i))
  }
  print(p)
  return(p)
}

#-------------

calc_confusion_matrix2 = function(dataset,prediction,truth) {
  classes = levels(dataset[,label_col(dataset)])
  res = matrix(rep(0,length(classes)**2),nrow=length(classes))
  rownames(res) = classes
  colnames(res) = classes
  for (r in 1:length(prediction)) {
    row = which(classes==truth[r])
    col = which(classes==prediction[r])
    res[row,col] = res[row,col]+1
  }
  return(res)
}

calc_confusion_matrix = function(dataset,prediction,rows) {
  classes = levels(dataset[,label_col(dataset)])
  res = matrix(rep(0,length(classes)**2),nrow=length(classes))
  rownames(res) = classes
  colnames(res) = classes
  for (r in rows) {
    row = which(classes==dataset[r,label_col(dataset)])
    col = which(classes==prediction[r])
    res[row,col] = res[row,col]+1
  }
  return(res)
}

print_confusion_matrix = function(text,mat) {
  bars = paste(rep("=",1+nchar(text)),collapse="")
  cat(paste("\n",bars,"\n",text,":\n",bars,"\n",sep=""))
  print(mat)
  cat(paste("\nNumber of instances:",sum(mat),"\n"))
  cat(paste("Number of errors   :",sum(mat)-sum(diag(mat)),"\n"))
}




# Part of software for the ECML'14 submission "Reliability maps: A tool to enhance probability estimates and improve classification accuracy"
# Author: Meelis Kull, meelis.kull@bristol.ac.uk, meelis.kull@gmail.com
# Written February-April 2014.
# All rights reserved. 

args = commandArgs(trailingOnly = TRUE)
seed = as.numeric(args[1])
train_size = as.numeric(args[2])
bwcal = as.numeric(args[3])
bwrel = as.numeric(args[4])
approxbinsize = as.numeric(args[5])
print(paste('seed',seed))

# seed = 1
# train_size = 400
# approxbinsize = 200
# bwcal = 0.01
# bwrel = 0.1

source("ecml14func.r")

test_size = 10000

fileprefix = paste("wang.tr",train_size,".bwc",bwcal,".bwr",bwrel,".bs",approxbinsize,".seed",seed,sep="")

sink(paste(fileprefix,".stdout.txt",sep=""), append=FALSE, split=TRUE)

reporting = function(test_fold,prefix,plot=TRUE) {
  dm = test_fold$datamatrix
  kl50 = median(dm$mat[test_fold$rows,paste(prefix,"kl",sep=".")])
  kl95 = as.numeric(quantile(dm$mat[test_fold$rows,paste(prefix,"kl",sep=".")],0.95))
  bene = mean(dm$mat[test_fold$rows,paste(prefix,"bene",sep=".")])
  nene = mean(dm$mat[test_fold$rows,paste(prefix,"nene",sep=".")])
  nenehcalcvar    = mean(dm$mat[test_fold$rows,paste(prefix,"nenehcalcvar",sep=".")])
  nenehcalcsvar    = mean(dm$mat[test_fold$rows,paste(prefix,"nenehcalcsvar",sep=".")])
  benehcalcrelvar = mean(dm$mat[test_fold$rows,paste(prefix,"benehcalcrelvar",sep=".")])
  nenehcalcrelvar = mean(dm$mat[test_fold$rows,paste(prefix,"nenehcalcrelvar",sep=".")])
  ss = strsplit(prefix,"\\.")[[1]]
  res = data.frame(seed=seed,trainsize=train_size,bwcal=bwcal,bwrel=bwrel,binsize=approxbinsize,name=prefix,modeltype=ss[2],cvtype=ss[3],code=ss[4],kl50=kl50,kl95=kl95,bene=bene,nene=nene, nenehcalcvar  = nenehcalcvar ,nenehcalcsvar  = nenehcalcsvar , benehcalcrelvar= benehcalcrelvar , nenehcalcrelvar= nenehcalcrelvar )
  if (plot) {
    nnn = paste(prefix,c("cal","low","upp","hcal","hlow","hupp"),sep=".")
    p = ggplot(dm$mat)
    p = p + geom_ribbon(aes_string(x=prefix,ymin=nnn[2],ymax=nnn[3]),fill="yellow")
    p = p + geom_line(  aes_string(x=prefix,y=nnn[1]))
    p = p + geom_line(  aes_string(x=prefix,y=nnn[2]))
    p = p + geom_line(  aes_string(x=prefix,y=nnn[3]))
    p = p + geom_point( aes_string(x=prefix,y=nnn[4]))
    p = p + geom_point( aes_string(x=prefix,y=nnn[5]))
    p = p + geom_point( aes_string(x=prefix,y=nnn[6]))
    p = p + ylim(c(0,1))
    print(p)
  }
  print(t(res))
  return(res)
}

set.seed(seed)
total_size =  test_size+train_size
uniwang = create_uniwang()
dm = create_datamatrix("uniwang")
print('generating data...')
dm$add_instances(rbind(uniwang$generate(test_size,stratified=TRUE),uniwang$generate(train_size,stratified=TRUE)))
task = create_task(dm,uniwang$features,uniwang$label)
ecoc_1vsR = create_ecoc(task,type="1vsRest")
ecoc_1vs1 = create_ecoc(task,type="1vs1")
tmp_mat = 2*abs(ecoc_1vs1$mat)-1
colnames(tmp_mat) = paste(c(12,13,14,15,23,24,25,34,35,45),"vsR",sep="")
ecoc_2vsR = create_ecoc(task,type="custom",custom_mat=tmp_mat)
ecoc_all = create_ecoc(task,type="custom",custom_mat=cbind(ecoc_1vsR$mat,ecoc_2vsR$mat))

print('calculating ecoc labels...')
dm$add_feature(task$labels,ecoc_all$calc_labels,ecoc_all$n_cols,paste(task$labels,colnames(ecoc_all$mat),sep="."))
post_1vR_cols = paste("post",colnames(ecoc_1vsR$mat),sep=".")
print('calculating posterior...')
dm$add_feature(task$features,uniwang$posterior,uniwang$n_classes,post_1vR_cols)
add_cols = function(cols) { return(function(x) {return(apply(x[,cols],1,sum))}) }
for (i in 1:(ecoc_all$n_classes-1)) {
  for (j in (i+1):(ecoc_all$n_classes)) {
    dm$add_feature(post_1vR_cols,add_cols(c(i,j)),1,paste("post.",i,j,"vsR",sep=""))
  }
}
dm$add_feature(post_1vR_cols,function(x) ecoc_all$classes_factor[apply(x,1,which.max)],1,"post.map")

print('creating folds')
folds = list(test=create_fold(dm,1:test_size,name="test"),train=create_fold(dm,(test_size+1):total_size,name="train"))
#folds = create_fold(dm)$split(4,proportions=c(0.4,2,10,10),stratified=TRUE,task=task,fold_names=c("tr400","tr2k","tr10k","test"))
test_fold = folds$test
test_posterior = dm$mat[test_fold$rows,post_1vR_cols]
test_labels = dm$mat[test_fold$rows,task$labels]

map_func = function(x) ecoc_all$classes_factor[apply(x,1,which.max)]

bri_func = function(x) {
  k = ncol(x)-1
  se = 0
  for (i in 1:k) {
    se = se+1/2*(x[,i]-1*(x[,k+1]==ecoc_all$classes_factor[i]))^2
  }
  return(se)
}

print('creating plan')
plan = data.frame()
#train_sizes = data.frame(names=c("tr400","tr2k","tr10k"),sizes=c(400,2000,10000))
#train_sizes$names = as.character(train_sizes$names)
#model_types = data.frame(names=c("svm")) #,learners=c(train_svm))
model_types = data.frame(names=c("lr")) #,learners=c(train_svm))
model_types$names = as.character(model_types$names)
#split_methods = paste(c(44,31,24,22))
split_methods = paste(c(44)) #,31,24,22))
modcal_splits = data.frame(t(matrix(
 c(1,1234,1234), nrow=3)))
# c(1,1234,1234,
#   2,123 ,4   ,
#   2,124 ,3   ,
#   2,134 ,2   ,
#   2,234 ,1   ,
#   3,12  ,1234,
#   3,34  ,1234,
#   4,12  ,34  ,
#   4,34  ,12  ),nrow=3)))
modcalrel_splits = cbind(modcal_splits,modcal_splits[,3])
colnames(modcalrel_splits) = c("method","mod","cal","rel")
all_splits = unique(as.numeric(as.matrix(modcalrel_splits[,2:4])))
rel_clust_learner = create_segmentwise_diana_learner(approxbinsize=approxbinsize,clustsize=10,dist_calc=eucl_distance_matrix_calculator)
ecoc_list = list(`1vsR`=ecoc_1vsR)
ecoc_list = list(all=ecoc_all)
ecoc_list = list(`1vsR`=ecoc_1vsR,`2vsR`=ecoc_2vsR,all=ecoc_all)
ecoc_cols_list = list(`1vsR`=1:5)
ecoc_cols_list = list(all=1:15)
ecoc_cols_list = list(`1vsR`=1:5,`2vsR`=6:15,all=1:15)
rmse_res = list()
mae_res = list()
err_res = list()
all_models_list = list()
report = data.frame()
report2 = data.frame()
for (i1 in 1:1) {
  name1 = "train"
  print(name1)
  cal_reg_learner = create_np_regression_learner(bwgiven=bwcal,ckertype='epanechnikov',bwtype='fixed',minvalue=1/1000,maxvalue=999/1000,regtype="ll",itmax=100)
  rel_reg_learner = create_np_regression_learner(bwgiven=bwrel,ckertype='epanechnikov',bwtype='fixed',minvalue=1/1000,maxvalue=999/1000,regtype="ll",itmax=100)
  train_fold = folds[[name1]]
  cv_folds = train_fold$split(4,proportions=c(1/4,1/4,1/4,1/4),stratified=TRUE,task=task,fold_names=1:4)
  split_folds = list()
  for (sp in all_splits) {
    rows = vector()
    for (cv_fold_name in strsplit(paste(sp),NULL)[[1]]) {
      rows = c(rows,cv_folds[[cv_fold_name]]$rows)
    }
    split_folds[[paste(sp)]] = create_fold(dm,rows=rows,name=paste(sp))
  }
  for (i2 in 1:nrow(model_types)) {
    name2 = paste(name1,model_types$names[i2],sep=".")
    print(name2)
    for (i3 in 1:nrow(modcalrel_splits)) {
      name3 = paste(name2,".s",i3,sep="")
      print(name3)
      mod_fold = split_folds[[paste(modcalrel_splits$mod[i3])]]
      cal_fold = split_folds[[paste(modcalrel_splits$cal[i3])]]
      rel_fold = split_folds[[paste(modcalrel_splits$rel[i3])]]
      if (model_types$names[i2]=='svm') {
        mod = train_svm(mod_fold,task,ecoc_all)
      } else if (model_types$names[i2]=='lr') {
        mod = train_lr(mod_fold,task,ecoc_all)
      } else { stop('unknown model class') }
      name4hcalvec = vector()
      name4hrelvec = vector()
      name4calvec = vector()
      name4relvec = vector()
      name4cvarvec = vector()
      name4csvarvec = vector()
      name4crelvarvec = vector()
      for (i4 in 1:length(mod)) {
	print(i4)
	code_name = colnames(ecoc_all$mat)[i4]
        name4 = paste(name3,code_name,sep=".")
	print(name4)
	name4hcal = paste(name4,"hcal",sep=".")
	name4hrel = paste(name4,"hrel",sep=".")
	name4hvar = paste(name4,"hvar",sep=".")
	name4hlow = paste(name4,"hlow",sep=".")
	name4hupp = paste(name4,"hupp",sep=".")
	name4hcalvec[i4] = name4hcal
	name4hrelvec[i4] = name4hrel
	name4cal = paste(name4,"cal",sep=".")
	name4rel = paste(name4,"rel",sep=".")
	name4var = paste(name4,"var",sep=".")
	name4low = paste(name4,"low",sep=".")
	name4upp = paste(name4,"upp",sep=".")
	name4calvec[i4] = name4cal
	name4relvec[i4] = name4rel
	name4cvar = paste(name4,"cvar",sep=".")
	name4csvar = paste(name4,"csvar",sep=".")
	name4crelvar = paste(name4,"crelvar",sep=".")
	name4cvarvec[i4] = name4cvar
	name4csvarvec[i4] = name4csvar
	name4crelvarvec[i4] = name4crelvar
	
	name4kl = paste(name4,"kl",sep=".")
	name4bene            = paste(name4,"bene",sep=".")
	name4nene            = paste(name4,"nene",sep=".")
	name4nenehcalcvar    = paste(name4,"nenehcalcvar",sep=".")
	name4nenehcalcsvar    = paste(name4,"nenehcalcsvar",sep=".")
	name4benehcalcrelvar = paste(name4,"benehcalcrelvar",sep=".")
	name4nenehcalcrelvar = paste(name4,"nenehcalcrelvar",sep=".")

        label_col = paste(task$labels,code_name,sep=".")

	print("calculating model output...")
	all_models_list[[name4]] = mod[[i4]]
        dm$add_feature(task$features,mod[[i4]]$calc,1,name4)

        print("calibration map from labels...")
	cal_map = calibration_map_learner(cal_fold,name4,label_col,cal_reg_learner)
	dm$add_feature(name4,cal_map,1,name4hcal)

        print("reliability map from labels...")
	rel_map = reliability_from_labels_learner(rel_fold,task$features,name4,name4hcal,label_col,rel_reg_learner,rel_clust_learner)
        dm$add_feature(name4,rel_map,1,name4hrel)
      
        w = mod[[i4]]$model$coefficients
        code = ecoc_all$mat[,code_name]
	print("true cal...")
        dm$add_feature(name4,create_true_cal_map(w,uniwang,code),1,name4cal)
	print("true rel...")
        dm$add_feature(name4,create_true_rel_map(w,uniwang,code),1,name4rel)
	print("confidence intervals...")
	dm$add_feature(c(name4cal,name4rel),function(x) (1-x[,2])*x[,1]*(1-x[,1]),1,name4var)
	dm$add_feature(c(name4cal,name4var),function(x) beta_lower(x[,1],sqrt(x[,2]),q=0.05),1,name4low)
	dm$add_feature(c(name4cal,name4var),function(x) beta_lower(x[,1],sqrt(x[,2]),q=0.95),1,name4upp)
	dm$add_feature(c(name4hcal,name4hrel),function(x) (1-x[,2])*x[,1]*(1-x[,1]),1,name4hvar)
	dm$add_feature(c(name4hcal,name4hvar),function(x) beta_lower(x[,1],sqrt(x[,2]),q=0.05),1,name4hlow)
	dm$add_feature(c(name4hcal,name4hvar),function(x) beta_lower(x[,1],sqrt(x[,2]),q=0.95),1,name4hupp)

	print("const var...")
	cvar = mean(dm$mat[rel_fold$rows,name4hvar])
	csvar = (mean(sqrt(dm$mat[rel_fold$rows,name4hvar])))^2
	crel = mean(dm$mat[rel_fold$rows,name4hrel])
	dm$add_feature(c(),function(x) cvar,1,name4cvar)
	dm$add_feature(c(),function(x) csvar,1,name4csvar)
	dm$add_feature(name4hcal,function(x) (1-crel)*x*(1-x),1,name4crelvar)
	
	print("kullback-leibler...")
	dm$add_feature(c(name4cal,name4var,name4hcal,name4hvar),function(x) kl_dist(x[,1],x[,2],x[,3],x[,4]),1,name4kl)
	print("beta energy...")
	dm$add_feature(c(name4cal,name4var,name4hcal,name4hvar),function(x) half_beta_energy(x[,1],x[,2],x[,3],x[,4],niter=1000),1,name4bene)
	print(mean(dm$mat[[name4bene]]))
	print("norm energy...")
	dm$add_feature(c(name4cal,name4var,name4hcal,name4hvar),function(x) half_norm_energy(x[,1],x[,2],x[,3],x[,4],niter=1000),1,name4nene)
	print(mean(dm$mat[[name4nene]]))
	print("norm energy hcal-cvar...")
	dm$add_feature(c(name4cal,name4var,name4hcal,name4cvar),function(x) half_norm_energy(x[,1],x[,2],x[,3],x[,4],niter=1000),1,name4nenehcalcvar)
	print(mean(dm$mat[[name4nenehcalcvar]]))
	print("norm energy hcal-csvar...")
	dm$add_feature(c(name4cal,name4var,name4hcal,name4csvar),function(x) half_norm_energy(x[,1],x[,2],x[,3],x[,4],niter=1000),1,name4nenehcalcsvar)
	print(mean(dm$mat[[name4nenehcalcsvar]]))
	print("beta energy hcal-crelvar...")
	dm$add_feature(c(name4cal,name4var,name4hcal,name4crelvar),function(x) half_beta_energy(x[,1],x[,2],x[,3],x[,4],niter=1000),1,name4benehcalcrelvar)
	print(mean(dm$mat[[name4benehcalcrelvar]]))
	print("norm energy hcal-crelvar...")
	dm$add_feature(c(name4cal,name4var,name4hcal,name4crelvar),function(x) half_norm_energy(x[,1],x[,2],x[,3],x[,4],niter=1000),1,name4nenehcalcrelvar)
	print(mean(dm$mat[[name4nenehcalcrelvar]]))
	if (any(is.na(dm$mat[,c(name4bene,name4nene,name4nenehcalcvar,name4nenehcalcsvar,name4benehcalcrelvar,name4nenehcalcrelvar)]))) {
	  stop("this error should not happen")
	}

	report = rbind(report,reporting(test_fold,name4,plot=FALSE))
      }
      for (code_name in names(ecoc_list)) {
        print(code_name)
	ecoc = ecoc_list[[code_name]]
        print("ls-ecoc with hcal...")
	name_prob = paste(name3,code_name,"lsecoc",ecoc$classes_factor,sep=".")
	name_map = paste(name3,code_name,"lsecoc","map",sep=".")
        dm$add_feature(name4hcalvec[ecoc_cols_list[[code_name]]],ecoc$ls_ecoc,5,name_prob)
        dm$add_feature(name_prob,map_func,1,name_map)

        print("ls-ecoc-v with cvar...")
	name_prob = paste(name3,code_name,"lsecoccvar",ecoc$classes_factor,sep=".")
	name_map = paste(name3,code_name,"lsecoccvar","map",sep=".")
        dm$add_feature(c(name4hcalvec[ecoc_cols_list[[code_name]]],name4cvarvec[ecoc_cols_list[[code_name]]]),function(x) ecoc$ls_ecoc_v(x[,1:(ncol(x)/2)],x[,(ncol(x)/2+1):(ncol(x))]),5,name_prob)
        dm$add_feature(name_prob,map_func,1,name_map)

        print("ls-ecoc-v with csvar...")
	name_prob = paste(name3,code_name,"lsecoccsvar",ecoc$classes_factor,sep=".")
	name_map = paste(name3,code_name,"lsecoccsvar","map",sep=".")
        dm$add_feature(c(name4hcalvec[ecoc_cols_list[[code_name]]],name4csvarvec[ecoc_cols_list[[code_name]]]),function(x) ecoc$ls_ecoc_v(x[,1:(ncol(x)/2)],x[,(ncol(x)/2+1):(ncol(x))]),5,name_prob)
        dm$add_feature(name_prob,map_func,1,name_map)

        print("ls-ecoc-v with crelvar...")
	name_prob = paste(name3,code_name,"lsecoccrelvar",ecoc$classes_factor,sep=".")
	name_map = paste(name3,code_name,"lsecoccrelvar","map",sep=".")
        dm$add_feature(c(name4hcalvec[ecoc_cols_list[[code_name]]],name4crelvarvec[ecoc_cols_list[[code_name]]]),function(x) ecoc$ls_ecoc_v(x[,1:(ncol(x)/2)],x[,(ncol(x)/2+1):(ncol(x))]),5,name_prob)
        dm$add_feature(name_prob,map_func,1,name_map)

        print("ls-ecoc-r with hrel...")
	name_prob = paste(name3,code_name,"lsecocr",ecoc$classes_factor,sep=".")
	name_map = paste(name3,code_name,"lsecocr","map",sep=".")
        dm$add_feature(c(name4hcalvec[ecoc_cols_list[[code_name]]],name4hrelvec[ecoc_cols_list[[code_name]]]),function(x) ecoc$ls_ecoc_r(x[,1:(ncol(x)/2)],x[,(ncol(x)/2+1):(ncol(x))]),5,name_prob)
        dm$add_feature(name_prob,map_func,1,name_map)

        print("ls-ecoc-r with rel...")
	name_prob = paste(name3,code_name,"lsecocri",ecoc$classes_factor,sep=".")
	name_map = paste(name3,code_name,"lsecocri","map",sep=".")
        dm$add_feature(c(name4calvec[ecoc_cols_list[[code_name]]],name4relvec[ecoc_cols_list[[code_name]]]),function(x) ecoc$ls_ecoc_r(x[,1:(ncol(x)/2)],x[,(ncol(x)/2+1):(ncol(x))]),5,name_prob)
        dm$add_feature(name_prob,map_func,1,name_map)

      }
    }
    for (i3 in 1:length(split_methods)) {
      name3 = paste(name2,".m",i3,sep="")
      print(name3)
      splvec = which(modcalrel_splits$method==i3)
      from3 = paste(name2,".s",splvec,sep="")
      for (code_name in names(ecoc_list)) {
        name4 = paste(name3,code_name,sep=".")
        from4 = paste(from3,code_name,sep=".")
	print(name4)
	ecoc = ecoc_list[[code_name]]
	for (n5 in c('lsecoc','lsecoccvar','lsecoccsvar','lsecoccrelvar','lsecocr','lsecocri')) {
          name5 = paste(name4,n5,sep=".")
          from5 = paste(from4,n5,sep=".")
	  print(name5)
	  name6vec = vector()
          for (i6 in 1:ecoc_all$n_classes) {
	    name6 = paste(name5,ecoc$classes_factor[i6],sep=".")
	    name6vec[i6] = name6
	    from6 = paste(from5,ecoc$classes_factor[i6],sep=".")
	    if (length(from6)>1) {
	      dm$add_feature(from6,function(x) apply(x,1,mean),1,name6)
	    } else {
	      dm$add_feature(from6,function(x) x,1,name6)
	    }
	  }
	  name5map = paste(name5,"map",sep=".")
	  dm$add_feature(name6vec,map_func,1,name5map)
	  mat_prob = dm$mat[test_fold$rows,name6vec]
	  mat_map = dm$mat[test_fold$rows,name5map]
	  name5mse = paste(name5,"mse",sep=".")
	  name5mae = paste(name5,"mae",sep=".")
	  name5err = paste(name5,"err",sep=".")
	  name5bri = paste(name5,"bri",sep=".")
	  name5bbri = paste(name5,"bbri",sep=".")
	  name5pea = paste(name5,"pea",sep=".")
	  name5berr = paste(name5,"berr",sep=".")
	  dm$add_feature(c(name6vec,post_1vR_cols),mse_func,1,name5mse)
	  dm$add_feature(c(name6vec,post_1vR_cols),mae_func,1,name5mae)
	  dm$add_feature(c(name6vec,post_1vR_cols),pea_func,1,name5pea)
	  dm$add_feature(c(name6vec,task$labels),bri_func,1,name5bri)
	  dm$add_feature(c(name6vec,"post.map"),bri_func,1,name5bbri)
	  dm$add_feature(c(name5map,task$labels),neq_func,1,name5err)
	  dm$add_feature(c(name5map,"post.map"),neq_func,1,name5berr)
          ss = strsplit(name5,"\\.")[[1]]
          report2 = rbind(report2,data.frame(seed=seed,trainsize=train_size,bwcal=bwcal,bwrel=bwrel,binsize=approxbinsize,name=name5,modeltype=ss[2],cvmethod=ss[3],codemat=ss[4],lsmethod=ss[5],
	    rmse=sqrt(mean(dm$mat[test_fold$rows,name5mse])),
	    mae=mean(dm$mat[test_fold$rows,name5mae]),
	    pea=mean(dm$mat[test_fold$rows,name5pea]),
	    bri=mean(dm$mat[test_fold$rows,name5bri]),
	    bbri=mean(dm$mat[test_fold$rows,name5bbri]),
	    err=mean(dm$mat[test_fold$rows,name5err]),
	    berr=mean(dm$mat[test_fold$rows,name5berr])
	    ))
	}
      }
    }
    print('---')
  }
}


write.csv(report,file=paste(fileprefix,".energy.csv",sep=""),quote=FALSE,row.names=FALSE)
write.csv(report2,file=paste(fileprefix,".ecoc.csv",sep=""),quote=FALSE,row.names=FALSE)

# save.image(paste("~/ecml14experiments/2014_04_14_08_21_25_wang.",fileprefix,".Rdata",sep=""))


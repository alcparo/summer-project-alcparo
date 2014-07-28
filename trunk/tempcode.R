args = commandArgs(trailingOnly = TRUE)
dataset = as.character(args[1])
bwcal = as.numeric(args[2])
bwrel = as.numeric(args[3])
approxbinsize = as.numeric(args[4])

setwd("research/scripts/myresearch");
dataset="diabetes_num";
bwcal = 0.2;
bwrel = 0.2;
approxbinsize = 200;

fileprefix = paste(dataset,".bwc",bwcal,".bwr",bwrel,".bs",approxbinsize,sep="")
polydegrees = list()
polydegrees[["page-blocks"]] = 1
polydegrees[["yeast"      ]] = 3
polydegrees[["segment"    ]] = 1
polydegrees[["sat"        ]] = 1
polydegrees[["vehicle"    ]] = 1
polydegrees[["shuttle"    ]] = 2
polydegrees[["data_banknote_authentication"    ]] = 1
polydegrees[["diabetes_num"    ]] = 1

source("binary_functions.r")

sink(paste(fileprefix,".stdout.txt",sep="."), append=FALSE, split=TRUE)

set.seed(1)
if (dataset=="sat") {
  arff = read.csv("sat.csv",sep=" ",header=FALSE)
  } else if (dataset=="shuttle") {
    arff = read.csv("shuttle.tst",sep=" ",header=FALSE)
    } else {
      arff = read.arff(file=paste(dataset,"arff",sep="."))
    }

    colnames(arff) = c(paste('f',1:(ncol(arff)-1),sep=""),'class') 
# arff$class = as.numeric(as.factor(arff$class)) #### as.factor changing the values from (0,1) to (1,2)

dm = create_datamatrix(dataset)
dm$add_instances(arff)
task = create_task(dm,colnames(arff)[1:(ncol(arff)-1)],colnames(arff)[ncol(arff)])

print(task$n_classes);
if (task$n_classes<2) {
  stop("does not work for less than 2 classes") ### new
  } else if (task$n_classes>15) {
    stop("does not work for more than 15 classes")
    } else if(task$n_classes==2){ ### new
     print('BINARY'); 
     ecoc = create_ecoc(task,type="binary")
     } else if (task$n_classes==3) {
      ecoc = create_ecoc(task,type="1vsR")
      } else {
        if (task$n_classes==4) {
          mat = read.csv("bch_n7_k4_t1_sum3_dist4_t.csv")
          } else {
            mat = read.csv("bch_n15_k5_t3_sum7_dist8_t.csv")
          }
          mat = 2*mat-1
          colnames(mat) = paste("code",1:ncol(mat),sep="")
          matrows = sample(1:task$n_classes,task$n_classes)
          ecoc = create_ecoc(task,type="custom",custom_mat=as.matrix(mat[matrows,]))
        }

        map_func = function(x) ecoc$classes_factor[apply(x,1,which.max)]

        print('calculating ecoc labels...')
        dm$add_feature(task$labels,ecoc$calc_labels,ecoc$n_cols,paste(task$labels,colnames(ecoc$mat),sep="."))

        print('creating plan')
        plan = data.frame()
# model_types = data.frame(names=c("svm")) #,learners=c(train_svm))
# model_types = data.frame(names=c("lr")) #,learners=c(train_svm))
# model_types = data.frame(names=c("svmpol","svmrad","lr")) #,learners=c(train_svm))
model_types = data.frame(names=c("svmpol","lr")) #,learners=c(train_svm)) ##### ADD THE LEARNING ALGORITHMS
modcal_splits = data.frame(t(matrix(c(1,1234,1234),nrow=3)))
#   2,123 ,4   ,
#   2,124 ,3   ,
#   2,134 ,2   ,
#   2,234 ,1   ),nrow=3)))
#   3,123 ,1234,
#   3,124 ,1234,
#   3,134 ,1234,
#   3,234 ,1234,
#   4,12  ,1234,
#   4,34  ,1234,
#   5,12  ,34  ,
#   5,34  ,12  ),nrow=3)))
modcalrel_splits = cbind(modcal_splits,modcal_splits[,3])
colnames(modcalrel_splits) = c("method","mod","cal","rel")
n_methods = max(modcalrel_splits$method)
all_splits = unique(as.numeric(as.matrix(modcalrel_splits[,2:4])))
cal_reg_learner = create_np_regression_learner(bwgiven=bwcal,ckertype='epanechnikov',bwtype='fixed',minvalue=1/1000,maxvalue=999/1000,regtype="ll",itmax=100)
rel_reg_learner = create_np_regression_learner(bwgiven=bwrel,ckertype='epanechnikov',bwtype='fixed',minvalue=1/1000,maxvalue=999/1000,regtype="ll",itmax=100)
rel_clust_learner = create_segmentwise_diana_learner(approxbinsize=approxbinsize,clustsize=10,dist_calc=eucl_distance_matrix_calculator)
err_res = list()

print('creating folds')
folds = create_fold(dm)$split(10,proportions=rep(1,10),stratified=TRUE,task=task,fold_names=paste(1:10))

n_cols_before_cv = dm$m


# cv=1

for (cv in 1:10) {

  name1 = paste("cv",cv,sep="")

  test_fold = folds[[cv]]
  train_fold = create_fold(dm,setdiff(1:dm$n,test_fold$rows),name="test")
  test_labels = dm$mat[test_fold$rows,task$labels]

  # features normalization from 0..1
  norm_features = paste(name1,task$features,sep=".")
  for (i in 1:length(task$features)) {
    ff = dm$mat[train_fold$rows,task$features[i]]
    mm = mean(ff)
    ss = sqrt(var(ff))
    if (ss==0) {
      ss = 1
    }
    dm$add_feature(task$features[i],function(x) (x-mm)/ss,1,norm_features[i])
  }
  modtask = create_task(dm,norm_features,task$labels)
  ecoc$task = modtask
  
  cv_folds = train_fold$split(4,proportions=c(1/4,1/4,1/4,1/4),stratified=TRUE,task=modtask,fold_names=1:4)
  split_folds = list()
  for (sp in all_splits) {  ## ????????????? just one split fold?
    rows = vector()
    for (cv_fold_name in strsplit(paste(sp),NULL)[[1]]) {
      rows = c(rows,cv_folds[[cv_fold_name]]$rows)
    }
    split_folds[[paste(sp)]] = create_fold(dm,rows=rows,name=paste(sp)) #paste("cv",fold_name,sep="")) 
  }

#ADD ALL learning algorithms
  for (i2 in 1:nrow(model_types)) { #model_types = data.frame(names=c("svmpol","lr")) 
    name2 = paste(name1,model_types$names[i2],sep=".")
    print(name2)

    for (i3 in 1:nrow(modcalrel_splits)) {
      name3 = paste(name2,".s",i3,sep="")
      print(name3)

      mod_fold = split_folds[[paste(modcalrel_splits$mod[i3])]]
      cal_fold = split_folds[[paste(modcalrel_splits$cal[i3])]]
      rel_fold = split_folds[[paste(modcalrel_splits$rel[i3])]]
      if (model_types$names[i2]=='svmrad') { ##### ALGORITHM LEARNING TYPE
        mod = train_svm(mod_fold,modtask,ecoc,kernel="radial")
      } else if (model_types$names[i2]=='svmpol') {
        mod = train_svm(mod_fold,modtask,ecoc,kernel="polynomial",degree=polydegrees[[dataset]])
      } else if (model_types$names[i2]=='lr') {
        mod = train_lr(mod_fold,modtask,ecoc)
      } else { stop('unknown model class') } ### ADD ANOTHER LEARNING MODEL

      name4calvec = vector()
      name4relvec = vector()
      name4cvarvec = vector()
      name4csvarvec = vector()
      name4crelvarvec = vector()
      
      for (i4 in 1:length(mod)) {
        code_name = colnames(ecoc$mat)[i4]
        name4 = paste(name3,code_name,sep=".")
        print(name4)

        name4cal = paste(name4,"hcal",sep=".")
        name4var = paste(name4,"hvar",sep=".")
        name4rel = paste(name4,"hrel",sep=".")

        name4calvec[i4] = name4cal
        name4relvec[i4] = name4rel
        name4cvar = paste(name4,"cvar",sep=".")
        name4csvar = paste(name4,"csvar",sep=".")
        name4crelvar = paste(name4,"crelvar",sep=".")
        name4cvarvec[i4] = name4cvar
        name4csvarvec[i4] = name4csvar
        name4crelvarvec[i4] = name4crelvar

        label_col = paste(modtask$labels,code_name,sep=".")

        print("calculating model output...")
        dm$add_feature(modtask$features,mod[[i4]]$calc,1,name4)

        print("calibration map from labels...")
        cal_map = calibration_map_learner(cal_fold,name4,label_col,cal_reg_learner)
        dm$add_feature(name4,cal_map,1,name4cal)

        print("reliability map from labels...")
        rel_map = reliability_from_labels_learner(rel_fold,modtask$features,name4,name4cal,label_col,rel_reg_learner,rel_clust_learner)
        dm$add_feature(name4,rel_map,1,name4rel)

        print("const var...")
        dm$add_feature(c(name4cal,name4rel),function(x) (1-x[,2])*x[,1]*(1-x[,1]),1,name4var)
        cvar = mean(dm$mat[rel_fold$rows,name4var])
        csvar = (mean(sqrt(dm$mat[rel_fold$rows,name4var])))^2
        crel = mean(dm$mat[rel_fold$rows,name4rel])
        dm$add_feature(c(),function(x) cvar,1,name4cvar)
        dm$add_feature(c(),function(x) csvar,1,name4csvar)
        dm$add_feature(name4cal,function(x) (1-crel)*x*(1-x),1,name4crelvar)
 

      } ### END i4

      print("plotting calibration and reliability maps")
      ggsave(paste(dataset,name3,"cal.pdf",sep="."),plot_calrel(dm,name3,ecoc,calrel=".hcal"),width=5,height=5)
      ggsave(paste(dataset,name3,"rel.pdf",sep="."),plot_calrel(dm,name3,ecoc,calrel=".hrel"),width=5,height=5)
      print("ls-ecoc with hcal...")
      name_prob = paste(name3,"lsecoc",ecoc$classes_factor,sep=".")
      name_map = paste(name3,"lsecoc","map",sep=".")
      dm$add_feature(name4calvec,ecoc$ls_ecoc,modtask$n_classes,name_prob)
      dm$add_feature(name_prob,map_func,1,name_map)

      print("ls-ecoc-v with cvar...")
      name_prob = paste(name3,"lsecoccvar",ecoc$classes_factor,sep=".")
      name_map = paste(name3,"lsecoccvar","map",sep=".")
      dm$add_feature(c(name4calvec,name4cvarvec),function(x) ecoc$ls_ecoc_v(x[,1:(ncol(x)/2)],x[,(ncol(x)/2+1):(ncol(x))]),modtask$n_classes,name_prob)
      dm$add_feature(name_prob,map_func,1,name_map)

      print("ls-ecoc-v with csvar...")
      name_prob = paste(name3,"lsecoccsvar",ecoc$classes_factor,sep=".")
      name_map = paste(name3,"lsecoccsvar","map",sep=".")
      dm$add_feature(c(name4calvec,name4csvarvec),function(x) ecoc$ls_ecoc_v(x[,1:(ncol(x)/2)],x[,(ncol(x)/2+1):(ncol(x))]),modtask$n_classes,name_prob)
      dm$add_feature(name_prob,map_func,1,name_map)

      print("ls-ecoc-v with crelvar...")
      name_prob = paste(name3,"lsecoccrelvar",ecoc$classes_factor,sep=".")
      name_map = paste(name3,"lsecoccrelvar","map",sep=".")
      dm$add_feature(c(name4calvec,name4crelvarvec),function(x) ecoc$ls_ecoc_v(x[,1:(ncol(x)/2)],x[,(ncol(x)/2+1):(ncol(x))]),modtask$n_classes,name_prob)
      dm$add_feature(name_prob,map_func,1,name_map)

      print("ls-ecoc-r with hrel...")
      name_prob = paste(name3,"lsecocr",ecoc$classes_factor,sep=".")
      name_map = paste(name3,"lsecocr","map",sep=".")
      dm$add_feature(c(name4calvec,name4relvec),function(x) ecoc$ls_ecoc_r(x[,1:(ncol(x)/2)],x[,(ncol(x)/2+1):(ncol(x))]),modtask$n_classes,name_prob)
      dm$add_feature(name_prob,map_func,1,name_map)
      
      
    } ### END i3

    for (i3 in 1:n_methods) {
      name3 = paste(name2,".m",i3,sep="")
      print(name3)
      splvec = which(modcalrel_splits$method==i3)
      from3 = paste(name2,".s",splvec,sep="")
      for (n4 in c('lsecoc','lsecoccvar','lsecoccsvar','lsecoccrelvar','lsecocr')) {
#      for (n4 in c('lsecoc')) {
        name4 = paste(name3,n4,sep=".")
        from4 = paste(from3,n4,sep=".")
        print(name4)
        name5vec = vector()
        for (i5 in 1:ecoc$n_classes) {
          name5 = paste(name4,ecoc$classes_factor[i5],sep=".")
          name5vec[i5] = name5
          from5 = paste(from4,ecoc$classes_factor[i5],sep=".")
          if (length(from5)>1) {
            dm$add_feature(from5,function(x) apply(x,1,mean),1,name5)
          } else {
            dm$add_feature(from5,function(x) x,1,name5)
          }
        }
        name4map = paste(name4,"map",sep=".")
        dm$add_feature(name5vec,map_func,1,name4map)
        mat_prob = dm$mat[test_fold$rows,name5vec]
        mat_map = dm$mat[test_fold$rows,name4map]
        err_res[[name4]] = err(mat_map,test_labels)
      }
    } ### END i3

    d = data.frame(method=names(err_res),err=as.numeric(err_res))
    d = d[order(d$err),]
    print(d)
    print('---')


  } ### END i2

  dm$delete_features((n_cols_before_cv+1):dm$m)
  write.csv(d,file=paste(fileprefix,".tmp.results.csv",sep=""),quote=FALSE,row.names=FALSE)
  
} ### END CV

write.csv(d,file=paste(fileprefix,".results.csv",sep=""),quote=FALSE,row.names=FALSE)

file_path =  gsub(paste(paste("~/ecml14experiments/",Sys.time(),sep=""),fileprefix,"rdata",sep="."), pattern=" ", replacement=".")
save.image(file_path);

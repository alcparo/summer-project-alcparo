args = commandArgs(trailingOnly = TRUE)
dataset = as.character(args[1])
bwcal = as.numeric(args[2])
bwrel = as.numeric(args[3])
approxbinsize = as.numeric(args[4])
p_degree = as.numeric(args[5])

# setwd("~/summer-project-alcparo/");
# dataset="diabetes";
# bwcal = 0.2;
# bwrel = 0.2;
# p_degree = 2;
approxbinsize = 200;

# NEW VARIABLES
ensemble_size = 200;

fileprefix = paste(dataset,".bwc",bwcal,".bwr",bwrel,".bs",approxbinsize,".pdeg",p_degree,sep="")
polydegrees = list()

polydegrees[["diabetes"    ]] = p_degree
polydegrees[["australian"    ]] = p_degree
polydegrees[["breast-cancer"    ]] = p_degree
polydegrees[["german.numer"    ]] = p_degree
polydegrees[["mushrooms"    ]] = p_degree
polydegrees[["svmguide1"    ]] = p_degree
polydegrees[["cod-rna"    ]] = p_degree
polydegrees[["splice"    ]] = p_degree
polydegrees[["liver-disorders"    ]] = p_degree


source("NEW_binary_functions.r")

sink(paste(fileprefix,".stdout.txt",sep="."), append=FALSE, split=TRUE)

set.seed(1)

arff = read.arff(file=paste("datasets/arff/",dataset,".arff",sep=""))
colnames(arff) = c(paste('f',1:(ncol(arff)-1),sep=""),'class') 

dm = create_datamatrix(dataset)
dm$add_instances(arff)
task = create_task(dm,colnames(arff)[1:(ncol(arff)-1)],colnames(arff)[ncol(arff)])

if(task$n_classes==2){ ### new
	print('BINARY'); 
	ecoc = create_ecoc(task,type="binary")
} else {
	stop("only work for datasets with binary-class")
}

map_func = function(x) ecoc$classes_factor[apply(x,1,which.max)]

print('calculating ecoc labels...')
dm$add_feature(task$labels,ecoc$calc_labels,ecoc$n_cols,paste(task$labels,colnames(ecoc$mat),sep="."))

print('creating plan')
plan = data.frame()
model_types = data.frame(names=c("svmpol","lr","rf")) #,learners=c(train_svm)) ##### ADD THE LEARNING ALGORITHMS


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
mse_res = list()

print('creating folds')
folds = create_fold(dm)$split(10,proportions=rep(1,10),stratified=TRUE,task=task,fold_names=paste(1:10))

n_cols_before_cv = dm$m

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

	# change the task to the normalized fold
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

  # train_rows = split_folds[['1234']]$rows;

	# LOOP

	for (i2 in 1:nrow(model_types)) { #### 

    # bootstrap_rows = sample(train_rows,replace=T);
    # split_folds[['1234']]$rows = bootstrap_rows;  

    name2 = paste(name1,model_types$names[i2],sep=".")
    print(name2)

    for (i3 in 1:nrow(modcalrel_splits)) {

      name3 = paste(name2,".s",i3,sep="")
      print(name3)

      mod_fold = split_folds[[paste(modcalrel_splits$mod[i3])]]
      cal_fold = split_folds[[paste(modcalrel_splits$cal[i3])]]
      rel_fold = split_folds[[paste(modcalrel_splits$rel[i3])]]

      if (model_types$names[i2]=='svmrad') {
        mod = train_svm(mod_fold,modtask,ecoc,kernel="radial")
      } else if (model_types$names[i2]=='svmpol') {
        mod = train_svm(mod_fold,modtask,ecoc,kernel="polynomial",degree=polydegrees[[dataset]])
      } else if (model_types$names[i2]=='lr') {
        mod = train_lr(mod_fold,modtask,ecoc)
      } else if (model_types$names[i2]=='rf') {
        mod = train_rf(mod_fold,modtask,ecoc)
      } else { stop('unknown model class') }

      name4vec = vector()
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

        name4vec[i4] = name4
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

        # print("const var...")
        # dm$add_feature(c(name4cal,name4rel),function(x) (1-x[,2])*x[,1]*(1-x[,1]),1,name4var)
        # cvar = mean(dm$mat[rel_fold$rows,name4var])
        # csvar = (mean(sqrt(dm$mat[rel_fold$rows,name4var])))^2
        # crel = mean(dm$mat[rel_fold$rows,name4rel])
        # dm$add_feature(c(),function(x) cvar,1,name4cvar)
        # dm$add_feature(c(),function(x) csvar,1,name4csvar)
        # dm$add_feature(name4cal,function(x) (1-crel)*x*(1-x),1,name4crelvar)
      } ## END i4 

      # print("plotting calibration and reliability maps")
      # ggsave(paste(dataset,name3,"cal.pdf",sep="."),plot_calrel(dm,name3,ecoc,calrel=".hcal"),width=5,height=5)
      # ggsave(paste(dataset,name3,"rel.pdf",sep="."),plot_calrel(dm,name3,ecoc,calrel=".hrel"),width=5,height=5)
      
      print("model output mean")
      name_prob = paste(name3,"prob",sep=".")
      name_map = paste(name3,"prob","map",sep=".")
      name_squared_error = paste(name3, "prob", "se", sep=".")
      dm$add_feature(name4vec,function(x) apply(x,1,mean),1, name_prob) 
      dm$add_feature(name_prob, function(x) ((x-as.numeric(as.character(dm$mat$class)))^2 + ((1-x)-ifelse(as.numeric(as.character(dm$mat$class))==0,1,0))^2)/2, 1, name_squared_error);
      dm$add_feature(name_prob,function(x) ifelse(x>0.5,1,0),1, name_map)

      #SE funciton
      #function(x) ((x-as.numeric(as.character(dm$mat$class)))^2 + ((1-x)-ifelse(as.numeric(as.character(dm$mat$class))==0,1,0))^2)/2     

      print("ls-ecoc with hcal...")
      name_prob = paste(name3,"lsecoc",ecoc$classes_factor,sep=".")
      name_map = paste(name3,"lsecoc","map",sep=".")
      name_squared_error = paste(name3, "lsecoc", "se", sep=".")
      dm$add_feature(name4calvec,ecoc$ls_ecoc,modtask$n_classes,name_prob)
      dm$add_feature(name_prob, function(x) ((x[1]-as.numeric(as.character(dm$mat$class)))^2 + (x[2]-ifelse(as.numeric(as.character(dm$mat$class))==0,1,0))^2)/2, 1, name_squared_error);
      dm$add_feature(name_prob,map_func,1,name_map)

      # print("ls-ecoc-v with cvar...")
      # name_prob = paste(name3,"lsecoccvar",ecoc$classes_factor,sep=".")
      # name_map = paste(name3,"lsecoccvar","map",sep=".")
      # dm$add_feature(c(name4calvec,name4cvarvec),function(x) ecoc$ls_ecoc_v(x[,1:(ncol(x)/2)],x[,(ncol(x)/2+1):(ncol(x))]),modtask$n_classes,name_prob)
      # dm$add_feature(name_prob,map_func,1,name_map)

      # print("ls-ecoc-v with csvar...")
      # name_prob = paste(name3,"lsecoccsvar",ecoc$classes_factor,sep=".")
      # name_map = paste(name3,"lsecoccsvar","map",sep=".")
      # dm$add_feature(c(name4calvec,name4csvarvec),function(x) ecoc$ls_ecoc_v(x[,1:(ncol(x)/2)],x[,(ncol(x)/2+1):(ncol(x))]),modtask$n_classes,name_prob)
      # dm$add_feature(name_prob,map_func,1,name_map)

      # print("ls-ecoc-v with crelvar...")
      # name_prob = paste(name3,"lsecoccrelvar",ecoc$classes_factor,sep=".")
      # name_map = paste(name3,"lsecoccrelvar","map",sep=".")
      # dm$add_feature(c(name4calvec,name4crelvarvec),function(x) ecoc$ls_ecoc_v(x[,1:(ncol(x)/2)],x[,(ncol(x)/2+1):(ncol(x))]),modtask$n_classes,name_prob)
      # dm$add_feature(name_prob,map_func,1,name_map)

      print("ls-ecoc-r with hrel...")
      name_prob = paste(name3,"lsecocr",ecoc$classes_factor,sep=".")
      name_map = paste(name3,"lsecocr","map",sep=".")
      name_squared_error = paste(name3, "lsecocr", "se", sep=".")
      dm$add_feature(c(name4calvec,name4relvec),function(x) ecoc$ls_ecoc_r(x[,1:(ncol(x)/2)],x[,(ncol(x)/2+1):(ncol(x))]),modtask$n_classes,name_prob)
      dm$add_feature(name_prob, function(x) ((x[1]-as.numeric(as.character(dm$mat$class)))^2 + (x[2]-ifelse(as.numeric(as.character(dm$mat$class))==0,1,0))^2)/2, 1, name_squared_error);
      dm$add_feature(name_prob,map_func,1,name_map)

    } ### END i3

    for (i3 in 1:n_methods) {
      name3 = paste(name2,".m",i3,sep="")
      print(name3)
      splvec = which(modcalrel_splits$method==i3)
      from3 = paste(name2,".s",splvec,sep="")
      # for (n4 in c('lsecoc','lsecoccvar','lsecoccsvar','lsecoccrelvar','lsecocr')) {
     for (n4 in c('prob','lsecoc', 'lsecocr')) {
      if(n4!='prob'){
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
        name4se = paste(name4,"se",sep=".")
        print(name4se)
        dm$add_feature(name5vec, function(x) ((x[1]-as.numeric(as.character(dm$mat$class)))^2 + (x[2]-ifelse(as.numeric(as.character(dm$mat$class))==0,1,0))^2)/2, 1, name4se);
      
        dm$add_feature(name5vec,map_func,1,name4map)
        mat_prob = dm$mat[test_fold$rows,name5vec]
        } else{
          name3aux = paste(name2,".s",i3,sep="")
          print(name3aux)
          name4 = paste(name3aux,n4,sep=".")
          print(name4)
          name4map = paste(name4,"map",sep=".")
          print(name4map)
          name4se = paste(name4,"se",sep=".")
          print(name4se)
        }
        
        mat_map = dm$mat[test_fold$rows,name4map]
        mat_se = dm$mat[test_fold$rows,name4se]
        err_res[[name4]] = err(mat_map,test_labels)
        mse_res[[name4]] = mse(mat_se)
      }
    }

    d = data.frame(method=names(err_res),err=as.numeric(err_res),mse=as.numeric(mse_res))
    d = d[order(d$err),]
    print(d)
    print('---')


	} ### END i2 model_types

  dm$delete_features((n_cols_before_cv+1):dm$m)
  write.csv(d,file=paste(fileprefix,".tmp.results.csv",sep=""),quote=FALSE,row.names=FALSE)
  


}

mean_list=list()
sd_list=list()
mean_mse_list=list()
sd_mse_list=list()

for(i in 1:nrow(model_types)) {
  for(i2 in c('prob','lsecoc','lsecocr')){
    if(i2=='prob'){
      method = c(paste('cv',1:10,".",model_types$names[i],".s1.",i2,sep=""))
    }else{
      method = c(paste('cv',1:10,".",model_types$names[i],".m1.",i2,sep=""))
    }
      index = unlist(lapply(method, function(x)(which(x==d$method))))
      m = mean(d$err[index])
      m_mse = mean(d$mse[index])
      sd = sd(d$err[index]);
      sd_mse = sd(d$mse[index]);
      mean_list[[paste(model_types$names[i],i2,sep=".")]] = m;
      sd_list[[paste(model_types$names[i],i2,sep=".")]] = sd;

      mean_mse_list[[paste(model_types$names[i],i2,sep=".")]] = m_mse;
      sd_mse_list[[paste(model_types$names[i],i2,sep=".")]] = sd_mse;
  }
}

mean_sd_df = data.frame(method=names(mean_list), err_rate_mean=as.numeric(mean_list), err_rate_sd=as.numeric(sd_list), mse_mean=as.numeric(mean_mse_list), mse_sd=as.numeric(sd_mse_list) );
mean_sd_df = mean_sd_df[order(mean_sd_df$method),];
print(mean_sd_df);
file_path =  gsub(paste(paste("~/summer-project-alcparo/experiments/",Sys.time(),sep=""),fileprefix,"rdata",sep="."), pattern=" ", replacement=".")
save.image(file_path);





mean_list=list()
sd_list=list()
mean_mse_list=list()
sd_mse_list=list()
teste=NULL;


data_plot = data.frame()
for(i in 1:nrow(model_types)) {
  for(i2 in c('prob','lsecoc','lsecocr')){
    if(i2=='prob'){
      method = c(paste('cv',1:10,".",model_types$names[i],".s1.",i2,sep=""))
    }else{
      method = c(paste('cv',1:10,".",model_types$names[i],".m1.",i2,sep=""))
    }
      index = unlist(lapply(method, function(x)(which(x==d$method))))
     # m = mean(d$err[index])
     # m_mse = mean(d$mse[index])
      data_plot[1:10,paste(i2,model_types[[1]][i],sep="_")] = data.frame(d$mse[index]);
      
      
  }

  
}

 data_melt = melt(data_plot);
 data_melt[1:30, "group"] = 1;
 data_melt[31:60, "group"] = 2;
 data_melt[61:90, "group"] = 3;

data_melt[c(1:10,31:40,61:70),'method']='raw';
data_melt[c(11:20,41:50,71:80),'method']='cal';
data_melt[c(21:30,51:60,81:90),'method']='rel';



 #aux_plot = qplot(x=variable, y=value, data=data_melt, aes(group=group));
#aux_plot = qplot(x=variable, y=value, data=data_melt);
aux_plot = ggplot(data=data_melt, aes(x=variable, y=value));
aux_plot + geom_boxplot(aes(fill=factor(group,labels=c("SVM","LR","RF")), interaction(group,method))) + labs(x="Methods", y="Values" ,fill="Classifier");
aux_plot + geom_boxplot(aes(fill=factor(group,labels=c("SVM","LR","RF")))) + labs(x="Methods", y="Values" ,fill="Classifier");
aux_plot + geom_boxplot(aes(fill=factor(group,labels=c("SVM","LR","RF")))) + labs(x="Methods", y="Values" ,fill="Classifier") + scale_x_discrete(labels=c('Raw', 'w/Cal', 'w/Rel', 'Raw', 'w/Cal', 'w/Rel', 'Raw', 'w/Cal', 'w/Rel'));

aux_plot + geom_boxplot(aes(fill=factor(group,labels=c("SVM","LR","RF")))) + labs(x="Methods", y="Values" ,fill="Classifier") + scale_x_discrete(limits=c('prob_svmpol', 'prob_lr', 'prob_rf', 'lsecoc_svmpol', 'lsecoc_lr', 'lsecoc_rf', 'lsecocr_svmpol', 'lsecocr_lr', 'lsecocr_rf'), labels=c('', 'Raw','','','w/Calibration','','','w/Reliability',''));




      limits <- aes(ymax = ymax_vec, ymin=ymin_vec);

      ymax_vec = as.vector(t(lim_mean)+t(lim_sd))
      ymin_vec = as.vector(t(lim_mean)-t(lim_sd))
      maxmin_frame = data.frame(ymax=ymax_vec, ymin=ymin_vec);

      aggregate(data_melt$value, by=list(data_melt$group), FUN=sd)[2]


      lim_sd = aggregate(data_melt$value, by=list(data_melt$group_ind), FUN=sd)[2];
      lim_mean = aggregate(data_melt$value, by=list(data_melt$group_ind), FUN=mean)[2];


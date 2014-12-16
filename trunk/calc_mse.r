calc_mse = function(){

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
		}

		 print("model output mean")
      name_prob = paste(name3,"prob",sep=".")
      name_map = paste(name3,"prob","map",sep=".")
      name_squared_error = paste(name3, "prob", "se", sep=".")
      dm$add_feature(name_prob, function(x) ((x-as.numeric(as.character(dm$mat$class)))^2 + ((1-x)-ifelse(as.numeric(as.character(dm$mat$class))==0,1,0))^2)/2, 1, name_squared_error);
      
		print("ls-ecoc with hcal...")
      name_prob = paste(name3,"lsecoc",ecoc$classes_factor,sep=".")
      name_map = paste(name3,"lsecoc","map",sep=".")
      name_squared_error = paste(name3, "lsecoc", "se", sep=".")
      dm$add_feature(name_prob, function(x) ((x[1]-as.numeric(as.character(dm$mat$class)))^2 + (x[2]-ifelse(as.numeric(as.character(dm$mat$class))==0,1,0))^2)/2, 1, name_squared_error);
      
       print("ls-ecoc-r with hrel...")
      name_prob = paste(name3,"lsecocr",ecoc$classes_factor,sep=".")
      name_map = paste(name3,"lsecocr","map",sep=".")
      name_squared_error = paste(name3, "lsecocr", "se", sep=".")
     dm$add_feature(name_prob, function(x) ((x[1]-as.numeric(as.character(dm$mat$class)))^2 + (x[2]-ifelse(as.numeric(as.character(dm$mat$class))==0,1,0))^2)/2, 1, name_squared_error);
      

		}

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


    }

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


		}
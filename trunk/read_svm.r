source("f_read.libsvm.r");
source("RWeka");

# d = read.libsvm( input_file, dimensionality )
# http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary.html

files = list.files("datasets/dat/");
files_lenght = length(files);
for (i in 1:files_lenght) {

	# file_name = files[i];
	file_name = strsplit(files[i],split=".dat",fixed=TRUE);
	# file_name = paste(name_split,".arff",sep="");
	file_dir_arff = paste("datasets/arff/",file_name,".arff",sep="");
	if (file_name=="diabetes" && !file.exists(file_dir_arff)) {
		dim = 8;
	} else if(file_name=="german.numer" && !file.exists(file_dir_arff)) {
		dim = 24;
	} else if(file_name=="svmguide1" && !file.exists(file_dir_arff)) {
		dim = 4;
	} else if(file_name=="mushrooms" && !file.exists(file_dir_arff)) {
		dim = 112;
	} else if(file_name=="breast-cancer" && !file.exists(file_dir_arff)) {
		dim = 10;
	} else if(file_name=="australian" && !file.exists(file_dir_arff)) {
		dim = 14;
	} else if(file_name=="cod-rna" && !file.exists(file_dir_arff)) {
		dim = 8;
	} else if(file_name=="splice" && !file.exists(file_dir_arff)) {
		dim = 60;
	} else if(file_name=="liver-disorders" && !file.exists(file_dir_arff)) {
		dim = 6;
	} else {
		next();
	}
	

	d = read.libsvm(paste("datasets/dat/",file_name,".dat",sep=""),dim);
	d_rows = nrow(d);
	d_cols = ncol(d);
	d_new = matrix(nrow=d_rows,ncol=d_cols);
	d_new[,d_cols] = d[,1];
	d_new[,-d_cols] = d[,-1];

	colnames(d_new) = c(paste('f',1:(d_cols-1),sep=""),'class');
	rownames(d_new) = c(1:d_rows);	
	d_new = data.frame(d_new);
	d_new[,d_cols] = factor(d_new[,d_cols], labels=c(1,0));
	
	

	print(d_new[1:5,]);
	write.arff(d_new,file_dir_arff);
	print(paste(file_dir_arff," created", sep=""));

}





mean_list=list()
for(i in 1:nrow(model_types)) {
	for(i2 in c('prob','lsecoc','lsecocr')){
		if(i2=='prob'){
			method = c(paste('cv',1:10,".",model_types$names[i],".s1.",i2,sep=""))
		}else{
			method = c(paste('cv',1:10,".",model_types$names[i],".m1.",i2,sep=""))
		}
			index = unlist(lapply(method, function(x)(which(x==d$method))))
			m = mean(d$err[index])
			mean_list[[paste(model_types$names[i],i2,sep=".")]] = m;
	}
}

mean_df = data.frame(method=names(mean_list), mean=as.numeric(mean_list));
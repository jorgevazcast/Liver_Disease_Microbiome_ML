set.seed(12345)
library(caret)
library(doParallel)
library(Boruta)
#library(xgboost)
library(pROC)  # For AUC calculation
library(microbiome)
library(phyloseq)
library(cvAUC)
library(coin)    
library(DMwR)
library(rstatix)

#XGboost_param_matrix <- function(Matrix = "xgbGrid_mid"){

	#### Tune matrix
#	xgbGrid_large <- expand.grid(
#	    nrounds = c(50, 100, 150),  # Test different numbers of boosting rounds
#	    eta = c(0.01, 0.05, 0.1, 0.2),  # More granular values for learning rate
#	    max_depth = c(3, 6, 9, 12),  # Include deeper trees
#	    gamma = c(0, 0.1, 0.2, 0.5),  # Explore different levels of regularization
#	    colsample_bytree = c(0.6, 0.8, 1.0),  # Vary the feature sampling rate
#	    min_child_weight = c(1, 5, 10),  # Different levels of minimum child weight
#	    subsample = c(0.6, 0.8, 1.0)  # Vary the sample size for training
#	)

#	xgbGrid_mid <- expand.grid(
#	    nrounds = c(50, 100,1000),        # Focus on fewer boosting rounds
#	    eta = c(0.05, 0.1, 0.2),          # Skip the smallest value (0.01) to prioritize faster learning
#	    max_depth = c(3, 6, 9),           # Drop the deepest trees (12) to reduce overfitting risk
#	    gamma = c(0, 0.2),                # Test fewer levels of regularization
#	    colsample_bytree = c(0.6, 0.8),   # Focus on the most likely useful ranges
#	    min_child_weight = c(1, 5),       # Reduce the range to test
#	    subsample = c(0.6, 0.8,1.0)          
#	)

#	xgbGrid_small <- expand.grid(
#	    nrounds = c(50, 100),       # Limit boosting rounds to two values
#	    eta = c(0.1, 0.2),          # Focus on faster learning rates
#	    max_depth = c(3, 6),        # Limit to shallower trees
#	    gamma = c(0, 0.1),          # Test minimal regularization
#	    colsample_bytree = c(0.8),  # Use a single, common value
#	    min_child_weight = c(1, 5), # A small range for regularization
#	    subsample = c(0.8)          # A single value to save time
#	)

#	xgbGrid_basic <- expand.grid(
#	    nrounds = c(50),              # Single value for boosting rounds
#	    eta = c(0.1),                 # Focus on one learning rate
#	    max_depth = c(3, 6),          # Test two levels of tree depth
#	    gamma = c(0),                 # Single value for minimal regularization
#	    colsample_bytree = c(0.8),    # Fixed value
#	    min_child_weight = c(1),      # Single value for regularization
#	    subsample = c(0.8)            # Fixed value
#	)


#	list_ret <- list(xgbGrid_large,xgbGrid_mid,xgbGrid_small,xgbGrid_basic)
#	names(list_ret) <- c("xgbGrid_large","xgbGrid_mid","xgbGrid_small","xgbGrid_basic")
	
#	return(list_ret[[Matrix]])

#}

# Function to calculate mode
calculate_mode <- function(x, return_all = FALSE) {
    x <- x[!is.na(x)]
    if(length(x) == 0) return(NA)
    
    ux <- unique(x)
    tab <- tabulate(match(x, ux))
    max_freq <- max(tab)
    
    # Valores con máxima frecuencia
    modes <- ux[tab == max_freq]
    
    if(return_all) {
        return(modes)  # Devuelve vector con todos los valores empatados
    } else {
        return(modes[1])  # Devuelve solo el primero
    }
}



Matt_Coef <- function (conf_matrix){
	TP <- conf_matrix$table[1,1]
	TN <- conf_matrix$table[2,2]
	FP <- conf_matrix$table[1,2]
	FN <- conf_matrix$table[2,1]

	mcc_num <- (TP*TN - FP*FN)
	mcc_den <- as.double((TP+FP))*as.double((TP+FN))*as.double((TN+FP))*as.double((TN+FN))

	mcc_final <- mcc_num/sqrt(mcc_den)
	return(mcc_final)
}

F1_Score <- function(conf_matrix) {

	TP <- conf_matrix$table[1, 1]  # True Positives
	FP <- conf_matrix$table[1, 2]  # False Positives
	FN <- conf_matrix$table[2, 1]  # False Negatives
	TN <- conf_matrix$table[2, 2]  # True Negatives

	# precision and recall
	precision <- TP / (TP + FP)
	recall <- TP / (TP + FN)
  
	# F1 Score
	F1 <- 2 * (precision * recall) / (precision + recall)
  
	return(F1)
}



norm_RMSE <- function(x=c(),range_x= MEspDF$Variable  ){
	ret <- (x/(max(range_x) - min(range_x))) * 100
	return(ret)
}


Stats_CV_function <- function(Stats_CV){
	# Remove the 'Fold' column for calculations
	metrics <- Stats_CV[, -1]  # Remove first column (Fold)

	# Calculate all statistics for each metric
	summary_stats <- data.frame(
	    Metric = colnames(metrics),
	    N = apply(metrics, 2, function(x) sum(!is.na(x))),  # Count non-missing values
	    N_missing = apply(metrics, 2, function(x) sum(is.na(x))),  # Count missing values
	    Mean = apply(metrics, 2, mean, na.rm = TRUE),
	    Median = apply(metrics, 2, median, na.rm = TRUE),
	    #Mode = apply(metrics, 2, calculate_mode),
	    SD = apply(metrics, 2, sd, na.rm = TRUE),
	    Min = apply(metrics, 2, min, na.rm = TRUE),
	    Max = apply(metrics, 2, max, na.rm = TRUE)
	)

	# Round to 4 decimal places for readability (except N and N_missing)
	summary_stats[, c("Mean", "Median", "SD", "Min", "Max")] <- 
	    round(summary_stats[, c("Mean", "Median" , "SD", "Min", "Max")], 4)

	# Add Range column
	summary_stats$Range <- round(summary_stats$Max - summary_stats$Min, 4)

	# Add Coefficient of Variation (CV = SD/Mean)
	summary_stats$CV_percent <- round((summary_stats$SD / summary_stats$Mean) * 100, 3)

	# Print the results
	print(summary_stats)

	return(summary_stats)
}


# Function to create boxplot for a single metric with N
create_metric_boxplot <- function(data, metric_name) {
    # Calculate N (non-missing values)
    n_valid <- sum(!is.na(data[[metric_name]]))
    
    # Create plot
    p <-  ggplot(data, aes(x = 1, y = .data[[metric_name]])) +
	geom_boxplot(fill = "white", outlier.shape = 23) +
	#geom_jitter() +
	geom_point() +
	theme_bw() +
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),plot.title = element_text(face = "bold", size = 11)) +
	labs(title = paste0(metric_name, " (N = ", n_valid, ")"),x = metric_name, y = metric_name) +
	stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "red")
    
    return(p)
}


stats_model_func_continuous <- function(FitModel = NA , Test_df = data.frame(), return_data.frame = T, df.complete = data.frame() ){
	Predicted <- predict(FitModel, newdata = Test_df)
	Stats <- postResample(pred = Predicted, obs = Test_df$Variable)
	Stats <- c(Stats,Cor=cor(Predicted,Test_df$Variable), Cor_Spearman=cor(Predicted,Test_df$Variable,method="spearman"))	
	RMSE_norm <- norm_RMSE( x = Stats["RMSE"] , range_x = df.complete$Variable )
	names(RMSE_norm) <- "RMSE_norm"
	Stats <- c(Stats,RMSE_norm)
	Mean_diff_abs <- mean(abs(Predicted -Test_df$Variable))
	names(Mean_diff_abs) <- "Mean_diff_abs"
	Stats <- c(Stats,Mean_diff_abs)
	if(return_data.frame  == T){ Stats <- data.frame(t(as(Stats,"matrix")))}
	return(Stats)

}

# Calculate classification metrics (AUC, accuracy, kappa, MCC, F1, sensitivity, specificity, precision)
# Two usage modes:
#   1) Pass FitModel + Test_df: predictions are generated internally
#   2) Pass PredProbs + PredClass + RealClass: use pre-computed predictions (set FitModel = NA)
stats_model_func_discrete <- function(FitModel = NA, Test_df = data.frame(), return_data.frame = T, PredProbs=NA, PredClass=NA, RealClass = NA ){
	
	if(all(is.na(FitModel)) & !all(is.na(PredProbs)) & !all(is.na(PredClass)) & !all(is.na(RealClass)) ){
		pred_probs <- PredProbs
		pred_class <- PredClass
								
	}else if( !all(is.na(FitModel)) & nrow(Test_df) != 0){
		pred_probs <- predict_ML(indf = Test_df, MLmodel = FitModel, prop = T )		
		pred_class <- predict_ML(indf = Test_df, MLmodel = FitModel, prop = F )		
		RealClass <- Test_df$Variable
	}
	
	if( class(RealClass) == "character" ){ RealClass <- factor(RealClass) }
	
	if( length(unique(RealClass)) == 2 ){
		roc_curve <- pROC::roc(RealClass, pred_probs[, 2], levels = levels(RealClass))
		auc_value <- as.numeric(pROC::auc(roc_curve))
	}else{
		classes <- levels(RealClass)
		roc_list <- list()
		auc_list <- numeric(length(classes))
		names(auc_list) <- classes
		for (i in seq_along(classes)) {
			this_class <- classes[i]
			binary_response <- ifelse(RealClass == this_class, 1, 0)
			this_prob <- pred_probs[, this_class]
			roc_obj <- roc(binary_response, this_prob)
			roc_list[[this_class]] <- roc_obj
			auc_list[this_class] <- auc(roc_obj)
		}
	}
	
	# Confusion Matrix
	conf_matrix <- confusionMatrix(pred_class, RealClass)
	accuracy <- conf_matrix$overall["Accuracy"]
	kappa <- conf_matrix$overall["Kappa"]
	MCC <- Matt_Coef(conf_matrix)
	F1 <- F1_Score(conf_matrix)
	
	if( length(unique(RealClass)) == 2 ){
		sensitivity <- conf_matrix$byClass["Sensitivity"]
		specificity <- conf_matrix$byClass["Specificity"]
		precision <- conf_matrix$byClass["Precision"]
		PPV <- conf_matrix$byClass["Pos Pred Value"]  # = Precision
		NPV <- conf_matrix$byClass["Neg Pred Value"]
	
		Stats <- c(auc_value, accuracy, kappa, MCC, F1, sensitivity, specificity, precision, PPV, NPV) 	
		names(Stats) <- c("AUC","accuracy","kappa","MCC","F1_score","sensitivity","specificity","precision","PPV","NPV")
		
	}else{
		sensitivity_list <- conf_matrix$byClass[,"Sensitivity"]
		specificity_list <- conf_matrix$byClass[,"Specificity"]
		precision_list <- conf_matrix$byClass[,"Precision"]
		PPV_list <- conf_matrix$byClass[,"Pos Pred Value"]
		NPV_list <- conf_matrix$byClass[,"Neg Pred Value"]
		
		auc_value <- median(auc_list, na.rm=T)
		sensitivity <- median(sensitivity_list, na.rm=T)
		specificity <- median(specificity_list, na.rm=T)
		precision <- median(precision_list, na.rm=T)
		PPV <- median(PPV_list, na.rm=T)
		NPV <- median(NPV_list, na.rm=T)
		
		names(auc_list) <- paste0("AUC.",names(auc_list))
		names(sensitivity_list) <- gsub("Class: ","sens.",names(sensitivity_list))
		names(specificity_list) <- gsub("Class: ","spec.",names(specificity_list))											
		names(precision_list) <- gsub("Class: ","prec.",names(precision_list))
		names(PPV_list) <- gsub("Class: ","PPV.",names(PPV_list))
		names(NPV_list) <- gsub("Class: ","NPV.",names(NPV_list))
																
		Stats <- c(auc_value, accuracy, kappa, MCC, F1, sensitivity, specificity, precision, PPV, NPV) 	
		names(Stats) <- c("AUC","accuracy","kappa","MCC","F1_score","sensitivity","specificity","precision","PPV","NPV")
		
		## Add values per category
		Stats <- c(Stats, c(auc_list, sensitivity_list, specificity_list, precision_list, PPV_list, NPV_list))		
	}	
		
	if(return_data.frame == T){ Stats <- data.frame(t(as(Stats,"matrix")))}
	return(Stats)
}



stats_model_func <- function(Fit.Model = NA, Test.df = data.frame(), RetDF = T, df_complete = data.frame() ){

	if(is.factor(Test.df$Variable)){
		RetStats <- stats_model_func_discrete(FitModel = Fit.Model, Test_df = Test.df , return_data.frame = RetDF )	
	}
	if(is.character(Test.df$Variable)){
		RetStats <- stats_model_func_discrete(FitModel = Fit.Model, Test_df = Test.df , return_data.frame = RetDF )	
	}
	if(is.numeric(Test.df$Variable)){
		RetStats <- stats_model_func_continuous(FitModel = Fit.Model , Test_df = Test.df , return_data.frame = RetDF , df.complete = df_complete )	
	}
	if(is.double(Test.df$Variable)){
		RetStats <- stats_model_func_continuous(FitModel = Fit.Model , Test_df = Test.df , return_data.frame = RetDF, df.complete = df_complete )	
	}
	return(RetStats)
}

##################################################
##### PLOT FUNCTIONS FOR binary classifiers  #####

plot_roc_simple_curve <- function(PredProbs, PredClass, RealClass, mcc_value) {
	
	if (class(RealClass) == "character") RealClass <- factor(RealClass)
	if (class(PredClass) == "character") PredClass <- factor(PredClass)
	n_classes <- length(levels(RealClass))
	
	if (n_classes == 2) {
		# Binary ROC
		roc_curve <- pROC::roc(RealClass, PredProbs[, 2], levels = levels(RealClass), quiet = TRUE)
		auc_value <- as.numeric(pROC::auc(roc_curve))
		
		p_roc <- pROC::ggroc(roc_curve, size = 1.2, color = "steelblue") +
			geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray50") +
			annotate("text", x = 0.2, y = 0.15, label = paste0("AUC = ", round(auc_value, 3)), size = 5) +
			annotate("text", x = 0.2, y = 0.10, label = paste0("MCC = ", round(mcc_value, 3)), size = 5) +
			theme_bw(base_size = 12) +
			labs(title = "ROC Curve", x = "Specificity", y = "Sensitivity") +
			theme(plot.title = element_text(hjust = 0.5, face = "bold"))
		
	} else {
		# Multiclass ROC (One-vs-Rest)
		classes <- levels(RealClass)
		roc_list <- list()
		auc_list <- c()
		mcc_list <- c()
		
		for (this_class in classes) {
			# ROC and AUC per class
			binary_response <- ifelse(RealClass == this_class, 1, 0)
			roc_list[[this_class]] <- pROC::roc(binary_response, PredProbs[, this_class], quiet = TRUE)
			auc_list[this_class] <- as.numeric(pROC::auc(roc_list[[this_class]]))
			
			# MCC per class (One-vs-Rest) with NA handling
			binary_pred <- ifelse(PredClass == this_class, 1, 0)
			TP <- sum(binary_response == 1 & binary_pred == 1)
			TN <- sum(binary_response == 0 & binary_pred == 0)
			FP <- sum(binary_response == 0 & binary_pred == 1)
			FN <- sum(binary_response == 1 & binary_pred == 0)
			
			# Check for edge cases
			denom_parts <- c((TP + FP), (TP + FN), (TN + FP), (TN + FN))
			if (any(denom_parts == 0)) {
				mcc_list[this_class] <- 0
			} else {
				denom <- sqrt(prod(denom_parts))
				mcc_list[this_class] <- ((TP * TN) - (FP * FN)) / denom
			}
		}
		
		# Median values (handle NAs)
		median_auc <- median(auc_list, na.rm = TRUE)
		median_mcc <- median(mcc_list, na.rm = TRUE)
		
		# Create legend labels with AUC and MCC per class
		legend_labels <- paste0(
			gsub("Enterotype\\.", "", classes),  # Shorter names
			" (AUC=", round(auc_list, 2), 
			", MCC=", round(mcc_list, 2), ")"
		)
		names(legend_labels) <- classes
		
		p_roc <- pROC::ggroc(roc_list, size = 1) +
			geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray50") +
			annotate("text", x = 0.3, y = 0.15, 
				label = paste0("Median AUC = ", round(median_auc, 3), "\nMedian MCC = ", round(median_mcc, 3)), 
				size = 4, hjust = 0) +
			theme_bw(base_size = 12) +
			labs(title = "ROC Curve (One-vs-Rest)", x = "Specificity", y = "Sensitivity", color = "Class") +
			theme(
				plot.title = element_text(hjust = 0.5, face = "bold"),
				legend.position = "right"
			) +
			scale_color_brewer(palette = "Set1", labels = legend_labels)
	}
	
	return(p_roc)
}

plot_auc_function <- function(list_ggroc_curves=list(), list_ggroc_auc = c() , Variable = "",  Categories_variable = c() ){

	list_sd_roc_figures <- list()
	sum_plots <- 1
	#### AUC mean and standar deviation #####
	if(length(Categories_variable) == 2){		
		mean_ROC_list <- mean_ROC_function(list_ggroc_curves, Variable)
		mean_ROC_df <-  mean_ROC_list$mean_ROC 

		auc_text<-paste0("AUC","\n", Variable,"=",mean_ROC_list$AUC, "/", mean_ROC_list$AUCsd)

		AUCMeanAUCsd <- plot_meanROC(mean_ROC_plot=mean_ROC_df, titleplot=paste("AUC",Variable), auc_text_vec = auc_text)
		AUCMeanAUCsd <- AUCMeanAUCsd + labs(title = paste("AUC",Variable), 
				subtitle = paste("Mean =", round(mean(list_ggroc_auc),digits=2) , "SD =", round(sd(list_ggroc_auc),digits=2) ) 
			  )
		ggsave(file = paste0(Variable,"_mean_sd_AUC.pdf"), AUCMeanAUCsd, width = 10, height = 7)
		list_sd_roc_figures[[sum_plots]] <- AUCMeanAUCsd
				
		#### AUC lines #####
		names(list_ggroc_curves) <- paste("CV",1:length(list_ggroc_curves))
		AUCplotLines <- ggroc(list_ggroc_curves)
		AUCplotLines <- AUCplotLines + geom_abline(slope = 1, intercept = 1, lty = 2, colour = 'black')+
			  labs(title = paste("AUC",Variable), 
				subtitle = paste("Mean =", round(mean(list_ggroc_auc),digits=2) , "SD =", round(sd(list_ggroc_auc),digits=2) ) ) + 
			  theme_bw()  # + scale_x_continuous(limits = c(0, 1))
		ggsave(file = paste0(Variable,"_AUC_CV_lines.pdf"), AUCplotLines, width = 13, height = 9)


	}else{
		### Obtain the AUC curve per categorie
		for(i in Categories_variable ){

			temp_list_ggroc_curves <- list()
			sum <- 1
			for(subggroc in list_ggroc_curves){
				temp_list_ggroc_curves[[sum]] <- subggroc[[i]]
				sum <- sum + 1	
			}

			temp_list_ggroc_auc <- c()
			for(sub_ggroc_auc in list_ggroc_auc){
				temp_list_ggroc_auc <- c(temp_list_ggroc_auc,sub_ggroc_auc[[i]])
			}
					
			
			mean_ROC_list <- mean_ROC_function(temp_list_ggroc_curves, i)
			mean_ROC_df <-  mean_ROC_list$mean_ROC 
			auc_text<-paste0("AUC","\n", i,"=",mean_ROC_list$AUC, "/", mean_ROC_list$AUCsd)	

			AUCMeanAUCsd <- plot_meanROC(mean_ROC_plot=mean_ROC_df, titleplot=paste("AUC",i), auc_text_vec = auc_text)		
			AUCMeanAUCsd <- AUCMeanAUCsd + labs(title = paste("AUC",i), 
					subtitle = paste("Mean =", round(mean(temp_list_ggroc_auc),digits=2) , 
					"SD =", round(sd(temp_list_ggroc_auc),digits=2) ) )
			
			ggsave(file = paste0(i,"_mean_sd_AUC.pdf"), AUCMeanAUCsd, width = 10, height = 7)
			
			sum_plots <- sum_plots + 1
			list_sd_roc_figures[[sum_plots]] <- AUCMeanAUCsd

			#### AUC lines #####
			names(temp_list_ggroc_curves) <- paste("CV",1:length(temp_list_ggroc_curves))
			AUCplotLines <- ggroc(temp_list_ggroc_curves)
			AUCplotLines <- AUCplotLines + geom_abline(slope = 1, intercept = 1, lty = 2, colour = 'black')+
				  labs(title = paste("AUC",i), 
					subtitle = paste("Mean =", round(mean(temp_list_ggroc_auc),digits=2) , 
					"SD =", round(sd(temp_list_ggroc_auc),digits=2) ) ) + 
					theme_bw()  # + scale_x_continuous(limits = c(0, 1))
			ggsave(file = paste0(i,"_AUC_CV_lines.pdf"), AUCplotLines, width = 10, height = 7)
			
			rm(sum,mean_ROC_list,temp_list_ggroc_curves,mean_ROC_df,auc_text,AUCMeanAUCsd,temp_list_ggroc_auc,AUCplotLines)		

		}


	}
	saveRDS(file="ROC_sd_curves.rds",list_sd_roc_figures)
}

mean_ROC_function <- function(list_ggroc_curves, Variable){
	list_sens <- lapply(X = list_ggroc_curves,FUN = function(x) { x$sensitivities} )
		          
	list_spec <- lapply(X = list_ggroc_curves, FUN = function(x) { specificities = x$specificities} )

	list_auc <- data.frame(sapply(X = list_ggroc_curves, FUN = function(x) { x$auc}))

	colnames(list_auc)= Variable
	auc1=colMeans(list_auc)
	auc1<-round(auc1,2) 

	sdAUC <- round( sd(list_auc[,1]) ,2)
	seAUC <- round( sdAUC / sqrt(length(list_auc[,1])) ,2)
	names(sdAUC) <- Variable
	names(seAUC) <- Variable

	biggerCV_partition<- max(unlist(lapply(list_sens,length)))
	sens_df <-   matrix(NA, length(list_sens) , biggerCV_partition)  
	spec_df <-   matrix(NA, length(list_sens), biggerCV_partition)  
	colnames(sens_df) <- c(1:biggerCV_partition)
	colnames(spec_df) <- c(1:biggerCV_partition)

	for(inds in 1:length(list_sens) ){
	  
	  # inds <- 19
	  #print(length(list_sens[[inds]]))
	  #print(length(list_spec[[inds]]))
	  
	  sens_vec  <- list_sens[[inds]]
	  names(sens_vec) <- 1:length(sens_vec)
	  
	  spec_vec  <- list_spec[[inds]]
	  names(spec_vec) <- 1:length(spec_vec)
	  
	  sens_df[inds,names(sens_vec)] <- sens_vec
	  spec_df[inds,names(spec_vec)] <- spec_vec
	  
	  #sens_df <- rbind(sens_df, list_sens[[inds]])
	  #spec_df <- rbind(spec_df, list_spec[[inds]])
	  
	}

	mean.sens_df <- apply(sens_df,2,function(x){mean(x,na.rm=T)})
	mean.spec_df  <- apply(spec_df,2,function(x){mean(x,na.rm=T)})

	se.spec_df  <- apply(spec_df,2,function(x){sd(x,na.rm=T)/sqrt(length(x))})
	se.sens_df  <- apply(sens_df,2,function(x){sd(x,na.rm=T)/sqrt(length(x))})

	sd.spec_df  <- apply(spec_df,2,function(x){ sd(x,na.rm=T) })
	sd.sens_df  <- apply(sens_df,2,function(x){ sd(x,na.rm=T) })
	
	  
	mean_ROC <- data.frame(sens = mean.sens_df, spec = mean.spec_df, se.spec_df,se.sens_df,sd.spec_df,sd.sens_df, Var= Variable)

	list_ret <- list(auc1, sdAUC, seAUC, mean_ROC)
	names(list_ret) <- c("AUC","AUCsd","AUCse","mean_ROC")
	
	return(list_ret)
}	


plot_meanROC <- function(mean_ROC_plot=data.frame(), titleplot="", auc_text_vec = c(), Deviation = "SE"){ # Deviation = "SE" or Deviation = "SD"
#	plotret <- ggplot(mean_ROC_plot, aes(x=1-spec, y=sens,color="red", group=Var)) +  
	plotret <- ggplot(mean_ROC_plot, aes(x=1-spec, y=sens)) +    
	  geom_path(size = 1, linetype = "solid", color = 'blue')+
	  geom_abline(slope = 1, intercept = 0, lty = 2, colour = 'red')+   
	  theme_minimal() + 
	  labs(x="FPR", y="TPR", title = titleplot)+
	  annotate("text", x=0.87, y= 0.18,label = auc_text_vec, size=4.5)+
	  theme(text = element_text(size = 13))+
	  theme(axis.text.x = element_text(colour="black", size=12),axis.ticks.x=element_blank())+
	  theme(axis.text.y = element_text(colour="black", size=12),axis.ticks.y=element_blank())+
	  theme(panel.border = element_rect(colour="black", fill=NA))+
	 # scale_color_manual(values=c( "#e6b729", "#e26731", "#a29a3b", "#6f9fbf"), guide="none")+
	 # scale_color_manual(values=c( "#e6b729", "#e26731", "#a29a3b", "#6f9fbf") )+
	  scale_x_continuous(limits = c(0, 1))
	  
	if(Deviation == "SE"){
		plotret <- plotret + geom_errorbar(aes(ymin=sens-se.sens_df, ymax=sens+se.sens_df), width=.008) 
	}else if(Deviation == "SD"){
		plotret <- plotret + geom_errorbar(aes(ymin=sens-sd.sens_df, ymax=sens+sd.sens_df), width=.008) 
	}
	return(plotret)
}


ROC_AUC_func <- function(Test_df,pred_probs){
	classes <- levels(Test_df$Variable)
	roc_list <- list()
	auc_list <- numeric(length(classes))
	names(auc_list) <- classes

	for (i in seq_along(classes)) {
		this_class <- classes[i]
	  
		# Create binary labels: 1 for this class, 0 for all others
		binary_response <- ifelse(Test_df$Variable == this_class, 1, 0)
		  
		# Get predicted probabilities for this class
		this_prob <- pred_probs[, this_class]  # assuming pred_probs has column names matching classes
		  
		# Compute ROC
		roc_obj <- pROC::roc(binary_response, this_prob)
		roc_list[[this_class]] <- roc_obj
		auc_list[this_class] <- auc(roc_obj)
	}

	list_ret <- list(roc_list,auc_list)
	names(list_ret) <- c("roc_list","auc_list")
	return(list_ret)
}
	


load_RData <- function(file) {
	# Load the data into a temporary environment
	temp_env <- new.env()
	load(file, envir = temp_env)
  
	# Get the first object (assuming there's only one)
	obj <- get(ls(temp_env)[1], envir = temp_env)
  
	# Return the renamed object
	return(obj)
}

data_format <- function(in.obj, Variable,PrevCutoff = 0.2){

	if(class(in.obj) == "phyloseq"){
		physeq_obj <- in.obj
		Metadata <- data.frame(sample_data(physeq_obj))
		Metadata$Variable <- Metadata[,colnames(Metadata) %in% Variable]
		Metadata <- Metadata[complete.cases(Metadata$Variable),]
	
		ret_phylo <- prune_samples( rownames(Metadata), physeq_obj  )

		Prevalence <- microbiome::prevalence(ret_phylo)
		in_phylo_filter <- prune_taxa( names(Prevalence[Prevalence >= PrevCutoff]) , ret_phylo)

		OTUs <- otu_table(in_phylo_filter)
		Abundance_data <-  t(OTUs[,match(rownames(Metadata), colnames(OTUs))])
			

		if(all(rownames(Metadata) == rownames(Abundance_data))){
			MEspDF<- data.frame(Variable = Metadata$Variable,Abundance_data)  
			rownames(MEspDF) <- rownames(Metadata)	
		}
		else{
			stop("Error: Not matching metadata and taxa table rownames")
		}

			
	}
	
	return(MEspDF)

}


cont_res_wt <- function(subDF = data.frame(), discrete = "", continuous = "", Paired = F, ID = NULL){
	
	# Check if ID is required for paired test
	if(Paired && is.null(ID)){
		stop("Paired = TRUE requires ID column to ensure correct sample pairing")
	}
	
	# Select columns based on whether ID is provided
	if(!is.null(ID)){
		subDF <- subDF[, match(c(discrete, continuous, ID), colnames(subDF))]
		colnames(subDF) <- c("discrete", "continuous", "ID")
	} else {
		subDF <- subDF[, match(c(discrete, continuous), colnames(subDF))]
		colnames(subDF) <- c("discrete", "continuous")
	}
	
	subDF <- subDF[complete.cases(subDF),]
	subDF$discrete <- factor(subDF$discrete)
	
	level1 = levels(subDF$discrete)[1]
	level2 = levels(subDF$discrete)[2]
	
	# Split by group
	df1 <- subDF[subDF$discrete == level1, ]
	df2 <- subDF[subDF$discrete == level2, ]
	
	# If paired, sort by ID to ensure correct pairing
	if(Paired){
		df1 <- df1[order(df1$ID), ]
		df2 <- df2[order(df2$ID), ]
		
		# Verify same IDs in both groups
		if(!identical(df1$ID, df2$ID)){
			stop("Paired test requires same IDs in both groups")
		}
	}
	
	group1 <- df1$continuous
	group2 <- df2$continuous
	
	WT <- wilcox.test(group1, group2, paired = Paired)
	
	mean.level1 <- mean(group1)
	median.level1 <- median(group1)
	sd.level1 <- sd(group1)
	N_level1 <- length(group1)
	
	mean.level2 <- mean(group2)
	median.level2 <- median(group2)
	sd.level2 <- sd(group2)
	N_level2 <- length(group2)
	
	#####   Estimate the effect size   #####
	Z <- statistic(independence_test(continuous ~ discrete, data = subDF), type = "standardized")	
	N <- nrow(subDF)
	EffectSize_r <- c(abs(Z) / sqrt(N))
	names(EffectSize_r) <- "EffectSize_r"
	
	if(median.level1 != median.level2){
		Dominant <- ifelse(median.level1 - median.level2 > 0, as.character(level1), as.character(level2))	
	} else if(mean.level1 != mean.level2){
		Dominant <- ifelse(mean.level1 - mean.level2 > 0, as.character(level1), as.character(level2))		
	} else {
		Dominant <- "Undetermined"
	}
	
	retDF <- data.frame(Discrete = discrete, Continuous = continuous, Level1 = level1, N_level1, Level2 = level2, N_level2,
		mean.level1, mean.level2, median.level1, median.level2, sd.level1, sd.level2,
		statistic = WT$statistic, EffectSize_r, p.value = WT$p.value, N = nrow(subDF), Dominant)
	return(retDF)
}

cont_res_kw <- function(subDF = data.frame(),discrete = "",continuous = ""){
	subDF <- subDF[,match(c(discrete,continuous), colnames(subDF))]
	colnames(subDF) <-c("discrete","continuous")
	subDF <- subDF[complete.cases(subDF),]
	subDF$discrete <- factor(subDF$discrete)
	
	KT <- kruskal.test(continuous ~ discrete, data = subDF)
	Mean <- c(by(subDF$continuous , subDF$discrete,mean))
	SD <- c(by(subDF$continuous , subDF$discrete,sd))
	N <- c(by(subDF$continuous , subDF$discrete,length))

	names(Mean) <- paste0("Mean.",names(Mean))
	names(SD) <- paste0("SD.",names(SD))
	names(N) <- paste0("N.",names(N))
	KT.statistic = KT$statistic
	KT.p.value = KT$p.value; names(KT.p.value) <- "p.value"
	kw_effsize <- kruskal_effsize(continuous ~ discrete, data = subDF )
	Effsize_eta2 <- kw_effsize$effsize
	names(Effsize_eta2) <- "Effsize_eta2"
	TotalN <- nrow(subDF)
	names(TotalN) <- "Total_N"
	ret <-c(Mean,SD,N, KT.statistic,Effsize_eta2,KT.p.value,TotalN)		
	

	ret.df <- data.frame(t(ret))
	retDF <- cbind(data.frame(Discrete=discrete,Continuous=continuous) , ret.df)
	
	return(retDF)
}

predict_ML <- function(phylo.obj = NULL, indf = data.frame(), MLmodel, prop = F){

	if( !is.null(phylo.obj) & nrow(indf) == 0 ){
		in.df <- data.frame(t(otu_table(phylo.obj)))
		if (!taxa_are_rows(phylo.obj)) { in.df <- data.frame(t(in.df)) }	
	} else if( is.null(phylo.obj) & nrow(indf) != 0 ){
		in.df <- indf
	} else {
		stop("Provide either phylo.obj OR indf, not both or neither")
	}

	#####################################
	### remove the cohort unique taxa ###
	cohort_specific_taxa <- setdiff(colnames(in.df), colnames(MLmodel$trainingData)) 
	if(length(cohort_specific_taxa) > 0){
	    in.df <- in.df[, !colnames(in.df) %in% cohort_specific_taxa, drop = FALSE]
	}

	############################
	### Add the missing taxa ###
	missing_taxa <- setdiff(colnames(MLmodel$trainingData),colnames(in.df))
	missing_taxa <- missing_taxa[!missing_taxa %in% c(".outcome")]
	AddMissingTaxa<- data.frame(matrix(0,nrow(in.df),length(missing_taxa)))
	rownames(AddMissingTaxa) <- rownames(in.df) 
	colnames(AddMissingTaxa) <- missing_taxa

	df2predict <- cbind(in.df,AddMissingTaxa)

	### Match the taxa ####
	Taxa2Match <- colnames(MLmodel$trainingData) 
	Taxa2Match <- Taxa2Match[!Taxa2Match %in% c(".outcome")]
	df2predict <- df2predict[, match(Taxa2Match, colnames(df2predict)), drop = FALSE]
	
	#Richness <- richness(phylo.obj)
	#Richness <- Richness[match(rownames(df2predict),rownames(Richness)),]
	#df2predict$Richness <- Richness$observed
	# predict(MLmodel,df2predict,type = "prob")
	if(prop == T){
		ret_pred <- predict(MLmodel,df2predict,type = "prob")
		rownames(ret_pred) <- rownames(df2predict)
	}else{
		ret_pred <- predict(MLmodel,df2predict)	
		names(ret_pred) <- rownames(df2predict)
	}

	return(ret_pred)
}

# Function to decide SMOTE configuration based on class imbalance
smart_smote <- function(data, threshold = 1/5,
                        perc.over.default = 200, perc.under.default = 200,
                        perc.over.adjust = 100, perc.under.adjust = 150) {

	# Extract response variable from formula
	tbl <- table(data$Variable)

	# Identify minority and majority class sizes
	min_class <- min(tbl)
	max_class <- max(tbl)
	ratio <- min_class / max_class

	message("Class ratio (minority / majority) = ", round(ratio, 3))

	# Decide which SMOTE parameters to use based on imbalance
	if (ratio < threshold) {
		message("⚠️ Strong imbalance detected → using adjusted SMOTE params")
		new_data <- SMOTE(Variable ~ ., data, perc.over = perc.over.adjust, perc.under = perc.under.adjust)
	} else {
		message("✅Mild imbalance → using default SMOTE params")
		new_data <- SMOTE(Variable ~ ., data, perc.over = perc.over.default, perc.under = perc.under.default)
	}

	return(new_data)
}

############################################
#####        Importance plot           #####
############################################

#Imp_plot_discrete_function <- function(importance_df = data.frame(), MeanImp=data.frame(), MEspDF = data.frame() ){
Imp_plot_discrete_function <- function(importance_df = data.frame(), MEspDF = data.frame()) {
	
	#### Imp Plot #####
	Mean <- c(by(importance_df$Overall, importance_df$Feature, mean))
	SD <- c(by(importance_df$Overall, importance_df$Feature, sd))
	Feature <- names(Mean)
	MeanImp <- data.frame(Feature, Mean, SD)
	
	###  Direction of the taxa abundance  ####
	Taxa2CompPhylo_df <- reshape2::melt(MEspDF, id.vars = "Variable")
	N_levels <- length(unique(MEspDF$Variable))
	
	rank_diff_df <- data.frame()
	for (i in MeanImp$Feature) {
		tempDF <- subset(Taxa2CompPhylo_df, variable == i)
		
		# Calculate dominant class (highest median abundance)
		medians <- aggregate(value ~ Variable, data = tempDF, FUN = median)
		dominant_class <- as.character(medians$Variable[which.max(medians$value)])
		
		# Statistical test
		if (N_levels == 2) {
			wt <- wilcox.test(value ~ Variable, data = tempDF)
			p_value <- wt$p.value
		} else {
			kt <- kruskal.test(value ~ Variable, data = tempDF)
			p_value <- kt$p.value
		}
		
		rank_diff_df <- rbind(rank_diff_df, data.frame(
			Feature = i,
			Dominant = dominant_class,
			p.value = p_value
		))
	}
	
	MeanImp$Increase <- rank_diff_df[match(MeanImp$Feature, rank_diff_df$Feature), "Dominant"]
	rank_diff_df$q.value <- p.adjust(rank_diff_df$p.value, method = "BH")
	SigFeatures <- ifelse(rank_diff_df$q.value < 0.1, "*", "")
	names(SigFeatures) <- rank_diff_df$Feature
	MeanImp$WTSignificant <- SigFeatures[match(MeanImp$Feature, names(SigFeatures))]
	
	#############################
	#### Write the results  #####
	write.table(rank_diff_df, "RankDiff_SelectedFeture_importance.tsv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	write.table(MeanImp, "MeanImportance.tsv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	
	#################
	#### Plots  #####
	MeanImp$Feature <- paste(MeanImp$WTSignificant, MeanImp$Feature)
	MeanImp <- MeanImp[order(MeanImp$Mean), ]
	MeanImp$Feature <- factor(as.character(MeanImp$Feature), as.character(MeanImp$Feature))
	
	p <- ggplot(MeanImp, aes(x = Feature, y = Mean, fill = Increase)) + 
		geom_bar(stat = "identity") +
		geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
		coord_flip() + 
		theme_bw() + 
		ylab("Mean Relative Importance") +
		ggtitle("Feature Importance")
	
	if (N_levels == 2) {
		p <- p + scale_fill_manual(values = c("blue", "red"))
	} else {
		p <- p + scale_fill_brewer(palette = "Set1")
	}
	
	# Adjust height
	n_features <- length(unique(MeanImp$Feature))
	plot_height <- ifelse(n_features > 160, 20, ifelse(n_features > 100, 17, 15))
	ggsave(file = "Imp.pdf", p, width = 10, height = plot_height)
	
	return(p)
}


############################################
##### Best Feature selection functions #####
############################################

feature_selection <- function( List_models = list(), cutoff = 0.5 ){
	all_features <- c()
	for(i in seq_along(List_models)) {  
		Model <- List_models[[i]]
		Features <- colnames(Model$trainingData)
		Features <- Features[!Features %in% ".outcome"]
		all_features <- c(all_features,Features)
	}
    
	Prev_CV <- table(all_features)  / length(List_models)
	return(Prev_CV[Prev_CV > cutoff])
}

feature_selection_adaptive <- function(List_models = list(), cutoffs = c(0.6, 0.5, 0.4, 0.3)) {
  # Adaptive feature selection with decreasing cutoff thresholds
  #
  # Attempts to select features starting from the highest cutoff.
  # If no features pass the threshold, it automatically decreases
 # to the next cutoff level.
  #
  # Args:
  #   List_models: list of trained caret models from CV folds
  #   cutoffs: numeric vector of thresholds to try (descending order)
  #
  # Returns:
  #   Named numeric vector with selected features and their prevalence
  
  # === TRY EACH CUTOFF IN ORDER ===
  for(cutoff in cutoffs) {
    Best_features <- feature_selection(List_models, cutoff = cutoff)
    
    if(length(Best_features) > 0) {
      # Success: features found at this cutoff
      cat("Features selected with cutoff =", cutoff, ":", length(Best_features), "\n")
      return(Best_features)
    } else {
      # No features passed: warn and try lower cutoff
      warning(paste0("Cutoff ", cutoff, " selected no features. Lowering threshold..."))
    }
  }
  
  # === FALLBACK: TOP 10 FEATURES ===
  # If no cutoff worked, return most frequent features
  warning("No cutoff worked. Returning top 10 features by prevalence.")
  
  # Collect all features across models
  all_features <- c()
  for(i in seq_along(List_models)) {  
    Model <- List_models[[i]]
    Features <- colnames(Model$trainingData)
    Features <- Features[!Features %in% ".outcome"]
    all_features <- c(all_features, Features)
  }
  
  # Calculate prevalence and return top 10
  Prev_CV <- table(all_features) / length(List_models)
  Prev_CV <- sort(Prev_CV, decreasing = TRUE)
  
  return(head(Prev_CV, 10))
}

####################################################
##### Best hyperparameters selection functions #####
####################################################

# Function to estimate the mean for the MCC and AUC for the best model selection
mean_auc_acc_param_function <- function(parameter_list=list(), performance_parameters=data.frame(), parameter2check=""  ){
	auc_mcc_stats <- data.frame()
	for(i in parameter_list){
		sub_df <- performance_parameters[performance_parameters[[parameter2check]] == i,]
		stats_vec <- apply( sub_df[,c("AUC","MCC")] , 2 , mean )
		stats_vec <- data.frame(t(data.frame(stats_vec)))
		rownames(stats_vec) <- NULL
		stats_vec$Parameter <- i
		auc_mcc_stats <- rbind( auc_mcc_stats, stats_vec )
	}
	return(auc_mcc_stats)
}

select_best_model_hyperparameter <- function(performance_summary, 
                              tie_threshold = 0.01,
                              verbose = TRUE) {
    # Select best model based on hierarchical criteria:
    # 1. Primary: MCC (higher is better, non-NA preferred over NA)
    # 2. Secondary: AUC (if all models have MCC=NA)
    # 3. Tie-break: Use AUC if MCC is tied
    #
    # Args:
    #   performance_summary: data.frame with "Parameter", "MCC", "AUC" columns
    #   tie_threshold: numeric, difference to consider a tie (default: 0.01)
    #   verbose: logical, print selection process (default: TRUE)
    #
    # Returns:
    #   list with optimal_parameter, optimal_model, selection_method
    
    # Validate input
    required_cols <- c("Parameter", "MCC", "AUC")
    if(!all(required_cols %in% names(performance_summary))) {
        stop("DataFrame must contain: Parameter, MCC, AUC")
    }
    
    if(nrow(performance_summary) == 0) {
        stop("DataFrame is empty")
    }
    
    # Remove rows where both MCC and AUC are NA
    valid_models <- performance_summary[
        !(is.na(performance_summary$MCC) & is.na(performance_summary$AUC)), 
    ]
    
    if(nrow(valid_models) == 0) {
        stop("No valid models: all have MCC=NA and AUC=NA")
    }
    
    if(verbose) {
        cat("\n=== MODEL SELECTION ===\n")
        cat("Total models:", nrow(performance_summary), "\n")
        cat("Valid models:", nrow(valid_models), "\n")
    }
    
    # Separate models with/without MCC
    models_with_mcc <- valid_models[!is.na(valid_models$MCC), ]
    models_without_mcc <- valid_models[is.na(valid_models$MCC), ]
    
    if(verbose) {
        cat("Models with MCC:", nrow(models_with_mcc), "\n")
        cat("Models with MCC=NA:", nrow(models_without_mcc), "\n\n")
    }
    
    # STEP 1: Primary criterion - MCC
    if(nrow(models_with_mcc) > 0) {
        # Use models with MCC (ignore those with MCC=NA)
        candidates <- models_with_mcc
        primary_metric <- "MCC"
        
        if(verbose) {
            cat("Primary criterion: MCC\n")
        }
        
        # Order by MCC (descending)
        candidates <- candidates[order(candidates$MCC, decreasing = TRUE), ]
        best_value <- candidates$MCC[1]
        
        # Check for ties in MCC
        tied <- candidates[abs(candidates$MCC - best_value) <= tie_threshold, ]
        
        if(nrow(tied) > 1) {
            # Tie in MCC → use AUC as tie-breaker
            if(verbose) {
                cat("Tie detected in MCC (n=", nrow(tied), ")\n", sep="")
                cat("Tie-break using AUC\n")
            }
            tied <- tied[order(tied$AUC, decreasing = TRUE, na.last = TRUE), ]
            optimal_model <- tied[1, ]
            selection_method <- "MCC (tie-break: AUC)"
        } else {
            # Clear winner by MCC
            optimal_model <- candidates[1, ]
            selection_method <- "MCC"
        }
        
    } else {
        # STEP 2: All models have MCC=NA → use AUC
        candidates <- models_without_mcc
        primary_metric <- "AUC"
        
        if(verbose) {
            cat("All models have MCC=NA\n")
            cat("Secondary criterion: AUC\n")
        }
        
        # Order by AUC (descending)
        candidates <- candidates[order(candidates$AUC, decreasing = TRUE), ]
        best_value <- candidates$AUC[1]
        
        # Check for ties in AUC
        tied <- candidates[abs(candidates$AUC - best_value) <= tie_threshold, ]
        
        if(nrow(tied) > 1) {
            if(verbose) {
                cat("Tie detected in AUC (n=", nrow(tied), ")\n", sep="")
                cat("Selecting first model\n")
            }
            optimal_model <- tied[1, ]
            selection_method <- "AUC (tie)"
        } else {
            optimal_model <- candidates[1, ]
            selection_method <- "AUC"
        }
    }
    
    # Print result
    if(verbose) {
        cat("\n=== SELECTED MODEL ===\n")
        cat("Parameter:", optimal_model$Parameter, "\n")
        cat("MCC:", ifelse(is.na(optimal_model$MCC), "NA", 
                          round(optimal_model$MCC, 4)), "\n")
        cat("AUC:", ifelse(is.na(optimal_model$AUC), "NA", 
                          round(optimal_model$AUC, 4)), "\n")
        cat("Method:", selection_method, "\n\n")
    }
    
    # Return result
    return(list(
        optimal_parameter = optimal_model$Parameter,
        optimal_model = optimal_model,
        selection_method = selection_method
    ))
}

best_hyperparameter_selection <- function(List_models = list(), 
                                         Performance = data.frame(),
                                         tie_threshold = 0.01,
                                         verbose = TRUE) {
    # Select optimal hyperparameters from nested CV models
    # Uses mode for each parameter, with MCC/AUC tie-breaking
    #
    # Args:
    #   List_models: list of trained caret models (one per CV fold)
    #   Performance: data.frame with Fold, AUC, MCC columns
    #   tie_threshold: threshold for tie-breaking (default: 0.01)
    #   verbose: print selection process (default: TRUE)
    #
    # Returns:
    #   data.frame with one row containing optimal parameter values
    
    # === INPUT VALIDATION ===
    if(length(List_models) == 0) {
        stop("List_models is empty")
    }
    
    if(nrow(Performance) == 0) {
        stop("Performance data.frame is empty")
    }
    
    if(!"Fold" %in% names(Performance)) {
        stop("Performance must have 'Fold' column")
    }
    
    if(verbose) {
        cat("\n=== HYPERPARAMETER SELECTION ===\n")
        cat("Models analyzed:", length(List_models), "\n")
    }
    
    # === EXTRACT HYPERPARAMETERS FROM ALL MODELS ===
    rownames(Performance) <- Performance$Fold
    
    Model_parameters_df <- data.frame()
    for(i in seq_along(List_models)) { 
        Model <- List_models[[i]]
        model_params <- Model$bestTune
        # REMOVED: features_list (was never used)
        Model_parameters_df <- rbind(Model_parameters_df, model_params)
    }
    rownames(Model_parameters_df) <- names(List_models)
    
    # === MATCH PERFORMANCE DATA TO MODELS ===
    Performance <- Performance[match(rownames(Model_parameters_df), rownames(Performance)), ]
    Parameters_Performance <- cbind(Performance, Model_parameters_df)
    
    # === SELECT OPTIMAL VALUE FOR EACH PARAMETER ===
    vec_parameters <- colnames(Model_parameters_df)
    return_parameters_list <- vector("list", length(vec_parameters))
    
    if(verbose) {
        cat("Parameters to optimize:", paste(vec_parameters, collapse = ", "), "\n\n")
    }
    
    for(i in seq_along(vec_parameters)) {  # IMPROVED: use seq_along
        param_name <- vec_parameters[i]
        
        if(verbose) {
            cat("--- Parameter:", param_name, "---\n")
        }
        
        # Calculate mode (most frequent value)
        Modes <- calculate_mode(Parameters_Performance[[param_name]], 
                               return_all = TRUE)
        
        if(verbose) {
            freq_table <- table(Parameters_Performance[[param_name]])
            cat("Distribution:", paste(names(freq_table), "=", freq_table, collapse=", "), "\n")
            cat("Mode(s):", paste(Modes, collapse=", "), "\n")
        }
        
        if(length(Modes) == 1) {
            # Single mode - use it directly
            Param <- Modes
            
            if(verbose) {
                cat("Selection: ", Param, " (clear mode)\n\n", sep="")
            }
            
        } else {
            # Multiple modes (tie) - break tie using MCC/AUC
            if(verbose) {
                cat("Tie detected - using MCC/AUC to break tie\n")
            }
            
            mean_auc_mcc <- mean_auc_acc_param_function(
                parameter_list = Modes, 
                performance_parameters = Parameters_Performance, 
                parameter2check = param_name
            )
            
            Param_list <- select_best_model_hyperparameter(
                performance_summary = mean_auc_mcc,
                tie_threshold = tie_threshold,
                verbose = verbose  # IMPROVED: pass verbose flag
            )
            
            Param <- Param_list$optimal_parameter
        }
        
        return_parameters_list[[i]] <- Param
    }
    
    names(return_parameters_list) <- vec_parameters
    
    # === RETURN AS DATA.FRAME (for direct use in tuneGrid) ===
    optimal_params_df <- as.data.frame(return_parameters_list)
    
    if(verbose) {
        cat("\n=== OPTIMAL HYPERPARAMETERS ===\n")
        print(optimal_params_df)
        cat("\n")
    }
    
    return(optimal_params_df)  # IMPROVED: return data.frame instead of list
}

##########################################################
########      FINAL MODEL TRAINING FUNCTIONS      ########
##########################################################

## Train final model using mode hyperparameters and selected features ####
train_final_model <- function(Best_features = c(), df_best_hyperparameter = data.frame(), 
		in_df = c(), ML_method = "rf", 
		cv_folds = 5, cv_repeats = 5,
		is_imbalanced = FALSE) { 
  
	# No CV - just train with best hyperparameters
	#ctrl <- trainControl(method = "none")
	fitControl <- trainControl( method = "repeatedcv",number = cv_folds,repeats = cv_repeats,allowParallel = T)

  
	# Select columns: target + selected features
	DF_TRAIN <- in_df[, match(c("Variable", Best_features), colnames(in_df))]
  
	# Handle class imbalance if needed
	if (isTRUE(is_imbalanced)) { 
		cat("\n Data balanced \n")
		Train_selected <- smart_smote(data = DF_TRAIN) 
		BEST_MODEL <- train(Variable ~ ., data = Train_selected, method = ML_method, trControl = fitControl, tuneGrid = df_best_hyperparameter)
	} else {
	BEST_MODEL <- train(Variable ~ ., data = DF_TRAIN, method = ML_method, trControl = fitControl, tuneGrid = df_best_hyperparameter)
	}
  
	return(BEST_MODEL)
}


## Select the best model ####
select_best_model <- function(stats_model = data.frame()) {
  
  # Sort models by primary metric (MCC) in descending order
  # If there is a tie, use AUC as a secondary criterion
  stats_model_ordered <- stats_model[
    order(-stats_model$MCC, -stats_model$AUC),
  ]
  
  # Return the top-ranked model
  return(stats_model_ordered[1, ])
}


calc_stats <- function(x) {
    x <- x[!is.na(x)]  # Remover NAs
    c(
        N = length(x),
        Mean = mean(x),
        Median = median(x),
        SD = sd(x),
        Min = min(x),
        Max = max(x),
        Range = max(x) - min(x)
    )
}

#### Slect the best model according the performance of the rest of the outer loop cv
best_performances_outer_loop_function <- function( Stats_CV_model_FIT = data.frame(), Task = NA, list_CV_models = list() ){

	if(Task == "classification"){
		# AUC
		stats_AUC <- t(sapply(split(Stats_CV_model_FIT$AUC, Stats_CV_model_FIT$Model_Fold), calc_stats))
		stats_AUC <- as.data.frame(stats_AUC)
		stats_AUC$Model_Fold <- rownames(stats_AUC)
		stats_AUC <- stats_AUC[order(stats_AUC$Median,decreasing=T),]
		
		# MCC
		stats_MCC <- t(sapply(split(Stats_CV_model_FIT$MCC, Stats_CV_model_FIT$Model_Fold), calc_stats))
		stats_MCC <- as.data.frame(stats_MCC)
		stats_MCC$Model_Fold <- rownames(stats_MCC)
		stats_MCC <- stats_MCC[order(stats_MCC$Median,decreasing=T),]

		# Rename columns to distinguish MCC and AUC metrics
		colnames(stats_MCC) <- paste0("MCC_", colnames(stats_MCC))
		colnames(stats_AUC) <- paste0("AUC_", colnames(stats_AUC))

		# Merge both dataframes by Model_Fold
		stats_combined <- merge(stats_MCC, stats_AUC, by.x = "MCC_Model_Fold", by.y = "AUC_Model_Fold")
		colnames(stats_combined)[1] <- "Model_Fold"

		# Sort by: 
		#   1) MCC Median (descending) - primary criterion
		#   2) AUC Median (descending) - tiebreaker
		stats_combined <- stats_combined[order(-stats_combined$MCC_Median, -stats_combined$AUC_Median), ]
		write.table(stats_combined,"stats_summary_MCC_AUC_outerloop.tsv",col.names=T,row.names = F,quote=FALSE,sep = "\t")
		
		# Display sorted results with key metrics
		print(stats_combined[, c("Model_Fold", "MCC_Median", "MCC_Mean", "MCC_SD", "AUC_Median", "AUC_Mean", "AUC_SD")])

		# Report best model based on ranking criteria
		cat("\nBest Model:", as.character(stats_combined$Model_Fold[1]), "\n")
		cat("MCC Median:", round(stats_combined$MCC_Median[1], 4), "\n")
		cat("AUC Median:", round(stats_combined$AUC_Median[1], 4), "\n")
			
		best_model.FIT <- list_CV_models[[ stats_combined[1,]$Model_Fold ]]

		median_cv_AUC.FIT <- stats_combined$AUC_Median
		names(median_cv_AUC.FIT) <- stats_combined$Model_Fold
		median_cv_MCC.FIT <- stats_combined$MCC_Median
		names(median_cv_MCC.FIT) <- stats_combined$Model_Fold		
		

		Stats_CV_model_FIT$Model_Fold <- factor(as.character(Stats_CV_model_FIT$Model_Fold) , names(sort(median_cv_AUC.FIT,decreasing=T)) )
		pAUC <- ggplot(Stats_CV_model_FIT, aes(x=Model_Fold, y=AUC)) + geom_boxplot() + geom_jitter() + theme_bw() + ggtitle("FIT Models AUC (10 CV)")

		Stats_CV_model_FIT$Model_Fold <- factor(as.character(Stats_CV_model_FIT$Model_Fold) , names(sort(median_cv_MCC.FIT,decreasing=T)) )
		pMCC <- ggplot(Stats_CV_model_FIT, aes(x=Model_Fold, y=MCC)) + geom_boxplot() + geom_jitter() + theme_bw() + ggtitle("FIT Models MCC (10 CV)")

		p3 <- ggarrange(pMCC,pAUC)
		ggsave("Outerloop_MCC_AUC_boxplots.pdf",p3, width = 12, height = 5 )
		
	}
	
	if(Task == "regression"){
		median_cv_RMSE.FIT   <- c( by(Stats_CV_model_FIT$RMSE, Stats_CV_model_FIT$Model_Fold ,median) )


		best_model.FIT <- sort(median_cv_RMSE.FIT,decreasing=F)[1]

		

		Stats_CV_model_FIT$Model_Fold <- factor(as.character(Stats_CV_model_FIT$Model_Fold) , rev(names(sort(median_cv_RMSE.FIT,decreasing=T))) )
		pFIT <- ggplot(Stats_CV_model_FIT, aes(x=Model_Fold, y=RMSE_norm)) + geom_boxplot() + geom_jitter() + theme_bw() + 
			ggtitle("FIT Models RMSE norm (10 CV)")
		
		
		Stats_CV_model_PREDICT$Model_Fold <- factor(as.character(Stats_CV_model_PREDICT$Model_Fold) , rev(names(sort(median_cv_RMSE.PREDICT,decreasing=T))) )
		pPREDICT <- ggplot(Stats_CV_model_PREDICT, aes(x=Model_Fold, y=RMSE_norm)) + geom_boxplot() + geom_jitter() + theme_bw() +
			 ggtitle("PREDICT RMSE norm (10 CV)")
			 
		p3 <- ggarrange(pFIT,pPREDICT)
		ggsave("RMSE_boxplots.pdf",p3, width = 12, height = 5 )
	
	}	
	return(best_model.FIT)
} 


#############################################
########      Balance functions      ########
#############################################
check_and_report_balance <- function(tbl, 
                                     output_file = "data_balance_report.txt",
                                     verbose = TRUE) {
    # Check balance
    is_imbalanced <- (min(tbl) / max(tbl)) < (1/5)
    ratio <- min(tbl) / max(tbl)
    balance_status <- ifelse(is_imbalanced, "IMBALANCED", "BALANCED")
    
    # Print to console (brief)
    if(verbose) {
        cat("\nData is", tolower(balance_status), "\n")
        cat("min/max ratio:", round(ratio, 4), "\n")
    }
    
    # Write detailed report to file
    report <- c(
        "=== DATA BALANCE REPORT ===",
        "",
        paste("Status:", balance_status),
        paste("min/max ratio:", round(ratio, 4)),
        paste("Threshold: 0.2 (1/5)"),
        "",
        "Class distribution:",
        paste("  ", names(tbl), ":", tbl, 
              paste0("(", round(100*tbl/sum(tbl), 1), "%)")),
        "",
        paste("Total samples:", sum(tbl)),
        paste("SMOTE applied:", is_imbalanced),
        paste("Date:", Sys.time())
    )
    
    writeLines(report, output_file)
    
    if(verbose) {
        cat("Report saved to:", output_file, "\n")
    }
    
    return(is_imbalanced)
}





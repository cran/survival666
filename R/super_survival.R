#' @title
#' Batch survival analysis
#' @description
#' The super_survival() is the frist function of package survival666.
#' You can use the function to perform batch operations on the gene expression
#' matrix to obtain the survival analysis of the entire transcriptome
#' @param exp Gene expression matrix .We need to use the count matrix here.Note
#' that symbol is required for the column name and sample for the row name.
#' @param time Time data for survival analysis.You need to make sure the sample
#' as the row name is consistent with the gene expression matrix
#' @param status living status of patient in survival analysis,You need to make sure the sample
#' as the row name is consistent with the gene expression matrix too.
#' @param pval Full in the minimum significance you can accpet here.It depend
#' on your need, but the number should not be too big,too broad results are meaningless.
#' @param title Give a name to the file to be output
#' @param path Set a save path for the output file
#' @export
#' @examples
#' \dontrun{super_survival(exp=survivaldata,
#' time = clinic[,3],
#' status = clinic[,2],
#' pval = 0.0001,
#' title = 'survival_rank',
#' path = '/yourdir/')}
#' @importFrom survminer surv_cutpoint
#' @importFrom survminer surv_categorize
#' @importFrom survminer surv_fit
#' @importFrom survminer surv_pvalue
#' @importFrom survival Surv
#' @importFrom survival coxph
#' @importFrom utils write.table
#' @importFrom utils read.table
#' @importFrom utils write.csv


super_survival<-function(exp,
                         time,
                         status,
                         pval=0.0001,
                         title='survival_rank',
                         path){
  if(!dir.exists(path))dir.create(path)
  for (i in 1:ncol(exp)) {
    tryCatch({
      data_forsurv<-data.frame(exp=as.numeric(exp[,i]),
                               time,
                               status)
      index<-which(!is.na(data_forsurv$status)&data_forsurv$time>0)
      data_forsurv = data_forsurv[index,]

      surv.cut<-survminer::surv_cutpoint(data_forsurv,time = 'time',event = 'status',variables = 'exp',minprop = 0.3)

      opticut = surv.cut$cutpoint$cutpoint

      data_for_surv_final = survminer::surv_categorize(surv.cut, labels = c(0, 1))


      fit = survminer::surv_fit(survival::Surv(time, status)~exp ,
                                data = data_for_surv_final)
      pval = survminer::surv_pvalue(fit)$pval
      cox = survival::coxph(survival::Surv(time, status)~exp, data = data_for_surv_final)
      cox_summary = summary(cox)
      cox_res = c(cox_summary$coefficients[,c(1:2,5)], cox_summary$conf.int[,c(3,4)])
      surv_res = matrix(c(as.numeric(table(data_for_surv_final$exp)), pval, cox_res), nrow = 1)
      colnames(surv_res) = c("Num_Low", "Num_High", "Log-rank p", "coef", "HR", "Cox p", "Lower95%", "Upper95%")

      if (surv_res[1,3]<pval) {
        utils::write.table(data.frame(Symbol = colnames(exp)[i],surv_res),
                           file = paste0(path,  "surv_res.txt"),
                           row.names = F, col.names = F,
                           fileEncoding = "GBK",append = T)
      }
    },
    error=function(e){})
  }
  test<-utils::read.table(paste(path,'surv_res.txt',sep = ''))
  colnames(test)<-c("gene","Num_Low", "Num_High", "Log-rank p", "coef", "HR", "Cox p", "Lower95%", "Upper95%")
  test<-test[order(test$`Log-rank p`),]
  utils::write.csv(test,file = paste(path,title,'.csv',sep = ''),
                   row.names = F)
  file.remove(paste(path,'surv_res.txt',sep = ''))
}

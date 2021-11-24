#' @title
#' check the survival analysis
#' @description
#' The checkout_surviva() is the second function of package survival666.
#' We use the function to filter the results of super_survival().
#' @param exp This is a transcriptome count matrix,Note that the matrix
#' requires symbol as column name and sample as row name
#' @param survivalrank This is the result calculated by the frist function
#' in the survival666 package.
#' However,that does not prevent you from getting results that you could get
#' using other methods,as long as three conditions are met.
#' Frist,you need to calculate the p value based on the expression of a
#' particular phenotype and the level of gene expresion
#' Second,your results must be large enough,preferably for the entire
#' transcriptome,otherwise filtering accuracy will be reduced
#' Finally,you need to eliminate genes that are not significant enough,and if
#' you are using a regular PC,it is recommended to have no more than 2000 lines
#' @param pcol Digtal vector,represents the column of P value in survivalrank
#' @param symbolcol  Digtal vector,represents the column of symbol in survivalrank
#' @param path Path to output file,do not omit '/',example:'d:/file/'.
#' @export
#' @examples
#' \dontrun{checkout_survival(exp=survivaldata,
#' survivalrank=survivalrank[1:300,],
#' pcol=4,
#' symbolcol=1,
#' path='/yourdir/')}
#' @seealso
#' checkout_survival  \code{\link{super_survival}}
#' @importFrom stats cor
#' @importFrom utils write.csv


checkout_survival<-function(exp,survivalrank,pcol=4,symbolcol=1,path){
  if(!dir.exists(path))dir.create(path)
  if(length(survivalrank[,symbolcol][grep('^[0-9]',survivalrank[,symbolcol])])!=0){
    survivalrank<-survivalrank[-grep('^[0-9]',survivalrank[,symbolcol]),]
    exp<-exp[,survivalrank[survivalrank[,symbolcol]%in%colnames(exp),symbolcol]]
  }
  rownames(survivalrank)<-survivalrank[,symbolcol]
  survivalrank<-data.frame(survivalrank,confidence=NA)

  for (i in 1:nrow(survivalrank)) {
    tryCatch({
      p_value<-survivalrank[i,pcol]
      target_gene<-survivalrank[i,symbolcol]
      index<-which(survivalrank[,pcol] <p_value)
      if(length(index)<1){
        survivalrank[i,'confidence']<-survivalrank[i,pcol]
      }else{
        symbol<-survivalrank[index,symbolcol]

        data_test<-data.frame(pearson_target_gene=stats::cor(exp[,c(symbol,target_gene)])[,target_gene],
                              p_value=survivalrank[names(stats::cor(exp[,c(symbol,target_gene)])[,target_gene]),pcol])
        data_test<-data.frame(data_test,additional_P_value=10^(log10(data_test[,2])*abs(data_test[,1])))
        data_test<-data.frame(data_test,Reciprocal=1/data_test[,3])


        if(sum(data_test[index,'Reciprocal'])/data_test[target_gene,'Reciprocal']>=1){
          survivalrank[i,'confidence']<-'Need to test'
        }else{
          survivalrank[i,'confidence']<-1/(data_test[target_gene,'Reciprocal']-sum(data_test[index,'Reciprocal']))

        }
      }
    }
    ,error=function(e){print(survivalrank[i,symbolcol],'make mistakes')})
  }
  utils::write.csv(survivalrank,file =  paste0(path,'checkout.csv',''))
}

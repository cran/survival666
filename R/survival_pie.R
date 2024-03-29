#' @title
#' Visualize the factors that influence survival analysis
#' @description
#' The checkout_survival() is the third function of survival666.
#' We use the function to visualize the co-expression interference of a specific
#' gene in survival analysis.
#'
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
#' @param target target gene symbols
#' @param pcol Digtal vector,represents the column of P value in survivalrank
#' @param symbolcol  Digtal vector,represents the column of symbol in survivalrank
#' @param path Path to output file,do not omit '/',example:'d:/file/'.
#' @export
#' @examples
#' \dontrun{survival_pie(exp=survivaldata,
#' survivalrank=survivalrank[1:300,],
#' target="EAF2",
#' pcol=4,
#' symbolcol=1,
#' path='/yourdir/')}
#' @seealso
#' survival_pie  \code{\link{super_survival}}
#' @importFrom stats cor
#' @importFrom utils write.csv
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 coord_polar
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggsave


survival_pie<-function(exp,survivalrank,target,pcol=4,symbolcol=1,path){
  if(!dir.exists(path))dir.create(path)
  if(length(survivalrank[,symbolcol][grep('^[0-9]',survivalrank[,symbolcol])])!=0){
    survivalrank<-survivalrank[-grep('^[0-9]',survivalrank[,symbolcol]),]
    exp<-exp[,survivalrank[survivalrank[,symbolcol]%in%colnames(exp),symbolcol]]
  }
  rownames(survivalrank)<-survivalrank[,symbolcol]
  survivalrank<-data.frame(survivalrank,confidence=NA)

  for (i in target) {
    tryCatch({
      p_value<-survivalrank[i,pcol]
      target_gene<-survivalrank[i,symbolcol]
      index<-which(survivalrank[,pcol] <p_value)
      if(length(index)<1){
        print(paste(i,'is the symbol of maximum significance, no need to check'),sep = ' ')
      }else{
        symbol<-survivalrank[index,symbolcol]

        data_test<-data.frame(pearson_target_gene=stats::cor(exp[,c(symbol,target_gene)])[,target_gene],
                              p_value=survivalrank[names(stats::cor(exp[,c(symbol,target_gene)])[,target_gene]),pcol])
        data_test<-data.frame(data_test,additional_P_value=10^(log10(data_test$p_value)*abs(data_test$pearson_target_gene)))
        data_test<-data.frame(data_test,Reciprocal=1/data_test$additional_P_value)



          data_test<-data_test[order(data_test$Reciprocal,decreasing = T ),]
          utils::write.csv(data_test,paste0(path,i,'_check_procedure.csv',sep=''))
          dt1<-data.frame(
            gene=rownames(data_test)[2:11],
            Reciprocal=data_test[2:11,'Reciprocal']
          )
          p1<-ggplot2::ggplot(dt1, ggplot2::aes(x = "", y = dt1[,'Reciprocal'],fill=dt1[,'gene'])) +
            ggplot2::geom_bar(stat = "identity", width = 1) +
            ggplot2::coord_polar(theta = "y") +
            ggplot2::labs(x = "", y = "", fill = 'Genes',title = "The top ten most interfering genes")

          ggplot2::ggsave(filename = paste0(path,paste0(i,'p1.pdf')))
          dt2<-data.frame(
            gene=c(paste0('gene ',i),'additional significance'),
            Reciprocal=c(data_test[1,'Reciprocal'],sum(data_test[-1,'Reciprocal']))
          )
          p2<-ggplot2::ggplot(dt2, ggplot2::aes(x = "", y = dt2[,'Reciprocal'],fill=dt2[,'gene'])) +
            ggplot2::geom_bar(stat = "identity", width = 1) +
            ggplot2::coord_polar(theta = "y") +
            ggplot2::labs(x = "", y = "", fill = 'Genes',title = "The interference with target gene")+

          ggplot2::ggsave(filename = paste0(path,paste0(i,'p2.pdf')))

      }
    }
    ,error=function(e){print(paste0(survivalrank[i,symbolcol],'make mistakes'))})
  }
}

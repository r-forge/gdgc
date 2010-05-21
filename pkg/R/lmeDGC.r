#linkage='average'
#k=9
#n=3

#calcularQ=function(k,n,linkage="average")
#{
#Q=c()
#trat=as.factor(rep(seq(1,k,1),n))
#for (i in (1:5000))
#{
#y=rnorm(k*n,0,1)
#datos=as.data.frame(cbind(y,trat))
#datos$trat=as.factor(datos$trat)
#modelo=lm(y~trat,data=datos)
#matrizM=calcularmatrizM(modelo,datos)
#MDE=MatrizDeDiferenciasDeMediasEstandarizadas(modelo,datos,"trat",matrizM)
#D=as.dist(MDE$distancias)
#H=hclust(D,method=linkage)
#plot(H)
#Q=c(Q,max(H$height))
#}
#quantile(Q,c(0.90,0.95,0.99),type=6)
#}
#
#
#source("lmeLSmeans.r")


#matrizM<-calcularmatrizM(modelo,datos)


#termino="Tamano"


MatrizDeDiferenciasDeMediasEstandarizadas<-function(modelo,datos,termino,matrizM)
{
  factores<-strsplit(termino, ":")

  if (class(modelo)=="lme") coeficientes<-modelo$coefficients$fixed
  if (class(modelo)=="gls") coeficientes<-modelo$coefficients
  if (class(modelo)=="lm")  coeficientes<-modelo$coefficients[complete.cases(modelo$coefficients)]

  f<-factores[[1]]
  lista<-paste(rep(f[1],nrow(datos)),datos[,f[1]],sep='')
  if (length(f)>1) for (i in (2:length(f)))
    {
    lista=paste(lista,paste(rep(f[i],nrow(datos)),datos[,f[i]],sep=''),sep=':')
    }
   n=as.vector(table(lista))
  lista<-unique(lista)
  tabla<-as.vector(combinacion.lineal.para.obtener.media(lista[1],matrizM,coeficientes))

  if (length(lista)>1)
  for (i in (2:length(lista))) tabla<-rbind(tabla,as.vector(combinacion.lineal.para.obtener.media(lista[i],matrizM,coeficientes)))
  rownames(tabla)=lista

  if (class(modelo)=="lme") covar<-tabla%*%modelo$varFix%*%t(tabla)
  if (class(modelo)=="gls") covar<-tabla%*%modelo$varBeta%*%t(tabla)
  if (class(modelo)=="lm") covar<-tabla%*%vcov(modelo)%*%t(tabla)

  medias<-tabla%*%coeficientes;
  D<-as.matrix(dist(medias))

  orden=sort(medias,decreasing=TRUE,index.return=TRUE)$ix
  medias<-medias[orden]
  covar<-covar[orden,orden]
  D<-D[orden,orden]

 # Calculate standard error of mean differences
 # For that we will calculate the matrix of linear combiantions that
 # produde de differences, and apply this matrix to covariance matrix of
 # means

  for (k in (1:nrow(D)))
  {
   M<-D*0
   for (i in (1:nrow(D)))
        {
        for (j in (1:nrow(D)))
            if (k==j) M[i,j]=1
            else if (i>k) M[i,j]<-(-1)*as.integer(((j-k)==(i-k)))
                 else     M[i,j]<-(-1)*as.integer(((k-j)==(k-i)))
         }
# calcular los errores estandares de las diferencias con respecto a la k-esima media

   S<-1/sqrt(diag(M%*%covar%*%t(M)))
   D[k,]<-D[k,]*S
   rownames(D)=c()
   cn=colnames(D)
   l=strsplit(termino,":")[[1]]
   for (i in 1:length(l)) cn=sub(l[i],"",cn)
   colnames(D)=cn;
   rm(cn);

  }
 result=list(medias=medias,n=n,distancias=D,covarmedias=covar)
 result
}

lmeDGC<-function(MDE,Q,showplot=TRUE,Title="")
{
  D=MDE$distancias
  medias=MDE$medias
  covar=MDE$covarmedias
  n=round(mean(MDE$n))
#  q=calcularQ(nrow(D),n,"average")
# Q=q[2]
  D<-as.dist(D)
  H<-hclust(D,method = "average")
  if (showplot==TRUE) {plot(H,main=Title,sub="",xlab="",ylab="Q",cex=0.8);abline(h=Q)}
  indices=cutree(H,h=Q)
 result<-as.data.frame(cbind(medias,diag(sqrt(covar)),indices))
 colnames(result)=c('medias','ee','inidices')
 result
}
#MDE=MatrizDeDiferenciasDeMediasEstandarizadas(modelo,datos,"Tamano:Episperma",matrizM)
#lmeDGC(MDE,2.3)

 
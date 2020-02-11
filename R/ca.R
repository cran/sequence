ca <- function (x,nfac=3,isup=0,jsup=0,histev=FALSE,grr=FALSE,grc=FALSE,grrc=FALSE,grlist = rbind(c(1,2),c(1,3),c(2,3)),prtm=FALSE,prtevr=FALSE,prtevc=FALSE,eps=10E-10) 
{
currentOptions<-options()
on.exit(options(currentOptions))
options(width=160)
testlsup<-isup[1]!=0
testcsup<-jsup[1]!=0
#------------------------------------------------------------------------------
#					Parametres
#------------------------------------------------------------------------------
#				x	    : data.frame - minimal dimension 3x3
#				isup 	    : illustrative rows. If none 0, else list or vector with their numbers.
#				jsup	    :llustrative columns. Same as isup. 
#				histev = TRUE : Plot the histogram of eigenvalues. Default=TRUE
#				grr = FALSE0: plot graph of rows. Default = FALSE 
#				grc = FALSE : plot graph of columns. Default = FALSE
#				grrc = F: Affichage des graphiques simultanes
#				grlist = rbind(c(1,2),c(1,3),v(2,3)). Choice of factorial plans. Default = 1-2, 1-3, 2-3
#				prtm = TRUE : Print the data frame. Default=TRUE 
# 				prtevr : Print the eigenvectors for the rows
#                               prtevc : Print the eigenvectors for columns.
#				eps = 10E-10 : Accuracy for calculation of eigenvalues
#------------------------------------------------------------------------------
# 	Function partLines. Splits the matrix into active and illustrative 
#	rows (resp columns)
#------------------------------------------------------------------------------

partLines<- function (x,liste,cols=FALSE)	
# x : matrix to split
# liste : list of the items to extract and move into matrix z 
{
if (!is.matrix(x)) stop (message="***** Error : x not a matrix ****")
if(cols) x<-t(x)
n<-dim(x)[1] # number of rows, or number of columns if cols is TRUE
Un <- c(1:n)/c(1:n)
for (i in liste) Un[i]<-0
row2<-length(liste)	# size of the illustrative matrix
row1<-n-row2		# size of the active matrix
ncol<-dim(x)[2]
m<-list(y<-matrix(0),z<-matrix(0))
m$y<-matrix(0,row1,ncol)
m$z<-matrix(0,row2,ncol)
j<-0
k<-0
for (i in 1:n) if (Un[i]==1) {j<-j+1
			      m$y[j,]<-x[i,]
			     }
			else {k<-k+1
			      m$z[k,]<-x[i,]
			     }
if(cols) { m$y<-t(m$y)
           m$z<-t(m$z)
         }
m
}

				
#------------------------------------------------------------------------------
# 		initial settings			
#------------------------------------------------------------------------------
# Check the class of x
if(!is.data.frame(x)) stop(message = "\n ***** ERROR **** x not a data.frame \n")
if(length(dim(x))!=2) stop(message = "\n ***** ERROR **** x must have two dimensions \n")
if (!is.matrix(grlist)) stop (message="***** Erreur : grlist is not a matrix. Use rbind for a single graph. Example : grlist=rbind(c(1,2)) ****")
# Extraction of the matrix from dthe data.fame given in argument.
# Its first column is the list of alphanumeric identifier.
ligne <-  "\n----------------------------------------------------------------------------"
TitreG<-  "\n                        Correspondence analysis"
Souligne<-"\n               ---------------------------------------\n"
message(ligne,TitreG,Souligne)

ngraph<-dim(grlist)[1]

message("\n\nParameters :\n")
message("------------\n")
message("nfac =",nfac,"\tisup =",isup,"\tjsup =",jsup,"\n")
message("histev =",histev, "\tgrr =",grr,"\tgrc =",grc,"\tgrrc =",grrc,"\n")
message("ngraph =",ngraph,"\tgrlist =")
for (i in 1:ngraph) message(grlist[i,],"\t")
message("\nprtevr =",prtevr,"\tprtevc =",prtevc,"\tprtm=",prtm)
message("\teps =",eps,"\n")
id<-as.character(x[,1])
idC<-as.character(names(x)[2:length(names(x))])



#------------------------------------------------------------------------------
# 		extraction of the matrix from the data.frame
#------------------------------------------------------------------------------
n<-dim(x)[1]
p<-dim(x)[2]-1
M<-as.matrix(x[1:n,2:(p+1)])
if(testlsup){
		m<-partLines(M,isup)
		M<-m$y
		Mlsup<-as.matrix(m$z)
                # print(Mlsup)
		idsup<-partLines(as.matrix(id),isup)$z
		id<-partLines(as.matrix(id),isup)$y
		rm(m)
	    }
if(testcsup){
		m<-partLines(M,jsup,cols=TRUE)
		M<-m$y
		Mcsup<-as.matrix(m$z)
		idCsup<-partLines(t(as.matrix(idC)),jsup,cols=TRUE)$z
		idC<-partLines(t(as.matrix(idC)),jsup,cols=TRUE)$y
		if(exists("Mlsup")) Mlsup<-as.matrix(partLines(Mlsup,jsup,cols=TRUE)$y)
		rm(m)
	    }
TitreG<-  "\n                       Identifiers of active rows \n"
message(TitreG)
print (t(as.data.frame(id)))
TitreG<-  "\n                       Identifiers of active columns \n"
message(TitreG)
print (idC)
TitreG<-  "\n                       x: data matrix:\n"
Souligne<-   "                      --------------\n"
if(prtm) {message(TitreG, Souligne)
          print (x)
	  if(testlsup)
	  {
	  message("\n Matrix of illustrative rows\n")
	  print(as.data.frame(cbind(idsup,Mlsup)))
          }
	  if(testcsup)
          {
	  message("\n Matrix of illustrative columns\n")
	  print(as.data.frame(rbind(idCsup,Mcsup)))
          }		
	}
	
# 				    management of illustrative elements
if(testlsup || testcsup){
			message(ligne,"\nIllustrative elements:\n")

			if(exists("idsup")) {
					   message("\n Illustrative rows:\n")
				           print (idsup)
					  }
			if(exists("idCsup")) {
					   message("\n Illustrative column:\n")
				           print (idCsup)
					  }
                        }
# Calculation of matrix dimensions.
n<-dim(M)[1]
p<-dim(M)[2]
message("\nMatrix dimensions:\nn = ",n," ; p = ",p,"\n")
if((n<3)|| (p<3))stop(message = "\n ***** ERROR **** x must be at least of dimension 3 x 3 *****\n")
total<-sum(M)
message("\nTotal sum of the matrix: ",total,"\n")
#------------------------------------------------------------------------------
# 			Division by total sum
#------------------------------------------------------------------------------

M<-M/total
rm(x)
#------------------------------------------------------------------------------
# 					Margins
#------------------------------------------------------------------------------
#
# 					margins p.j

Unn<-c(1:n)/c(1:n)
Pj<-t(matrix(Unn)) %*% M
PDj<-diag(c(Pj))
inPj<-diag(c(1/Pj))
# 					pi. maargins
#                                    	Unp is the vector "one" of dimension p
Unp<-c(1:p)/c(1:p)
Pi<- M %*% matrix(Unp)

PDi<-diag(c(Pi))
inPi<-diag(c(1/Pi))

#------------------------------------------------------------------------------
# 				calculation of the matrix to diagonalize
#------------------------------------------------------------------------------

M2<-M %*% sqrt(inPj)

M3<-t(M2) %*% inPi %*% M2
rm(M2)


# 								Diagonalisation

L<-eigen(M3)
Lambda<-c(L$values)

Vec<-as.matrix(L$vectors)

message(ligne,"\n Eigenvalues \n")
message("\n\n\n Trivial eigenvalue verification: Lambda[1] = ",Lambda[1])
# print (Lambda)
#				eigenvalues cleanig and various calculations
pourcent<-c(1:p)*0
cumul<-pourcent
Inert<-sum(Lambda[2:p])

for (i in 2:p) { 
                 if(Lambda[i]<eps) Lambda[i]<-0
		 pourcent[i]<-(Lambda[i]/Inert)*100
		 if(i==2) cumul[i]<-pourcent[i] else cumul[i]<-pourcent[i]+cumul[i-1]
		}
#  conformance of the number of extractible factors
test<-(Lambda!=0)
nval<-sum(test)-1
if(nval<nfac) {
		message("\n *** WARNING ***\nYou asked ",nfac,
		    " factors, but there are only ",nval, 
		    "non zero eigenvalues\nnfac reduced to ",nval,"\n")
		nfac=nval
              }

if(nval<=1) stop(message="One single eigenvalue. CA cancelled \n")
#------------------------------------------------------------------------------
# 			Eigenvalues edition
#------------------------------------------------------------------------------

message("\n\n\n Statistics on eigenvalues:")
    message("\n --------------------------\n\n")
u<-as.data.frame(cbind(Lambda,pourcent,cumul))[2:p,]
if(histev){
	nev=dim(u)[1]
        dev.new(2)
        barplot(u[, 2],names.arg=as.character(2:(nev+1)), xlab = "rank", ylab = "Inertia (%)")
	}

print(format(u,digits=4))
rm(u)               
#
#					null eigenvalues are set to 1
#					to permit the following calculations
Lambda<-Lambda+(Lambda==0)
if(prtevr)
 	{
	message("\n    Rows eigenvectors\n")
	  message("    -----------------\n")
	u<-data.frame(Vec[,2:p])
	a<-as.character(c(1:(p-1)))
	a<-paste("V",a,sep="")
	names(u)<-a
	print (u)
	rm(u)
	}
#		Calculation of column eigenvectors
ILambdaD<-diag(1/Lambda)
# print(ILambdaD)
VecC<-sqrt(inPi) %*% M %*% sqrt(inPj) %*% Vec %*% sqrt(ILambdaD<-diag(1/Lambda))
if(prtevc)
 	{
	message("\n    Column eigenvectors\n")
	  message("    -------------------\n")
	u<-data.frame(VecC[,2:p])
	a<-as.character(c(1:(p-1)))
	a<-paste("V",a,sep="")
	names(u)<-a
	print (u)
	rm(u)
	}
#------------------------------------------------------------------------------
#			Calculation of factorial coordinates
#------------------------------------------------------------------------------

# 						Rows factorial coordinates

LFacL<- inPi %*% M %*% sqrt(inPj) %*% Vec
# 						Columns factorial corrdinates

LFacC<-inPj %*% t(M) %*% sqrt(inPi) %*%  VecC
#						Illustrative rows coordinates

if(testlsup)	{
		 nlsup<-length(isup)
		 Unnsup <- as.matrix(rep(1,p))
		 PDisup <- Mlsup %*% Unnsup
		 InPisup<- diag(c(1/PDisup),nrow=nlsup)
		 PDisup<-diag(c(PDisup),nrow=nlsup)
		 Mlsup<- InPisup %*% Mlsup %*% sqrt(inPj)
		 LFacLsup<- Mlsup %*% Vec
		} 


#						Illustrative columns coordinates
if(testcsup)	{
		 ncsup<-length(jsup)
		 Uncsup<-as.matrix(rep(1,n))
		 PDjsup<-t(Uncsup) %*% Mcsup
		 InPjsup<-diag(c(1/PDjsup),nrow=ncsup)
                 PDjsup<-diag(c(PDjsup),nrow=ncsup)
		 Mcsup<-sqrt(inPi) %*% Mcsup %*% InPjsup
		 LFaccsup<- t(Mcsup) %*% VecC   
		}

#------------------------------------------------------------------------------
#						Row graphs
#------------------------------------------------------------------------------
if(grr)
	{
	if(exists("nlsup")) fac<-rbind(LFacL,LFacLsup) else fac<-LFacL
	for (i in 1:ngraph)
		{
		dev.new()
		index<-ngraph+1-i
		f1<-grlist[index,1]+1
		f2<-grlist[index,2]+1
		xl<-paste("Factor",as.character(f1-1))
		yl<-paste("Factor",as.character(f2-1))
		plot(fac[,f1],fac[,f2],xlab=xl,ylab=yl,type="n")
		abline(0,0)
		abline(v=0)
		titre<-paste("Rows graph\n Axes",as.character(f1-1),"-",as.character(f2-1))
        	title(titre)
		for (k in 1:n) legend(LFacL[k,f1],LFacL[k,f2],as.character(id[k]),bty="n",xjust=0.5,yjust=0.5,pch=20)
		if(exists("nlsup")){
				for (k in 1:nlsup) 
				legend(LFacLsup[k,f1],LFacLsup[k,f2],
				as.character(idsup[k]),bty="n",xjust=0.5,yjust=0.5,pch=19,text.col="blue")
				}
		}
	rm(fac)
	}
#------------------------------------------------------------------------------
#				Columns graph
#------------------------------------------------------------------------------
if(grc)
	{
	if(exists("ncsup")) fac<-rbind(LFacC,LFaccsup) else fac<-LFacC
	for (i in 1:ngraph)
		{
		dev.new()
		index<-ngraph+1-i
		f1<-grlist[index,1]+1
		f2<-grlist[index,2]+1
		xl<-paste("Factor",as.character(f1-1))
		yl<-paste("Factor",as.character(f2-1))
		plot(fac[,f1],fac[,f2],xlab=xl,ylab=yl,type="n")
		abline(0,0)
		abline(v=0)
		titre<-paste("Columns graph \n Axes",
				as.character(f1-1),"-",as.character(f2-1))
        	title(titre)
		for (k in 1:p) legend(LFacC[k,f1],LFacC[k,f2],
				as.character(idC[k]),bty="n",xjust=0.5,yjust=0.5,pch=20)
		if(exists("ncsup")){
				    for (k in 1:ncsup) legend(LFaccsup[k,f1],LFaccsup[k,f2],
				    as.character(idCsup[k]),bty="n",xjust=0.5,yjust=0.5,pch=19,text.col="green")
				   }
		}
	rm(fac)
	}
#------------------------------------------------------------------------------
#				Combined graphs
#------------------------------------------------------------------------------
if(grrc)
	{
	for (i in 1:ngraph)
		{
		dev.new()
		index<-ngraph+1-i
		f1<-grlist[index,1]+1
		f2<-grlist[index,2]+1
		xl<-paste("Factor",as.character(f1-1))
		yl<-paste("Factor",as.character(f2-1))
		xm=min(LFacL[, f1])*1.1
	        xM=max(LFacL[, f1])*1.1
	        ym=min(LFacL[, f2])*1.1
	        yM=max(LFacL[, f2])*1.1
                plot(LFacL[, f1], LFacL[, f2], xlab = xl, ylab = yl, type = "n",xlim=c(xm,xM),ylim=c(ym,yM))
		abline(0,0)
		abline(v=0)
		titre<-paste("Combined graph \n Axes",as.character(f1-1),"-",as.character(f2-1))
        	title(titre)
		for (k in 1:n) legend(LFacL[k,f1],LFacL[k,f2],as.character(id[k]),bty="n",xjust=0.5,yjust=0.5)
		for (k in 1:p) legend(LFacC[k,f1],LFacC[k,f2],as.character(idC[k]),bty="o",bg="red",xjust=0.5,yjust=0.5)
		if(exists("nlsup")){
				    for (k in 1:nlsup) 
				    legend(LFacLsup[k,f1],LFacLsup[k,f2],
				    as.character(idsup[k]),bty="n",xjust=0.5,yjust=0.5,pch=19,text.col="blue")
				   }
		if(exists("ncsup")){
				    for (k in 1:ncsup) legend(LFaccsup[k,f1],LFaccsup[k,f2],
				    as.character(idCsup[k]),bty="o",bg="pink",xjust=0.5,yjust=0.5,pch=19,text.col="green")
				   }
		}
	}
#------------------------------------------------------------------------------
#				Calculation of rows loadings
#------------------------------------------------------------------------------

#			Factor loadings
# 			inertia of each row divided by the total inertia of the axis
#                       Sum to one for the set of rows

corL<- PDi %*% LFacL^2 %*% ILambdaD
#			Cosine squared
# 		 	association row-axis
#                       Sum to one for all axes
#
# 
#			Matrix of rows barycenters
B<-matrix(1,n,p) %*% PDj
B<-sqrt(B)
#
#			Matrix of scaled coordinates
#
R<-inPi %*% M %*% sqrt(inPj)
#
#			Matrix of z-scored coordinates
#
CR<-R-B
rm(R)
rm(B)
#
#			Squares of z-scored coordinates
#
CR<-CR^2
#
#			Calculation of rows loading
#
cos2L <- t(diag(c(1/(CR %*% Unp)))) %*% LFacL^2
#------------------------------------------------------------------------------
#				Calculation of illustrative rows loadings
#------------------------------------------------------------------------------
#					loadings
if(testlsup) {
	corLsup<- (PDisup/total) %*% LFacLsup^2 %*% ILambdaD
#					squared cosines
	B<-matrix(1,nlsup,p) %*% PDj
	B<-sqrt(B)
	R<-Mlsup 
	CR<-R-B
	rm(R)
	rm(B)
	CR<-CR^2
	cos2Lsup <- t(diag(c(1/(CR %*% Unp)),nrow=nlsup)) %*% LFacLsup^2
	}
#------------------------------------------------------------------------------
#				Calculation of columns loadings
#------------------------------------------------------------------------------

corC<- PDj %*% LFacC^2 %*% ILambdaD
B<-matrix(1,p,n) %*% PDi
B<-sqrt(B)
R<-inPj %*% t(M) %*% sqrt(inPi)
CR<-R-B
rm(R)
rm(B)
CR<-CR^2
cos2C <- diag(c(1/(CR %*% Unn))) %*% LFacC^2
#------------------------------------------------------------------------------
#				Calculation of illustrative columns loadings
#------------------------------------------------------------------------------
#						loadings
if(testcsup)
	{						
	corCsup<- (PDjsup/total) %*% LFaccsup^2 %*% ILambdaD
	B<-matrix(1,ncsup,n) %*% PDi
	B<-sqrt(B)
	R<-t(Mcsup)
	CR<-R-B
	rm(R)
	rm(B)
	CR<-CR^2
	cos2Csup <- diag(c(1/(CR %*% as.matrix(Unn))),nrow=ncsup) %*% LFaccsup^2
	}
#------------------------------------------------------------------------------
#				Summary of rows elements
#------------------------------------------------------------------------------
FL<-as.character(c(1:n)*0)
for (i in 1:nfac) 
		{
		j<-nfac+2-i
		FL<-cbind("|",round(LFacL[,j],4),round(corL[,j],4),round(cos2L[,j],4),FL)
                }
FL<-cbind(as.character(id),FL)
ncoledit<-1+4*nfac
u<-as.character(1:ncoledit)
u[]<-"======="
for(i in 1:nfac) u[(i-1)*4+2]<-"+"
v<-u
v[]<-"       "
for(i in 1:nfac) v[(i-1)*4+4]<-paste("f",i,sep="")
FL<-as.data.frame(rbind(as.character(u),FL[,1:ncoledit],as.character(u)))
message(ligne,"\n\n\n\n                  ROWS STATISTICS\n")
message("                  ---------------\n\n\n\n")

kchar<-NULL
lst<-c("|","FAC","LOAD","COS2")
for (i in 1:nfac) kchar<-cbind(kchar,t(lst))
kchar<-cbind("ident",kchar)
names(FL)<-as.character(kchar)
message(v,"\n")
print(format(FL,nsmall=6,digits=4))
rm(u)
rm(v)
rm(FL)
#------------------------------------------------------------------------------
#				Summary of illustrative rows elements
#------------------------------------------------------------------------------
if(testlsup)
{
FLsup<-as.character(c(1:nlsup)*0)
for (i in 1:nfac) 
		{
		j<-nfac+2-i
		FLsup<-cbind("|",round(LFacLsup[,j],4),round(corLsup[,j],4)
		,round(cos2Lsup[,j],4),FLsup)
                }
FLsup<-cbind(as.character(idsup),FLsup)
ncoledit<-1+4*nfac
u<-as.character(1:ncoledit)
u[]<-"======="
for(i in 1:nfac) u[(i-1)*4+2]<-"+"
v<-u
v[]<-"       "
for(i in 1:nfac) v[(i-1)*4+4]<-paste("f",i,sep="")
FLsup<-as.data.frame(rbind(as.character(u),FLsup[,1:ncoledit],as.character(u)))
message("\n                  ILLUSTRATIVE ROWS\n")
message(  "                  -----------------\n\n\n\n")

kchar<-NULL
lst<-c("|","FAC","LOAD","COS2")
for (i in 1:nfac) kchar<-cbind(kchar,t(lst))
kchar<-cbind("ident",kchar)
names(FLsup)<-as.character(kchar)
message(v,"\n")
print(format(FLsup,nsmall=6,digits=4))
rm(u)
rm(v)
rm(FLsup)
}
#------------------------------------------------------------------------------
#			Summary of column elements
#------------------------------------------------------------------------------
FC<-as.character(c(1:p)*0)
for (i in 1:nfac) 
		{
		j<-nfac+2-i
		FC<-cbind("|",round(LFacC[,j],4),round(corC[,j],4),round(cos2C[,j],4),FC)
                }
FC<-cbind(as.character(idC[1:p]),FC)
ncoledit<-1+4*nfac
u<-as.character(1:ncoledit)
u[]<-"======="
for(i in 1:nfac) u[(i-1)*4+2]<-"+"
v<-u
v[]<-"       "
for(i in 1:nfac) v[(i-1)*4+4]<-paste("f",i,sep="")
FC<-as.data.frame(rbind(as.character(u),FC[,1:ncoledit],as.character(u)))
message(ligne,"\n\n\n\n                  COLUMNS STATISTICS\n")
message("                  ------------------\n\n\n\n")
kchar<-NULL
lst<-c("|","FAC","LOAD","COS2")
for (i in 1:nfac) kchar<-cbind(kchar,t(lst))
kchar<-cbind("ident",kchar)
names(FC)<-as.character(kchar)
message(v,"\n")
print(format(FC,nsmall=6,digits=4))
rm(u)
rm(v)
rm(FC)
#------------------------------------------------------------------------------
#			Summary of illustrative columns
#------------------------------------------------------------------------------
#%
if(testcsup)
   {
	FCsup<-as.character(c(1:ncsup)*0)
	for (i in 1:nfac) 
			{
			j<-nfac+2-i
			FCsup<-cbind("|",round(LFaccsup[,j],4),round(corCsup[,j],4)
			,round(cos2Csup[,j],4),FCsup)
	                }
	idCsup<-as.matrix(idCsup)
	FCsup<-cbind(as.character(idCsup),FCsup)
	ncoledit<-1+4*nfac
	u<-as.character(1:ncoledit)
	u[]<-"======="
	for(i in 1:nfac) u[(i-1)*4+2]<-"+"
	v<-u
	v[]<-"       "
	for(i in 1:nfac) v[(i-1)*4+4]<-paste("f",i,sep="")
	FCsup<-as.data.frame(rbind(as.character(u),FCsup[,1:ncoledit],as.character(u)))
	message("\n\n\n\n                  ILLUSTRATIVE COLUMNS\n")
	message("                  --------------------\n")
	kchar<-NULL
	lst<-c("|","FAC","LOAD","COS2")
	for (i in 1:nfac) kchar<-cbind(kchar,t(lst))
	kchar<-cbind("ident",kchar)
	names(FCsup)<-as.character(kchar)
	message(v,"\n")
	print(format(FCsup,nsmall=6,digits=4))
	rm(u)
	rm(v)
	rm(FCsup)
   }
####################################################################################
#                                                                                  #
#    Rieturn of results in an object of class ca                                   #
#                                                                                  #
####################################################################################

if (any(isup !=0)) {Pisup <- PDisup %*% rep(1,length(isup))/(total+sum(PDisup))}
if (any(jsup !=0)) {Pjsup <- PDjsup %*% rep(1,length(jsup))/(total+sum(PDjsup))}

    res <- NULL
NoLsup<- all(isup==0)
PrincRows<-data.frame(id,Pi,LFacL[, 2:(nfac + 1)],rep("pri",n))
SupRows<-if(NoLsup) NULL else data.frame(idsup,Pisup,matrix(LFacLsup[,2:(nfac+1)],nlsup,nfac),rep("ill",nlsup))
if(!NoLsup) names(SupRows)<-names(PrincRows)
    res$fr <- rbind(PrincRows,SupRows)
    res$fr<-as.data.frame(res$fr)
    names(res$fr)<-c("id","w",paste(rep("f",nfac),as.character(c(1:nfac)),sep=""),"type")

NoCsup<- all(jsup == 0)
PrinCols<-data.frame(as.vector(idC),as.vector(Pj),LFacC[, 2:(nfac + 1)],rep("pri",p))
SupCols<-if(NoCsup) NULL else data.frame(as.vector(idCsup),as.vector(Pjsup),matrix(LFaccsup[, 2:(nfac + 1)],ncsup,nfac),rep("ill",ncsup ) )
if(!NoCsup) names(SupCols)<-names(PrinCols)
    res$fc <- rbind(PrinCols,SupCols)
    res$fc<-as.data.frame(res$fc)
    names(res$fc)<-c("id","w",paste(rep("f",nfac),as.character(c(1:nfac)),sep=""),"type")

    class(res) <- "ca"
    invisible(res)
}

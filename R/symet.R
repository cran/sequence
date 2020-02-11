symet <- function (x, sym=TRUE,charge=TRUE,ident=FALSE) 

{
# Verification of x status
if(!is.data.frame(x)) stop(message = "\n ***** ERROR **** x not a data.frame \n")
if(length(dim(x))!=2) stop(message = "\n ***** ERROR **** x must have two dimensions\n")
n<-dim(x)[1]
if(ident) p<-dim(x)[2]-1 else p<-dim(x)[2]
if(n != p) stop("*** ERROR *** this function works only on square matrices\n")
if(n <=3) stop ("*** ERROR *** dimension must be greater or equal to 3\n")
if (ident) M<-as.matrix(x[1:n,2:(p+1)]) else M<-as.matrix(x)
# symmetrization
if(sym) M <- M + t(M)
# loading
if(charge)  diag(M) <- M %*% rep(1,n)
y<-x
if(ident) y[,2:(n+1)]<-M else y<-M
invisible(y)
}

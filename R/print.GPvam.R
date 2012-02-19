print.GPvam <-
function (x, ...) 
{
    cat("Object of class 'GPvam'. Use functions 'plot' and 'summary'.", "\n")
    cat("Number of iterations: ",x$iter,"\n",sep="")
    cat("Log-likelihood: ",x$loglik,"\n",sep="")
    cat("Parameter estimates:\n")
    print(x$parameters)
    cat("\n")
    
}

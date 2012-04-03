plot.GPvam <-
function (x, ..., alpha=.1) 
{
    devAskNewPage(ask = TRUE)
   c.level <- qnorm(1 - alpha/2)
    for (i in 1:length(unique(x$teach.effects$teacher_year))) {
        for (j in i:length(unique(x$teach.effects$teacher_year))) {
            temp.df <- x$teach.effects[x$teach.effects$teacher_year == 
                i & x$teach.effects$effect_year == j, ]
            temp.df <- temp.df[order(temp.df$EBLUP), ]
            plotCI(temp.df$EBLUP, uiw = c.level * temp.df$std_error, 
                barcol = 2, xlab = "Ranked Teachers", ylab = "Teacher Effect")
            title(paste("Year ", i, " Teachers' Year ", j," Effect\nwith ",(1-alpha)*100,"% Confidence Intervals", sep = ""))
            abline(h = 0)
        }
    }
    qqnorm(x$cresid, main = "Normal Q-Q Plot\n for raw conditional residuals")
    qqline(x$cresid)
    qqnorm(x$sresid, main = "Normal Q-Q Plot\n for scaled conditional residuals")
    qqline(x$sresid)
    plot(x$yhat.s, x$sresid, main = "Scaled Conditional Residuals\n(by inverse Cholesky root\nof conditional error matrix)", 
        xlab = "Predicted", ylab = "Residuals")
    plot(x$yhat, x$cresid, main = "Raw conditional residuals", 
        xlab = "Predicted", ylab = "Residuals")
    plot(x$yhat.m, x$mresid, main = "Raw marginal residuals", 
        xlab = "Predicted", ylab = "Residuals")
    devAskNewPage(ask = FALSE)
    invisible(x)
}

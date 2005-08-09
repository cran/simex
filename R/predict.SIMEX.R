"predict.SIMEX" <-
function(object, newdata, ...)
{
new.object <- object$model
new.object$coefficients <- object$coefficients
predict(new.object, newdata= data.frame(newdata), ...)
}


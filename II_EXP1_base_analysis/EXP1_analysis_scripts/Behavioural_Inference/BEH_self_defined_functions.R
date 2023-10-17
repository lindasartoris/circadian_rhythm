pi_wrap <- function ( angle){ ####function that force any angle value to lie between -pi and + pi
  while(angle>pi){
    angle <- angle - 2*pi
  }
  while(angle<=-pi){
    angle <- angle + 2*pi
  }
  return(angle) 
}

####inclination_angle = function that converts an oriented angle difference between two vectors to an inclination angle (i.e. acute angle between the inclinations of the two lines on which the vectors lie)
#### values will range from 0 (parallel lines) to pi/2 (orthogonal lines)
inclination_angle <- function(oriented_angle){
  if (is.na(oriented_angle)){return(NA)}else{
    abs_angle <- abs(pi_wrap(oriented_angle))
    return (min (abs_angle, pi-abs_angle))
  }
}

## custom arrows function 
## quiver plot for vectors
arrows.az <- function(x, y, azimuth, rho, HeadWidth, ..., units=c("degrees", "radians"), Kol, Lwd)
{
  units <- match.arg(units)
  az.rad <- switch(units,
                   degrees=azimuth*pi/180,
                   radians=azimuth)
  arrows(x, y, x+cos(az.rad)*rho, y+sin(az.rad)*rho, length=HeadWidth,
         ..., col=Kol, lwd=Lwd)
}

# prepare for the display of classification performance (not used in the output)
display <- function(prediction, reference, Hitclass) {
  cm <- caret::confusionMatrix(data=prediction, reference=reference,
                               mode = "sens_spec", positive=Hitclass)
  print(cm$table)
  data.frame(Accuracy=cm$overall[1], Sensitivity=cm$byClass[1], Specificity=cm$byClass[2],
             row.names=NULL)
}

#output a text file from a list of objects (for SIRUS Rules output)
fnlist <- function(x, fil){ z <- deparse(substitute(x))
cat(z, "\n", file=fil)
nams=names(x) 
for (i in seq_along(x) ){ cat(nams[i], "\t",  x[[i]], "\n", 
                              file=fil, append=TRUE) }
}

trim_short_interactions <- function (interaction_table,trim_length_sec,duration_column_name){
  return (interaction_table[which(interaction_table[,duration_column_name]>trim_length_sec),])
}

# get R2 from regression
# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

lm_eqn_coef_p <- function (df)
{
  #### calculate the slope and p-value for simple linear regression ####
  # use example: project/scPipeline/fancy_analysis/mouse_human_comparison_ctrl.ipynb
  lmp <- function (modelobject) {
    # get p value of a lm model
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  #
  df$y <- df$value
  df$x <- df$percentage
  m = lm(y ~ x, df)
  p.raw <- lmp(m)
  if (is.na(lmp(m)))
    p.raw <- 1
  p <- signif(p.raw, 3)
  #     if (p.raw < 0.05) {
  #         p <- "p-value < 0.05 **"
  #     }
  #     else {
  #         p1 <- format(as.double(p.raw), digits = 2)
  #         p <- paste("p-value =", p1)
  #     }
  eq <- paste("Slope = ", format(as.double(coef(m)[2]), digits = 2),
              "\n", "P-value = ", p, sep = "")
  # as.character(eq)
  return(c(as.double(coef(m)[2]), p.raw, as.character(eq)))
}





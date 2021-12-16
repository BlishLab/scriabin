`%notin%` <- Negate(`%in%`)
resample <- function(x, ...) x[sample.int(length(x), ...)]

# .noGenerics <- TRUE
# .conflicts.OK <- TRUE

.onLoad <- .First.lib <- function(lib, pkg)
{
    library.dynam("pdfCluster", pkg, lib)
#library.dynam("func_c_spdep",pkg, lib)
}


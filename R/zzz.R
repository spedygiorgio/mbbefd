# File overview: package hooks called by R when the package is attached
# (library(mbbefd)) or unloaded.
# TODO: Add comment
# 
# Author: Giorgio
###############################################################################

#adding a startup message
#for future version, should test : if(verbose <- getOption("verbose")) 

# Print a start-up message with the package name, version, date (if any) and
# the URL where bugs can be reported. packageStartupMessage() lets users silence
# it with suppressPackageStartupMessages().
.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  # 'Date' is an optional DESCRIPTION field: do not print an empty line when absent
  datemsg <- if(!is.null(desc$Date)) paste0('Date:     ', desc$Date, '\n') else ''
  packageStartupMessage('Package:  ', desc$Package, '\n',
                        'Version:  ', desc$Version, '\n', 
                        datemsg,
                        'BugReport: ', desc$BugReports, '\n\n')
}




# for unloading dynamic libraries

# Unload the compiled C++ code (shared library) when the package is unloaded.
.onUnload <- function (libpath) {
  library.dynam.unload("mbbefd", libpath)
}
# TODO: Add comment
# 
# Author: Giorgio
###############################################################################

#adding a startup message
#for future version, should test : if(verbose <- getOption("verbose")) 

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

.onUnload <- function (libpath) {
  library.dynam.unload("mbbefd", libpath)
}
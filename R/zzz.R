# TODO: Add comment
# 
# Author: Giorgio
###############################################################################

#adding a startup message
#for future version, should test : if(verbose <- getOption("verbose")) 

.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage('Package:  ', desc$Package, '\n',
                        'Version:  ', desc$Version, '\n', 
                        'Date:     ', desc$Date, '\n',
                        'BugReport: ', desc$BugReports, '\n\n',
                        "Please note <d,q,p,r>mbbefd functions have been split into:\n", 
                        "<d,q,p,r>mbbefd and <d,q,p,r>MBBEFD, depending whether (a,b) or (g,b) are used.\n\n",
                        "Please update your code where necessary.")
}




# for unloading dynamic libraries

.onUnload <- function (libpath) {
  library.dynam.unload("mbbefd", libpath)
}
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
                        'BugReport: ', desc$BugReports, '\n',
                        "Please note <d,q,p,r>mbbefd functions have been split into  <d,q,p,r>mbbefd and  <d,q,p,r>MBBEFD",'\n',
                        "depending respectively whether (a,b) or (g,b) are used. Update your code whether necessary")
}




# for unloading dynamic libraries

.onUnload <- function (libpath) {
  library.dynam.unload("mbbefd", libpath)
}
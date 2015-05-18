# TODO: Add comment
# 
# Author: Giorgio
###############################################################################

#adding a startup message

.onLoad <- function(libname, pkgname) {
packageStartupMessage('This is mbbefd package. ','\n',
                        "Please note <d,q,p,r>mbbefd functions have been split into  <d,q,p,r>mbbefd and  <d,q,p,r>MBBEFD",'\n',
                        "depending respectively whether (a,b) or (g,b) are used. Update your code whether necessary"
  )
}


# for unloading dynamic libraries

.onUnload <- function (libpath) {
  library.dynam.unload("mbbefd", libpath)
}
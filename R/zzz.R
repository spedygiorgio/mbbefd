# TODO: Add comment
# 
# Author: Giorgio
###############################################################################


# display version number and date when the package is loaded
# .onAttach <- function(libname, pkgname) {
# 	desc  <- packageDescription(pkgname, libname)
# 	packageStartupMessage(
# 			'Version:  ', desc$Version, '\n', 
# 			'Date:     ', desc$Date, '\n',
# 			'Author:   ', 'Giorgio Alfredo Spedicato Ph.D, C.Stat ACAS'
# 	)
# }

#adding a startup message

.onLoad <- function(libname, pkgname) {
  
  desc  <- packageDescription(pkgname, libname)
  
  packageStartupMessage('This is package: ',desc$Package, '\n',
    'Version:  ', desc$Version, '\n', 
    'Date:     ', desc$Date, '\n'
  )
  message("Please note <d,q,p,r>mbbefd functions have been split into  <d,q,p,r>mbbefd and  <d,q,p,r>MBBEFD")
  message("depending the group of parameters being used. Update your code whether necessary")
}


# for unloading dynamic libraries

.onUnload <- function (libpath) {
  library.dynam.unload("mbbefd", libpath)
}
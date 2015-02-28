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

# for unloading dynamic libraries

.onUnload <- function (libpath) {
  library.dynam.unload("mbbefd", libpath)
}
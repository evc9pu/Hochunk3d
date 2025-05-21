AC_INIT(release.2009,version-1.0)

#
# the resulting Makefile (from running configure) is either serial or parallel.
# the same file cannot do both (need to rename and rerun, etc.)
#

AC_ARG_ENABLE(profiling, [  --enable-profiling      Enables the profiling flags]
)
AC_ARG_ENABLE(debug, [  --enable-debug          Enables the full debugging flags
])
AC_ARG_ENABLE(static, [  --enable-static         Enables static compilation])

# Static flags
if test "$enable_static" == "yes"
then  static_ifort='--static'
  static_nag=''  static_gfortran=''
  static_g95=''
  static_pgi='-Bstatic'
else
  static_ifort=''
  static_nag=''
  static_gfortran=''
  static_g95=''
  static_pgi=''
fi
# Profiling flags

if test "$enable_profiling" == "yes"
then
  profiling_ifort='-g -save-temps'
  profiling_nag=''
  profiling_gfortran=''
  profiling_g95=''
  profiling_pgi=''
else
  profiling_ifort=''
  profiling_nag=''
  profiling_gfortran=''
  profiling_g95=''
  profiling_pgi=''
fi

# Debugging flags

if test "$enable_debug" == "yes"
then
  debug_nag='-C=all -gline -nan'
  debug_ifort='-check all -warn all -warn nodec -warn interfaces -gen_interfaces -traceback -g'
  debug_gfortran='-Wall -fbounds-check -g'
  debug_g95='-Wall -fbounds-check -g'
  debug_pgi='-Mbounds -traceback -g'
else
  debug_nag=''
  debug_ifort=''
  debug_gfortran=''
  debug_g95=''
  debug_pgi=''
fi

extra_ifort="-heap-arrays -error-limit 1 $debug_ifort $profiling_ifort $static_ifort"
#extra_gfortran="-g -ffree-line-length-none $debug_gfortran $profiling_gfortran $static_gfortran"
extra_gfortran=" -ffree-line-length-none $debug_gfortran $profiling_gfortran $static_gfortran"
#extra_g95="-g -ffree-line-length-huge $debug_g95 $profiling_g95 $static_g95"
extra_g95="-ffree-line-length-huge $debug_g95 $profiling_g95 $static_g95"
extra_pgi="$debug_pgi $profiling_pgi $static_pgi"
extra_f95="-maxcontin=50 -f2003 -colour -ieee=full $debug_f95 $profiling_f95 $static_f95"

# define variables for either SERIAL or PARALLEL version:

AC_ARG_ENABLE(mpich,[  --enable-mpich          Use  MPICH mpif90 compiler])
if test "$enable_mpich" != "yes"

############################################################
# Serial version

then

AC_PROG_FC([ifort f95 gfortran xlf g95])
if test "$FC" == f95
then
  AC_SUBST(extra_serial, "$extra_nag")
elif test "$FC" == nagfor
then
  AC_SUBST(extra_serial, "$extra_nag")
elif test "$FC" == pgf90
then
  AC_SUBST(extra_serial, "$extra_pgi")
elif test "$FC" == pgf95
then
  AC_SUBST(extra_serial, "$extra_pgi")
elif test "$FC" == pgfortran
then
  AC_SUBST(extra_serial, "$extra_pgi")
elif test "$FC" == ifort
then
  AC_SUBST(extra_serial, "$extra_ifort")
elif test "$FC" == gfortran
then
  AC_SUBST(extra_serial, "$extra_gfortran")
elif test "$FC" == g95
then
  AC_SUBST(extra_serial, "$extra_g95")
else
  AC_MSG_ERROR([Not enabled for Fortran compiler: $FC])
fi
AC_SUBST(ttsre,"ttsre")
AC_SUBST(ttsre_message,"Serial version compiled: ttsre")

############################################################
#
# PARALLEL VERSION, if desired
#

else
   AC_MSG_NOTICE([Attempting to configure for parallel version])
   AC_CHECK_PROG(parallel_compiler,mpif90,true,false)
   if test "$parallel_compiler"
   then
      mpi_compiler=`mpif90 -show | awk {'print $1'}`
      AC_SUBST(mpi_fc, 'mpif90 -DMPI')
   else
      AC_ERROR(["mpif90 needs to be present"])
   fi

# Remove any leading path from the compiler
   mpi_compiler=`basename $mpi_compiler`

# gfortran sometimes includes version number, so check whether this is the case
   mpi_compiler_short=`echo $mpi_compiler | cut -c1-8`
   if test "$mpi_compiler_short" == gfortran
   then
     mpi_compiler=gfortran
   fi

   if test "$mpi_compiler" == f95
   then
     AC_SUBST(extra_mpi, "$extra_nag")
   elif test "$mpi_compiler" == nagfor
   then
     AC_SUBST(extra_mpi, "$extra_nag")
   elif test "$mpi_compiler" == pgf90
   then
     AC_SUBST(extra_mpi, "$extra_pgi")
   elif test "$mpi_compiler" == pgf95
   then
     AC_SUBST(extra_mpi, "$extra_pgi")
   elif test "$mpi_compiler" == pgfortran
   then
     AC_SUBST(extra_mpi, "$extra_pgi")
   elif test "$mpi_compiler" == ifort
   then
     AC_SUBST(extra_mpi, "$extra_ifort")
   elif test "$mpi_compiler" == gfortran
   then
     AC_SUBST(extra_mpi, "$extra_gfortran")
   elif test "$mpi_compiler" == g95
   then
     AC_SUBST(extra_mpi, "$extra_g95")
   else
     AC_MSG_ERROR([Could not determine MPI compiler: $mpi_compiler])
   fi

# i don't think that this is necessary for our application, but let's run
# the underlying fortran compiler through some tests
   echo "underyling fortran compiler is $mpi_compiler...\c"
   echo "Run some tests."
   AC_PROG_FC([$mpi_compiler])

# force compiler variables to be for a parallel compilation
   AC_SUBST(ttsre,"p_ttsre")
   AC_SUBST(ttsre_message,"Parallel version compiled: p_ttsre")
   AC_SUBST(FC, "$mpi_fc")
   AC_SUBST(extra_serial,"$extra_mpi")
fi


############################################################
# Check for FITSIO if requested

AC_ARG_ENABLE(fits,[  --enable-fits           Replaces ASCII output with FITS output])
AC_ARG_ENABLE(legacy_cfitsio,[  --enable-legacy-cfitsio Required to compile with older versions of cfitsio])

if test "$enable_fits" == "yes"
then

	AC_MSG_NOTICE([Using FITS output])

	echo 'checking for cfitsio... \c'

	AC_ARG_WITH(cfitsio,
		[  --with-cfitsio=<path>   path containing libcfitsio.a],
		[CFITSIO_lib_check="$with_cfitsio"],
		[CFITSIO_lib_check="/usr/local/lib/ /opt/local/lib/ /usr/lib/ /sw/lib/"])

	CFITSIO_libdir=
	for dir in $CFITSIO_lib_check
	do
		if test -d "$dir" && ( test -f "$dir/libcfitsio.a" )
		then
			CFITSIO_libdir=$dir
			break
		fi
	done

	if test -z "$CFITSIO_libdir"
	then
		echo 'no'
		AC_MSG_ERROR([If you have libcfitsio.a installed, please use the --with-cfitsio=path option to specify its parent directory])
	fi
		
	echo "$CFITSIO_libdir"/libcfitsio.a
	
	AC_SUBST(output,'output_FITS.o')
	AC_SUBST(nested_image,'type_nested_image.o')
	AC_SUBST(CFITSIO_libdir,"-L$CFITSIO_libdir -lcfitsio")
	if test "$enable_legacy_cfitsio" == "yes"
	then
		AC_SUBST(fits,'cfitsio_legacy.o lib_cfitsio.o')
	else
		AC_SUBST(fits,'lib_cfitsio.o')
	fi
else
	AC_MSG_NOTICE([Using ASCII output (specify --enable-fits to enable FITS file output)])
	AC_SUBST(output,'output.o')
	AC_SUBST(nested_image,'type_nested_image_nofits.o')
	AC_SUBST(CFITSIO_libdir,'')
	AC_SUBST(fits,'')
fi

# Output Makefile

AC_OUTPUT(Makefile)

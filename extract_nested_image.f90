program extract_nested_image
  
  use type_nested_image
  use constants
  use lib_array
  
  implicit none
  
  character(len=100) :: filename,filename_out,cwav,clevel
  real :: wav,nu
  integer :: level,inu,n_nu,unit
  type(nested_image) :: image
  real,allocatable :: array(:,:)
  
  call getarg(1,filename)
  call getarg(2,cwav)
  call getarg(3,clevel)
  call getarg(4,filename_out)
  
  read(cwav,*) wav
  read(clevel,*) level

  call nested_image_read(image,filename)
  
  nu = c / (wav * 1.e-4)
    
  inu = ipos(image%nu_min,image%nu_max,log10(nu),image%n_nu)
  
  allocate(array(image%imsize,image%imsize))
  
  call get_single_image(image,level,inu,array)
  
  call fits_open_new(unit,filename_out,confirm=.true.)
  call fits_write_primary_header(unit,-32,(/image%imsize,image%imsize/),.false.)
  call fits_write_array(unit,array)
  call fits_close(unit)
  
end program extract_nested_image
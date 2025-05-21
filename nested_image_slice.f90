program nested_slice

  use type_nested_image
  implicit none

  character(len=200) :: filename,cslice,output
  integer :: slice

  type(nested_image) :: image,image_new

  real :: nu_step

  call getarg(1,filename)
  call getarg(2,cslice) ; read(cslice,*) slice
  call getarg(3,output)

  call nested_image_read(image,filename)

  image_new = image
  deallocate(image_new%image)
  allocate(image_new%image(image%n_levels,1))
  image_new%image(:,1) = image%image(:,slice)

  image_new%n_nu = 1

  nu_step = (image%nu_max - image%nu_min) / real(image%n_nu)

  image_new%nu_min = image%nu_min + nu_step * real(slice-1)
  image_new%nu_max = image%nu_min + nu_step * real(slice)

  call nested_image_write(image_new,output)

end program nested_slice



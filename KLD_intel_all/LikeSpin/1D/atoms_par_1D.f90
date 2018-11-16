MODULE atoms_par_1D

IMPLICIT NONE

integer, parameter, private :: nSymbols = 32
character(len=2), private :: symbols(nSymbols)
double precision, private :: bsr(nSymbols)
integer, private :: n_bc(nSymbols)

data symbols / 'H ','Li','Be','B ','C ','N ','O ','F ',                  &
               'Na','Mg','Al','Si','P ','S ','Cl','K ',                  &
               'Ca','Sc','Ti','V ','Cr','Mn','Fe','Co',                  &
               'Ni','Cu','Zn','Ga','Ge','As','Se','Br'                   /

data bsr     / 0.661404096d0,1.370051342d0,0.992106144d0,0.803133545d0,  &
               0.661404096d0,0.614160946d0,0.566917797d0,0.472431497d0,  &
               1.700753390d0,1.417294491d0,1.181078743d0,1.039349294d0,  &
               0.944862994d0,0.944862994d0,0.944862994d0,2.078698587d0,  &
               1.700753390d0,3.023568000d0,2.645633000d0,2.551135500d0,  &
               2.645633000d0,2.645633000d0,2.645633000d0,2.551135500d0,  &
               2.551135500d0,2.551135500d0,2.551135500d0,1.228321893d0,  &
               1.138559908d0,1.086592443d0,1.086592443d0,1.086592443d0   /

data n_bc    / 20,25,25,25,25,25,25,25,30,30,30,30,30,30,30,35,35,35,35, &
               35,35,35,35,35,35,35,35,35,35,35,35,35                    /

CONTAINS

DOUBLE PRECISION FUNCTION get_bsr(atomic_symbol)

IMPLICIT NONE

character(len=2), intent(in) :: atomic_symbol

integer :: i

  get_bsr = 0.d0

  do i = 1, nSymbols
    if (symbols(i) == atomic_symbol) then
      get_bsr = bsr(i)
      exit
    end if
  end do

END FUNCTION get_bsr

INTEGER FUNCTION get_nbc(atomic_symbol)

IMPLICIT NONE

character(len=2), intent(in) :: atomic_symbol

integer :: i

  get_nbc = 0.d0

  do i = 1, nSymbols
    if (symbols(i) == atomic_symbol) then
      get_nbc = n_bc(i)
      exit
    end if
  end do

END FUNCTION get_nbc

END MODULE atoms_par_1D

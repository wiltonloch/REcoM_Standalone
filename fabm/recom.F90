#include "fabm_driver.h"

module awi_recom

    use fabm_types, only: type_base_model

    implicit none

    private

    type, extends(type_base_model), public :: type_awi_recom
    contains
        procedure :: initialize
    end type type_awi_recom

contains

    subroutine initialize(self, configunit)
        class(type_awi_recom), intent(inout), target :: self
        integer, intent(in) :: configunit

    end subroutine initialize

end module awi_recom

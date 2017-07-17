!> This module contains variables that describe the connectivity of the mesh and
!> also data structures for connectivity-related communication.
module mod_connectivity

   use mod_indices, only: ngridshi
   implicit none

   integer, dimension(2,-1:1,-1:1,ngridshi), save :: neighbor
   integer, dimension(2,0:3,0:3,ngridshi), save :: neighbor_child
   integer, dimension(-1:1,-1:1,ngridshi), save :: neighbor_type
   logical, dimension(-1:1,-1:1,ngridshi), save :: neighbor_active
   logical, dimension(-1:1,-1:1), save :: leveljump
   

   integer, dimension(ngridshi), save :: igrids
   integer, save :: igridstail

   integer, dimension(ngridshi), save :: igrids_active
   integer, save :: igridstail_active
   integer, dimension(ngridshi), save :: igrids_passive
   integer, save :: igridstail_passive

   integer, dimension(2), save :: nrecv_fc, nsend_fc

   integer, save :: nrecv_bc_srl, nsend_bc_srl, nrecv_bc_r, nsend_bc_r,&
       nrecv_bc_p, nsend_bc_p
!$OMP THREADPRIVATE(leveljump)

end module mod_connectivity

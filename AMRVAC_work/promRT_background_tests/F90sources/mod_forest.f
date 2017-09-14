module mod_forest
   implicit none

   type tree_node_ptr
      type(tree_node), pointer :: node
   end type tree_node_ptr

   type tree_node
      integer :: ig1,ig2, level, igrid, ipe
      logical :: leaf, active
      type(tree_node_ptr) :: parent, child(2,2), neighbor(2,2), next, prev
   end type tree_node

   type(tree_node_ptr), dimension(:,:), allocatable, save :: tree_root
   type(tree_node_ptr), dimension(:,:), allocatable, save :: igrid_to_node
   type(tree_node_ptr), dimension(:), allocatable, save :: level_head,&
       level_tail
   integer, dimension(:,:), allocatable, save :: sfc, sfc_iglevel1
   integer, dimension(:,:), allocatable, save :: iglevel1_sfc
   integer, dimension(:), allocatable, save :: sfc_to_igrid, igrid_to_sfc
   integer, dimension(:), allocatable, save :: sfc_phybound
   integer, dimension(:), allocatable, save :: Morton_start, Morton_stop

   integer, dimension(:), allocatable, save :: Morton_sub_start,&
       Morton_sub_stop

   logical, dimension(:,:), allocatable, save :: coarsen, refine, buffer,&
       igrid_inuse

   integer, save :: nleafs, nleafs_active, nglev1
   integer, dimension(:), allocatable, save :: nleafs_level

end module mod_forest

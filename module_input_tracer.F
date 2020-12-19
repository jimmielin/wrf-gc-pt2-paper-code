!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               __          _______  ______       _____  _____                 !
!               \ \        / /  __ \|  ____|     / ____|/ ____|                !
!                \ \  /\  / /| |__) | |__ ______| |  __| |                     !
!                 \ \/  \/ / |  _  /|  __|______| | |_ | |                     !
!                  \  /\  /  | | \ \| |         | |__| | |____                 !
!                   \/  \/   |_|  \_\_|          \_____|\_____|                !
!                                                                              !
!--------- v2.0 FOR PAPER SUBMITTION - PERMANENT ARCHIVAL (20201219) ----------!
!
! WRF-GC: GEOS-Chem-powered Chemistry Add-On for WRF Model
! 
! Version 2 developed by Xu Feng <fengx7@pku.edu.cn>, Haipeng Lin <hplin@g.harvard.edu>
!    April 2020, SUSTech, Atmospheric Chemistry and Climate Group
! 
! Originally developed by Haipeng Lin <hplin@g.harvard.edu>, Xu Feng <fengx7@pku.edu.cn>
!    January 2018, Peking University, Dept of Atmospheric and Oceanic Sciences
!    Correspondence to: Tzung-May Fu <fuzm@sustech.edu.cn>
!
! VERSION INFORMATION:
!    This copy of WRF-GC is permanently archived for submission of the WRF-GC, Pt. 2
!    two-way with ARI, ACI (chemical feedback) paper.
!
!    THIS VERSION IS PERMANENTLY ARCHIVED FOR SCIENTIFIC REPLICABILITY.
!    DO NOT USE THIS VERSION IN PRODUCTION, AS IT MAY NOT BE THE LATEST VERSION.
!
!    It may not be the latest version. We always recommend using the latest version
!    of WRF-GC, as important fixes may have been incorporated.
!    Please visit http://wrf.geos-chem.org for information.
!
! COPYRIGHT STATEMENT:
!    Permission is hereby granted, free of charge, to any person obtaining a copy
!   of this software and associated documentation files (the "Software"), to 
!   use, copy, modify the Software, and to permit persons to whom the Software is
!   furnished to do so, subject to the following conditions:
!
!   - The above copyright notice and this permission notice shall be included in all
!   copies or substantial portions of the Software.
! 
!   - The Software, modified in part or in full may not be redistributed without
!   express permission from the copyright holder(s).
! 
!   Except as contained in this notice or in attribution, the name of the WRF-GC model
!   shall not be used as an endorsement for distributing modified copies of the
!   Software without prior written permission from the copyright holder(s).
! 
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!   SOFTWARE.
! 
!  WRF and the GEOS-Chem model, GCHP are (c) their original authors.
!
!------------------------------------------------------------------------------- 

!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: module_input_tracer.F
!
! !DESCRIPTION: Module Input Tracer inputs tracer quantities from WRF input/bdy data
!  into the chem data structure in WRF.
!\\
!\\
! !INTERFACE:
module module_input_tracer
!
! !USES:
!
  use module_state_description, only: param_first_scalar
!
! !REMARKS:
!  This module provides interfaces to existing WRF endpoints, adapted to work with (GI)GC.
!
! !REVISION HISTORY:
!  29 May 2018 - H.P. Lin  - First crack, based on WRF-Chem implementation.
!------------------------------------------------------------------------------
!BOC
contains
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initialize_tracer
!
! !DESCRIPTION: Subroutine initializes tracer quantities in WRF-GC. If chem_in_opt
!  enabled, data is read from initial conditions based on input_chem_data, and as
!  such the routine does nothing. Tracer options in original WRF-Chem are unsupported.
!\\
!\\
! !INTERFACE:
!
   subroutine initialize_tracer(chem, chem_in_opt, &
                                tracer_opt, num_chem, &
                                ids, ide, jds, jde, kds, kde, & ! domain dims
                                ims, ime, jms, jme, kms, kme, & ! memory dims
                                ips, ipe, jps, jpe, kps, kpe, & ! patch  dims
                                its, ite, jts, jte, kts, kte)
!
! !INPUT PARAMETERS:
!
      integer, intent(in)    :: chem_in_opt, tracer_opt, num_chem
      integer, intent(in)    :: ids, ide, jds, jde, kds, kde
      integer, intent(in)    :: ims, ime, jms, jme, kms, kme
      integer, intent(in)    :: ips, ipe, jps, jpe, kps, kpe
      integer, intent(in)    :: its, ite, jts, jte, kts, kte
      real, dimension(ims:ime, kms:kme, jms:jme, num_chem), intent(inout) :: chem
      if (chem_in_opt == 1) return
   end subroutine initialize_tracer
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: flow_dep_bdy_tracer
!
! !DESCRIPTION: This sets zero gradient conditions for outflow and a set profile
!  for inflow in the boundary specified region. Note that the field must be un-
!  staggered. Velocities u, v are used to check sign. spec_zone is width of outer
!  specified bcs set here. (Based on WRF-Chem code, JD Aug 2000)
!\\
!\\
! !INTERFACE:
!
   subroutine flow_dep_bdy_tracer(chem, &
                                  chem_bxs, chem_btxs, &
                                  chem_bxe, chem_btxe, &
                                  chem_bys, chem_btys, &
                                  chem_bye, chem_btye, &
                                  dt, &
                                  spec_bdy_width, z, &
                                  have_bcs_chem, &
                                  u, v, tracer_opt, alt, &
                                  t, pb, p, t0, p1000mb, rcp, ph, phb, g, &
                                  spec_zone, ic, &
                                  ids, ide, jds, jde, kds, kde, & ! domain dims
                                  ims, ime, jms, jme, kms, kme, & ! memory dims
                                  ips, ipe, jps, jpe, kps, kpe, & ! patch  dims
                                  its, ite, jts, jte, kts, kte)

      implicit none
!
! !INPUT PARAMETERS:
!
      integer, intent(in)    :: tracer_opt
      integer, intent(in)    :: ids, ide, jds, jde, kds, kde
      integer, intent(in)    :: ims, ime, jms, jme, kms, kme
      integer, intent(in)    :: ips, ipe, jps, jpe, kps, kpe
      integer, intent(in)    :: its, ite, jts, jte, kts, kte
      integer, intent(in)    :: spec_zone, spec_bdy_width, ic
      real, intent(in)       :: dt

      real, dimension(jms:jme, kds:kde, spec_bdy_width), intent(in) :: chem_bxs, chem_bxe, chem_btxs, chem_btxe
      real, dimension(ims:ime, kds:kde, spec_bdy_width), intent(in) :: chem_bys, chem_bye, chem_btys, chem_btye
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: z
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: alt
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: u
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: v
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: ph, phb, t, pb, p
      real, intent(in) :: g, rcp, t0, p1000mb
!
! !INPUT/OUTPUT PARAMETERS:
!
      real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: chem
!
! !LOCAL VARIABLES:
!
      integer    :: i, j, k, numgas
      integer    :: ibs, ibe, jbs, jbe, itf, jtf, ktf
      integer    :: i_inner, j_inner
      integer    :: b_dist
      integer    :: i_bdy_method
      real       :: tempfac, convfac
      real       :: tracer_bv_def
      logical, optional :: have_bcs_chem
!
! !REMARKS:
!  This is based upon the WRF-Chem original code (c) original authors.
!
! !REVISION HISTORY:
!  29 May 2018 - H.P. Lin   - First crack, based upon original WRF-Chem code.
!EOP
!------------------------------------------------------------------------------
!BOC
      tracer_bv_def = 1.e-30
      ibs = ids
      ibe = ide - 1
      itf = min(ite, ide - 1)
      jbs = jds
      jbe = jde - 1
      jtf = min(jte, jde - 1)
      ktf = kde - 1

      i_bdy_method = 0
      if (have_bcs_chem) i_bdy_method = 6
      if (ic .lt. param_first_scalar) i_bdy_method = 0

      if (jts - jbs .lt. spec_zone) then
      ! y-start boundary
         do j = jts, min(jtf, jbs + spec_zone - 1)
            b_dist = j - jbs
            do k = kts, ktf
               do i = max(its, b_dist + ibs), min(itf, ibe - b_dist)
                  i_inner = max(i, ibs + spec_zone)
                  i_inner = min(i_inner, ibe - spec_zone)
                  if (v(i, k, j) .lt. 0.) then
                     chem(i, k, j) = chem(i_inner, k, jbs + spec_zone)
                  else
                     if (i_bdy_method .eq. 0) then
                        chem(i, k, j) = tracer_bv_def
                     else if (i_bdy_method .eq. 6) then
                        call bdy_tracer_value(chem(i, k, j), chem_bys(i, k, 1), chem_btys(i, k, 1), dt, ic)
                     else
                        chem(i, k, j) = tracer_bv_def
                     endif
                  endif
               enddo
            enddo
         enddo
      endif
      if (jbe - jtf .lt. spec_zone) then
      ! y-end boundary
         do j = max(jts, jbe - spec_zone + 1), jtf
            b_dist = jbe - j
            do k = kts, ktf
               do i = max(its, b_dist + ibs), min(itf, ibe - b_dist)
                  i_inner = max(i, ibs + spec_zone)
                  i_inner = min(i_inner, ibe - spec_zone)
                  if (v(i, k, j + 1) .gt. 0.) then
                     chem(i, k, j) = chem(i_inner, k, jbe - spec_zone)
                  else
                     if (i_bdy_method .eq. 0) then
                        chem(i, k, j) = tracer_bv_def
                     else if (i_bdy_method .eq. 6) then
                        call bdy_tracer_value(chem(i, k, j), chem_bye(i, k, 1), chem_btye(i, k, 1), dt, ic)
                     else
                        chem(i, k, j) = tracer_bv_def
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

      if (its - ibs .lt. spec_zone) then
      ! x-start boundary
         do i = its, min(itf, ibs + spec_zone - 1)
            b_dist = i - ibs
            do k = kts, ktf
               do j = max(jts, b_dist + jbs + 1), min(jtf, jbe - b_dist - 1)
                  j_inner = max(j, jbs + spec_zone)
                  j_inner = min(j_inner, jbe - spec_zone)
                  if (u(i, k, j) .lt. 0.) then
                     chem(i, k, j) = chem(ibs + spec_zone, k, j_inner)
                  else
                     if (i_bdy_method .eq. 0) then
                        chem(i, k, j) = tracer_bv_def
                     else if (i_bdy_method .eq. 6) then
                        call bdy_tracer_value(chem(i, k, j), chem_bxs(j, k, 1), chem_btxs(j, k, 1), dt, ic)
                     else
                        chem(i, k, j) = tracer_bv_def
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

      if (ibe - itf .lt. spec_zone) then
      ! x-end boundary
         do i = max(its, ibe - spec_zone + 1), itf
            b_dist = ibe - i
            do k = kts, ktf
               do j = max(jts, b_dist + jbs + 1), min(jtf, jbe - b_dist - 1)
                  j_inner = max(j, jbs + spec_zone)
                  j_inner = min(j_inner, jbe - spec_zone)
                  if (u(i + 1, k, j) .gt. 0.) then
                     chem(i, k, j) = chem(ibe - spec_zone, k, j_inner)
                  else
                     if (i_bdy_method .eq. 0) then
                        chem(i, k, j) = tracer_bv_def
                     else if (i_bdy_method .eq. 6) then
                        call bdy_tracer_value(chem(i, k, j), chem_bxe(j, k, 1), chem_btxe(j, k, 1), dt, ic)
                     else
                        chem(i, k, j) = tracer_bv_def
                     endif
                  endif
               enddo
            enddo
         enddo
      endif
   end subroutine flow_dep_bdy_tracer
!EOC

   subroutine set_tracer(dtstep, ktau, pbl_h, tracer, t, tracer_opt, num_tracer, &
                         z, ht, ids, ide, jds, jde, kds, kde, &
                         ims, ime, jms, jme, kms, kme, &
                         its, ite, jts, jte, kts, kte)
      integer, intent(in)    :: ktau, tracer_opt, num_tracer
      integer, intent(in)    :: ids, ide, jds, jde, kds, kde
      integer, intent(in)    :: ims, ime, jms, jme, kms, kme
      integer, intent(in)    :: its, ite, jts, jte, kts, kte
      real, dimension(ims:ime, kms:kme, jms:jme, num_tracer), intent(inout) :: tracer
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: t, z
      real, dimension(ims:ime, jms:jme), intent(in) :: PBL_H, HT
      real, intent(in) :: dtstep
      integer :: count_trop, count_pbl

   end subroutine set_tracer

   subroutine bdy_tracer_value(trac, trac_b, trac_bt, dt, ic)
      implicit none

      real, intent(out)   :: trac
      real, intent(in)    :: trac_b, trac_bt
      real, intent(in)    :: dt
      integer, intent(in) :: ic

      real :: epsilc = 1.e-12

      trac = max(epsilc, trac_b + trac_bt*dt)

      return
   end subroutine bdy_tracer_value

END MODULE module_input_tracer

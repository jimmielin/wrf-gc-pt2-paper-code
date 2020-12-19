!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GC_Stateful_Mod.F
!
! !DESCRIPTION: Module GIGC\_Stateful\_Mod stores derived type objects related
!  to Grid-Independent GEOS-Chem (GIGC)'s State structures for a coupled multi-domain
!  model. (hplin, 6/4/18)
!\\
!\\
! !INTERFACE:
module GC_Stateful_Mod
!
! !USES:
!
    use ErrCode_Mod
    use Input_Opt_Mod
    use State_Met_Mod
    use State_Chm_Mod
    use State_Diag_Mod
    use State_Grid_Mod
    use DiagList_Mod, only: DgnList
    use HCO_TYPES_MOD, only: ConfigObj
    use HCO_State_Mod, only: HCO_State
    use HCOX_State_Mod, only: Ext_State

    implicit none
    private

!
! !PUBLIC MEMBER FUNCTIONS:
!
    public :: GIGC_State_Boot

    public :: GIGC_State_Init
    public :: GIGC_State_Final

    public :: GIGC_State_Get_Status
    public :: GIGC_State_Get_Opt
    public :: GIGC_State_Get_DiagList
    public :: GIGC_State_Get_Met
    public :: GIGC_State_Get_Chm
    public :: GIGC_State_Get_Diag
    public :: GIGC_State_Get_Grid
    public :: GIGC_State_Get_HCO
    public :: GIGC_State_Get_HCOX

    public :: GIGC_State_Set_Opt
    public :: GIGC_State_Set_Met
    public :: GIGC_State_Set_Chm
    public :: GIGC_State_Set_Diag
    public :: GIGC_State_Set_Grid

    ! HCO / HCOX are passed as POINTERs. For updates they can be performed in-place,
    ! but after initialization they need to be saved, hence setter methods need
    ! to be available. (hplin, 2/26/19)
    public :: GIGC_State_Set_HCO
    public :: GIGC_State_Set_HCOX

!
! !PRIVATE TYPES:
!
    ! Derived type for domain properties
    type GIGC_Stateful_Object
        integer                                :: ID = -999
        logical                                :: Init = .false.
        type(MetState)                         :: State_Met
        type(ChmState)                         :: State_Chm
        type(DgnState)                         :: State_Diag
        type(GrdState)                         :: State_Grid
        type(HCO_State), pointer               :: HcoState
        type(Ext_State), pointer               :: ExtState
    end type GIGC_Stateful_Object

    ! Global options objects
    type(OptInput)                             :: Global_Input_Opt
    type(ConfigObj), pointer                   :: Global_HcoConfig => NULL()
    type(DgnList)                              :: Global_DiagList

    ! Stateful objects
#if defined ( EXTERNAL_GRID ) || defined( MODEL_ )
    !-----------------------------------------------------------------
    !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
    !
    ! The max number of domain configurations can be flexible when coupling
    ! with an external model. (hplin, 6/7/18)
    !-----------------------------------------------------------------
    integer                                    :: EXTERNAL_MAX_DOM = 8
    type(GIGC_Stateful_Object), dimension(1:8) :: States
#else
    !-----------------------------------------------------------------
    !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
    !
    ! For GEOS-Chem "Classic", the world is seen as a contiguous set of
    ! columns with no domain extensions possible.
    !-----------------------------------------------------------------
    integer                                    :: EXTERNAL_MAX_DOM = 1
    type(GIGC_Stateful_Object), dimension(1:1) :: States
#endif
    logical                                    :: Init = .false.

!
! !REMARKS:
!  This module is intended to be the ONLY stateful module in WRF-GC, meaning that all other
!  modules should not be statefully programmed (in an ideal world, of course) & only rely on
!  this module to retrieve/set derived type objects with storage within CPUs.
!
!  Note: Fortran passes derived type objects by intrinsic copy, not by reference.
!  This means that there are some maneuverings required in your coding, specifically relating to
!  how to work with Input_Opt, a global variable. Namely, the workflow in external models would be
!  (EM - External Model, CH - GIGC_CHUNK_MOD, ST - GC_Stateful_Mod)
!
!  > On Initialization <
!   EM - check if initialized GIGC_State, if not
!         ST - call GIGC_State_Boot, initializes Global_Input_Opt, Global_HcoConfig
!      - call GIGC_Chunk_Init, do all the chunk initialization routines
!         CH - operate on datas and return
!      - call GIGC_State_Init, ... to store data into appropriate containers by ID - opens a new container
!      - if needing other operations, call GIGC_State_Get_... to get initialized datas or _Set_ to set.
!
!  > On Column-Code Running <
!   EM - call GIGC_State_Get_... to get initialized datas to pass into GIGC_CHUNK_MOD
!         CH - operate on datas and return
!   EM - call GIGC_State_Set_Opt, Met, Chm ... to store data into appropriate containers by ID
!      - operate as necessary
!      - call GIGC_State_Set_Opt, Met, Chm ... to store data if operated on, otherwise skip
!
!  Note the drastic workflow change in GIGC_Chunk_Init: This means that the part where options are read
!  is NO LONGER in GIGC_Chunk_Init but part of an upper model. This facilitates future I/O changes, as
!  I/O should be really handled by a parent model or a gridded component anyway.
!
! !REVISION HISTORY:
!  04 Jun 2018 - H.P. Lin  - First crack.
!  07 Jun 2018 - H.P. Lin  - Added workflow description and clarified code.
!  28 Dec 2019 - H.P. Lin  - Update to GEOS-Chem 12.6.3
!------------------------------------------------------------------------------
!BOC
contains
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Boot
!
! !DESCRIPTION: Subroutine GIGC\_State\_Boot initializes the stateful module and
!  read configuration variables to store in Input_Opt, HcoConfig. (hplin, 6/7/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Boot(am_I_Root,  MPI_COMM, &
                               State_Grid, RC)
!
! !USES:
! 
        USE GC_Environment_Mod, only: GC_Allocate_All
        USE INPUT_MOD,          only: Read_Input_File
        USE HCO_CONFIG_MOD,     only: Config_Readfile
        USE HCO_DIAGN_MOD,      only: DiagnFileOpen
        USE LINOZ_MOD,          only: Linoz_Read
        USE DiagList_Mod,       only: Init_DiagList, Print_DiagList

!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: MPI_COMM           ! MPI Communicator #
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(GrdState), intent(inout) :: State_Grid         ! Grid State object
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  This code was written with the intent to supercede GIGC_Get_Options, with
!  the role of reading configuration files & initializing global vars such as
!  Input_Opt, HcoConfig, ... into Stateful_Mod. Note that in this design, you
!  cannot change "input.geos"/"HEMCO_Config.rc" across domains, so they are still
!  not independent.
!
!  It is OK to call this method more than once - it will not crash and instead work
!  gracefully.
!
! !REVISION HISTORY:
!  07 Jun 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer :: HCO_DIAGN_LUN

        ! Assume success
        RC = GC_SUCCESS

        ! If we are already initialized, return
#if defined ( LINUX_GFORTRAN )
        if(Init .eqv. .true.) then
#else
        if(Init .eq. .true.) then
#endif
            return
        endif

        ! Initialize Input_Opt fields to zeros or equivalent (v12)
        call Set_Input_Opt( am_I_Root, Global_Input_Opt, RC )

        ! Some necessary set-up for GEOS-Chem HP
#if defined ( EXTERNAL_GRID ) || defined( MODEL_ )
        !-----------------------------------------------------------------
        !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
        !
        ! We have to set up some HPC and Root CPU parameters into Input_Opt
        ! for GEOS-Chem HP / GIGC. This must be done after GC_Allocate_All
        ! so it is not overwritten by GC_Allocate_All. (hplin, 6/13/18)
        !-----------------------------------------------------------------
        Global_Input_Opt%isMPI   = .true.
        Global_Input_Opt%AmIRoot = am_I_Root
        Global_Input_Opt%MPIComm = MPI_COMM
        ! Global_Input_Opt%NPES  = ... ! To set!

        ! Set some DEFAULT time-steps which will be overwritten by GIGC_Chunk_Mod later on
        Global_Input_Opt%TS_CHEM = 10   ! Chemistry timestep [min]
        Global_Input_Opt%TS_EMIS = 10   ! Chemistry timestep [min]
        Global_Input_Opt%TS_DYN  = 20   ! Dynamic   timestep [min]
        Global_Input_Opt%TS_CONV = 20   ! Dynamic   timestep [min]
#endif

        ! Read input.geos, now done on all threads since GC(HP) 12.2.0
        ! to remove GIGC MPI wrapper dependency
        ! Read input.geos at very beginning of simulation on every thread
        CALL Read_Input_File( am_I_Root, Global_Input_Opt, State_Grid, RC )
        if(RC /= GC_SUCCESS) then 
            write(6, *) "STOP GC_Stateful_Mod :: Return Code /= GC_SUCCESS (Read_Input_File)"
            return
        endif
        ! Echo info
        write(6, *) '### GC_Stateful_Mod: Root CPU, after READ_INPUT_FILE'

        ! In the ESMF/MPI environment, we can get the total overhead ozone
        ! either from the met fields (GIGCsa) or from the Import State (GEOS-5)
        Global_Input_Opt%USE_O3_FROM_MET = .TRUE.

        ! Read LINOZ climatology
        ! This is here because it is read into Input_Opt (hplin, 12/28/19)
        IF ( Global_Input_Opt%LLINOZ ) THEN
            CALL Linoz_Read( am_I_Root, Global_Input_Opt, RC ) 
            if(RC /= GC_SUCCESS) then 
                write(6, *) "STOP GC_Stateful_Mod :: Return Code /= GC_SUCCESS (Linoz_Read)"
                return
            endif

            ! Echo info
            write(6, *) '### GC_Stateful_Mod: Root CPU, after LINOZ_READ'
        ENDIF

        ! Read HEMCO Configuration File
        ! No longer required as EMISSIONS_INIT will handle it (for now) and pass back later.
        ! (hplin, 8/6/18)
        ! call Config_ReadFile(am_I_Root, Global_HcoConfig, Global_Input_Opt%HcoConfigFile, 0, RC)

        ! Read HEMCO Diagnostic File
        ! call DiagnFileOpen(am_I_Root, Global_HcoConfig, HCO_DIAGN_LUN, RC)

        ! Read HISTORY.rc & initialize diagnostics list object
        ! Ported from gigc_historyexports_mod, v12
        CALL Init_DiagList(am_I_Root, "HISTORY.rc", Global_DiagList, RC)
        if(RC /= GC_SUCCESS) then 
            write(6, *) "STOP GC_Stateful_Mod :: Return Code /= GC_SUCCESS (Init_DiagList)"
            return
        endif

        CALL Print_DiagList(am_I_Root, Global_DiagList, RC)

        !
        ! Allocate GEOS-Chem module arrays,
        !
        ! GC_Allocate_All's task is to initialize the CMN_ module arrays.
        ! Input_Opt is no longer initialized as of GEOS-Chem v12, so this should be called
        ! AFTER Input_Opt is read, as some options depend on Input_Opt%LUCX, for example.
        ! (hplin, 8/7/18)
        !
        ! Generic sizings are passed in -- they will be switched by GIGC_Switch_Dims at runtime.
        !
        ! Allocate all lat/lon arrays
        CALL GC_Allocate_All( am_I_Root, Global_Input_Opt, State_Grid, RC )
        if(RC /= GC_SUCCESS) then 
            write(6, *) "STOP GC_Stateful_Mod :: Return Code /= GC_SUCCESS (GC_Allocate_All)"
            return
        endif

        ! We are finished!
        if(RC /= GC_SUCCESS) then 
            write(6, *) "STOP GC_Stateful_Mod :: Return Code /= GC_SUCCESS"
        endif
        
        ! We are now initialized!
        Init = .true.
    
    end subroutine GIGC_State_Boot
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Init
!
! !DESCRIPTION: Subroutine GIGC\_State\_Init initializes a new stateful array object
!  and stores the given data inside the container. (hplin, 6/7/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Init(am_I_Root, &
                               ID, State_Met, State_Chm, State_Diag, State_Grid, &
                               HcoState, ExtState, &
                               RC)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
        type(MetState), intent(in)    :: State_Met
        type(ChmState), intent(in)    :: State_Chm
        type(DgnState), intent(in)    :: State_Diag
        type(GrdState), intent(in)    :: State_Grid
        type(HCO_State), pointer      :: HcoState
        type(Ext_State), pointer      :: ExtState
!
! !INPUT/OUTPUT PARAMETERS:
!
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs.
!
! !REVISION HISTORY:
!  07 Jun 2018 - H.P. Lin   - Initial version.
!  28 Dec 2019 - H.P. Lin   - Adapt to State_Grid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Assume success
        RC = GC_SUCCESS

        ! Check if this domain ID has been initialized already...
        ! Note: To simply the code flow we have repeated the loop twice for different uses.
        ! This makes the code clearer. Note that we don't need to worry about slowness here
        ! since this routine doesn't run more than EXTERNAL_MAX_DOM times, also
        ! it should be quite small (In WRF, it is at most 24, we are using 8 for GCHP)
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        if(Found_Index .ne. -1) then
            write(6, *) "%%% GC_Stateful_Mod: Detected duplicate initialization of ID", ID, "%%%"
            write(6, *) "Check your code for GIGC_State_Init duplicate calls on the same ID."

            write(6, *) "You are attempting to initialize ID", ID
            write(6, *) "But domains already stored are at States%ID", States%ID
            write(6, *) "Of which Found_Index =", Found_Index, "matches ID..."

            ! FIXME hplin: Should we fail silently or store State_Met etc. for this call?
            stop
        endif

        ! Check for a free slot.
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(.not. States(D)%Init) then
                Found_Index = D
                exit
            endif
        enddo

        if(Found_Index .eq. -1) then
            write(6, *) "%%% GC_Stateful_Mod: We have no more space (max number of ext. domains is", EXTERNAL_MAX_DOM
            write(6, *) "You might need to change the headers in GC_Stateful_Mod.f90"

            write(6, *) "You are attempting to initialize ID", ID
            write(6, *) "But domains already stored are at States%ID", States%ID
            stop
        endif

        ! Write into this slot.
        States(D)%Init       = .true.
        States(D)%ID         = ID
        States(D)%State_Met  = State_Met
        States(D)%State_Chm  = State_Chm
        States(D)%State_Diag = State_Diag
        States(D)%State_Grid = State_Grid
        States(D)%HcoState   => HcoState
        States(D)%ExtState   => ExtState
    end subroutine GIGC_State_Init
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Final
!
! !DESCRIPTION: Subroutine GIGC\_State\_Final deallocates the stateful array object
!  from the container. (hplin, 6/10/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Final(am_I_Root, &
                                ID, RC)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
!
! !INPUT/OUTPUT PARAMETERS:
!
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs.
!
! !REVISION HISTORY:
!  10 Jun 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Assume success
        RC = GC_SUCCESS

        ! Find for this domain ID
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        ! If this domain ID is unavailable, then throw an error
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GC_Stateful_Mod: The ID requested to be finalized in GIGC_State_Final could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif

        ! Deallocate state structures regarding this ID... This is from GIGC_Chunk_Final
        ! FIXME (hplin): Is not finished yet

        ! Finalize GEOS-Chem
        ! Deallocate fields of the Chemistry State object
        CALL Cleanup_State_Chm(am_I_Root, States(Found_Index)%State_Chm, RC)
        IF (am_I_Root) write(6, *) 'GC_Stateful_Mod State_Chm Finalize ID =', ID

        ! Deallocate fields of the Meteorology State object
        CALL Cleanup_State_Met(am_I_Root, States(Found_Index)%State_Met, RC)
        IF (am_I_Root) write(6, *) 'GC_Stateful_Mod State_Met Finalize ID =', ID

        ! Deallocate fields of the Diagnostics State object
        CALL Cleanup_State_Diag(am_I_Root, States(Found_Index)%State_Diag, RC)
        IF (am_I_Root) write(6, *) 'GC_Stateful_Mod State_Diag Finalize ID =', ID
    end subroutine GIGC_State_Final
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Get_Status
!
! !DESCRIPTION: Subroutine GIGC\_State\_Get\_Status verifies a domain ID's status
!  inside the container. (hplin, 6/10/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Get_Status(am_I_Root, ID, Status)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
!
! !INPUT/OUTPUT PARAMETERS:
!
        logical, intent(inout)        :: Status             ! Have we found it?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs, which is why this function exists.
!
! !REVISION HISTORY:
!  10 Jun 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Find for this domain ID
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        Status = Found_Index .ne. -1
    end subroutine GIGC_State_Get_Status
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Get_Opt
!
! !DESCRIPTION: Subroutine GIGC\_State\_Get_Opt gets Input_Opt for GLOBAL.
!  (hplin, 6/10/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Get_Opt(am_I_Root, Input_Opt, HcoConfig)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)                :: am_I_Root          ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(OptInput), intent(inout)      :: Input_Opt          ! Input_Opt object.
        type(ConfigObj), pointer, optional :: HcoConfig          ! HEMCO Configuration object (optional)
!
! !REMARKS:
!  Note that the Input Options in GIGC are global across GIGC_Stateful_Object.
!
! !REVISION HISTORY:
!  10 Jun 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
        ! Are we initialized? If not, we cannot proceed.
        if(.not. Init) then
            write(6, *) "%%% GC_Stateful_Mod: Not initialized but requested GIGC_State_Get_Opt. Stop."
            stop
        endif

        Input_Opt = Global_Input_Opt

        ! Do we need to return the HEMCO Configuration object?
        if(present(HcoConfig)) then
            HcoConfig => Global_HcoConfig
        endif
    end subroutine GIGC_State_Get_Opt
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Get_DiagList
!
! !DESCRIPTION: Subroutine GIGC\_State\_Get_DiagList gets DiagList for GLOBAL.
!  (hplin, 6/10/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Get_DiagList(am_I_Root, DiagList)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)                :: am_I_Root          ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(DgnList), intent(inout)      :: DiagList           ! DiagList object.
!
! !REMARKS:
!  Note that the diagnostics list in GIGC is global across GIGC_Stateful_Object.
!
! !REVISION HISTORY:
!  06 Aug 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
        ! Are we initialized? If not, we cannot proceed.
        if(.not. Init) then
            write(6, *) "%%% GC_Stateful_Mod: Not initialized but requested GIGC_State_Get_DiagList. Stop."
            stop
        endif

        DiagList = Global_DiagList
    end subroutine GIGC_State_Get_DiagList
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Get_Met
!
! !DESCRIPTION: Subroutine GIGC\_State\_Get\_Met gets the meterology state from
!  the container. (hplin, 6/10/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Get_Met(am_I_Root, ID, State_Met, RC)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(MetState), intent(inout) :: State_Met          ! Met state.
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs.
!
! !REVISION HISTORY:
!  10 Jun 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Assume success
        RC = GC_SUCCESS

        ! Find for this domain ID
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        ! If this domain ID is unavailable, then throw an error
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif

        State_Met = States(Found_Index)%State_Met
    end subroutine GIGC_State_Get_Met
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Get_Chm
!
! !DESCRIPTION: Subroutine GIGC\_State\_Get\_Chm gets the chemistry state from
!  the container. (hplin, 6/10/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Get_Chm(am_I_Root, ID, State_Chm, RC)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(ChmState), intent(inout) :: State_Chm          ! Chemistry state.
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs.
!
! !REVISION HISTORY:
!  10 Jun 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Assume success
        RC = GC_SUCCESS

        ! Find for this domain ID
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        ! If this domain ID is unavailable, then throw an error
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif

        State_Chm = States(Found_Index)%State_Chm
    end subroutine GIGC_State_Get_Chm
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Get_Diag
!
! !DESCRIPTION: Subroutine GIGC\_State\_Get\_Diag gets the diagnostics state from
!  the container. (hplin, 6/10/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Get_Diag(am_I_Root, ID, State_Diag, RC)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(DgnState), intent(inout) :: State_Diag         ! Diagnostic state.
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs.
!
! !REVISION HISTORY:
!  10 Jun 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Assume success
        RC = GC_SUCCESS

        ! Find for this domain ID
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        ! If this domain ID is unavailable, then throw an error
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif

        State_Diag = States(Found_Index)%State_Diag
    end subroutine GIGC_State_Get_Diag
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Get_Chm
!
! !DESCRIPTION: Subroutine GIGC\_State\_Get\_Grid gets the grid state from
!  the container. (hplin, 6/10/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Get_Grid(am_I_Root, ID, State_Grid, RC)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(GrdState), intent(inout) :: State_Grid         ! Grid state.
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs.
!
! !REVISION HISTORY:
!  28 Dec 2019 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Assume success
        RC = GC_SUCCESS

        ! Find for this domain ID
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        ! If this domain ID is unavailable, then throw an error
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif

        State_Grid = States(Found_Index)%State_Grid
    end subroutine GIGC_State_Get_Grid
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Get_HCO 
!
! !DESCRIPTION: Subroutine GIGC\_State\_Get\_HCO gets the HEMCO state from
!  the container to be wrapped into HCO_INTERFACE_MOD. (hplin, 6/10/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Get_HCO(am_I_Root, ID, HcoState, RC)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(HCO_State), pointer      :: HcoState           ! HEMCO state.
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs. Also, HcoState and ExtStates are pointers.
!
! !REVISION HISTORY:
!  10 Jun 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Assume success
        RC = GC_SUCCESS

        ! Find for this domain ID
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        ! If this domain ID is unavailable, then throw an error
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif

        HcoState => States(Found_Index)%HcoState
    end subroutine GIGC_State_Get_HCO
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Get_HCOX
!
! !DESCRIPTION: Subroutine GIGC\_State\_Get\_HCOX gets the HEMCO extensions state
!  from the container to be wrapped into HCO_INTERFACE_MOD. (hplin, 6/10/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Get_HCOX(am_I_Root, ID, ExtState, RC)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(Ext_State), pointer      :: ExtState           ! HEMCO extensions state.
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs. Also, HcoState and ExtStates are pointers.
!
! !REVISION HISTORY:
!  10 Jun 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Assume success
        RC = GC_SUCCESS

        ! Find for this domain ID
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        ! If this domain ID is unavailable, then throw an error
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif

        ExtState => States(Found_Index)%ExtState
    end subroutine GIGC_State_Get_HCOX
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Set_Opt
!
! !DESCRIPTION: Subroutine GIGC\_State\_Set_Opt sets Input_Opt for GLOBAL.
!  (hplin, 6/11/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Set_Opt(am_I_Root, Input_Opt, HcoConfig)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)                :: am_I_Root          ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(OptInput), intent(inout)      :: Input_Opt          ! Input_Opt object.
        type(ConfigObj), pointer, optional :: HcoConfig          ! HEMCO Configuration object (optional)
!
! !REMARKS:
!  Note that the Input Options in GIGC are global across GIGC_Stateful_Object.
!
! !REVISION HISTORY:
!  11 Jun 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
        ! Are we initialized? If not, we cannot proceed.
        if(.not. Init) then
            write(6, *) "%%% GC_Stateful_Mod: Not initialized but requested GIGC_State_Set_Opt. Stop."
            write(6, *) "(Developer Note: Are you trying to set Input_Opt manually instead? This option"
            write(6, *) " is currently unsupported, but you can edit GC_Stateful_Mod manually.)"
            stop
        endif

        Global_Input_Opt = Input_Opt

        ! If we have HcoConfig, point it
        if(present(HcoConfig)) then
            Global_HcoConfig => HcoConfig
        endif
    end subroutine GIGC_State_Set_Opt
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Set_Met
!
! !DESCRIPTION: Subroutine GIGC\_State\_Set\_Met sets the meterology state from
!  the container. (hplin, 6/10/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Set_Met(am_I_Root, ID, State_Met, RC)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(MetState), intent(inout) :: State_Met          ! Met state.
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs.
!
! !REVISION HISTORY:
!  10 Jun 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Assume success
        RC = GC_SUCCESS

        ! Find for this domain ID
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        ! If this domain ID is unavailable, then throw an error
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif

        States(Found_Index)%State_Met = State_Met
    end subroutine GIGC_State_Set_Met
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Set_Chm
!
! !DESCRIPTION: Subroutine GIGC\_State\_Set\_Chm sets the chemistry state from
!  the container. (hplin, 6/10/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Set_Chm(am_I_Root, ID, State_Chm, RC)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(ChmState), intent(inout) :: State_Chm          ! Chemistry state.
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs.
!
! !REVISION HISTORY:
!  10 Jun 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Assume success
        RC = GC_SUCCESS

        ! Find for this domain ID
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        ! If this domain ID is unavailable, then throw an error
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif

        States(Found_Index)%State_Chm = State_Chm
    end subroutine GIGC_State_Set_Chm
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Set_Diag
!
! !DESCRIPTION: Subroutine GIGC\_State\_Set\_Diag sets the diagnostics state from
!  the container. (hplin, 6/13/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Set_Diag(am_I_Root, ID, State_Diag, RC)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(DgnState), intent(inout) :: State_Diag         ! Diagnostics state.
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs.
!
! !REVISION HISTORY:
!  13 Jun 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Assume success
        RC = GC_SUCCESS

        ! Find for this domain ID
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        ! If this domain ID is unavailable, then throw an error
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif

        States(Found_Index)%State_Diag = State_Diag
    end subroutine GIGC_State_Set_Diag
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Set_Grid
!
! !DESCRIPTION: Subroutine GIGC\_State\_Set\_Grid sets the grid state from
!  the container. (hplin, 6/13/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Set_Grid(am_I_Root, ID, State_Grid, RC)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(GrdState), intent(inout) :: State_Grid         ! Grid state.
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs.
!
! !REVISION HISTORY:
!  28 Dec 2019 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Assume success
        RC = GC_SUCCESS

        ! Find for this domain ID
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        ! If this domain ID is unavailable, then throw an error
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif

        States(Found_Index)%State_Grid = State_Grid
    end subroutine GIGC_State_Set_Grid
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Set_HCO 
!
! !DESCRIPTION: Subroutine GIGC\_State\_Set\_HCO sets the HEMCO state from
!  the container to be wrapped into HCO_INTERFACE_MOD. (hplin, 2/26/19)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Set_HCO(am_I_Root, ID, HcoState, RC)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(HCO_State), pointer      :: HcoState           ! HEMCO state.
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs. Also, HcoState and ExtStates are pointers.
!
! !REVISION HISTORY:
!  26 Feb 2019 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Assume success
        RC = GC_SUCCESS

        ! Find for this domain ID
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        ! If this domain ID is unavailable, then throw an error
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif

        States(Found_Index)%HcoState => HcoState
    end subroutine GIGC_State_Set_HCO
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Set_HCOX
!
! !DESCRIPTION: Subroutine GIGC\_State\_Set\_HCOX sets the HEMCO extensions state
!  from the container to be wrapped into HCO_INTERFACE_MOD. (hplin, 2/26/19)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Set_HCOX(am_I_Root, ID, ExtState, RC)
!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(Ext_State), pointer      :: ExtState           ! HEMCO extensions state.
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  Note that the array indeces in GC_Stateful_Mod::States might not necessarily
!  correspond to their domain IDs. Also, HcoState and ExtStates are pointers.
!
! !REVISION HISTORY:
!  26 Feb 2019 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer                       :: D
        integer                       :: Found_Index = -1

        ! Assume success
        RC = GC_SUCCESS

        ! Find for this domain ID
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo

        ! If this domain ID is unavailable, then throw an error
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif

        States(Found_Index)%ExtState => ExtState
    end subroutine GIGC_State_Set_HCOX
!EOC
end module GC_Stateful_Mod
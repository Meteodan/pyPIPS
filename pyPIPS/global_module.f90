!########################################################################
!########################################################################
!#########                                                      #########
!#########                  module rsa_table                    #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

MODULE rsa_table

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine read in the precalculated scattering amplitude
! from the tables, which are calcualted using T-matrix method.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 4/17/2008
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  SAVE

!-----------------------------------------------------------------------
! PARAMETER
! dsr : rain drop size,  rsa: scattering amplitude for rain drop
! dss : snow drop size,  ssa: scattering amplitude for snow aggregate
! dsh : hail drop size,  hsa: scattering amplitude for hailstone
! dsg : grpl drop size,  gsa: scattering amplitude for graupel  
!-----------------------------------------------------------------------
  INTEGER, PARAMETER :: nd = 112, ns = 8, nfw = 21, nvarden = 19
  INTEGER, PARAMETER :: nd_r = 100 ! 100 
  INTEGER, PARAMETER :: nd_s = 112
  INTEGER, PARAMETER :: nd_g = 875 ! 875 
  INTEGER, PARAMETER :: nd_h = 875 ! 875 
  REAL, PRIVATE, ALLOCATABLE :: rsa(:,:)
  REAL, PRIVATE, ALLOCATABLE :: ssa(:,:,:)
  REAL, PRIVATE, ALLOCATABLE :: hsa(:,:,:,:)
  REAL, PRIVATE, ALLOCATABLE :: gsa(:,:,:,:)
  REAL, ALLOCATABLE :: dsr(:), dss(:), dsh(:), dsg(:)

  COMPLEX :: far_b(nd_r), fbr_b(nd_r), far_f(nd_r), fbr_f(nd_r)
  COMPLEX :: fas_b(nd_s,nfw), fbs_b(nd_s,nfw), fas_f(nd_s,nfw), fbs_f(nd_s,nfw)
  COMPLEX :: fah_b(nd_h,nfw,nvarden), fbh_b(nd_h,nfw,nvarden), fah_f(nd_h,nfw,nvarden), fbh_f(nd_h,nfw,nvarden)
  COMPLEX :: fag_b(nd_g,nfw,nvarden), fbg_b(nd_g,nfw,nvarden), fag_f(nd_g,nfw,nvarden), fbg_f(nd_g,nfw,nvarden)

  CONTAINS

  SUBROUTINE read_table (rsafndir,vardendir,ndr,nds,ndg,ndh,bin_opt)
!-----------------------------------------------------------------------
! Read radar scattering amplitudes tables calculated using T-matrix
! method.
!-----------------------------------------------------------------------

    INTEGER :: ndr,nds,ndg,ndh ! separate bin numbers for different categories
    INTEGER :: bin_opt ! = 0: original bulk, 1: Takahashi bins, 2: PARSIVEL bins (not yet implemented)
    INTEGER :: istatus, i, j, k, i_rhog, i_rhoh
    INTEGER, PARAMETER :: nfw = 21, varden = 19
    CHARACTER (LEN=256) :: rsafndir,vardendir
    CHARACTER (LEN=256) :: rsafn
    CHARACTER (LEN=3), DIMENSION(nfw) :: extn = (/'000','005','010',    &
           '015','020','025','030','035','040','045','050','055','060', &
           '065','070','075','080','085','090','095','100'/)
    CHARACTER (LEN=3), DIMENSION(varden) :: denextn = (/'050','100','150', &
           '200','250','300','350','400','450','500','550','600','650','700',&
           '750','800','850','900','950'/)
    CHARACTER (LEN=256) :: head


!   IF(bin_opt == 0) THEN ! Original bins for bulk schemes
!     ndr = nd
!     nds = nd
!     ndg = nd
!     ndh = nd
!   ENDIF

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------
    ALLOCATE(dsr      (ndr),stat=istatus)
    !CALL check_alloc_status(istatus, "dsr")
    ALLOCATE(dss      (nds),stat=istatus)
    !CALL check_alloc_status(istatus, "dss")
    ALLOCATE(dsh      (ndh),stat=istatus)
    !CALL check_alloc_status(istatus, "dsh")
    ALLOCATE(dsg      (ndg),stat=istatus)
    !CALL check_alloc_status(istatus, "dsg")
    ALLOCATE(rsa   (ns,ndr),stat=istatus)
    !CALL check_alloc_status(istatus, "rsa")
    ALLOCATE(ssa(ns,nds,nfw),stat=istatus)
    !CALL check_alloc_status(istatus, "ssa")
    ALLOCATE(hsa(ns,ndh,nfw,varden),stat=istatus)
    !CALL check_alloc_status(istatus, "hsa")
    ALLOCATE(gsa(ns,ndg,nfw,varden),stat=istatus)
    !CALL check_alloc_status(istatus, "gsa")

!-----------------------------------------------------------------------
!  Read rain
!-----------------------------------------------------------------------
!   rsafn = TRIM(rsafndir)//'/SCTT_RAIN_fw100.dat'
    write (rsafn, '(A, A)') trim(rsafndir), '/SCTT_RAIN_fw100.dat'
    OPEN(UNIT=51,FILE=TRIM(rsafn),STATUS='old',FORM='formatted')
    READ(51,*) head

    DO j=1,ndr
      IF(bin_opt == 0) THEN
        READ(51,'(f5.2,8e13.5)') dsr(j), (rsa(i,j),i=1,8)
      ELSEIF(bin_opt == 1) THEN
        READ(51,'(f10.6,8e13.5)') dsr(j), (rsa(i,j),i=1,8)
      ENDIF
    ENDDO
    CLOSE(51)

    far_b = CMPLX(rsa(1,:),rsa(2,:))
    fbr_b = CMPLX(rsa(3,:),rsa(4,:))
    far_f = CMPLX(rsa(5,:),rsa(6,:))
    fbr_f = CMPLX(rsa(7,:),rsa(8,:))

!-----------------------------------------------------------------------
!  Read snow
!-----------------------------------------------------------------------
    IF(bin_opt == 0) THEN ! Only read snow table for original bulk schemes for now
                          ! Eventually will add Takahashi to this
    DO k=1,nfw
!     rsafn = TRIM(rsafndir)//'/SCTT_SNOW_fw'//extn(k)//'.dat'
      write (rsafn, '(A, A, A, A)') trim(rsafndir), '/SCTT_SNOW_fw', extn(k), '.dat'
      OPEN(UNIT=51,FILE=TRIM(rsafn),STATUS='old',FORM='formatted')
      READ(51,*) head

      DO j=1,nds
        IF(bin_opt == 0) THEN
          READ(51,'(f5.2,8e13.5)') dss(j), (ssa(i,j,k),i=1,8)
        ELSEIF(bin_opt == 1) THEN
          READ(51,'(f10.6,8e13.5)') dss(j), (ssa(i,j,k),i=1,8)
        ENDIF
      ENDDO
      CLOSE(51)
    ENDDO

    fas_b = CMPLX(ssa(1,:,:),ssa(2,:,:))
    fbs_b = CMPLX(ssa(3,:,:),ssa(4,:,:))
    fas_f = CMPLX(ssa(5,:,:),ssa(6,:,:))
    fbs_f = CMPLX(ssa(7,:,:),ssa(8,:,:))
    
    ENDIF

!-----------------------------------------------------------------------
!  Read hail
!-----------------------------------------------------------------------
    DO i_rhoh = 1, varden
      DO k=1, nfw
!       rsafn = TRIM(rsafndir)//'/SCTT_HAIL_fw'//extn(k)//'.dat'
        write (rsafn, '(A, A, A, A, A, A)') trim(vardendir), '/icedenrimed', denextn(i_rhoh), '/SCTT_GRPL_fw', extn(k), '.dat'
        OPEN(UNIT=51,FILE=TRIM(rsafn),STATUS='old',FORM='formatted')
        READ(51,*) head

        DO j=1,ndh
          IF(bin_opt == 0) THEN
            READ(51,'(f5.2,8e13.5)') dsh(j), (hsa(i,j,k,i_rhoh),i=1,8)
          ELSEIF(bin_opt == 1) THEN
            READ(51,'(f10.6,8e13.5)') dsh(j), (hsa(i,j,k,i_rhoh),i=1,8)
          ENDIF
        ENDDO
        CLOSE(51)
      ENDDO
    ENDDO

    fah_b = CMPLX(hsa(1,:,:,:),hsa(2,:,:,:))
    fbh_b = CMPLX(hsa(3,:,:,:),hsa(4,:,:,:))
    fah_f = CMPLX(hsa(5,:,:,:),hsa(6,:,:,:))
    fbh_f = CMPLX(hsa(7,:,:,:),hsa(8,:,:,:))

!-----------------------------------------------------------------------
!  Read graupel
!-----------------------------------------------------------------------
    DO i_rhog = 1, varden
      DO k=1, nfw
!       rsafn = TRIM(rsafndir)//'/SCTT_GRPL_fw'//extn(k)//'.dat'
        write (rsafn, '(A, A, A, A, A, A)') trim(vardendir), '/icedenrimed', denextn(i_rhog), '/SCTT_GRPL_fw', extn(k), '.dat'
        OPEN(UNIT=51,FILE=TRIM(rsafn),STATUS='old',FORM='formatted')
        READ(51,*) head

        DO j=1,ndg
          IF(bin_opt == 0) THEN
            READ(51,'(f5.2,8e13.5)') dsg(j), (gsa(i,j,k,i_rhog),i=1,8)
          ELSEIF(bin_opt == 1) THEN
            READ(51,'(f10.6,8e13.5)') dsg(j), (gsa(i,j,k,i_rhog),i=1,8)
          ENDIF
        ENDDO
        CLOSE(51)
      ENDDO
    ENDDO

    fag_b = CMPLX(gsa(1,:,:,:),gsa(2,:,:,:))
    fbg_b = CMPLX(gsa(3,:,:,:),gsa(4,:,:,:))
    fag_f = CMPLX(gsa(5,:,:,:),gsa(6,:,:,:))
    fbg_f = CMPLX(gsa(7,:,:,:),gsa(8,:,:,:))

    deallocate(rsa,ssa,hsa,gsa)

  END SUBROUTINE read_table

END MODULE rsa_table



 Module MPI_P     !Parallel                                                               ! MPI parallel calculation
      include 'mpif.h'
      integer, parameter :: sumtype = 17
      integer, parameter :: sizetype = 10
      integer, parameter :: masternode = 0
      integer            :: Process             ! 0 .. Numnodes-1
      integer            :: NumNodes
      integer            :: MPIerror
      integer            :: MPIstatus(3)
      real(8)            :: Pt0,Ptf
 contains


   subroutine P_start                                                          ! initialization of parallel calculation
        call MPI_INIT( MPIerror )                                        
          if(MPIerror/=0) print *,' P_start: error in MPI_INIT'
        call MPI_COMM_RANK( MPI_COMM_WORLD, Process, MPIerror )          
          if(MPIerror/=0) print *,' P_start: error in MPI_COMM_RANK'
        call MPI_COMM_SIZE( MPI_COMM_WORLD, Numnodes, MPIerror )         
          if(MPIerror/=0) print *,' P_start: error in MPI_COMM_SIZE'
        call cpu_time(Pt0)
        call MSG('Start parallel calculation on '//Pfstr(Numnodes)//' processors    ')
   end subroutine P_start


   subroutine P_sendI(Int)                                                     ! send data to all parallel processors
        integer     :: Int
        call MPI_BCAST(Int,1,MPI_INTEGER,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendI: error in BCAST'
   end subroutine P_sendI


   subroutine P_sendL(Int)                                                     ! send data to all parallel processors
        logical     :: Int
        call MPI_BCAST(Int,1,MPI_LOGICAL,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendL: error in BCAST'
   end subroutine P_sendL


   subroutine P_sendI1(Int,N)                                                     ! send data to all parallel processors
        integer     :: Int(N)
        integer     :: N
        call MPI_BCAST(Int,N,MPI_INTEGER,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendI: error in BCAST'
   end subroutine P_sendI1


   subroutine P_sendI2(Int,N1,N2)                                                     ! send data to all parallel processors
        integer     :: Int(N1,N2)
        integer     :: N1,N2
        call MPI_BCAST(Int(1,1),N1*N2,MPI_INTEGER,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendI: error in BCAST'
   end subroutine P_sendI2


   subroutine P_sendR(R)
        real*8      :: R
        call MPI_BCAST(R,1,MPI_DOUBLE_PRECISION,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendR: error in BCAST'
   end subroutine P_sendR


   subroutine P_sendR1(R,N)
        integer     :: N
        real*8      :: R(N)
        call MPI_BCAST(R,N,MPI_DOUBLE_PRECISION,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendR: error in BCAST'
   end subroutine P_sendR1


   subroutine P_sendR2(R,N1,N2)
        integer     :: N1,N2
        real*8      :: R(N1,N2)
        call MPI_BCAST(R(1,1),N1*N2,MPI_DOUBLE_PRECISION,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendR: error in BCAST'
   end subroutine P_sendR2


   subroutine P_sendR3(R,N1,N2,N3)
        integer     :: N1,N2,N3
        real*8      :: R(N1,N2,N3)
        call MPI_BCAST(R(1,1,1),N1*N2*N3,MPI_DOUBLE_PRECISION,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendR: error in BCAST'
   end subroutine P_sendR3


   subroutine P_sendA(A)
        character(*)     :: A
        integer          :: l
        l = len(A)
        call MPI_BCAST(A,l,MPI_CHARACTER,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendA: error in BCAST'
   end subroutine P_sendA


   subroutine P_receiveR(R,Nn,ist,ngr)                                         ! collect array R from parallel processors
        integer    :: Nn,ist,ngr
        real*8     :: R(Nn)
        call MPI_Barrier(MPI_COMM_WORLD, MPIerror)                             ! wait all processors
          if(MPIerror/=0) print *,' P_receiveR: error in MPI_BARRIER'
        call MPI_GATHER(R(ist),ngr,MPI_DOUBLE_PRECISION,R,ngr,MPI_DOUBLE_PRECISION,masternode,MPI_COMM_WORLD,MPIerror)   ! Masternode receives data from others cores 
          if(MPIerror/=0) print *,' P_receiveR: error in MPI_GATHER'
   end subroutine P_receiveR


   subroutine P_receiveZ(R,Nn,ist,ngr)                                         ! collect array R from parallel processors
        integer    :: Nn,ist,ngr
        complex(8) :: R(Nn)
        call MPI_Barrier(MPI_COMM_WORLD, MPIerror)                             ! wait all processors
          if(MPIerror/=0) print *,' P_receiveR: error in MPI_BARRIER'
        call MPI_GATHER(R(ist),ngr,MPI_DOUBLE_COMPLEX,R,ngr,MPI_DOUBLE_COMPLEX,masternode,MPI_COMM_WORLD,MPIerror)   ! Masternode receives data from others cores 
          if(MPIerror/=0) print *,' P_receiveR: error in MPI_GATHER'
   end subroutine P_receiveZ


   subroutine P_wait                                         ! wait all processes
        call MPI_Barrier(MPI_COMM_WORLD, MPIerror)                             ! wait all processors
          if(MPIerror/=0) print *,' P_wait: error in MPI_BARRIER'
   end subroutine P_wait


   subroutine P_calc_group(Nn,ist,ngr)              ! calculate groups for calculation
    integer :: ist,ngr,Nn                         ! parallel calculation from ist+1 to ist+ngr
    integer :: ifin
    ngr = Nn/NumNodes                             ! number of group of points for calculation
    if(mod(Nn,NumNodes) /= 0) ngr = ngr + 1
    ist = Process*ngr                             ! starting point
    ifin = ist + ngr                              ! finil point
    if(ifin > Nn) ngr = Nn - ist                  ! correct number of group for last process
   end subroutine P_calc_group


   subroutine MSG(comment)
    character(*)          :: comment
 !   integer               :: l1
 !   l1 = len(trim(adjustl(comment)))
    if(Process==0) print 1,trim(adjustl(comment))
 1  format(A)
   end subroutine MSG


   subroutine MSGT(comment)
    character(*)          :: comment
    integer               :: l1     !,l2
    character(len=60)     :: cp
    call cpu_time(Ptf)
    cp = ' '
    l1 = len(trim(adjustl(comment)))
    if(l1>=60) then
     if(Process==0) print 1,trim(adjustl(comment)),(Ptf-Pt0)
    else
  !   l2 = 60 - l1
     cp(1:l1) = trim(adjustl(comment))
     if(Process==0) print 2,cp,(Ptf-Pt0)
    endif
    Pt0 = Ptf
 1  format(A,' (',F10.3,' s)')
 2  format(A60,' (',F10.3,' s)')
   end subroutine MSGT


   subroutine P_stop                                                ! finalizing parallel calculation
        call MPI_FINALIZE(MPIerror)
          if(MPIerror/=0) print *,' P_stop: error in MPI_FINALIZE'
   end subroutine P_stop


     character(len=5) function Pfstr(k)                             ! Convert an integer to character(5)
      integer, intent(in) :: k
      write(Pfstr,'(I5)') k
     end function Pfstr


 end module MPI_P



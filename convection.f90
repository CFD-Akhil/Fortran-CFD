
!     Convection example problem solver with zero gradient exit conditions
!     compile using g95 compiler with:
!     g95 -O2 -o convection.exe convection.f90
!     -O2 sets an optimization switch, and the -o sets the output file name

      implicit none
!     nx and ny are number of interior finite volumes in each direction
      integer, parameter :: nx=50, ny=50
      real,dimension(0:nx+1) :: x
      real,dimension(0:ny+1) :: y
      real,dimension(0:nx+1,0:ny+1)::phi,AP,APT,AE,AET,AW,AWT,AN,ANT,AS,AST, &
                                     SU,SUT,SP,SPT
      real :: z,dx,dy,beta,inletvelocity,rho,phiwest,phisouth, &
              mdote,mdotw,mdotn,mdots
      integer :: i,j,iter
!**************************************************************************
!     set cell size
      dx=1./nx;   dy=1./ny

!     set blending factor
      write(*,*)'enter blending factor where 0<=beta<=1'
      read(*,*)beta

!     set u and v inlet velocity
      inletvelocity=2.0  !each component
!     set density
      rho=1.
!     set phi values for inlet boundary conditions
      phiwest=1.0; phisouth=0.0

!     compute mdot  (same for all faces and all cells)
      mdote=rho*inletvelocity*dy;    mdotw=rho*inletvelocity*dy
      mdotn=rho*inletvelocity*dx;    mdots=rho*inletvelocity*dx
!**************************************************************************
!     initialize solution to zero
      phi=0

!     initialize su,sut,sp,spt
      SU=0; SP=0; SUT=0; SPT=0

!     set AP,AE, etc. over interior cells
!          note AE, AW, etc. are Upwinding Coefficients
!          AET, AWT, etc. are Central Differencing Coefficients
      AE=0; AW=mdotw; AN=0;  AS=mdots
      AET=-mdote/2.;  AWT=mdotw/2.; ANT=-mdotn/2.; AST=mdots/2.

!     set (overwrite) boundary conditions where necessary

!     west side
!     cut links
      AW(1,:)=0.0;   AWT(1,:)=0.0
!     set boundary values
      SU(1,:)=mdotw*phiwest;   SUT(1,:)=mdotw*phiwest
      SP(1,:)=-mdote;   SPT(1,:)=-mdote

!     east side
      AE(nx,:)=0.0;   AET(nx,:)=0.0
      SP(nx,:)=-(mdote-mdotw);  SPT(nx,:)=-(mdote-mdotw)

!     south side
      AS(:,1)=0.0;   AST(:,1)=0.0
      SU(:,1)=mdots*phisouth;   SUT(:,1)=mdots*phisouth
      SP(:,1)=-mdotn;  SPT(:,1)=-mdotn

!     north side
      AN(:,ny)=0.0;   ANT(:,ny)=0.0
      SP(:,ny)=-(mdotn-mdots);   SPT(:,ny)=-(mdotn-mdots)

!     southwest corner
      AS(1,1)=0.0;   AW(1,1)=0.0
      AST(1,1)=0.0;  AWT(1,1)=0.0
!     next four lines determine west wall bc on southwest cell (i.e., phiwest=0 or 1)
!     to set phiwest to 0 on west side of southwest cell:
         SU(1,1)=mdotw*0 + mdots*phisouth
         SUT(1,1)=mdotw*0 + mdots*phisouth
!     OR to let let phiwest equal 1 on the west side of southwest cell, uncomment:
!        SU(1,1)=mdotw*phiwest + mdots*phisouth
!        SUT(1,1)=mdotw*phiwest + mdots*phisouth
      SP(1,1)=-mdote-mdotn
      SPT(1,1)=-mdote-mdotn

!     northwest corner
      AN(1,ny)=0.0;   AW(1,ny)=0.0
      ANT(1,ny)=0.0;  AWT(1,ny)=0.0
      SU(1,ny)=mdotw*phiwest;   SUT(1,ny)=mdotw*phiwest
      SP(1,ny)=-mdote-(mdotn-mdots);   SPT(1,ny)=-mdote-(mdotn-mdots) 

!     southeast corner
      AS(nx,1)=0.0;   AE(nx,1)=0.0
      AST(nx,1)=0.0;  AET(nx,1)=0.0
      SP(nx,1)=-(mdote-mdotw)-mdotn;  SPT(nx,1)=-(mdote-mdotw)-mdotn

!     northeast corner
      AN(nx,ny)=0.0;   AE(nx,ny)=0.0
      ANT(nx,ny)=0.0;  AET(nx,ny)=0.0
      SP(nx,ny)=-(mdote-mdotw)-(mdotn-mdots) 
      SPT(nx,ny)=-(mdote-mdotw)-(mdotn-mdots)

!     compute AP and APT
      AP=AE+AW+AN+AS-SP
      APT=AET+AWT+ANT+AST-SPT

!     iterative solution using deferred correction scheme
      do iter=1,100000   !ensure enough iterations for convergence (needs to be a large number)
      do i=1,nx
      do j=1,ny
      phi(i,j)=(AE(i,j)*phi(i+1,j)+AW(i,j)*phi(i-1,j)+AN(i,j)*phi(i,j+1)+AS(i,j)*phi(i,j-1)+SU(i,j))/AP(i,j)- &
               beta/AP(i,j)*( (APT(i,j)*phi(i,j)-AET(i,j)*phi(i+1,j)-AWT(i,j)*phi(i-1,j)-ANT(i,j)*phi(i,j+1)- &
               AST(i,j)*phi(i,j-1)-SUT(i,j)) - (AP(i,j)*phi(i,j)-AE(i,j)*phi(i+1,j)-AW(i,j)*phi(i-1,j) - &
               AN(i,j)*phi(i,j+1)-AS(i,j)*phi(i,j-1) - SU(i,j)) )
50    continue
      end do; end do; end do

!*********************************************************************************************

!     WRITE comma separated variable DATA FILE FOR PARAVIEW
!     only writing interior cell values, interior points only, no boundary values written
      open(unit=10,file='results.csv')

!     compute x and y values (at cell center locations)
      y(1)=0.5*dy;  x(1)=0.5*dx
      do j=2,ny
        y(j)=y(j-1)+dy
      end do
      do i=2,nx
        x(i)=x(i-1)+dx
      end do

!     interior points only
      write(10,*)'x coord, y coord, z coord, phi'
      do i=1,nx
      do j=1,ny
      z=0.
        write(10,*)x(i),',',y(j),',',z,',',phi(i,j)
      end do; end do
      
      end

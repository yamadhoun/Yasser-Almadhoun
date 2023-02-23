#include "fintrf.h"
!
!   intersections in Fortran
!
! Copyright notice
!  --------------------------------------------------------------------
!  Copyright (C) 2022 IHE Delft & Deltares
!
!      Johan Reyns
!      j.reyns@un-ihe.org
!      Westvest 7
!      2611AX Delft
!
!  This library is free software: you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation, either
!  version 2.1 of the License, or (at your option) any later version.
!
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library. If not, see <http://www.gnu.org/licenses
!  --------------------------------------------------------------------     
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)

      implicit none 

      integer*8 :: plhs(*), prhs(*)
                
      integer*8 :: nlhs, nrhs
                
      integer*8 :: mxcreatedoublematrix, mxgetpr, mxCreateNumericMatrix
      integer*8 :: mxClassIDFromClassName
      integer*8 :: mxgetm, mxgetn
                
      integer*8 :: np1, np2, np0, m, n
                
      integer*8 :: xp_pr,yp_pr,x1_pr,y1_pr,x2_pr,y2_pr
      integer*8 :: indi_pr, indj_pr

      
      double precision, dimension(:),     allocatable    :: x1  ! poly1, row vector
      double precision, dimension(:),     allocatable    :: y1
      double precision, dimension(:),     allocatable    :: x2  ! poly2, row vector
      double precision, dimension(:),     allocatable    :: y2
      double precision, dimension(:),     pointer        :: xp  ! intersection, row vector
      double precision, dimension(:),     pointer        :: yp  ! intersection, row vector
      double precision, dimension(:),     pointer        :: indi
      double precision, dimension(:),     pointer        :: indj

      interface
         subroutine intersect(xi,yi,xj,yj,xp,yp,indi,indj)
            real*8, dimension(:),           intent(in)      :: xi,yi  ! coordinate of first polyline
            real*8, dimension(:),           intent(in)      :: xj,yj  ! coordinate of second polyline
            double precision, dimension(:), pointer         :: xp  ! coordinates of intersection points
            double precision, dimension(:), pointer         :: yp  ! coordinates of intersection points
            double precision, dimension(:), pointer         :: indi ! sum of segment number +fraction
            double precision, dimension(:), pointer         :: indj 
         end subroutine

         subroutine selfintersect(xi,yi,xp,yp,indi,indj)
            real*8, dimension(:),           intent(in)      :: xi,yi  ! coordinate of first polyline
            double precision, dimension(:), pointer         :: xp     ! coordinates of intersection points
            double precision, dimension(:), pointer         :: yp     ! coordinates of intersection points
            double precision, dimension(:), pointer         :: indi   ! sum of segment number +fraction
            double precision, dimension(:), pointer         :: indj 
         end subroutine
      end interface

!     handle poly1
      m = mxgetm(prhs(1))
      n = mxgetn(prhs(1))
      np1=max(m,n)

!     allocate
      allocate(x1(1:np1))
      allocate(y1(1:np1))

!     get pointers
      x1_pr = mxgetpr(prhs(1))  
      y1_pr = mxgetpr(prhs(2))

!     copy to fortran arrays
      call mxcopyptrtoreal8(x1_pr,x1,np1)
      call mxcopyptrtoreal8(y1_pr,y1,np1)

      if (nrhs==4) then
          
!         length poly2
          m = mxgetm(prhs(3))
          n = mxgetn(prhs(3))
          np2=max(m,n)
          
!         allocate     
          allocate(x2(1:np2))
          allocate(y2(1:np2))
          
!         Get pointers to the input arrays
          x2_pr = mxgetpr(prhs(3))  
          y2_pr = mxgetpr(prhs(4))
           
!         load the data into fortran arrays.
          call mxcopyptrtoreal8(x2_pr,x2,np2)
          call mxcopyptrtoreal8(y2_pr,y2,np2)
          
!         call the computational subroutine
          call intersect(x1,y1,x2,y2,xp,yp,indi,indj)
            
      elseif (nrhs==2) then
          ! self intersection
          call selfintersect(x1,y1,xp,yp,indi,indj)

      else
          call mexErrMsgIdAndTxt("mexinterx:ErrorIdIn",
     .           "Invalid number of inputs to mexinterx function.");
      endif

!     load the output into a matlab array.
      np0=size(xp)
      plhs(1) = mxcreatedoublematrix(1,np0,0)
      plhs(2) = mxcreatedoublematrix(1,np0,0)
      plhs(3) = mxcreatedoublematrix(1,np0,0)
      plhs(4) = mxcreatedoublematrix(1,np0,0)

      xp_pr    = mxgetpr(plhs(1))  
      yp_pr    = mxgetpr(plhs(2))  
      indi_pr  = mxgetpr(plhs(3))
      indj_pr  = mxgetpr(plhs(4))

      call mxcopyreal8toptr(xp,xp_pr,np0)
      call mxcopyreal8toptr(yp,yp_pr,np0)
      call mxcopyreal8toptr(indi,indi_pr,np0)
      call mxcopyreal8toptr(indj,indj_pr,np0)
    
      deallocate(x1)
      deallocate(y1)
      if (nrhs==4) then
         deallocate(x2)
         deallocate(y2)
      endif
      deallocate(indi)
      deallocate(indj)
      deallocate(xp)
      deallocate(yp)

      end subroutine mexFunction

        ! program testintersect
        ! 
        !   interface
        !      subroutine intersect(xi,yi,xj,yj,xp,yp,indxp,indyp)
        !         real*8, dimension(:),              intent(in)   :: xi,yi  ! coordinate of first polyline
        !         real*8, dimension(:),              intent(in)   :: xj,yj  ! coordinate of second polyline
        !         double precision, dimension(:), pointer         :: xp  ! coordinates of intersection points
        !         double precision, dimension(:), pointer         :: yp  ! coordinates of intersection points
        !         real*8, dimension(:), pointer           :: indxp,indyp ! indices of intersection points such that
        !       end subroutine
        !      
        !       subroutine selfintersect(xi,yi,xp,yp,indi,indj)
        !          real*8, dimension(:),           intent(in)      :: xi,yi  ! coordinate of first polyline
        !          double precision, dimension(:), pointer         :: xp     ! coordinates of intersection points
        !          double precision, dimension(:), pointer         :: yp     ! coordinates of intersection points
        !          double precision, dimension(:), pointer         :: indi   ! sum of segment number +fraction
        !          double precision, dimension(:), pointer         :: indj 
        !       end subroutine
        !   end interface
        ! 
        !    double precision, allocatable, dimension(:) :: x1,y1,x2,y2
        !    double precision, dimension(:), pointer :: xp,yp
        !    double precision, dimension(:), pointer :: indxp, indyp 
        ! 
        !    allocate(x1(6),y1(6),x2(3),y2(3))
        ! 
        !    x1=[real*8::0,1,2, 0.5, 1.5, 2.5]
        !    y1=[real*8::0,1,2, 0,1,2]
        !    x2=[2,1,0]
        !    y2=[0,1,2]
        ! 
        !    call intersect(x1,y1,x2,y2,xp,yp,indxp,indyp)
        !    write(*,*) 'Intersect 2 pols'
        !    write(*,*) xp
        !    write(*,*) yp
        !    write(*,*) indxp
        !    write(*,*) indyp
        ! 
        !    write(*,*) 'Selfintersect'
        !    deallocate(x1,y1)
        !    allocate(x1(11),y1(11))
        !    x1=[real*8:: 500,
        !.                490.043823074421,
        !.                461.500101955640,
        !.                415.434027279884,
        !.                353.553390593274,
        !.                277.592593469123,
        !.                191.159601266592,
        !.                97.4686125760648,
        !.                -2.30058506182734e-13,
        !.                -96.4867241639488,
        !.                -195.866845175306,
        !.                -284.400791182354,
        !.                -353.194070843977,
        !.                -412.377476426763,
        !.                -455.305010070887,
        !.                -481.474602664877,
        !.                -490.050568644990,
        !.                -481.474602664877,
        !.                -455.305010070887,
        !.                -412.377476426763,
        !.                -353.194070843976,
        !.                -284.400791182354,
        !.                -195.866845175305,
        !.                -96.4867241639485,
        !.                4.35186014496447e-15,
        !.                97.4686125760658,
        !.                191.159601266592,
        !.                277.592593469124,
        !.                353.553390593274,
        !.                415.434027279884,
        !.                461.500101955640,
        !.                490.043823074421,
        !.                500]
        ! 
        !    y1=[real*8 :: 0,
        !.                 -97.4686125760651,
        !.                 -191.159601266592,
        !.                 -277.592593469124,
        !.                 -353.553390593274,
        !.                 -415.434027279884,
        !.                 -461.500101955640,
        !.                 -490.043823074421,
        !.                 -504.974715677505,
        !.                 -485.168553089752,
        !.                 -472.864394041454,
        !.                 -425.676201508103,
        !.                 -353.194070843976,
        !.                 -275.560837118870,
        !.                 -188.593510187780,
        !.                 -95.7427557009659,
        !.                 2.18453545796162e-13,
        !.                 95.7427557009664,
        !.                 188.593510187780,
        !.                 275.560837118870,
        !.                 353.194070843976,
        !.                 425.676201508103,
        !.                 472.864394041454,
        !.                 485.168553089752,
        !.                 504.974715677505,
        !.                 490.043823074421,
        !.                 461.500101955640,
        !.                 415.434027279883,
        !.                 353.553390593274,
        !.                 277.592593469123,
        !.                 191.159601266592,
        !.                 97.4686125760649,
        !.                 0]
        !    call selfintersect(x1,y1,xp,yp,indxp,indyp)
        !    write(*,*) xp
        !    write(*,*) yp
        !    write(*,*) indxp
        !    write(*,*) indyp
        !    pause(10)
        ! 
        ! end program testintersect

      subroutine intersect(xi,yi,xj,yj,xp,yp,indxp,indyp)

         implicit none 

         real*8, dimension(:),              intent(in)   :: xi,yi  ! coordinate of first polyline
         real*8, dimension(:),              intent(in)   :: xj,yj  ! coordinate of second polyline
         real*8, dimension(:), pointer         :: xp  ! coordinates of intersection points
         real*8, dimension(:), pointer         :: yp  ! coordinates of intersection points
         real*8, dimension(:), pointer         :: indxp,indyp ! indices of intersection points such that
                                                                      ! they are between indi and indi+1
         logical, dimension(:), allocatable     :: remove, in_range
         integer*4                              :: ni, nj, i, j, ip, cnt
         integer*4                              :: np, ii, jj, k, n
         integer, allocatable, dimension(:)   :: indi, indj
         integer, allocatable, dimension(:,:)   :: pos
         real*8, allocatable, dimension(:)    :: indxpt, indypt
         real*8, allocatable, dimension(:,:)  :: dxmini, dxmaxi
         real*8, allocatable, dimension(:,:)  :: dymini, dymaxi
         real*8, allocatable, dimension(:,:)  :: dxminj, dxmaxj
         real*8, allocatable, dimension(:,:)  :: dyminj, dymaxj
         real*8, allocatable, dimension(:,:)  :: temp
         logical, allocatable, dimension(:,:) :: ij
         real*8, allocatable, dimension(:,:)  :: xy1, xy2, T, B
         real*8, allocatable, dimension(:,:)  :: dxy1, dxy2
         real*8, allocatable, dimension(:)    :: xpt,ypt
         real*8, allocatable, dimension(:,:,:):: AA


         interface
            subroutine find_2d_true(ij,pos1,pos2)              
               logical, dimension(:,:), intent(in)  :: ij
               integer, dimension(:), intent(out)   :: pos1
               integer, dimension(:), intent(out)   :: pos2
            end subroutine
 
         end interface

         ni=size(xi)-1
         nj=size(xj)-1
         np=max(ni+1,nj+1)

         allocate(xy1(ni+1,2),xy2(nj+1,2))
         allocate(dxy1(ni,2),dxy2(nj,2))

         xy1(:,1)=xi; xy1(:,2)=yi
         xy2(:,1)=xj; xy2(:,2)=yj

         dxy1 = xy1(2:ni+1,:)-xy1(1:ni,:)
         dxy2 = xy2(2:nj+1,:)-xy2(1:nj,:)

         indi=0; indj=0
         cnt=0
         allocate(dxmini(ni,nj), dxmaxi(ni,nj),dymini(ni,nj), 
     .            dymaxi(ni,nj))
         allocate(dxminj(ni,nj), dxmaxj(ni,nj),dyminj(ni,nj), 
     .            dymaxj(ni,nj))
         allocate(ij(ni,nj))

         dxmini=spread(min(xi(1:ni),xi(2:ni+1)),2,nj)
         dxmaxi=spread(max(xi(1:ni),xi(2:ni+1)),2,nj)
         dymini=spread(min(yi(1:ni),yi(2:ni+1)),2,nj)
         dymaxi=spread(max(yi(1:ni),yi(2:ni+1)),2,nj)

         dxminj=spread(min(xj(1:nj),xj(2:nj+1)),1,ni)
         dxmaxj=spread(max(xj(1:nj),xj(2:nj+1)),1,ni)
         dyminj=spread(min(yj(1:nj),yj(2:nj+1)),1,ni)
         dymaxj=spread(max(yj(1:nj),yj(2:nj+1)),1,ni)

         ij=(dxmini<=dxmaxj) .and. (dxmaxi>=dxminj) .and.
     .      (dymini<=dymaxj) .and. (dymaxi>=dyminj)

         np = count( ij == .true. )
         allocate(indi(np),indj(np))
         call find_2d_true(ij,indi,indj)

         allocate(remove(cnt))
         remove = .not.(isnan(sum(dxy1(indi,:) + dxy2(indj,:),2)))
         indi=pack(indi,remove)
         indj=pack(indj,remove)

         n = size(indi)
         allocate(T(4,n))
         T = 0.0
         allocate(AA(4,4,n))
         AA = 0.0
         AA([1, 2],3,:) = -1
         AA([3, 4],4,:) = -1
         AA([1, 3],1,:) = transpose(dxy1(indi,:))
         AA([2, 4],2,:) = transpose(dxy2(indj,:))
         allocate(B(4,n))

         B(1,:) = -xi(indi)
         B(2,:) = -xj(indj)
         B(3,:) = -yi(indi)
         B(4,:) = -yj(indj)
      
         allocate(temp(4,4))
         do k=1,n
            call matinv4(AA(:,:,k),temp)
            T(:,k) = matmul(temp,B(:,k))
         enddo
	
         allocate(in_range(n))
         in_range=.false.
	   in_range = T(1,:) >= 0 .and. T(2,:) >= 0 .and. T(1,:) < 1 .and. 
     .              T(2,:) < 1
	   xpt = pack(T(3,:),in_range)
	   ypt = pack(T(4,:),in_range)
	
		indxpt = pack(indi+T(1,:),in_range)
		indypt = pack(indj+T(2,:),in_range)

         np=size(xpt)
         allocate(xp(np),yp(np),indxp(np),indyp(np))
         xp=xpt
         yp=ypt
         indxp=indxpt
         indyp=indypt

         deallocate(indi,indj)
         deallocate(xy1,xy2)
         deallocate(dxy1,dxy2)
         deallocate(remove)
         deallocate(T)
         deallocate(AA)
         deallocate(B)
         deallocate(temp)
         deallocate(in_range)
  
      end subroutine intersect



      subroutine selfintersect(xi,yi,xp,yp,indxp,indyp)

         implicit none 

         real*8, dimension(:),    intent(in)   :: xi,yi  ! coordinate of first polyline
         real*8, dimension(:), pointer         :: xp  ! coordinates of intersection points
         real*8, dimension(:), pointer         :: yp  ! coordinates of intersection points
         real*8, dimension(:), pointer         :: indxp,indyp ! indices of intersection points such that
                                                                      ! they are between indi and indi+1
         logical, dimension(:), allocatable   :: remove, in_range
         integer*4                            :: ni, nj, i, j, ip, cnt
         integer*4                            :: np, ii, jj, k, n
         integer*4                            :: shp(2)
         integer, allocatable, dimension(:)   :: indi, indj
         integer, allocatable, dimension(:,:) :: pos
         real*8                               :: eps
         real*8, allocatable, dimension(:)    :: indxpt, indypt
         real*8, allocatable, dimension(:,:)  :: dxmini, dxmaxi
         real*8, allocatable, dimension(:,:)  :: dymini, dymaxi
         real*8, allocatable, dimension(:,:)  :: dxminj, dxmaxj
         real*8, allocatable, dimension(:,:)  :: dyminj, dymaxj
         real*8, allocatable, dimension(:,:)  :: temp
         logical, allocatable, dimension(:,:) :: ij
         real*8, allocatable, dimension(:,:)  :: xy1, xy2, T, B
         real*8, allocatable, dimension(:,:)  :: dxy1, dxy2
         real*8, allocatable, dimension(:)    :: xt,yt
         real*8, allocatable, dimension(:)    :: xj, yj
         real*8, allocatable, dimension(:)    :: xpt,ypt
         real*8, allocatable, dimension(:,:,:):: AA


         interface
            subroutine find_2d_true(ij,pos1,pos2)              
               logical, dimension(:,:), intent(in)  :: ij
               integer, dimension(:), intent(out)   :: pos1
               integer, dimension(:), intent(out)   :: pos2
            end subroutine
 
         end interface

         eps=1d-4
         np=size(xi)
         if (abs(xi(1)-xi(np))<eps .and. abs(yi(1)-yi(np))<eps) then
            allocate(temp(np,2))
            temp(:,1)=xi
            temp(:,2)=yi
            np=np-1
            allocate(xt(np),yt(np))
            xt=temp(1:np,1)
            yt=temp(1:np,2)
         else
            allocate(xt(np),yt(np))
            xt=xi
            yt=yi
         endif
         
         allocate(xj(np),yj(np))
         xj=xt
         yj=yt


         ni=size(xt)-1
         nj=size(xj)-1
         np=max(ni+1,nj+1)

         allocate(xy1(ni+1,2),xy2(nj+1,2))
         allocate(dxy1(ni,2),dxy2(nj,2))

         xy1(:,1)=xt; xy1(:,2)=yt
         xy2(:,1)=xj; xy2(:,2)=yj

         dxy1 = xy1(2:ni+1,:)-xy1(1:ni,:)
         dxy2 = xy2(2:nj+1,:)-xy2(1:nj,:)

         indi=0; indj=0
         cnt=0
         allocate(dxmini(ni,nj), dxmaxi(ni,nj),dymini(ni,nj), 
     .            dymaxi(ni,nj))
         allocate(dxminj(ni,nj), dxmaxj(ni,nj),dyminj(ni,nj), 
     .            dymaxj(ni,nj))
         allocate(ij(ni,nj))

         dxmini=spread(min(xt(1:ni),xt(2:ni+1)),2,nj)
         dxmaxi=spread(max(xt(1:ni),xt(2:ni+1)),2,nj)
         dymini=spread(min(yt(1:ni),yt(2:ni+1)),2,nj)
         dymaxi=spread(max(yt(1:ni),yt(2:ni+1)),2,nj)

         dxminj=spread(min(xj(1:nj),xj(2:nj+1)),1,ni)
         dxmaxj=spread(max(xj(1:nj),xj(2:nj+1)),1,ni)
         dyminj=spread(min(yj(1:nj),yj(2:nj+1)),1,ni)
         dymaxj=spread(max(yj(1:nj),yj(2:nj+1)),1,ni)

         ij=(dxmini<=dxmaxj) .and. (dxmaxi>=dxminj) .and.
     .      (dymini<=dymaxj) .and. (dymaxi>=dyminj)

         np = count( ij == .true. )
         allocate(indi(np),indj(np))
         call find_2d_true(ij,indi,indj)

         allocate(remove(cnt))
         remove = .not.(isnan(sum(dxy1(indi,:) + dxy2(indj,:),2)) 
     .                  .or. ( indj <= indi + 1))
         indi=pack(indi,remove)
         indj=pack(indj,remove)

         n = size(indi)
         allocate(T(4,n))
         T = 0.0
         allocate(AA(4,4,n))
         AA = 0.0
         AA([1, 2],3,:) = -1
         AA([3, 4],4,:) = -1
         AA([1, 3],1,:) = transpose(dxy1(indi,:))
         AA([2, 4],2,:) = transpose(dxy2(indj,:))
         allocate(B(4,n))

         B(1,:) = -xt(indi)
         B(2,:) = -xj(indj)
         B(3,:) = -yt(indi)
         B(4,:) = -yj(indj)
      
         if (allocated(temp)) then
            deallocate(temp)
         endif
         allocate(temp(4,4))
         do k=1,n
            call matinv4(AA(:,:,k),temp)
            T(:,k) = matmul(temp,B(:,k))
         enddo
	
         allocate(in_range(n))
         in_range=.false.
	   in_range = T(1,:) >= 0 .and. T(2,:) >= 0 .and. T(1,:) < 1 .and. 
     .              T(2,:) < 1
	   xpt = pack(T(3,:),in_range)
	   ypt = pack(T(4,:),in_range)
	
		indxpt = pack(indi+T(1,:),in_range)
		indypt = pack(indj+T(2,:),in_range)

         np=size(xpt)
         allocate(xp(np),yp(np),indxp(np),indyp(np))
         xp=xpt
         yp=ypt
         indxp=indxpt
         indyp=indypt

         deallocate(indi,indj)
         deallocate(xy1,xy2)
         deallocate(dxy1,dxy2)
         deallocate(remove)
         deallocate(T)
         deallocate(AA)
         deallocate(B)
         deallocate(temp)
         deallocate(in_range)
         deallocate(xt,yt)
  
      end subroutine selfintersect






      subroutine matinv4(A,B)
        !! Performs a direct calculation of the inverse of a 4×4 matrix.
        double precision, intent(in) :: A(4,4)   !! Matrix
        double precision, intent(out):: B(4,4)   !! Inverse matrix
        double precision             :: detinv
      
        ! Calculate the inverse determinant of the matrix
        detinv = 
     .   1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)
     .   *A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
     .   - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)
     .   *A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
     .   + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)
     .   *A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
     .   - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)
     .   *A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))
      
        ! Calculate the inverse of the matrix
        B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)
     .  *(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)
     .  *(A(3,2)*A(4,3)-A(3,3)*A(4,2)))

        B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)
     .  *(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)
     .  *(A(3,3)*A(4,1)-A(3,1)*A(4,3)))

        B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)
     .  *(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)
     .  *A(4,1)))

        B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*
     .  (A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)
     .  *A(4,2)))

        B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*
     .  (A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*
     .  A(4,3)))

        B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*
     .  (A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*
     .  A(4,1)))

        B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)
     .  *(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)
     .  *A(4,2)))

        B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)
     .  *(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)
     .  *A(4,1)))

        B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)
     .  *(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*
     .  A(4,2)))

        B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)
     .  *(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)
     .  *A(4,3)))

        B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*
     . (A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)
     . *A(4,1)))

        B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*
     .  (A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)
     . *A(4,2)))

        B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)
     .  *(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)
     .  *A(3,3)))

        B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*
     .  (A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)
     .  *A(3,1)))

        B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)
     .  *(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)
     .  *A(3,2)))

        B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)
     .  *(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)
     .  *A(3,1)))
      end subroutine


      subroutine find_2d_true(ij,pos1,pos2)
         implicit none

         logical, dimension(:,:), intent(in)  :: ij
         integer, dimension(:), intent(out)   :: pos1
         integer, dimension(:), intent(out)   :: pos2

         integer, dimension(:,:), allocatable :: ai, aj
         integer                 :: i,j,shp(2)

         shp=shape(ij)
         allocate(ai(shp(1),shp(2)),aj(shp(1),shp(2)))
         do j = 1,size(ij,2)
            do i = 1,size(ij,1)
               ai(i,j) = i
               aj(i,j) = j
            enddo
         enddo

         pos1 = pack( ai, ij == .true. )
         pos2 = pack( aj, ij == .true.)

      end subroutine


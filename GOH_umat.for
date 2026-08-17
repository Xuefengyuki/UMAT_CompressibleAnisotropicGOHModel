C     ABAQUS Subroutine for anisotropic GOH model
C
C     More Infos can be found at https://github.com/Xuefengyuki/UMAT_CompressibleAnisotropicGOHModel.git


C---------------------------------------------------------------------------------------
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
C---------------------------------------------------------------------------------------
C End of Base code
C---------------------------------------------------------------------------------------

      
      
      
C--------------------------------------------------------------------------------------- 
C Start of USER code
C---------------------------------------------------------------------------------------
C     Material parameters
      real*8 C10,D
      real*8 K1(2),K2(2),KAPPA(2)
      real*8 NV(3,2),Mfiber(3,3,2)
C
C     Deformation tensors
      real*8 C(3,3),Cinv(3,3)
      real*8 B(3,3),Binv(3,3)
      real*8 eye(3,3)
C
C     Determinants and invariants
      real*8 J_F,detF,detC,detB
      real*8 I1,I2,I3
      real*8 I1bar,Jm23
      real*8 I4(2),I4bar(2)
      real*8 Ealpha(2)
C
C     First derivatives
      real*8 dJdC(3,3)
      real*8 dI1bardC(3,3)
      real*8 dI4bardC(3,3,2)
      real*8 dEalphadC(3,3,2)
C
C     Second derivatives
      real*8 dCinvdC(3,3,3,3)
      real*8 d2JdC2(3,3,3,3)
      real*8 d2I1bardC2(3,3,3,3)
      real*8 d2I4bardC2(3,3,3,3,2)
      real*8 d2EalphadC2(3,3,3,3,2)
C
C     Fiber activation and coefficients
      real*8 Halpha(2)
      real*8 fiberCoeff(2)
      real*8 fiberTangent1(2)
      real*8 fiberTangent2(2)
C
C     Stress tensors
      real*8 S(3,3)
      real*8 sigma(3,3)
C
C     Material and spatial tangent tensors
      real*8 C0ijkl(3,3,3,3)
      real*8 Cvolijkl(3,3,3,3)
      real*8 Cijkl(3,3,3,3)
      real*8 Cadd(3,3,3,3)
      real*8 Caba(3,3,3,3)
C
C     Angle of symmetric fibers family
      real*8 beta
C
C     Integer variables
      integer i,j,k,l,m,n,p,q,a
      integer ind1(6),ind2(6)
C
      data ind1/1,2,3,1,1,2/
      data ind2/1,2,3,2,3,3/
C
C=======================================================================
C     GOH MODEL: TWO FIBER FAMILIES
C=======================================================================
      C10     =PROPS(1)
      D       =PROPS(2)
      K1(1)   =PROPS(3)
      K2(1)   =PROPS(4)
      KAPPA(1)=PROPS(5)
      K1(2)   =PROPS(6)
      K2(2)   =PROPS(7)
      KAPPA(2)=PROPS(8)
C
C     Fiber directions
      beta=90.d0/180.d0*3.141592653589793d0

      NV(1,1)=dcos(beta)
      NV(2,1)=dsin(beta)
      NV(3,1)=0.d0

      NV(1,2)=dcos(beta)
      NV(2,2)=-dsin(beta)
      NV(3,2)=0.d0
C
C     Identity tensor
      eye=0.d0
      do i=1,3
          eye(i,i)=1.d0
      enddo
C
C     Deformation gradient, C, B, and invariants
      call detmat33(DFGRD1,detF)
      call rightCauchy(DFGRD1,C,detC,Cinv,I1,I2,I3)
      call LeftCauchy(DFGRD1,B,detB,Binv)
      J_F=detF
C
C     Structural tensor M_alpha=N_alpha tensor N_alpha
      Mfiber=0.d0
      do a=1,2
          do i=1,3
              do j=1,3
                  Mfiber(i,j,a)=NV(i,a)*NV(j,a)
              enddo
          enddo
      enddo
C
C     dCinv/dC
      dCinvdC=0.d0
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      dCinvdC(i,j,k,l)
     #                    =-0.5d0*(Cinv(i,k)*Cinv(j,l)
     #                    +Cinv(i,l)*Cinv(j,k))
                  enddo
              enddo
          enddo
      enddo
C
C     dJ_F/dC
      dJdC=0.d0
      do i=1,3
          do j=1,3
              dJdC(i,j)=0.5d0*J_F*Cinv(i,j)
          enddo
      enddo
C
C     d2J_F/dC2
      d2JdC2=0.d0
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      d2JdC2(i,j,k,l)
     #                    =0.25d0*J_F*Cinv(i,j)*Cinv(k,l)
     #                    +0.5d0*J_F*dCinvdC(i,j,k,l)
                  enddo
              enddo
          enddo
      enddo
C
C     Modified first invariant
      Jm23=J_F**(-2.d0/3.d0)
      I1bar=Jm23*I1
C
C     dI1bar/dC
      dI1bardC=0.d0
      do i=1,3
          do j=1,3
              dI1bardC(i,j)=Jm23*(eye(i,j)
     #                    -(1.d0/3.d0)*I1*Cinv(i,j))
          enddo
      enddo
C
C     d2I1bar/dC2
      d2I1bardC2=0.d0
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      d2I1bardC2(i,j,k,l)
     #                    =Jm23*(-(1.d0/3.d0)*
     #                    (Cinv(k,l)*eye(i,j)
     #                    +Cinv(i,j)*eye(k,l))
     #                    +(1.d0/9.d0)*I1*
     #                    Cinv(i,j)*Cinv(k,l)
     #                    -(1.d0/3.d0)*I1*
     #                    dCinvdC(i,j,k,l))
                  enddo
              enddo
          enddo
      enddo
C
C     Fiber invariants and derivatives
      I4=0.d0
      I4bar=0.d0
      dI4bardC=0.d0
      d2I4bardC2=0.d0
      Ealpha=0.d0
      dEalphadC=0.d0
      d2EalphadC2=0.d0
      Halpha=0.d0
C
      do a=1,2
          do i=1,3
              do j=1,3
                  I4(a)=I4(a)+Mfiber(i,j,a)*C(i,j)
              enddo
          enddo
          I4bar(a)=Jm23*I4(a)
C
C         dI4bar/dC
          do i=1,3
              do j=1,3
                  dI4bardC(i,j,a)=Jm23*(Mfiber(i,j,a)
     #                    -(1.d0/3.d0)*I4(a)*Cinv(i,j))
              enddo
          enddo
C
C         d2I4bar/dC2
          do i=1,3
              do j=1,3
                  do k=1,3
                      do l=1,3
                          d2I4bardC2(i,j,k,l,a)
     #                        =Jm23*(-(1.d0/3.d0)*
     #                        (Cinv(k,l)*Mfiber(i,j,a)
     #                        +Mfiber(k,l,a)*Cinv(i,j))
     #                        +(1.d0/9.d0)*I4(a)*
     #                        Cinv(i,j)*Cinv(k,l)
     #                        -(1.d0/3.d0)*I4(a)*
     #                        dCinvdC(i,j,k,l))
                      enddo
                  enddo
              enddo
          enddo
C
C         Fiber strain measure E_alpha
          Ealpha(a)=KAPPA(a)*(I1bar-3.d0)
     #              +(1.d0-3.d0*KAPPA(a))*
     #              (I4bar(a)-1.d0)
C
C         dE_alpha/dC
          do i=1,3
              do j=1,3
                  dEalphadC(i,j,a)
     #                =KAPPA(a)*dI1bardC(i,j)
     #                +(1.d0-3.d0*KAPPA(a))*
     #                dI4bardC(i,j,a)
              enddo
          enddo
C
C         d2E_alpha/dC2
          do i=1,3
              do j=1,3
                  do k=1,3
                      do l=1,3
                          d2EalphadC2(i,j,k,l,a)
     #                        =KAPPA(a)*d2I1bardC2(i,j,k,l)
     #                        +(1.d0-3.d0*KAPPA(a))*
     #                        d2I4bardC2(i,j,k,l,a)
                      enddo
                  enddo
              enddo
          enddo
C
C         Fiber activation
          if (Ealpha(a).gt.0.d0) then
              Halpha(a)=1.d0
          else
              Halpha(a)=0.d0
          endif
      enddo
C
C     Second Piola-Kirchhoff stress
      S=0.d0
C
C     Isotropic matrix contribution
      do i=1,3
          do j=1,3
              S(i,j)=S(i,j)+2.d0*C10*dI1bardC(i,j)
          enddo
      enddo
C
C     Anisotropic fiber contribution
      do a=1,2
          fiberCoeff(a)=Halpha(a)*2.d0*K1(a)*Ealpha(a)
     #                 *dexp(K2(a)*Ealpha(a)*Ealpha(a))
          do i=1,3
              do j=1,3
                  S(i,j)=S(i,j)+fiberCoeff(a)
     #                    *dEalphadC(i,j,a)
              enddo
          enddo
      enddo
C
C     Volumetric contribution
      do i=1,3
          do j=1,3
              S(i,j)=S(i,j)+(1.d0/D)*(J_F*J_F-1.d0)
     #                    *Cinv(i,j)
          enddo
      enddo
C
C     Cauchy stress
      sigma=0.d0
      do i=1,3
          do j=1,3
              do m=1,3
                  do n=1,3
                      sigma(i,j)=sigma(i,j)+1.d0/J_F
     #                    *DFGRD1(i,m)*S(m,n)
     #                    *DFGRD1(j,n)
                  enddo
              enddo
          enddo
      enddo
C
C     Material tangent: C0=2*dS/dC
      C0ijkl=0.d0
C
C     Isotropic contribution
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      C0ijkl(i,j,k,l)=4.d0*C10
     #                    *d2I1bardC2(i,j,k,l)
                  enddo
              enddo
          enddo
      enddo
C
C     Anisotropic contribution
      do a=1,2
          fiberTangent1(a)=Halpha(a)*4.d0*K1(a)
     #        *dexp(K2(a)*Ealpha(a)*Ealpha(a))
     #        *(1.d0+2.d0*K2(a)*Ealpha(a)*Ealpha(a))
          fiberTangent2(a)=Halpha(a)*4.d0*K1(a)
     #        *Ealpha(a)
     #        *dexp(K2(a)*Ealpha(a)*Ealpha(a))
          do i=1,3
              do j=1,3
                  do k=1,3
                      do l=1,3
                          C0ijkl(i,j,k,l)
     #                        =C0ijkl(i,j,k,l)
     #                        +fiberTangent1(a)
     #                        *dEalphadC(i,j,a)
     #                        *dEalphadC(k,l,a)
     #                        +fiberTangent2(a)
     #                        *d2EalphadC2(i,j,k,l,a)
                      enddo
                  enddo
              enddo
          enddo
      enddo
C
C     Volumetric contribution
      Cvolijkl=0.d0
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      Cvolijkl(i,j,k,l)=4.d0/D*
     #                    ((1.d0+J_F**(-2.d0))
     #                    *dJdC(i,j)*dJdC(k,l)
     #                    +(J_F-1.d0/J_F)
     #                    *d2JdC2(i,j,k,l))
                      C0ijkl(i,j,k,l)=C0ijkl(i,j,k,l)
     #                    +Cvolijkl(i,j,k,l)
                  enddo
              enddo
          enddo
      enddo
C
C     Push forward material tangent
      Cijkl=0.d0
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      do m=1,3
                          do n=1,3
                              do p=1,3
                                  do q=1,3
                                      Cijkl(i,j,k,l)
     #                                    =Cijkl(i,j,k,l)+1.d0/J_F
     #                                    *C0ijkl(m,n,p,q)
     #                                    *DFGRD1(i,m)
     #                                    *DFGRD1(j,n)
     #                                    *DFGRD1(k,p)
     #                                    *DFGRD1(l,q)
                                  enddo
                              enddo
                          enddo
                      enddo
                  enddo
              enddo
          enddo
      enddo
C
C     Geometric contribution for ABAQUS spatial tangent
      Cadd=0.d0
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      Cadd(i,j,k,l)=0.5d0*
     #                    (eye(i,k)*sigma(j,l)
     #                    +eye(i,l)*sigma(j,k)
     #                    +eye(j,k)*sigma(i,l)
     #                    +eye(j,l)*sigma(i,k))
                  enddo
              enddo
          enddo
      enddo
C
C     Total spatial tangent
      Caba=0.d0
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      Caba(i,j,k,l)=Cijkl(i,j,k,l)
     #                    +Cadd(i,j,k,l)
                  enddo
              enddo
          enddo
      enddo
C
C     Rearrange stress and tangent into ABAQUS notation
      do i=1,6
          k=ind1(i)
          l=ind2(i)
          STRESS(i)=sigma(k,l)
          do j=1,6
              m=ind1(j)
              n=ind2(j)
              DDSDDE(i,j)=Caba(k,l,m,n)
          enddo
      enddo
C=======================================================================
C     END OF GOH MODEL
C=======================================================================
      
C End of USER code
C---------------------------------------------------------------------------------------
      RETURN
      END
      
c***********************************************************************

      subroutine rightCauchy(F,C,detC,Cinv,I1,I2,I3)
      implicit none
      integer  nint,i,j,k
      real*8   F(3,3),C(3,3),Cinv(3,3),C2(3,3)
      real*8   detC,I1,I2,I3  

C     Calculate the right Cauchy Green tensor
      do i=1,3
          do j=1,3
                 C(i,j)=0.0d0
              do k=1,3
                 C(i,j)=C(i,j)+F(k,i)*F(k,j)
              enddo
          enddo
      enddo
c     C square and invariants
      C2=matmul(C,C)
      I1=C(1,1)+C(2,2)+C(3,3)
      I2=1.d0/2.d0*(I1*I1-(C2(1,1)+C2(2,2)+C2(3,3)))
       
      call detmat33(C,detC)
      
      I3=detC
      
      call invertmat33(C,detC,Cinv)
      return
      end     


c***********************************************************************
      subroutine LeftCauchy(F,B,detB,Binv)
      implicit none
      integer  nint,i,j,k
      real*8   F(3,3),B(3,3),Binv(3,3),B2(3,3)
      real*8   detB,I1,I2,I3  

C     Calculate the right Cauchy Green tensor
      B=matmul(F,transpose(F))
      
      call detmat33(B,detB)
      
      call invertmat33(B,detB,Binv)
      return
      end     

c***********************************************************************
      subroutine detmat33(a,deta)
      implicit none
      real*8   a(3,3)
      real*8   deta

      deta=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
     #    -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
     #    +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

      return
      end
c***********************************************************************
c***********************************************************************
      subroutine invertmat33(a,deta,ainv)
      implicit none
      real*8   a(3,3),deta
      real*8   detainv
      real*8   ainv(3,3)

      detainv=1.0d0/deta
      ainv(1,1)=(+a(2,2)*a(3,3)-a(2,3)*a(3,2))*detainv
      ainv(1,2)=(-a(1,2)*a(3,3)+a(1,3)*a(3,2))*detainv
      ainv(1,3)=(+a(1,2)*a(2,3)-a(1,3)*a(2,2))*detainv
      ainv(2,1)=(-a(2,1)*a(3,3)+a(2,3)*a(3,1))*detainv
      ainv(2,2)=(+a(1,1)*a(3,3)-a(1,3)*a(3,1))*detainv
      ainv(2,3)=(-a(1,1)*a(2,3)+a(1,3)*a(2,1))*detainv
      ainv(3,1)=(+a(2,1)*a(3,2)-a(2,2)*a(3,1))*detainv
      ainv(3,2)=(-a(1,1)*a(3,2)+a(1,2)*a(3,1))*detainv
      ainv(3,3)=(+a(1,1)*a(2,2)-a(1,2)*a(2,1))*detainv

      return
      end
c***********************************************************************

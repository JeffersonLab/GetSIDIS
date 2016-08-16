      Program Cteq6Demo

C     Simple test program for using the publicly released CTEQ6 PDFs.
C     See header file of Cteq6Pdf-2004.f for listing of the CTEQ6 PDFS
C         and detailed instructions.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      Character Ans*1
      
      Dimension Z(10), XX(10)
      
      Data Iset, Iptn / 1, 1 /
      Data x, Q / 0.D0, 0.0D0 /
      Data (Z(I), I=1,10) /0.05, .1, .2, .3, .4, .5, .6, .7, .8, .9/
      Save
C                 This array samples 10^-4 < x < 0.9 in a more "uniform" way then either linear or log.
      Do I = 1, 10
      XX(I) = Z(I) **3
      EndDo

 5    Print *, 
     >'Select PDF set: Iset, Iparton; press / to keep defaults'
     >      , Iset, Iptn
      Read *, Iset, Iptn
      Call SetCtq6(Iset)
      
 3    Print *, 'Input x, Q; or 0,0 for default X array & Q'
      Print '(A, F10.5, F10.2)', 'Press / to keep defaults', x, Q
      read *, x, Q
 
      Print '(A, 3I6, F10.2)', 'Iset, Ihadron, Iparton, ' 
     >                , Iset, Ihad, Iptn

      
      Print '(A)', '    Z           Q           X           PD '
      If (abs(x).LT.1D-6 .and. abs(Q).LT.1D-6) Then
         Q = 5D0
         Do I = 1, 10
           PD = Ctq6pdf (Iptn, xx(I), Q)
           Print '(F10.3, 3(1pE13.4))', Z(I), Q, XX(I), PD
         EndDo
      Else
         PD = Ctq6pdf (Iptn, x, Q)
         ZZ = X**(1./3)
         Print '(F10.3, 3(1pE13.4))', ZZ, Q, X, PD
      EndIf

      Print *, 'Do another x, Q ?'
      Read *, Ans
      If (Ans .eq. 'y') Goto 3
      
      print *, 'Do another Iset, Iptn ?'
      Read *, Ans
      If (Ans .eq. 'y') Goto 5

      Stop
C                        ****************************
      END

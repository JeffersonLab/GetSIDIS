	program szar
	implicit none

c-----  Declaring variables
	integer iprtn1, iprtn2, i,n
        double precision x1,q1,x2,q2, gluon, Ctq6Pdf

c-----  Plotting the gluon distribution of 41 pdfs in one loop ----- 
	iprtn1 = 0 
	iprtn2 = 1 

	q1 = 1000.0d0
	q2 = 100.0d0

	x1 = 0.3456d0
	x2 = 0.567d0

c	do i=1,41
	Call SetCtq6(400)	  
	do n=1,10000000
	    gluon = Ctq6Pdf(iprtn1, X1, Q1)
	    gluon = Ctq6Pdf(iprtn2, X2, Q2)
	end do  
c	end do
	
	stop
	end



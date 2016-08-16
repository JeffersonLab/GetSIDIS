	program szar
	implicit none

c-----  Declaring the external functions -----
	integer cteq_pdf_alloc_id,cteq_pdf_alloc_name
        double precision cteq_pdf_evolvepdf,cteq_pdf_evolveas

c-----  Declaring variables
	integer pdf(41), iprtn, i,n
        double precision x, q, gluon

c-----  Create the pdf ojects ----
c	do i=1,41
c	  pdf(i) = cteq_pdf_alloc_id(99+i)
c	end do
	pdf(1) = cteq_pdf_alloc_id(100)

c-----  Plotting the gluon distribution of 41 pdfs in one loop ----- 
	iprtn = 0 
	q = 1000.0d0
	x = 0.3456d0
	i=1

c	do i=1,41
	do n=1,10000000
	    gluon = cteq_pdf_evolvepdf(pdf(i), iprtn, x, q)
	  end do
c	end do  
	
c-----  Destroying the objects -----
	call cteq_pdf_free(pdf(i))
c	do i=1,41
c	  call cteq_pdf_free(pdf(i))
c	end do
	
	stop
	end



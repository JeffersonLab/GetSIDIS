	program szar
	implicit none

c-----  Declaring the external functions -----
	integer cteq_pdf_alloc_id,cteq_pdf_alloc_name
        double precision cteq_pdf_evolvepdf,cteq_pdf_evolveas

c-----  Declaring variables
	integer pdf(41), iprtn, i
        double precision x, q

c-----  Create the pdf ojects ----
	do i=1,41
	  pdf(i) = cteq_pdf_alloc_id(99+i)
	end do

c-----  Plotting the gluon distribution of 41 pdfs in one loop ----- 
	iprtn = 0 
	q = 1000.0d0
	x = 0.1d0

 10	write(*,*) x,(cteq_pdf_evolvepdf(pdf(i), iprtn, x, q), i=1,41)
	x = x + 0.1d0
	if(x < 1.0d0) goto 10
	
c-----  Destroying the objects -----
	do i=1,41
	  call cteq_pdf_free(pdf(i))
	end do
	
	stop
	end



#include <cteqpdf.h>
#include <stdio.h>
#include <stdlib.h>


extern double ctq6pdf_(int *, double *, double *);
extern int setctq6_(int *);


int main()
{
  /*   We compare CTEQ65.00  pdf.  */
  int mode = 400;
  
  /*   Initialize the old fortran interface.   */
  setctq6_(&mode);
  
  /*   Creates a pdf object and initialize with the CTEQ66.00 table   */
  cteq_pdf_t *pdf = cteq_pdf_alloc_id(mode);
  
  
  /*   we plot the old and new pdf */
  int iprtn = 0;
  double x, q = 1100.0, f, c;

  printf("#  q = %f   iprtn = %d   as(q)/2pi = %e\n", q, iprtn, cteq_pdf_evolveas(pdf, q));
  printf("#  x        old f77 interface       C interface          error\n");
  
  for(x = 1e-8; x <= 1e-7; x+=1e-8) {
	f = ctq6pdf_(&iprtn, &x, &q);
	c = cteq_pdf_evolvepdf(pdf, iprtn, x, q);
	printf("%e      %e         %e         %e\n", x, f, c, f-c);
  }
  
  iprtn = 1;
  q = 2030.45;
  
  printf("#  q = %f   iprtn = %d   as(q)/2pi = %e\n", q, iprtn, cteq_pdf_evolveas(pdf, q));
  printf("#  x        old f77 interface       C interface          error\n");
  
  for(x = 1.0; x > 0.999; x-=0.0001) {
	f = ctq6pdf_(&iprtn, &x, &q);
	c = cteq_pdf_evolvepdf(pdf, iprtn, x, q);
	printf("%f      %e         %e         %e\n", x, f, c, f-c);
  }
  
  iprtn = -3;
  q = 100.0;
  
  printf("#  q = %f   iprtn = %d   as(q)/2pi = %e\n", q, iprtn, cteq_pdf_evolveas(pdf, q));
  printf("#  x        old f77 interface       C interface          error\n");
  
  for(x = 1.0; x > 0.01; x-=0.1) {
	f = ctq6pdf_(&iprtn, &x, &q);
	c = cteq_pdf_evolvepdf(pdf, iprtn, x, q);
	printf("%f      %e         %e         %e\n", x, f, c, f-c);
  }
  
  /*   Before exit destroy the odf object  */
  cteq_pdf_free(pdf);
  
  return 0;
}


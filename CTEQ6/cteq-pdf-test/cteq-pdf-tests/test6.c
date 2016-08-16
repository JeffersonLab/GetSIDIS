#include <cteqpdf.h>
#include <stdio.h>
#include <stdlib.h>




int main()
{
  /*   Creating and initializing the pdf object   */
  unsigned int i;
  cteq_pdf_t *pdf = cteq_pdf_alloc_id(400);
  
  
  /*   we plot all the 41 pdfs */
  int iprtn1 = 0, iprtn2 = 1;
  double x1 = 0.3456, q1 = 1000.0, gluon;
  double x2 = 0.567, q2 = 100.0;

  for(i = 0; i < 10000000; i++) {
	gluon = cteq_pdf_evolvepdf(pdf, iprtn1, x1, q1);
	gluon = cteq_pdf_evolvepdf(pdf, iprtn2, x2, q2);
  }
  cteq_pdf_free(pdf);
  
  return 0;
}


#include <cteqpdf.h>
#include <stdio.h>
#include <stdlib.h>




int main()
{
  /*   Creating and initializing the pdf object   */
  unsigned int i;
  cteq_pdf_t *pdf[41];
  
  for(i = 0; i < 41; i++) 
	pdf[i] = cteq_pdf_alloc_id(100+i);
  
  
  /*   we plot all the 41 pdfs */
  int iprtn = 0;
  double x, q = 1000.0;

  printf("#  q = %f   iprtn = %d\n", q, iprtn);
  printf("#  x     and the result of the 41 pdfs \n");
  
  for(x = 1.0; x > 0.01; x-=0.1) {
	printf("%f  ", x);
	for(i = 0; i < 41; i++)
	  printf("  %e", cteq_pdf_evolvepdf(pdf[i], iprtn, x, q));
	printf("\n");
  }

  /*   Before exit destroy the pdf objects  */
  for(i = 0; i < 41; i++) 
	cteq_pdf_free(pdf[i]);
  
  return 0;
}


/*
 * CheckForIota.cpp
 *
 *  Created on: Apr 2, 2018
 *      Author: asgibso
 */

#include <Teuchos_config.h>

int main(int argc, char** argv)
{
#ifndef HAVE_TEUCHOS_COMPLEX
#error " Teuchos COMPLEX is not enabled! To enable it, set Teuchos_ENABLE_COMPLEX=ON."
#endif
  return 0;
}



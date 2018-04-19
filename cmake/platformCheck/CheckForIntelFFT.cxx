//
// Test program to check for Intel FFT.
//

#include "mkl_dfti.h"

int main(int argc, char** argv)
{
	float _Complex x[32];
	float y[34];
	DFTI_DESCRIPTOR_HANDLE my_desc1_handle;
	DFTI_DESCRIPTOR_HANDLE my_desc2_handle;
	MKL_LONG status;
	//...put input data into x[0],...,x[31]; y[0],...,y[31]
	for(int i = 0; i < 32; ++i) {
		x[i] = 1.0 + static_cast<float>(i);
	}

	status = DftiCreateDescriptor( &my_desc1_handle, DFTI_SINGLE, DFTI_COMPLEX, 1, 32);
	status = DftiCommitDescriptor( my_desc1_handle );
	status = DftiComputeForward( my_desc1_handle, x);
	status = DftiFreeDescriptor(&my_desc1_handle);
	/* result is x[0], ..., x[31]*/
	status = DftiCreateDescriptor( &my_desc2_handle, DFTI_SINGLE, DFTI_REAL, 1, 32);
	status = DftiCommitDescriptor( my_desc2_handle);
	status = DftiComputeForward( my_desc2_handle, y);
	status = DftiFreeDescriptor(&my_desc2_handle);

	return 0;
}

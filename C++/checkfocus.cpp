#include "mex.h"
#include <windows.h>

bool matlab_has_focus(int matlab_hdl)
{
	return (int(GetForegroundWindow()) == matlab_hdl);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
	bool runfunc = true;

	if (nrhs != 1)
	{
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","One input required.");
		runfunc = false;
	}

	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfElements(prhs[0]) != 1)
	{
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input multiplier must be a scalar.");
		runfunc = false;
	}
	
	if (runfunc)
	{
		plhs[0] = mxCreateLogicalScalar(matlab_has_focus(mxGetScalar(prhs[0])));
	}

}
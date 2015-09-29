
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string.h>
#include <cmath>
#include <limits>
#include <assert.h>

#include "src/kovtun.h"
#include "src/energy.h"

#include "mex.h"

#define MATLAB_ASSERT(expr,msg) if (!(expr)) { mexErrMsgTxt(msg);}

#if !defined(MX_API_VER) || MX_API_VER < 0x07030000
typedef int mwSize;
typedef int mwIndex;
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MATLAB_ASSERT(nrhs == 2, "Wrong number of input arguments, expected 2" );
	MATLAB_ASSERT(nlhs <= 1, "Too many output arguments, expected 1"); 	
	
	//Fix input parameter order:
	const mxArray *uInPtr = (nrhs > 0) ? prhs[0] : NULL; //unary
	const mxArray *pInPtr = (nrhs > 1) ? prhs[1] : NULL; //pairwise
	
	//Fix output parameter order:
	mxArray **sOutPtr = (nlhs > 0) ? &plhs[0] : NULL; //solution
	
	// get unary potentials
	MATLAB_ASSERT(mxGetNumberOfDimensions(uInPtr) == 2, "Unary term array is not 2-dimensional");
	MATLAB_ASSERT(mxGetPi(uInPtr) == NULL, "Unary potentials should not be complex");
	
	mwSize numNodes = mxGetN(uInPtr);
	mwSize numLabels = mxGetM(uInPtr);

	MATLAB_ASSERT(numNodes >= 1, "The number of nodes is not positive");
	MATLAB_ASSERT(numLabels >= 1, "The number of labels is not positive");
	MATLAB_ASSERT(mxGetClassID(uInPtr) == mxDOUBLE_CLASS, "Expected mxDOUBLE_CLASS for input unary term argument");
	double* termW = (double*)mxGetData(uInPtr);

	//get pairwise potentials
	MATLAB_ASSERT(mxIsSparse(pInPtr), "Expected sparse array for neighbours");
	MATLAB_ASSERT(mxGetN(pInPtr) == numNodes && mxGetM(pInPtr) == numNodes,
	              "Neighbours array must be NumNodes x NumNodes in size");
	MATLAB_ASSERT(mxGetClassID(pInPtr) == mxDOUBLE_CLASS, "Expected mxDOUBLE_CLASS for neighbours array");
	MATLAB_ASSERT(mxGetPi(pInPtr) == NULL, "Pairwise potentials should not be complex");

	mwIndex colNum = (mwIndex)mxGetN(pInPtr);
	const mwIndex* ir = mxGetIr(pInPtr);
	const mwIndex* jc = mxGetJc(pInPtr);
	double*        pr = mxGetPr(pInPtr);

	//check pairwise terms
	mwSize numEdges = 0;
	for (mwIndex c = 0; c < colNum; ++c) {
		mwIndex rowStart = jc[c]; 
		mwIndex rowEnd   = jc[c+1]; 
		for (mwIndex ri = rowStart; ri < rowEnd; ++ri)  {
			mwIndex r = ir[ri];

			double dw = pr[ri];
			if( r < c) {
				++numEdges;
				MATLAB_ASSERT(dw >= 0, "Potts potentials should be positive");
			}
		}
	}

	
	//create MRF object
	Energy* energy = NULL;

	energy =  new Energy( numLabels );
	energy -> nvar = numNodes;
	energy -> npair = numEdges; 
	energy -> allocateMemory();

	
	// construct energy
	// add unary terms
	for(int iNode = 0; iNode < numNodes; ++iNode){
		for(int iLabel = 0; iLabel < numLabels; ++iLabel) {
			energy -> unaryCost[iNode][iLabel] = termW[iNode * numLabels + iLabel];
		}
	}

	//add pairwise terms
	int iEdge = 0;
	for (mwIndex c = 0; c < colNum; ++c) {
		mwIndex rowStart = jc[c]; 
		mwIndex rowEnd   = jc[c + 1]; 
		for (mwIndex ri = rowStart; ri < rowEnd; ++ri)  {
			mwIndex r = ir[ri];
			double dw = pr[ri];

			if( r < c) {
				energy -> pairIndex[iEdge][0] = c;	
				energy -> pairIndex[iEdge][1] = r;
				energy -> pairCost[iEdge] = dw;
				++iEdge;
			}
		}
	}	

	// run the method
	Grapht** graph = new Grapht* [numLabels];
	for(int  iLabel =0; iLabel < numLabels; ++iLabel)
		graph[ iLabel ] = new Grapht(energy -> nvar, energy -> npair);
	Kovtun* solver = new Kovtun(energy, graph, 1); // 1 => use multiple graphs and the Reuse method (see Alahari et al., 2008)

	// Projecting the energy function (Reduce method, see Alahari et al., 2008) 
	for(int iLabel = 0; iLabel < numLabels; ++iLabel)
	{
		solver -> findPersistent(iLabel);
		// energy -> Project( solver -> multiSolution );
	}

	// output the partial labeling
	if(sOutPtr != NULL)	{
		*sOutPtr = mxCreateNumericMatrix(numNodes, 1, mxDOUBLE_CLASS, mxREAL);
		double* segment = (double*)mxGetData(*sOutPtr);
	
		// Compute only the partially optimal solution 
		for(int iNode = 0; iNode < numNodes; ++iNode) {
			int lbl = solver->multiSolution[ iNode ];
			
			if(lbl != NOLABEL)
				segment[iNode] = lbl + 1;
			else
				segment[iNode] = 0;
		}
	}

	// done
	for(int iLabel = 0; iLabel < numLabels; ++iLabel)
		delete graph[ iLabel ];
	delete [] graph;
	delete solver;
	delete energy;
	
}


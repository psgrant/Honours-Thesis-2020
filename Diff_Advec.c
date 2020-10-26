#include <mex.h>
#include <matrix.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/*Initialise functions*/
int linearIndex();


/*Calculates flux of a diffusionadvectioj equation in 2D*/



/*Gateway Function*/
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {


	/*Ensure correct inputs and outputs*/
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "There is exactly 1 output");
	}

	if (nrhs != 9) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "This function requires 11 inputs");
	}






	double* phi, * phiOld, * SCVFacelengths, * SCVNormalsx, * SCVNormalsy, * CVAreas, * parameters, * nodeData, * elementData, * out;     /*Preallocate double array pointers*/

	/*Load in all the array data*/
	phi = mxGetPr(prhs[0]);
	phiOld = mxGetPr(prhs[1]);
	nodeData = mxGetPr(prhs[2]);
	SCVFacelengths = mxGetPr(prhs[3]);
	SCVNormalsx = mxGetPr(prhs[4]);
	SCVNormalsy = mxGetPr(prhs[5]);
	CVAreas = mxGetPr(prhs[6]);
	elementData = mxGetPr(prhs[7]);
	parameters = mxGetPr(prhs[8]);


	/* Load in all the input scalar variables*/
	double Dx = parameters[0];
	double Dy = parameters[1];
	double deltat = parameters[2];
	double theta = parameters[3];
	double xc = parameters[4];
	double yc = parameters[5];
	double curTime = parameters[6];
	double vx = parameters[7];
	double vy = parameters[8];

	/*Get number of nodes and elements*/
	const int NUMBER_OF_NODES = mxGetM(prhs[0]);
	const int NUMBER_OF_ELEMENTS = mxGetM(prhs[3]);
	double* F, * flux;

	int nodeOrder[3] = { 1, 2, 0 };

	double 	shapeFuncA, shapeFuncXPhi, shapeFuncYPhi, shapeFuncXPhiOld, shapeFuncYPhiOld, advectionPhi, advectionPhiOld, nDotv, explicit, implicit, SCVNx, SCVNy, faceLength, FVMConst, dirichletCondition, dirichletX, dirichletY;
	int nodeNumber[3];
	double nodePhi[3], nodePhiOld[3], nodeXLoc[3], nodeYLoc[3];



	F = (double*)mxMalloc(NUMBER_OF_NODES * sizeof(double));
	flux = (double*)mxCalloc(NUMBER_OF_NODES, sizeof(double));

	//F = mxGetPr(F_Out);
	/*Loop through each element*/
	for (int element = 0; element < NUMBER_OF_ELEMENTS; element++) {



		for (int i = 0; i < 3; i++) {

			/* Node numbers of the nodes in a given element*/
			nodeNumber[i] = elementData[linearIndex(element, i, NUMBER_OF_ELEMENTS)] - 1;

			/* Isolate x and y data*/
			nodeXLoc[i] = nodeData[linearIndex(nodeNumber[i], 1, NUMBER_OF_NODES)];
			nodeYLoc[i] = nodeData[linearIndex(nodeNumber[i], 2, NUMBER_OF_NODES)];

			/*Isolate phi and phiOld data*/
			nodePhi[i] = phi[nodeNumber[i]];
			nodePhiOld[i] = phiOld[nodeNumber[i]];

		}

		/*Calcualting the shape functions of each element for phi and phiold
		Remebering that it is constant throughout the whole element*/

		/*Get the divisor of cramers rule for the shapefunction of the element*/
		shapeFuncA = nodeXLoc[0] * (nodeYLoc[1] - nodeYLoc[2]) + nodeXLoc[1] * (nodeYLoc[2] - nodeYLoc[0]) + nodeXLoc[2] * (nodeYLoc[0] - nodeYLoc[1]);

		/*Get the numerator of the x term of the phi shape function*/
		shapeFuncXPhi = nodePhi[0] * (nodeYLoc[1] - nodeYLoc[2]) + nodePhi[1] * (nodeYLoc[2] - nodeYLoc[0]) + nodePhi[2] * (nodeYLoc[0] - nodeYLoc[1]);

		shapeFuncXPhi = (double)shapeFuncXPhi / shapeFuncA;

		/* Get the numerator for the y term of the phi shape function*/
		shapeFuncYPhi = nodePhi[0] * (nodeXLoc[2] - nodeXLoc[1]) + nodePhi[1] * (nodeXLoc[0] - nodeXLoc[2]) + nodePhi[2] * (nodeXLoc[1] - nodeXLoc[0]);

		shapeFuncYPhi = (double)shapeFuncYPhi / shapeFuncA;


		/*Get the numerator of the x term of the phiOld shape function*/
		shapeFuncXPhiOld = nodePhiOld[0] * (nodeYLoc[1] - nodeYLoc[2]) + nodePhiOld[1] * (nodeYLoc[2] - nodeYLoc[0]) + nodePhiOld[2] * (nodeYLoc[0] - nodeYLoc[1]);

		shapeFuncXPhiOld = (double)shapeFuncXPhiOld / shapeFuncA;

		/* Get the numerator for the y term of the phiOld shape function*/
		shapeFuncYPhiOld = nodePhiOld[0] * (nodeXLoc[2] - nodeXLoc[1]) + nodePhiOld[1] * (nodeXLoc[0] - nodeXLoc[2]) + nodePhiOld[2] * (nodeXLoc[1] - nodeXLoc[0]);

		shapeFuncYPhiOld = (double)shapeFuncYPhiOld / shapeFuncA;


		/* Loop over each node in the element (3 nodes)*/
		for (int node = 0; node < 3; node++) {

			SCVNx = SCVNormalsx[linearIndex(element, node, NUMBER_OF_ELEMENTS)];
			SCVNy = SCVNormalsy[linearIndex(element, node, NUMBER_OF_ELEMENTS)];
			faceLength = SCVFacelengths[linearIndex(element, node, NUMBER_OF_ELEMENTS)];
	
			// Apply averaging
			advectionPhi = 0.5 * (nodePhi[node] + nodePhi[nodeOrder[node]]);
			advectionPhiOld = 0.5 * (nodePhiOld[node] + nodePhiOld[nodeOrder[node]]);
			



			/*Calculate the implicit (theta = 1) and explicit (theta = 0) parts of the scheme*/
			implicit = theta * (SCVNx * faceLength * (Dx * shapeFuncXPhi - vx * advectionPhi) + SCVNy * faceLength * (Dy * shapeFuncYPhi - vy * advectionPhi));

			explicit = (1 - theta) * (SCVNx * faceLength * (Dx * shapeFuncXPhiOld - vx * advectionPhiOld) + SCVNy * faceLength * (Dy * shapeFuncYPhiOld - vy * advectionPhiOld));


			/* Add on to the running flux array for each node,
			also update the negative of that for the adjacent node*/
			flux[nodeNumber[node]] = flux[nodeNumber[node]] + (implicit + explicit);
			flux[nodeNumber[nodeOrder[node]]] = flux[nodeNumber[nodeOrder[node]]] - (implicit + explicit);
		}
	}

	/* Loop through every node and check to apply the "phi[node] - phiOld[node]* part for each node. This couldn't have been done above as the F array is updated at least twice for every node (once for each element). Also apply the FVM constant to the flux that was calculated above*/

	for (int node = 0; node < NUMBER_OF_NODES; node++) {


		// if the node is an internal node
		if (nodeData[linearIndex(node, 0, NUMBER_OF_NODES)] == 3) {

			FVMConst = deltat / CVAreas[node];
			F[node] = phi[node] - phiOld[node] - FVMConst * flux[node];
		}
		else {
			dirichletX = nodeData[linearIndex(node, 1, NUMBER_OF_NODES)];
			dirichletY = nodeData[linearIndex(node, 2, NUMBER_OF_NODES)];


			dirichletCondition = (double)1 / (4 * curTime + 1) * exp(-(pow(dirichletX - vx * curTime - xc, 2)) / (Dx * (4 * curTime + 1)) - (pow(dirichletY - vy * curTime - yc, 2)) / (Dy * (4 * curTime + 1)));

			F[node] = phi[node] - dirichletCondition;


		}

	}

	plhs[0] = mxCreateDoubleMatrix(NUMBER_OF_NODES, 1, mxREAL);
	out = mxGetPr(plhs[0]);

	for (int k = 0; k < NUMBER_OF_NODES; k++) {
		out[k] = F[k];
	}

	mxFree(F);
	mxFree(flux);

	return;
}

// Function that returns the linear index of a 2D array.
int linearIndex(int row, int col, int nrows) {
	int index = row + nrows * col;
	return(index);
}
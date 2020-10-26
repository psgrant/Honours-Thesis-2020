#include <mex.h>
#include <matrix.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/*Initialise functions*/
int linearIndex();
void generatekPAndPsiP();
void generatePsiP();
void generateKP();

/*Calculates flux of a diffusionadvectioj equation in 2D*/



/*Gateway Function*/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

	/*Ensure correct inputs and outputs*/
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "There is exactly 1 output");
	}

	if (nrhs != 23) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "This function requires 23 inputs");
	}

	double* h, * hOld, * H, * SCVFacelengths, * SCVNormalsx, * SCVNormalsy, * CVAreas, * Kxx, * Kzz, * aPar, * nPar, * mPar, * psiRes, * psiSat, * SCVAreas, * parameters, * nodeBLengths, * nodeCLengths, * HOld, * nodeData, * shapeFuncOld;     /*Preallocate double array pointers*/
	double * elementData, * nodeClass; /*Prealocate int array pointers*/

	int nodeOrder[3] = { 1, 2, 0 }, localNodeClass, nodeNumber[3];

	double 	shapeFuncA, shapeFuncXH, shapeFuncYH, shapeFuncXHOld, shapeFuncYHOld, faceH, faceHOld, nDotv, explicit, implicit, SCVNx, SCVNy, faceLength, FVMConst, dirichletCondition, dirichletX, dirichletY, localKxx, localKzz, facekP, facekPOld, localpsiP, localpsiPOld, riverFlux, localH, tempEvapo;
	double* psiP, * psiPOld, * kP, * kPOld, * out;
	double nodeH[3], nodeHOld[3], nodeXLoc[3], nodeYLoc[3];


	/*Load in all the array data*/
	h = mxGetPr(prhs[0]);
	hOld = mxGetPr(prhs[1]);
	nodeData = mxGetPr(prhs[2]);
	nodeClass = mxGetPr(prhs[3]);
	SCVFacelengths = mxGetPr(prhs[4]);
	SCVNormalsx = mxGetPr(prhs[5]);
	SCVNormalsy = mxGetPr(prhs[6]);
	CVAreas = mxGetPr(prhs[7]);
	elementData = mxGetPr(prhs[8]);
	Kxx = mxGetPr(prhs[9]);
	Kzz = mxGetPr(prhs[10]);
	aPar = mxGetPr(prhs[11]);
	nPar = mxGetPr(prhs[12]);;
	mPar = mxGetPr(prhs[13]);;
	psiRes = mxGetPr(prhs[14]);
	psiSat = mxGetPr(prhs[15]);
	SCVAreas = mxGetPr(prhs[16]);
	nodeBLengths = mxGetPr(prhs[17]);
	nodeCLengths = mxGetPr(prhs[18]);
	shapeFuncOld = mxGetPr(prhs[19]);
	psiPOld = mxGetPr(prhs[20]);
	kPOld = mxGetPr(prhs[21]);
	parameters = mxGetPr(prhs[22]);

	/* Load in all the input scalar variables*/
	double upwinding = parameters[0];
	double fluxLimiting = parameters[1];
	double deltat = parameters[2];
	double theta = parameters[3];
	double curTime = parameters[4];
	double rainfall = parameters[5];
	double Rb = parameters[6];
	double Hr = parameters[7];
	double HxR = parameters[8];
	double Kr = parameters[9];
	double pumpingRatesPerNode = parameters[10];
	double pumpingTurnedOn = parameters[11];

	/*Get number of nodes and elements*/
	const int NUMBER_OF_NODES = mxGetM(prhs[1]);
	const int NUMBER_OF_ELEMENTS = mxGetM(prhs[4]);
	const int NUMBER_OF_CORNER_NODES = mxGetM(prhs[19]);
	double* F, * flux;

	kP = (double*) mxCalloc(NUMBER_OF_NODES , sizeof(double));
	//kPOld = (double*) mxCalloc(NUMBER_OF_NODES , sizeof(double));
	psiP = (double*) mxCalloc(NUMBER_OF_NODES , sizeof(double));
	//psiPOld = (double*) mxCalloc(NUMBER_OF_NODES , sizeof(double));

	H = (double*) mxCalloc(NUMBER_OF_NODES, sizeof(double));
	HOld = (double*) mxCalloc(NUMBER_OF_NODES, sizeof(double));

	for (int i = 0; i < NUMBER_OF_NODES; i++) {
		H[i] = h[i] + nodeData[linearIndex(i, 1, NUMBER_OF_NODES)];
		HOld[i] = hOld[i] + nodeData[linearIndex(i, 1, NUMBER_OF_NODES)];
	}
	
	// get pointers to kp and psiP arrays
	/*generateKP(h, nPar, elementData, mPar, aPar, SCVAreas, CVAreas, NUMBER_OF_ELEMENTS, NUMBER_OF_NODES, kP);
	generateKP(hOld, nPar, elementData, mPar, aPar, SCVAreas, CVAreas, NUMBER_OF_ELEMENTS, NUMBER_OF_NODES, kPOld);
	generatePsiP(h, nPar, elementData, mPar, aPar, SCVAreas, CVAreas, psiRes, psiSat, NUMBER_OF_ELEMENTS, NUMBER_OF_NODES, psiP);
	generatePsiP(hOld, nPar, elementData, mPar, aPar, SCVAreas, CVAreas, psiRes, psiSat, NUMBER_OF_ELEMENTS, NUMBER_OF_NODES, psiPOld);*/
	//(double* h, double* nPar, double* elementData, double* mPar, double* aPar, double* SCVAreas, double* CVAreas, double* psiRes, double* psiSat, int NUMBER_OF_ELEMENTS, int NUMBER_OF_NODES, double* psiP)

	generatekPAndPsiP(h, nPar, elementData, mPar, aPar, SCVAreas, CVAreas, psiRes, psiSat, NUMBER_OF_ELEMENTS, NUMBER_OF_NODES, kP, psiP);
	//generatekPAndPsiP(hOld, nPar, elementData, mPar, aPar, SCVAreas, CVAreas, psiRes, psiSat, NUMBER_OF_ELEMENTS, NUMBER_OF_NODES, kPOld, psiPOld);

	F = (double*) mxCalloc(NUMBER_OF_NODES , sizeof(double));
	flux = (double*) mxCalloc(NUMBER_OF_NODES , sizeof(double));

	//F = mxGetPr(F_Out);
	/*Loop through each element*/
	for (int element = 0; element < NUMBER_OF_ELEMENTS; element++) {

		localKxx = Kxx[element];
		localKzz = Kzz[element];

		for (int i = 0; i < 3; i++) {

			/* Node numbers of the nodes in a given element*/
			nodeNumber[i] = elementData[linearIndex(element, i, NUMBER_OF_ELEMENTS)] - 1;

			/* Isolate x and y data*/

			nodeXLoc[i] = nodeData[linearIndex(nodeNumber[i], 0, NUMBER_OF_NODES)];
			nodeYLoc[i] = nodeData[linearIndex(nodeNumber[i], 1, NUMBER_OF_NODES)];

			/*Isolate h and hOld data*/
			nodeH[i] = H[nodeNumber[i]];
			nodeHOld[i] = HOld[nodeNumber[i]];

		}

		/*Calcualting the shape functions of each element for h and phiold
		Remebering that it is constant throughout the whole element*/

		/*Get the divisor of cramers rule for the shapefunction of the element*/
		shapeFuncA = nodeXLoc[0] * (nodeYLoc[1] - nodeYLoc[2]) + nodeXLoc[1] * (nodeYLoc[2] - nodeYLoc[0]) + nodeXLoc[2] * (nodeYLoc[0] - nodeYLoc[1]);

		/*Get the numerator of the x term of the h shape function*/
		shapeFuncXH = nodeH[0] * (nodeYLoc[1] - nodeYLoc[2]) + nodeH[1] * (nodeYLoc[2] - nodeYLoc[0]) + nodeH[2] * (nodeYLoc[0] - nodeYLoc[1]);

		shapeFuncXH = (double)shapeFuncXH / shapeFuncA;

		/* Get the numerator for the y term of the h shape function*/
		shapeFuncYH = nodeH[0] * (nodeXLoc[2] - nodeXLoc[1]) + nodeH[1] * (nodeXLoc[0] - nodeXLoc[2]) + nodeH[2] * (nodeXLoc[1] - nodeXLoc[0]);

		shapeFuncYH = (double)shapeFuncYH / shapeFuncA;

		shapeFuncXHOld = shapeFuncOld[linearIndex(element, 0, NUMBER_OF_ELEMENTS)];
		shapeFuncYHOld = shapeFuncOld[linearIndex(element, 1, NUMBER_OF_ELEMENTS)];

		/*Get the numerator of the x term of the hOld shape function*/
		//shapeFuncXHOld = nodeHOld[0] * (nodeYLoc[1] - nodeYLoc[2]) + nodeHOld[1] * (nodeYLoc[2] - nodeYLoc[0]) + nodeHOld[2] * (nodeYLoc[0] - nodeYLoc[1]);

		//shapeFuncXHOld = (double)shapeFuncXH / shapeFuncA;

		///* Get the numerator for the y term of the hOld shape function*/
		//shapeFuncYHOld = nodeHOld[0] * (nodeXLoc[2] - nodeXLoc[1]) + nodeHOld[1] * (nodeXLoc[0] - nodeXLoc[2]) + nodeHOld[2] * (nodeXLoc[1] - nodeXLoc[0]);

		//shapeFuncYHOld = (double)shapeFuncYH / shapeFuncA;


		/* Loop over each node in the element (3 nodes)*/
		for (int node = 0; node < 3; node++) {

			SCVNx = SCVNormalsx[linearIndex(element, node, NUMBER_OF_ELEMENTS)];
			SCVNy = SCVNormalsy[linearIndex(element, node, NUMBER_OF_ELEMENTS)];
			faceLength = SCVFacelengths[linearIndex(element, node, NUMBER_OF_ELEMENTS)];
			
			/*Apply upwinding, set by upwinding input parameter
			  Defaults to averaging if set to zero.*/
			if (upwinding == 1.0) {


				/*Decide which node's h and hOld data to use*/
				if (nodeH[node] >= nodeH[nodeOrder[node]]) {

					// Flow is in the direction of the SCV normal face 
					facekP = kP[nodeNumber[node]];
					facekPOld = kPOld[nodeNumber[node]];


					// Flow is not in the direction of th SCV normal face
				}
				else {

					facekP = kP[nodeNumber[nodeOrder[node]]];
					facekPOld = kPOld[nodeNumber[nodeOrder[node]]];
				}

				// Apply averaging
			}
			else {

				//facekP = (kP[nodeNumber[node]] * kP[nodeNumber[nodeOrder[node]]]) / (kP[nodeNumber[node]] + kP[nodeNumber[nodeOrder[node]]]);
				//facekPOld = (kPOld[nodeNumber[node]] * kPOld[nodeNumber[nodeOrder[node]]]) / (kPOld[nodeNumber[node]] + kPOld[nodeNumber[nodeOrder[node]]]);


				facekP = 0.5 * (kP[nodeNumber[node]] + kP[nodeNumber[nodeOrder[node]]]);
				facekPOld = 0.5 * (kPOld[nodeNumber[node]] + kPOld[nodeNumber[nodeOrder[node]]]);
			}
			/*Calculate the implicit (theta = 1) and explicit (theta = 0) parts of the scheme*/
			implicit = theta * faceLength *((SCVNx * facekP * localKxx * shapeFuncXH) + (SCVNy * facekP * localKzz * shapeFuncYH));

			explicit = (1 - theta) * faceLength  * ((SCVNx * facekPOld * localKxx * shapeFuncXHOld) + (SCVNy * facekPOld * localKzz * shapeFuncYHOld));


			/* Add on to the running flux array for each node,
			also update the negative of that for the adjacent node*/
			flux[nodeNumber[node]] = flux[nodeNumber[node]] + (implicit + explicit);
			flux[nodeNumber[nodeOrder[node]]] = flux[nodeNumber[nodeOrder[node]]] - (implicit + explicit);
		}
	}

	/* Loop through every node and check to apply the "phi[node] - phiOld[node]* part for each node. This couldn't have been done above as the F array is updated at least twice for every node (once for each element). Also apply the FVM constant to the flux that was calculated above*/
	
	for (int node = 0; node < NUMBER_OF_NODES; node++) {

		FVMConst = deltat / CVAreas[node];
		localNodeClass = nodeClass[node];
		// Interior nodes and all no flux boundaries
		F[node] = psiP[node] - psiPOld[node] - FVMConst * flux[node];
		localH = H[node];

		if (localH < 100) {

		switch (localNodeClass) {

		case 4: // Top Boundary
			F[node] = F[node] - FVMConst * nodeBLengths[node] * rainfall;
			break;

		case 5: // River boundary nodes (Rb < z < 100)
			riverFlux = 0;
			if (localH <= 100) {
				riverFlux = Kr * (Hr - H[node]) / HxR;
				if (localH < Rb) {
					riverFlux = Kr * (Hr - Rb) / HxR;
				}
			}
			F[node] = F[node] - FVMConst * nodeBLengths[node] * riverFlux;
			break;

		case 9: // Top right corner
			F[node] = F[node] - FVMConst * nodeCLengths[3] * rainfall;
			//if (psiP[node] > 0.025) {
			//	F[node] = F[node] + deltat * 0.035 * rainfall * pow(nodeData[linearIndex(node, 1, NUMBER_OF_NODES)] - 90, 2) / 100;
			//}
			break;

		case 10: // Top left corner (Rainfall and River)
			riverFlux = 0;
			if (localH <= 100) {
				riverFlux = Kr * (Hr - H[node]) / HxR;
				if (localH < Rb) {
					riverFlux = Kr * (Hr - Rb) / HxR;
				}
			}
			F[node] = F[node] - FVMConst * (nodeCLengths[13] * rainfall + nodeCLengths[5] * riverFlux);
			break;

		case 11: // Bottom river boundary (z = Rb)
			riverFlux = 0;
			if (localH <= 100) {
				riverFlux = Kr * (Hr - H[node]) / HxR;
				if (localH < Rb) {
					riverFlux = Kr * (Hr - Rb) / HxR;
				}
			}
			F[node] = F[node] - FVMConst *  nodeCLengths[14] * riverFlux;
			break;

		case 12: // Town pumping
			if (curTime > (365 * pumpingTurnedOn)) {
				F[node] = F[node] + deltat * pumpingRatesPerNode;
			}
			
			break;

		case 13: // Reparing vegetation zone
			
			if (psiP[node] > 0.165) {
					F[node] = F[node] + deltat * 0.05 * rainfall * pow(nodeData[linearIndex(node, 1, NUMBER_OF_NODES)] - 85 ,2) / 225;
			}

			if (nodeData[linearIndex(node, 1, NUMBER_OF_NODES)] == 100) {
				F[node] = F[node] - FVMConst * nodeBLengths[node] * rainfall;
			}

			break;

		case 14: // Crop zone


			if (psiP[node] > 0.165) {
				tempEvapo = nodeData[linearIndex(node, 1, NUMBER_OF_NODES)] - 95;
				tempEvapo = pow(tempEvapo, 2);
				F[node] = F[node] + deltat * 0.025 * rainfall *  tempEvapo / 25;
			}

			if (nodeData[linearIndex(node, 1, NUMBER_OF_NODES)] == 100) {
				F[node] = F[node] - FVMConst * nodeBLengths[node] * rainfall;
			}

			break;

		case 15: // plantation zone


			if (psiP[node] > 0.025) {
				tempEvapo = nodeData[linearIndex(node, 1, NUMBER_OF_NODES)] - 90;
				tempEvapo = pow(tempEvapo, 2);
				F[node] = F[node] + deltat * 0.0035 * rainfall * tempEvapo / 100;
			}

			if (nodeData[linearIndex(node, 1, NUMBER_OF_NODES)] == 100) {
				F[node] = F[node] - FVMConst * nodeBLengths[node] * rainfall;
			}

			break;

		default: // otherwise (no addition flux or sources)
			
			
			break;
		}
		}

	}


	// Output to MATLAB
	plhs[0] = mxCreateDoubleMatrix(NUMBER_OF_NODES, 1, mxREAL);
	out = mxGetPr(plhs[0]);

	for (int k = 0; k < NUMBER_OF_NODES; k++) {
		out[k] = F[k];
	}

	// Deallocate memory from calloc and mxNumericArray
	mxFree(psiP);
	//mxFree(psiPOld);
	mxFree(kP);
	//mxFree(kPOld);
	mxFree(F);
	mxFree(flux);
	mxFree(H);
	mxFree(HOld);

	return;
}

// Function that returns the linear index of a 2D array.
int linearIndex(int row, int col, int nrows) {
	int index = row + nrows * col;
	return(index);
}



// generate psip and kp
void generatekPAndPsiP(double* h, double* nPar, double* elementData, double* mPar, double* aPar, double* SCVAreas, double* CVAreas, double* psiRes, double* psiSat, int NUMBER_OF_ELEMENTS, int NUMBER_OF_NODES, double* kP, double* psiP) {

	// preallocate variables
	int nodeNumber;
	double nodeh, S, SCVA, CVA;
	double elementAPar;
	double elementNPar;
	double elementMPar;
	double temp1, temp2, temp3;
	double elementpsiRes, elementpsiSat;

	for (int element = 0; element < NUMBER_OF_ELEMENTS; element++) {

		// extract element constanta
		elementAPar = aPar[element];
		elementNPar = nPar[element];
		elementMPar = mPar[element];
		elementpsiRes = psiRes[element];
		elementpsiSat = psiSat[element];


		// loop over the three nodes
		for (int node = 0; node < 3; node++) {

			// Extract data about the given node
			nodeNumber = elementData[linearIndex(element, node, NUMBER_OF_ELEMENTS)] - 1;
			SCVA = SCVAreas[linearIndex(element, node, NUMBER_OF_ELEMENTS)];
			nodeh = h[nodeNumber];
			CVA = CVAreas[nodeNumber];

			// calcualte kP
			if (nodeh <= 0) {

				S = pow((1 + pow((-elementAPar * nodeh), elementNPar)), -elementMPar);

				
				temp1 = pow((S), 1 / elementMPar);
				temp2 = pow((1 - temp1), elementMPar);
				temp3 = pow(1 - temp2, 2);
				kP[nodeNumber] = kP[nodeNumber] + pow(S, 0.5) * temp3 * SCVA / CVA;
				psiP[nodeNumber] = psiP[nodeNumber] + (elementpsiRes + S * (elementpsiSat - elementpsiRes)) * SCVA / CVA;

			}
			else {
				kP[nodeNumber] = kP[nodeNumber] + SCVA / CVA;
				psiP[nodeNumber] = psiP[nodeNumber] + elementpsiSat * SCVA / CVA;
			}

			if (S > 1) {
				kP[nodeNumber] = 1;
			}


		}

	}
}
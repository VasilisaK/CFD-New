#include "methods.h"


// === CALCULATE METRIC METHOD ===
void B_CalcMetric(int NI, int NJ, std::vector<std::vector<double>> &X, std::vector<std::vector<double>> &Y,
	std::vector<std::vector<std::vector<double>>> &CellCenter, std::vector<std::vector<double>> &CellVolume,
	std::vector<std::vector<std::vector<double>>> &IFaceCenter, std::vector<std::vector<std::vector<double>>> &IFaceVector,
	std::vector<std::vector<std::vector<double>>> &JFaceCenter, std::vector<std::vector<std::vector<double>>> &JFaceVector)
{

	std::vector<double> r(2);
	int i, j, k;


	//=== FACE CENTERS AND FACE VECTORS ===
	//I - DIRECTION

	for (j = 0; j < NJ - 1; ++j)
	{
		for (i = 0; i < NI; ++i)
		{
			r[0] = X[i][j + 1] - X[i][j];  //r = vector from one node to another
			r[1] = Y[i][j + 1] - Y[i][j];
			IFaceVector[i][j][0] = r[1]; //IFaceVector = r rotated on 90 degree
			IFaceVector[i][j][1] = -r[0]; //IFaceVector directed to increasing I - index
			IFaceCenter[i][j][0] = 0.5*(X[i][j] + X[i][j + 1]);
			IFaceCenter[i][j][1] = 0.5*(Y[i][j] + Y[i][j + 1]);

		}
	}

	//J - DIRECTION

	for (j = 0; j < NJ; ++j)
	{
		for (i = 0; i < NI - 1; ++i)
		{
			r[0] = X[i + 1][j] - X[i][j];  //r = vector from one node to another
			r[1] = Y[i + 1][j] - Y[i][j];
			JFaceVector[i][j][0] = -r[1]; //JFaceVector = r rotated on - 90 degree
			JFaceVector[i][j][1] = r[0]; //JFaceVector directed to increasing J - index
			JFaceCenter[i][j][0] = 0.5*(X[i][j] + X[i + 1][j]);
			JFaceCenter[i][j][1] = 0.5*(Y[i][j] + Y[i + 1][j]);
		}
	}


	//=== CELL VOLUMES ===
	double termCV1, termCV2;
	for (j = 0; j < NJ - 1; ++j)
	{
		for (i = 0; i < NI - 1; ++i)
		{

			r[0] = X[i + 1][j + 1] - X[i][j];  //r = vector from one node to another
			r[1] = Y[i + 1][j + 1] - Y[i][j];
			termCV1 = 0.5 * (IFaceVector[i][j][0] * r[0] + IFaceVector[i][j][1] * r[1]);
			termCV2 = 0.5 * (JFaceVector[i][j][0] * r[0] + JFaceVector[i][j][1] * r[1]);
			CellVolume[i][j] = termCV1 + termCV2;
		}
	}

	//=== CELL CENTERS ===
	//FOR INNER CELLS : CENTER OF CONTOUR(sum of FaceCenter*FaceLength / Perimeter)

	double termCC1, termCC2, termCC3, termCC4, termCC5, termCC6;
	for (j = 1; j < NJ; ++j)
	{
		for (i = 1; i < NI; ++i)
		{

			for (k = 0; k < 2; ++k)
			{

				termCC1 = sqrt(IFaceVector[i - 1][j - 1][0] * IFaceVector[i - 1][j - 1][0] + IFaceVector[i - 1][j - 1][1] * IFaceVector[i - 1][j - 1][1]);
				termCC2 = sqrt(IFaceVector[i][j - 1][0] * IFaceVector[i][j - 1][0] + IFaceVector[i][j - 1][1] * IFaceVector[i][j - 1][1]);
				termCC3 = sqrt(JFaceVector[i - 1][j - 1][0] * JFaceVector[i - 1][j - 1][0] + JFaceVector[i - 1][j - 1][1] * JFaceVector[i - 1][j - 1][1]);
				termCC4 = sqrt(JFaceVector[i - 1][j][0] * JFaceVector[i - 1][j][0] + JFaceVector[i - 1][j][1] * JFaceVector[i - 1][j][1]);

				termCC5 = IFaceCenter[i - 1][j - 1][k] * termCC1 + IFaceCenter[i][j - 1][k] * termCC2 + JFaceCenter[i - 1][j - 1][k] * termCC3 + JFaceCenter[i - 1][j][k] * termCC4;
				termCC6 = termCC1 + termCC2 + termCC3 + termCC4;

				CellCenter[i][j][k] = termCC5 / termCC6;

			}
		}

	}


	//FOR DUMMY CELLS ON BOUNDARIES : CELL CENTER = FACE CENTER
	//I - BOUNDARIES---------------------------------------------------- -
	int NBOUND, IBOUND, IOUT;
	for (NBOUND = 0; NBOUND < 2; ++NBOUND)
	{
		if (NBOUND == 0)
		{
			IBOUND = 0; IOUT = 0;
		}
		else
		{
			IBOUND = NI - 1; IOUT = NI;
		}
		for (j = 0; j < NJ - 1; ++j)
			for (k = 0; k < 2; ++k)
			{
				CellCenter[IOUT][j + 1][k] = IFaceCenter[IBOUND][j][k];
			}
	}


	//J - BOUNDARIES---------------------------------------------------- -
	int JBOUND, JOUT;
	for (NBOUND = 0; NBOUND < 2; ++NBOUND)
	{
		if (NBOUND == 0)
		{
			JBOUND = 0; JOUT = 0;
		}
		else
		{
			JBOUND = NJ - 1; JOUT = NJ;
		}
		for (i = 0; i < NI - 1; ++i)
			for (k = 0; k < 2; ++k)
				CellCenter[i + 1][JOUT][k] = JFaceCenter[i][JBOUND][k];
	}

}


// === CALCULATE GRADIENT ===
// non-iterative method 
void B_CalcGradient(int NI, int NJ, std::vector<std::vector<double>> &P, std::vector<std::vector<std::vector<double>>> &Grad,
	std::vector<std::vector<std::vector<double>>> &CellCenter, std::vector<std::vector<double>> &CellVolume,
	std::vector<std::vector<std::vector<double>>> &IFaceCenter, std::vector<std::vector<std::vector<double>>> &IFaceVector,
	std::vector<std::vector<std::vector<double>>> &JFaceCenter, std::vector<std::vector<std::vector<double>>> &JFaceVector)
{

	int i, j, k;

	std::vector<std::vector<int>> NCELL(4);
	NCELL.resize(4);
	for (i = 0; i < NCELL.size(); ++i)
		NCELL[i].resize(2);

	double VOL;

	std::vector<double> RC(2);
	RC.resize(2);

	std::vector<std::vector<double>> RF(4);
	RF.resize(4);
	for (i = 0; i < RF.size(); ++i)
		RF[i].resize(2);

	std::vector<std::vector<double>> SF(4);
	SF.resize(4);
	for (i = 0; i < SF.size(); ++i)
		SF[i].resize(2);

	std::vector<double> RN(2);
	RN.resize(2);

	double	DC, DN, PF;

	int IFace, IN, JN;

	//Temporary arrays

	std::vector<std::vector<double>> tempRFRN(4);
	tempRFRN.resize(4);
	for (i = 0; i < tempRFRN.size(); ++i)
		tempRFRN[i].resize(2);

	std::vector<std::vector<double>> tempRFRC(4);
	tempRFRC.resize(4);
	for (i = 0; i < tempRFRC.size(); ++i)
		tempRFRC[i].resize(2);


	for (i = 1; i < NI; ++i)
	{
		for (j = 1; j < NJ; ++j)
		{

			NCELL[0][0] = i - 1;
			NCELL[0][1] = j;
			NCELL[1][0] = i + 1;
			NCELL[1][1] = j;
			NCELL[2][0] = i;
			NCELL[2][1] = j - 1;
			NCELL[3][0] = i;
			NCELL[3][1] = j + 1;


			for (k = 0; k < 2; ++k)
			{
				RF[0][k] = IFaceCenter[i - 1][j - 1][k];
				RF[1][k] = IFaceCenter[i][j - 1][k];
				RF[2][k] = JFaceCenter[i - 1][j - 1][k];
				RF[3][k] = JFaceCenter[i - 1][j][k];

				SF[0][k] = -IFaceVector[i - 1][j - 1][k];
				SF[1][k] = IFaceVector[i][j - 1][k];
				SF[2][k] = -JFaceVector[i - 1][j - 1][k];
				SF[3][k] = JFaceVector[i - 1][j][k];

				RC[k] = CellCenter[i][j][k];

			}

			VOL = CellVolume[i - 1][j - 1];


			for (IFace = 0; IFace < 4; ++IFace)
			{
				IN = NCELL[IFace][0];
				JN = NCELL[IFace][1];

				for (k = 0; k < 2; ++k)
				{
					RN[k] = CellCenter[IN][JN][k];
					tempRFRC[IFace][k] = RF[IFace][k] - RC[k];
					tempRFRN[IFace][k] = RF[IFace][k] - RN[k];
				}

				DC = sqrt(tempRFRC[IFace][0] * tempRFRC[IFace][0] + tempRFRC[IFace][1] * tempRFRC[IFace][1]);
				DN = sqrt(tempRFRN[IFace][0] * tempRFRN[IFace][0] + tempRFRN[IFace][1] * tempRFRN[IFace][1]);
				PF = (P[i][j] * DN + P[IN][JN] * DC) / (DC + DN);
				Grad[i][j][0] = Grad[i][j][0] + PF * SF[IFace][0];
				Grad[i][j][1] = Grad[i][j][1] + PF * SF[IFace][1];
			}

			Grad[i][j][0] = Grad[i][j][0] / VOL;
			Grad[i][j][1] = Grad[i][j][1] / VOL;

		}
	}


}

// === CALCULATE GRADIENT ===
// iterative method (actions during one iteration)
void B_CalcGradient_Iter(int NI, int NJ, std::vector<std::vector<double>> &P, std::vector<std::vector<std::vector<double>>> &Grad,
	std::vector<std::vector<std::vector<double>>> &CellCenter, std::vector<std::vector<double>> &CellVolume,
	std::vector<std::vector<std::vector<double>>> &IFaceCenter, std::vector<std::vector<std::vector<double>>> &IFaceVector,
	std::vector<std::vector<std::vector<double>>> &JFaceCenter, std::vector<std::vector<std::vector<double>>> &JFaceVector)
{

	int i, j, k;

	std::vector<std::vector<int>> NCELL(4);
	NCELL.resize(4);
	for (i = 0; i < NCELL.size(); ++i)
		NCELL[i].resize(2);

	double VOL;

	std::vector<double> RC(2);
	RC.resize(2);

	std::vector<std::vector<double>> RF(4);
	RF.resize(4);
	for (i = 0; i < RF.size(); ++i)
		RF[i].resize(2);

	std::vector<std::vector<double>> SF(4);
	SF.resize(4);
	for (i = 0; i < SF.size(); ++i)
		SF[i].resize(2);

	std::vector<double> RN(2);
	RN.resize(2);

	std::vector<double> RM(2);
	RM.resize(2);

	std::vector<double> GRM(2);
	GRM.resize(2);

	std::vector<double> GP(2);
	GP.resize(2);

	double	DC, DN, PF, PM;

	int IFace, IN, JN;

	//Temporary arrays

	std::vector<std::vector<double>> tempRFRN(4);
	tempRFRN.resize(4);
	for (i = 0; i < tempRFRN.size(); ++i)
		tempRFRN[i].resize(2);

	std::vector<std::vector<double>> tempRFRC(4);
	tempRFRC.resize(4);
	for (i = 0; i < tempRFRC.size(); ++i)
		tempRFRC[i].resize(2);

	std::vector<std::vector<double>> tempRFRM(4);
	tempRFRM.resize(4);
	for (i = 0; i < tempRFRM.size(); ++i)
		tempRFRM[i].resize(2);


	for (i = 1; i < NI; ++i)
	{
		for (j = 1; j < NJ; ++j)
		{

			NCELL[0][0] = i - 1;
			NCELL[0][1] = j;
			NCELL[1][0] = i + 1;
			NCELL[1][1] = j;
			NCELL[2][0] = i;
			NCELL[2][1] = j - 1;
			NCELL[3][0] = i;
			NCELL[3][1] = j + 1;


			for (k = 0; k < 2; ++k)
			{
				RF[0][k] = IFaceCenter[i - 1][j - 1][k];
				RF[1][k] = IFaceCenter[i][j - 1][k];
				RF[2][k] = JFaceCenter[i - 1][j - 1][k];
				RF[3][k] = JFaceCenter[i - 1][j][k];

				SF[0][k] = -IFaceVector[i - 1][j - 1][k];
				SF[1][k] = IFaceVector[i][j - 1][k];
				SF[2][k] = -JFaceVector[i - 1][j - 1][k];
				SF[3][k] = JFaceVector[i - 1][j][k];

				RC[k] = CellCenter[i][j][k];

			}

			VOL = CellVolume[i - 1][j - 1];

			GP[0] = 0.0;
			GP[1] = 0.0;

			for (IFace = 0; IFace < 4; ++IFace)
			{
				IN = NCELL[IFace][0];
				JN = NCELL[IFace][1];

				for (k = 0; k < 2; ++k)
				{
					RN[k] = CellCenter[IN][JN][k];
					tempRFRC[IFace][k] = RF[IFace][k] - RC[k];
					tempRFRN[IFace][k] = RF[IFace][k] - RN[k];
				}

				DC = sqrt(tempRFRC[IFace][0] * tempRFRC[IFace][0] + tempRFRC[IFace][1] * tempRFRC[IFace][1]);
				DN = sqrt(tempRFRN[IFace][0] * tempRFRN[IFace][0] + tempRFRN[IFace][1] * tempRFRN[IFace][1]);
				PM = (P[i][j] * DN + P[IN][JN] * DC) / (DC + DN);

				RM[0] = (CellCenter[i][j][0] * DN + CellCenter[IN][JN][0] * DC) / (DC + DN);
				RM[1] = (CellCenter[i][j][1] * DN + CellCenter[IN][JN][1] * DC) / (DC + DN);
				GRM[0] = (Grad[i][j][0] * DN + Grad[IN][JN][0] * DC) / (DC + DN);
				GRM[1] = (Grad[i][j][1] * DN + Grad[IN][JN][1] * DC) / (DC + DN);

				for (k = 0; k < 2; ++k)
				{
					tempRFRM[IFace][k] = RF[IFace][k] - RM[k];
				}
				PF = PM + tempRFRM[IFace][0] * GRM[0] + tempRFRM[IFace][1] * GRM[1];

				GP[0] = GP[0] + PF * SF[IFace][0];
				GP[1] = GP[1] + PF * SF[IFace][1];
			}

			Grad[i][j][0] = GP[0] / VOL;
			Grad[i][j][1] = GP[1] / VOL;

		}
	}


}



void B_CalcDiv(int NI, int NJ, std::vector<std::vector<std::vector<double>>> &V, std::vector<std::vector<double>> &Div,
	std::vector<std::vector<std::vector<double>>> &CellCenter, std::vector<std::vector<double>> &CellVolume,
	std::vector<std::vector<std::vector<double>>> &IFaceCenter, std::vector<std::vector<std::vector<double>>> &IFaceVector,
	std::vector<std::vector<std::vector<double>>> &JFaceCenter, std::vector<std::vector<std::vector<double>>> &JFaceVector,
	std::vector<std::vector<double>> &P, std::vector<std::vector<std::vector<double>>> &Grad, int mode)
{

	int i, j, k;

	std::vector<std::vector<int>> NCELL(4);
	NCELL.resize(4);
	for (i = 0; i < NCELL.size(); ++i)
		NCELL[i].resize(2);

	double VOL;

	std::vector<double> RC(2);
	RC.resize(2);

	std::vector<std::vector<double>> RF(4);
	RF.resize(4);
	for (i = 0; i < RF.size(); ++i)
		RF[i].resize(2);

	std::vector<std::vector<double>> SF(4);
	SF.resize(4);
	for (i = 0; i < SF.size(); ++i)
		SF[i].resize(2);

	std::vector<double> RN(2);
	RN.resize(2);

	std::vector<double> VF(2);
	VF.resize(2);

	double	DC, DN, PF, GC, GB, GN, PN;

	int IFace, IN, JN;

	//Temporary arrays

	std::vector<std::vector<double>> tempRFRN(4);
	tempRFRN.resize(4);
	for (i = 0; i < tempRFRN.size(); ++i)
		tempRFRN[i].resize(2);

	std::vector<std::vector<double>> tempRFRC(4);
	tempRFRC.resize(4);
	for (i = 0; i < tempRFRC.size(); ++i)
		tempRFRC[i].resize(2);


	for (i = 1; i < NI; ++i)
	{
		for (j = 1; j < NJ; ++j)
		{

			NCELL[0][0] = i - 1;
			NCELL[0][1] = j;
			NCELL[1][0] = i + 1;
			NCELL[1][1] = j;
			NCELL[2][0] = i;
			NCELL[2][1] = j - 1;
			NCELL[3][0] = i;
			NCELL[3][1] = j + 1;


			for (k = 0; k < 2; ++k)
			{
				RF[0][k] = IFaceCenter[i - 1][j - 1][k];
				RF[1][k] = IFaceCenter[i][j - 1][k];
				RF[2][k] = JFaceCenter[i - 1][j - 1][k];
				RF[3][k] = JFaceCenter[i - 1][j][k];

				SF[0][k] = -IFaceVector[i - 1][j - 1][k];
				SF[1][k] = IFaceVector[i][j - 1][k];
				SF[2][k] = -JFaceVector[i - 1][j - 1][k];
				SF[3][k] = JFaceVector[i - 1][j][k];

				RC[k] = CellCenter[i][j][k];

			}

			VOL = CellVolume[i - 1][j - 1];

			for (IFace = 0; IFace < 4; ++IFace)
			{
				IN = NCELL[IFace][0];
				JN = NCELL[IFace][1];

				for (k = 0; k < 2; ++k)
				{
					RN[k] = CellCenter[IN][JN][k];
					tempRFRC[IFace][k] = RF[IFace][k] - RC[k];
					tempRFRN[IFace][k] = RF[IFace][k] - RN[k];
				}

				DC = sqrt(tempRFRC[IFace][0] * tempRFRC[IFace][0] + tempRFRC[IFace][1] * tempRFRC[IFace][1]);
				DN = sqrt(tempRFRN[IFace][0] * tempRFRN[IFace][0] + tempRFRN[IFace][1] * tempRFRN[IFace][1]);
				VF[0] = (V[i][j][0] * DN + V[IN][JN][0] * DC) / (DC + DN);
				VF[1] = (V[i][j][1] * DN + V[IN][JN][1] * DC) / (DC + DN);

				if (mode == 0)
					Div[i][j] = Div[i][j] + VF[0] * SF[IFace][0] + VF[1] * SF[IFace][1];

				if (mode == 1)
				{
					PF = (P[i][j] * DN + P[IN][JN] * DC) / (DC + DN);
					Div[i][j] = Div[i][j] + PF * SF[IFace][0] * VF[0] + PF * SF[IFace][1] * VF[1];
				}

				if (mode == 2)
				{
					if ((SF[IFace][0] * VF[0] + SF[IFace][1] * VF[1]) > 0.0)
						PF = P[i][j];
					else
					{
						PF = P[IN][JN];
						if (DN <= 1e-6) PF = 2 * P[IN][JN] - P[i][j];
					}
					Div[i][j] = Div[i][j] + PF * SF[IFace][0] * VF[0] + PF * SF[IFace][1] * VF[1];
				}

				if (mode == 3)
				{
					if ((SF[IFace][0] * VF[0] + SF[IFace][1] * VF[1]) > 0.0)
						PF = P[i][j] + (RF[IFace][0] - CellCenter[i][j][0]) * Grad[i][j][0] + (RF[IFace][1] - CellCenter[i][j][1]) * Grad[i][j][1];
					else
					{
						PF = P[IN][JN] + (RF[IFace][0] - CellCenter[IN][JN][0]) * Grad[IN][JN][0] + (RF[IFace][1] - CellCenter[IN][JN][1]) * Grad[IN][JN][1];
						if (DN <= 1e-6)
						{
							PN = 2 * P[IN][JN] - P[i][j];
							GC = Grad[i][j][0] * (CellCenter[i][j][0] - RF[IFace][0]) + Grad[i][j][1] * (CellCenter[i][j][1] - RF[IFace][1]);
							GB = P[i][j] - P[IN][JN];
							GN = 4 * GB - 3 * GC;
							PF = PN + GN;
						}
					}
					Div[i][j] = Div[i][j] + PF * SF[IFace][0] * VF[0] + PF * SF[IFace][1] * VF[1];
				}


			}

			Div[i][j] = Div[i][j] / VOL;


		}
	}



}


void B_CalcRotZ(int NI, int NJ, std::vector<std::vector<std::vector<double>>> &V, std::vector<std::vector<double>> &RotZ,
	std::vector<std::vector<std::vector<double>>> &CellCenter, std::vector<std::vector<double>> &CellVolume,
	std::vector<std::vector<std::vector<double>>> &IFaceCenter, std::vector<std::vector<std::vector<double>>> &IFaceVector,
	std::vector<std::vector<std::vector<double>>> &JFaceCenter, std::vector<std::vector<std::vector<double>>> &JFaceVector)
{

	int i, j, k;

	std::vector<std::vector<int>> NCELL(4);
	NCELL.resize(4);
	for (i = 0; i < NCELL.size(); ++i)
		NCELL[i].resize(2);

	double VOL;

	std::vector<double> RC(2);
	RC.resize(2);

	std::vector<std::vector<double>> RF(4);
	RF.resize(4);
	for (i = 0; i < RF.size(); ++i)
		RF[i].resize(2);

	std::vector<std::vector<double>> SF(4);
	SF.resize(4);
	for (i = 0; i < SF.size(); ++i)
		SF[i].resize(2);

	std::vector<double> RN(2);
	RN.resize(2);

	std::vector<double> VF(2);
	VF.resize(2);

	double	DC, DN;

	int IFace, IN, JN;

	//Temporary arrays

	std::vector<std::vector<double>> tempRFRN(4);
	tempRFRN.resize(4);
	for (i = 0; i < tempRFRN.size(); ++i)
		tempRFRN[i].resize(2);

	std::vector<std::vector<double>> tempRFRC(4);
	tempRFRC.resize(4);
	for (i = 0; i < tempRFRC.size(); ++i)
		tempRFRC[i].resize(2);


	for (i = 1; i < NI; ++i)
	{
		for (j = 1; j < NJ; ++j)
		{

			NCELL[0][0] = i - 1;
			NCELL[0][1] = j;
			NCELL[1][0] = i + 1;
			NCELL[1][1] = j;
			NCELL[2][0] = i;
			NCELL[2][1] = j - 1;
			NCELL[3][0] = i;
			NCELL[3][1] = j + 1;


			for (k = 0; k < 2; ++k)
			{
				RF[0][k] = IFaceCenter[i - 1][j - 1][k];
				RF[1][k] = IFaceCenter[i][j - 1][k];
				RF[2][k] = JFaceCenter[i - 1][j - 1][k];
				RF[3][k] = JFaceCenter[i - 1][j][k];

				SF[0][k] = -IFaceVector[i - 1][j - 1][k];
				SF[1][k] = IFaceVector[i][j - 1][k];
				SF[2][k] = -JFaceVector[i - 1][j - 1][k];
				SF[3][k] = JFaceVector[i - 1][j][k];

				RC[k] = CellCenter[i][j][k];

			}

			VOL = CellVolume[i - 1][j - 1];

			for (IFace = 0; IFace < 4; ++IFace)
			{
				IN = NCELL[IFace][0];
				JN = NCELL[IFace][1];

				for (k = 0; k < 2; ++k)
				{
					RN[k] = CellCenter[IN][JN][k];
					tempRFRC[IFace][k] = RF[IFace][k] - RC[k];
					tempRFRN[IFace][k] = RF[IFace][k] - RN[k];
				}

				DC = sqrt(tempRFRC[IFace][0] * tempRFRC[IFace][0] + tempRFRC[IFace][1] * tempRFRC[IFace][1]);
				DN = sqrt(tempRFRN[IFace][0] * tempRFRN[IFace][0] + tempRFRN[IFace][1] * tempRFRN[IFace][1]);
				VF[0] = (V[i][j][0] * DN + V[IN][JN][0] * DC) / (DC + DN);
				VF[1] = (V[i][j][1] * DN + V[IN][JN][1] * DC) / (DC + DN);
				RotZ[i][j] = RotZ[i][j] + (VF[1] * SF[IFace][0] - VF[0] * SF[IFace][1]);
			}

			RotZ[i][j] = RotZ[i][j] / VOL;


		}
	}


}


void B_CalcLap(int NI, int NJ, std::vector<std::vector<std::vector<double>>> &V, std::vector<std::vector<double>> &Lap,
	std::vector<std::vector<std::vector<double>>> &CellCenter, std::vector<std::vector<double>> &CellVolume,
	std::vector<std::vector<std::vector<double>>> &IFaceCenter, std::vector<std::vector<std::vector<double>>> &IFaceVector,
	std::vector<std::vector<std::vector<double>>> &JFaceCenter, std::vector<std::vector<std::vector<double>>> &JFaceVector,
	std::vector<std::vector<double>> &P, std::vector<std::vector<std::vector<double>>> &Grad)
{

	int i, j, k;

	std::vector<std::vector<int>> NCELL(4);
	NCELL.resize(4);
	for (i = 0; i < NCELL.size(); ++i)
		NCELL[i].resize(2);

	double VOL;

	std::vector<double> RC(2);
	RC.resize(2);

	std::vector<std::vector<double>> RF(4);
	RF.resize(4);
	for (i = 0; i < RF.size(); ++i)
		RF[i].resize(2);

	std::vector<std::vector<double>> SF(4);
	SF.resize(4);
	for (i = 0; i < SF.size(); ++i)
		SF[i].resize(2);

	std::vector<double> RN(2);
	RN.resize(2);

	std::vector<double> NF(2);
	NF.resize(2);

	std::vector<double> RNC(2);
	RNC.resize(2);

	std::vector<double> GF(2);
	GF.resize(2);

	double	DC, DN, DNC, dpdn, dpdn_c;

	int IFace, IN, JN;

	//Temporary arrays

	std::vector<std::vector<double>> tempRFRN(4);
	tempRFRN.resize(4);
	for (i = 0; i < tempRFRN.size(); ++i)
		tempRFRN[i].resize(2);

	std::vector<std::vector<double>> tempRFRC(4);
	tempRFRC.resize(4);
	for (i = 0; i < tempRFRC.size(); ++i)
		tempRFRC[i].resize(2);


	for (i = 1; i < NI; ++i)
	{
		for (j = 1; j < NJ; ++j)
		{

			NCELL[0][0] = i - 1;
			NCELL[0][1] = j;
			NCELL[1][0] = i + 1;
			NCELL[1][1] = j;
			NCELL[2][0] = i;
			NCELL[2][1] = j - 1;
			NCELL[3][0] = i;
			NCELL[3][1] = j + 1;


			for (k = 0; k < 2; ++k)
			{
				RF[0][k] = IFaceCenter[i - 1][j - 1][k];
				RF[1][k] = IFaceCenter[i][j - 1][k];
				RF[2][k] = JFaceCenter[i - 1][j - 1][k];
				RF[3][k] = JFaceCenter[i - 1][j][k];

				SF[0][k] = -IFaceVector[i - 1][j - 1][k];
				SF[1][k] = IFaceVector[i][j - 1][k];
				SF[2][k] = -JFaceVector[i - 1][j - 1][k];
				SF[3][k] = JFaceVector[i - 1][j][k];

				RC[k] = CellCenter[i][j][k];

			}

			VOL = CellVolume[i - 1][j - 1];

			for (IFace = 0; IFace < 4; ++IFace)
			{
				IN = NCELL[IFace][0];
				JN = NCELL[IFace][1];

				for (k = 0; k < 2; ++k)
				{
					RN[k] = CellCenter[IN][JN][k];
					tempRFRC[IFace][k] = RF[IFace][k] - RC[k];
					tempRFRN[IFace][k] = RF[IFace][k] - RN[k];
				}

				DC = sqrt(tempRFRC[IFace][0] * tempRFRC[IFace][0] + tempRFRC[IFace][1] * tempRFRC[IFace][1]);
				DN = sqrt(tempRFRN[IFace][0] * tempRFRN[IFace][0] + tempRFRN[IFace][1] * tempRFRN[IFace][1]);
				DNC = sqrt((RN[0] - RC[0])*(RN[0] - RC[0]) + (RN[1] - RC[1])*(RN[1] - RC[1]));
				NF[0] = SF[IFace][0] / sqrt(SF[IFace][0] * SF[IFace][0] + SF[IFace][1] * SF[IFace][1]);
				NF[1] = SF[IFace][1] / sqrt(SF[IFace][0] * SF[IFace][0] + SF[IFace][1] * SF[IFace][1]);

				//for skew correction
				RNC[0] = (RN[0] - RC[0]) / DNC;
				RNC[1] = (RN[1] - RC[1]) / DNC;
				GF[0] = (Grad[i][j][0] * DN + Grad[IN][JN][0] * DC) / (DC + DN);
				GF[1] = (Grad[i][j][1] * DN + Grad[IN][JN][1] * DC) / (DC + DN);

				dpdn = (P[IN][JN] - P[i][j]) / DNC;

				if (DN < 1e-5)
				{
					dpdn_c = Grad[i][j][0] * NF[0] + Grad[i][j][1] * NF[1];
					dpdn = 5.0 / 3.0 * dpdn - 2.0 / 3.0 * dpdn_c; // 2 order
					GF[0] = Grad[i][j][0];
					GF[1] = Grad[i][j][1];
				}

				//skew correction
				dpdn = dpdn + (NF[0] - RNC[0]) * GF[0] + (NF[1] - RNC[1]) * GF[1];

				Lap[i][j] = Lap[i][j] + dpdn * sqrt(SF[IFace][0] * SF[IFace][0] + SF[IFace][1] * SF[IFace][1]);
			}

			Lap[i][j] = Lap[i][j] / VOL;


		}
	}


}


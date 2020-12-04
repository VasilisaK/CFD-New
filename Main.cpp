#include "methods.h"

void IsOpen(std::ifstream& oi);
void IsOpen(std::ofstream& of);
void B_OutputFields(std::ostream &OutputFile, int NI, int NJ, const std::vector<std::vector<double>> &X, const std::vector<std::vector<double>> &Y,
	const std::vector<std::vector<std::vector<double>>> &V, const std::vector<std::vector<double>> &P, const std::vector<std::vector<double>> &T,
	const std::vector<std::vector<std::vector<double>>> &Grad, 
	const std::vector<std::vector<double>> &Div, 
	const std::vector<std::vector<double>> &RotZ, 
	const std::vector<std::vector<double>> &Lap,
	const std::vector<std::vector<double>> &Terr);


int main()
{
	std::ifstream InputFile; //input unit
	std::ofstream OutputFile; //output unit
	char MeshFileName[30]; //name of file with computational mesh
	char InputFileName[] = "input.txt";//names of input and output files
	char OutputFileName[] = "data.plt";

	std::ifstream SolFile; //solution input unit
	char SolFileName[30];

	//===  READ INPUT FILE ===
	std::cout << "Read input file: " << InputFileName << "\n";
	InputFile.open(InputFileName);
	IsOpen(InputFile);

	InputFile >> MeshFileName;  //read name of file with computational mesh
	InputFile >> SolFileName;
	InputFile.close();
	//std::cout << MeshFileName << "\n";

	//===   READ NODES NUMBER(NI, NJ) FROM FILE WITH MESH ===
	std::cout << "Read nodes number from file: " << MeshFileName << "\n";
	std::ifstream MeshFile;
	MeshFile.open(MeshFileName);
	IsOpen(MeshFile);

//	int NI, NJ; //number of nodes
//	MeshFile >> NI >> NJ;
//	std::cout << "NI, NJ = " << NI << " " << NJ << "\n";

	int NI, NJ, NK; //number of nodes
	MeshFile >> NI >> NJ >> NK;
	std::cout << "NI, NJ = " << NI << " " << NJ << "\n";

	//=== ALLOCATE ALL ARRAYS ===
	std::cout << "Allocate arrays\n";
	int i;
	int j;

	//=== scalar arrays ===

	std::vector<std::vector<double>> X(NI);//mesh nodes X - coordinates
	X.resize(NI);
	for (i = 0; i < X.size(); ++i)
		X[i].resize(NJ);

	std::vector<std::vector<double>> Y(NI);//mesh nodes Y - coordinates
	Y.resize(NI);
	for (i = 0; i < Y.size(); ++i)
		Y[i].resize(NJ);

	std::vector<std::vector<double>> CellVolume(NI - 1);//Cell Volumes
	CellVolume.resize(NI - 1);
	for (i = 0; i < CellVolume.size(); ++i)
		CellVolume[i].resize(NJ - 1);

	std::vector<std::vector<double>> P(NI + 1);//Pressure
	P.resize(NI + 1);
	for (i = 0; i < P.size(); ++i)
		P[i].resize(NJ + 1);

	std::vector<std::vector<double>> T(NI + 1);//Temperature
	T.resize(NI + 1);
	for (i = 0; i < T.size(); ++i)
		T[i].resize(NJ + 1);

	std::vector<std::vector<double>> Tnew(NI + 1);//Temperature new field
	Tnew.resize(NI + 1);
	for (i = 0; i < Tnew.size(); ++i)
		Tnew[i].resize(NJ + 1);

	std::vector<std::vector<double>> Terr(NI + 1);//Temperature error
	Terr.resize(NI + 1);
	for (i = 0; i < Terr.size(); ++i)
		Terr[i].resize(NJ + 1);

	std::vector<std::vector<double>> Div(NI + 1);//Divergence
	Div.resize(NI + 1);
	for (i = 0; i < Div.size(); ++i)
		Div[i].resize(NJ + 1);

	std::vector<std::vector<double>> DivExact(NI + 1);//Divergence Exact
	DivExact.resize(NI + 1);
	for (i = 0; i < DivExact.size(); ++i)
		DivExact[i].resize(NJ + 1);

	std::vector<std::vector<double>> DivError(NI + 1);//Divergence Error
	DivError.resize(NI + 1);
	for (i = 0; i < DivError.size(); ++i)
		DivError[i].resize(NJ + 1);

	std::vector<std::vector<double>> RotZ(NI + 1);//Rotor-Z
	RotZ.resize(NI + 1);
	for (i = 0; i < RotZ.size(); ++i)
		RotZ[i].resize(NJ + 1);

	std::vector<std::vector<double>> RotZExact(NI + 1);//Rotor-Z Exact
	RotZExact.resize(NI + 1);
	for (i = 0; i < RotZExact.size(); ++i)
		RotZExact[i].resize(NJ + 1);

	std::vector<std::vector<double>> RotZError(NI + 1);//Rotor-Z Error
	RotZError.resize(NI + 1);
	for (i = 0; i < RotZError.size(); ++i)
		RotZError[i].resize(NJ + 1);

	std::vector<std::vector<double>> Lap(NI + 1);//Laplacian
	Lap.resize(NI + 1);
	for (i = 0; i < Lap.size(); ++i)
		Lap[i].resize(NJ + 1);

	std::vector<std::vector<double>> LapExact(NI + 1);//Laplacian Exact
	LapExact.resize(NI + 1);
	for (i = 0; i < LapExact.size(); ++i)
		LapExact[i].resize(NJ + 1);

	std::vector<std::vector<double>> LapError(NI + 1);//Laplacian Error
	LapError.resize(NI + 1);
	for (i = 0; i < LapError.size(); ++i)
		LapError[i].resize(NJ + 1);

	//=== vector arrays ===

	std::vector<std::vector<std::vector<double>>> CellCenter(NI + 1);//Cell Centers
	CellCenter.resize(NI + 1);
	for (i = 0; i < CellCenter.size(); ++i)
		CellCenter[i].resize(NJ + 1);
	for (i = 0; i < CellCenter.size(); ++i)
		for (j = 0; j < CellCenter[0].size(); ++j)
			CellCenter[i][j].resize(2);

	std::vector<std::vector<std::vector<double>>> IFaceCenter(NI);//Face Centers for I - faces
	IFaceCenter.resize(NI);
	for (i = 0; i < IFaceCenter.size(); ++i)
		IFaceCenter[i].resize(NJ - 1);
	for (i = 0; i < IFaceCenter.size(); ++i)
		for (j = 0; j < IFaceCenter[0].size(); ++j)
			IFaceCenter[i][j].resize(2);

	std::vector<std::vector<std::vector<double>>> IFaceVector(NI);//Face Vectors for I - faces
	IFaceVector.resize(NI);
	for (i = 0; i < IFaceVector.size(); ++i)
		IFaceVector[i].resize(NJ - 1);
	for (i = 0; i < IFaceVector.size(); ++i)
		for (j = 0; j < IFaceVector[0].size(); ++j)
			IFaceVector[i][j].resize(2);

	std::vector<std::vector<std::vector<double>>> JFaceCenter(NI - 1);//Face Centers for J - faces
	JFaceCenter.resize(NI - 1);
	for (i = 0; i < JFaceCenter.size(); ++i)
		JFaceCenter[i].resize(NJ);
	for (i = 0; i < JFaceCenter.size(); ++i)
		for (j = 0; j < JFaceCenter[0].size(); ++j)
			JFaceCenter[i][j].resize(2);

	std::vector<std::vector<std::vector<double>>> JFaceVector(NI - 1);//Face Vectors for I - faces
	JFaceVector.resize(NI - 1);
	for (i = 0; i < JFaceVector.size(); ++i)
		JFaceVector[i].resize(NJ);
	for (i = 0; i < JFaceVector.size(); ++i)
		for (j = 0; j < JFaceVector[0].size(); ++j)
			JFaceVector[i][j].resize(2);

	std::vector<std::vector<std::vector<double>>> V(NI + 1);//Velocity, 2 components
	V.resize(NI + 1);
	for (i = 0; i < V.size(); ++i)
		V[i].resize(NJ + 1);
	for (i = 0; i < V.size(); ++i)
		for (j = 0; j < V[0].size(); ++j)
			V[i][j].resize(2);

	std::vector<std::vector<std::vector<double>>> Grad(NI + 1);//Pressure Gradient, 2 components
	Grad.resize(NI + 1);
	for (i = 0; i < Grad.size(); ++i)
		Grad[i].resize(NJ + 1);
	for (i = 0; i < Grad.size(); ++i)
		for (j = 0; j < Grad[0].size(); ++j)
			Grad[i][j].resize(2);

	std::vector<std::vector<std::vector<double>>> GradError(NI + 1);//Pressure Gradient Error, 2 components
	GradError.resize(NI + 1);
	for (i = 0; i < GradError.size(); ++i)
		GradError[i].resize(NJ + 1);
	for (i = 0; i < GradError.size(); ++i)
		for (j = 0; j < GradError[0].size(); ++j)
			GradError[i][j].resize(2);

	std::vector<std::vector<std::vector<double>>> GradExact(NI + 1);//Pressure Gradient Error, 2 components
	GradExact.resize(NI + 1);
	for (i = 0; i < GradExact.size(); ++i)
		GradExact[i].resize(NJ + 1);
	for (i = 0; i < GradExact.size(); ++i)
		for (j = 0; j < GradExact[0].size(); ++j)
			GradExact[i][j].resize(2);



	//===  READ GRID ===
	std::cout << "Read mesh from file: " << MeshFileName << "\n";
	double tempZ;
	
	for (j = 0; j < NJ; ++j)
		for (i = 0; i < NI; ++i)		
		{
		//	MeshFile >> X[i][j] >> Y[i][j];
			MeshFile >> X[i][j] >> Y[i][j] >> tempZ;
		}
	MeshFile.close();

	//=== CALCULATE METRIC ===
	std::cout << "Calculate metric\n";
	B_CalcMetric(NI, NJ, X, Y, CellCenter, CellVolume, IFaceCenter, IFaceVector, JFaceCenter, JFaceVector);


	//=== INITIATE FIELDS ===
	std::cout << "Initiate fields\n";
	int k = 1;// maximum number of components in the 3rd direction

//	std::ifstream SolFile; //solution input unit
//	char SolFileName[] = "Re100-41.plt";
	std::string xxx;
	double Dump;

	SolFile.open(SolFileName);
	IsOpen(SolFile);

	std::getline(SolFile, xxx);
	std::getline(SolFile, xxx);
	
	for (j = 0; j < NJ + 1; ++j)
		for (i = 0; i < NI + 1; ++i)

		{
			SolFile >> Dump >> Dump >> V[i][j][0] >> V[i][j][1] >> Dump >> P[i][j] >> Dump >> Dump;
		}
	SolFile.close();



	//=== SOLVING EQUATION ===

	double CFL; //Courant number
	double dt; //pseudo-time step
	double min; //minimum mesh step size (mesh 41*41)
	int NT; // number of time steps
	int t; //counter
	double Re, Pr; //Re = U*L/nyu, U = 0.01 m/s, L = 0.01 m, nyu = 1e-6 m*m/s
					// Pr is for water at 20*C
	CFL = 0.1;
	min = 0.00025;
	dt = CFL * min;
	Re = 100.0;
	Pr = 6.2;
	NT = 10;

		int NIter = 20; //for gradient calculation
		int it;
		double max1, max2; //maximum X&Y error for gradient

	//initial field
	for (j = 0; j < NJ + 1; ++j)
	for (i = 0; i < NI + 1; ++i)
	{
		T[i][j] = 0.0;	
	}


	//pseudo-time cycle

	for (t = 0; t < NT; ++t)
	{

		//gradient calculation for other operators

		//Number of iterations    
		//int NIter = 1;  //for non-iterative method

		if (NIter == 1)
		{
			B_CalcGradient(NI, NJ, T, Grad, CellCenter, CellVolume, IFaceCenter, IFaceVector, JFaceCenter, JFaceVector);
		}
		else
		{
			for (it = 0; it < NIter; it++)
			{
				B_CalcGradient_Iter(NI, NJ, T, Grad, CellCenter, CellVolume, IFaceCenter, IFaceVector, JFaceCenter, JFaceVector);
			}

		}

		//divergence calculation for convection term
		//divergence calculation method (0, 1, 2 or 3)
		int mode = 3;
		B_CalcDiv(NI, NJ, V, Div, CellCenter, CellVolume, IFaceCenter, IFaceVector, JFaceCenter, JFaceVector, T, Grad, mode);

		//laplacian calculation for diffusion term
		B_CalcLap(NI, NJ, V, Lap, CellCenter, CellVolume, IFaceCenter, IFaceVector, JFaceCenter, JFaceVector, T, Grad);

		//solving the equation
		for (j = 1; j < NJ; ++j)
		{
			for (i = 1; i < NI; ++i)
			{
				Tnew[i][j] = T[i][j] - dt * ( Div[i][j] - (1.0/(Re*Pr)) * Lap[i][j] ); //new temperature field calculation
			}
		}

		//boundary conditions
		for (i = 0; i < NI + 1; ++i)
		{
			Tnew[i][0] = T[i][1]; //lower
			Tnew[i][NJ] = T[i][NJ - 1]; //upper
		}

		for (j = 0; j < NJ + 1; ++j)
		{
			Tnew[0][j] = 1.0; //left
			Tnew[NI][j] = 0.0; //right
		}

		//residual calculation
		for (j = 1; j < NJ; ++j)
		{
			for (i = 1; i < NI; ++i)
			{
				Terr[i][j] = (Tnew[i][j] - T[i][j]) / dt;
			}
		}

		//maximum error calculation
		max1 = 0.0;
		for (j = 1; j < NJ; ++j)
		{
			for (i = 1; i < NI; ++i)
			{
				if (Terr[i][j] > max1)
					max1 = Terr[i][j];
			}
		}
		
		//temperature field renewal
		for (j = 0; j < NJ + 1; ++j)
		{
			for (i = 0; i < NI + 1; ++i)
			{
				T[i][j] = Tnew[i][j];
			}
		}


		std::cout << "Iter = " << t + 1 << ", Max Error = " << max1 << "\n";
	}


	//=== OUTPUT FIELDS ===
	std::cout << "Output fields to file:" << OutputFileName << "\n";
	OutputFile.open(OutputFileName);
	IsOpen(OutputFile);

	B_OutputFields(OutputFile, NI, NJ, X, Y, V, P, T, Grad, Div, RotZ, Lap, Terr);
	OutputFile.close();

	system("pause");
	return 0;
}

double input(std::ifstream& inFile) {
	double x;
	inFile >> x;
	while (inFile.get() != '\n' && !inFile.eof()) {
		continue;
	}
	return x;
}

void IsOpen(std::ifstream& oi)
{
	if (!oi.is_open())
	{
		exit(EXIT_FAILURE);
	}
}

void IsOpen(std::ofstream& of)
{
	if (!of.is_open())
	{
		exit(EXIT_FAILURE);
	}
}

// === WRITING OUTPUT FIELDS ===
void B_OutputFields(std::ostream &OutputFile, int NI, int NJ, const std::vector<std::vector<double>> &X, const std::vector<std::vector<double>> &Y,
	const std::vector<std::vector<std::vector<double>>> &V, const std::vector<std::vector<double>> &P, const std::vector<std::vector<double>> &T,
	const std::vector<std::vector<std::vector<double>>> &Grad, 
	const std::vector<std::vector<double>> &Div, 
	const std::vector<std::vector<double>> &RotZ, 
	const std::vector<std::vector<double>> &Lap,
	const std::vector<std::vector<double>> &Terr)
{
	OutputFile << "VARIABLES = \"X\", \"Y\", \"Velocity-X\", \"Velocity-Y\", \"P\", \"T\", \"Grad-X\", \"Grad-Y\", \"DivP\", \"ZRot\", \"Lap\", \"Terr\"\n";
	OutputFile << "ZONE I=" << NI << " J=" << NJ << " DATAPACKING=BLOCK, VARLOCATION=([3-12]=CELLCENTERED)\n";
	int i, j;
	//writing X and Y
	for (j = 0; j < NJ; ++j)
		for (i = 0; i < NI; ++i)	
		{
			OutputFile << std::setprecision(20) << std::fixed << X[i][j] << "\n";
		}
	
	for (j = 0; j < NJ; ++j)
		for (i = 0; i < NI; ++i)	
		{
			OutputFile << std::setprecision(20) << std::fixed << Y[i][j] << "\n";
		}
	
	//writing Vecocity components
	for (j = 1; j < NJ; ++j)
		for (i = 1; i < NI; ++i)
		{
			OutputFile << std::setprecision(20) << std::fixed << V[i][j][0] << "\n";
		}
	
	for (j = 1; j < NJ; ++j)
		for (i = 1; i < NI; ++i)
		{
			OutputFile << std::setprecision(20) << std::fixed << V[i][j][1] << "\n";
		}
	
	//writing pressure
	for (j = 1; j < NJ; ++j)
		for (i = 1; i < NI; ++i)
		{
			OutputFile << std::setprecision(20) << std::fixed << P[i][j] << "\n";
		}

	//writing temperature
	for (j = 1; j < NJ; ++j)
		for (i = 1; i < NI; ++i)
		{
			OutputFile << std::setprecision(20) << std::fixed << T[i][j] << "\n";
		}

	//writing gradient	
	for (j = 1; j < NJ; ++j)
		for (i = 1; i < NI; ++i)
		{
			OutputFile << std::setprecision(20) << std::fixed << Grad[i][j][0] << "\n";
		}
		
	for (j = 1; j < NJ; ++j)
		for (i = 1; i < NI; ++i)
		{
			OutputFile << std::setprecision(20) << std::fixed << Grad[i][j][1] << "\n";
		}
	

	//writing divergence	
	for (j = 1; j < NJ; ++j)
		for (i = 1; i < NI; ++i)
		{
			OutputFile << std::setprecision(20) << std::fixed << Div[i][j] << "\n";
		}
	
	//writing rotor Z-component	
	for (j = 1; j < NJ; ++j)
		for (i = 1; i < NI; ++i)
		{
			OutputFile << std::setprecision(20) << std::fixed << RotZ[i][j] << "\n";
		}

	//writing laplacian
	for (j = 1; j < NJ; ++j)
		for (i = 1; i < NI; ++i)
		{
			OutputFile << std::setprecision(20) << std::fixed << Lap[i][j] << "\n";
		}

	for (j = 1; j < NJ; ++j)
		for (i = 1; i < NI; ++i)
		{
			OutputFile << std::setprecision(20) << std::fixed << Terr[i][j] << "\n";
		}
	
}


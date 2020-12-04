#pragma once
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <iomanip>
#include <numeric>
#include <cmath>
#include <string>

void B_CalcMetric(int NI, int NJ, std::vector<std::vector<double>> &X, std::vector<std::vector<double>> &Y,
	std::vector<std::vector<std::vector<double>>> &CellCenter, std::vector<std::vector<double>> &CellVolume,
	std::vector<std::vector<std::vector<double>>> &IFaceCenter, std::vector<std::vector<std::vector<double>>> &IFaceVector,
	std::vector<std::vector<std::vector<double>>> &JFaceCenter, std::vector<std::vector<std::vector<double>>> &JFaceVector);

void B_CalcGradient(int NI, int NJ, std::vector<std::vector<double>> &P, std::vector<std::vector<std::vector<double>>> &Grad,
	std::vector<std::vector<std::vector<double>>> &CellCenter, std::vector<std::vector<double>> &CellVolume,
	std::vector<std::vector<std::vector<double>>> &IFaceCenter, std::vector<std::vector<std::vector<double>>> &IFaceVector,
	std::vector<std::vector<std::vector<double>>> &JFaceCenter, std::vector<std::vector<std::vector<double>>> &JFaceVector);

void B_CalcGradient_Iter(int NI, int NJ, std::vector<std::vector<double>> &P, std::vector<std::vector<std::vector<double>>> &Grad,
	std::vector<std::vector<std::vector<double>>> &CellCenter, std::vector<std::vector<double>> &CellVolume,
	std::vector<std::vector<std::vector<double>>> &IFaceCenter, std::vector<std::vector<std::vector<double>>> &IFaceVector,
	std::vector<std::vector<std::vector<double>>> &JFaceCenter, std::vector<std::vector<std::vector<double>>> &JFaceVector);

void B_CalcDiv(int NI, int NJ, std::vector<std::vector<std::vector<double>>> &V, std::vector<std::vector<double>> &Div,
	std::vector<std::vector<std::vector<double>>> &CellCenter, std::vector<std::vector<double>> &CellVolume,
	std::vector<std::vector<std::vector<double>>> &IFaceCenter, std::vector<std::vector<std::vector<double>>> &IFaceVector,
	std::vector<std::vector<std::vector<double>>> &JFaceCenter, std::vector<std::vector<std::vector<double>>> &JFaceVector,
	std::vector<std::vector<double>> &P, std::vector<std::vector<std::vector<double>>> &Grad, int mode);

void B_CalcRotZ(int NI, int NJ, std::vector<std::vector<std::vector<double>>> &V, std::vector<std::vector<double>> &RotZ,
	std::vector<std::vector<std::vector<double>>> &CellCenter, std::vector<std::vector<double>> &CellVolume,
	std::vector<std::vector<std::vector<double>>> &IFaceCenter, std::vector<std::vector<std::vector<double>>> &IFaceVector,
	std::vector<std::vector<std::vector<double>>> &JFaceCenter, std::vector<std::vector<std::vector<double>>> &JFaceVector);

void B_CalcLap(int NI, int NJ, std::vector<std::vector<std::vector<double>>> &V, std::vector<std::vector<double>> &Lap,
	std::vector<std::vector<std::vector<double>>> &CellCenter, std::vector<std::vector<double>> &CellVolume,
	std::vector<std::vector<std::vector<double>>> &IFaceCenter, std::vector<std::vector<std::vector<double>>> &IFaceVector,
	std::vector<std::vector<std::vector<double>>> &JFaceCenter, std::vector<std::vector<std::vector<double>>> &JFaceVector,
	std::vector<std::vector<double>> &P, std::vector<std::vector<std::vector<double>>> &Grad);

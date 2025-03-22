#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <math.h>

using namespace std;

void split(const string& s, char delim, vector<string>& elems) {
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
}

vector<string> split(const string& s, char delim) {
	vector<string> elems;
	split(s, delim, elems);
	return elems;
}

struct Node {
	double x, y;
	bool BC;
};

struct Jakobian {
	double J[2][2];
	double J1[2][2];
	double detJ;
};

struct Element {
	int ID[4] = { 0 };
	Jakobian* Jakobian;
	double** H;
	double** C;
	double** Hbc;
	double* P;
};

struct Grid {
	int nN = 0;
	int nE = 0;
	Element* element;
	Node* node;
};

struct GlobalData {
	int SimulationTime = 0;
	int SimulationStepTime = 0;
	int Conductivity = 0;
	int Alfa = 0;
	int Tot = 0;
	int InitialTemp = 0;
	int Density = 0;
	int SpecificHeat = 0;
	int nN = 0;
	int nE = 0;
};

struct ElemUniv {
	double** dNdKsi;
	double** dNdEta;
	double** N;

	struct Surface {
		double** N;
	};

	Surface surface[4];
};

struct Solver {
	double** HG;
	double** CG;
	double* PG;
};

void readData(string file, GlobalData& globaldata, Grid& grid)
{
	fstream data;
	data.open(file, ios::in);
	if (!data.good()) {
		cout << "Cannot found file '" << file << "'" << endl;
		exit(0);
	}

	string linia;
	string curr;
	int curr_node = 0;
	int curr_elem = 0;

	Node node;
	Element element;

	while (getline(data, linia))
	{
		vector<string> slowa = split(linia, ' ');

		if (slowa[0] == "SimulationTime")
		{
			globaldata.SimulationTime = stoi(slowa[1]);
		}
		else if (slowa[0] == "SimulationStepTime")
		{
			globaldata.SimulationStepTime = stoi(slowa[1]);
		}
		else if (slowa[0] == "Conductivity")
		{
			globaldata.Conductivity = stoi(slowa[1]);
		}
		else if (slowa[0] == "Alfa")
		{
			globaldata.Alfa = stoi(slowa[1]);
		}
		else if (slowa[0] == "Tot")
		{
			globaldata.Tot = stoi(slowa[1]);
		}
		else if (slowa[0] == "InitialTemp")
		{
			globaldata.InitialTemp = stoi(slowa[1]);

		}
		else if (slowa[0] == "Density")
		{
			globaldata.Density = stoi(slowa[1]);
		}
		else if (slowa[0] == "SpecificHeat")
		{
			globaldata.SpecificHeat = stoi(slowa[1]);
		}
		else if (slowa[0] == "Nodes")
		{
			globaldata.nN = stoi(slowa[2]);
			grid.nN = stoi(slowa[2]);
			grid.node = new Node[grid.nN];
		}
		else if (slowa[0] == "Elements")
		{
			globaldata.nE = stoi(slowa[2]);
			grid.nE = stoi(slowa[2]);
			grid.element = new Element[grid.nE];
		}
		else if (slowa[0] == "*Node")
		{
			curr = "nodes";
		}
		else if (slowa[0] == "*Element,")
		{
			curr = "element";
		}
		else if (slowa[0] == "*BC") {
			curr = "bc";
		}
		else
		{
			if (curr == "nodes")
			{
				int i = 0;
				while (slowa[i] == "") i++;
				while (slowa[i + 1] == "") i++;
				node.x = stod(slowa[i + 1]);
				while (slowa[i + 1] == "") i++;
				node.y = stod(slowa[i + 2]);
				grid.node[curr_node] = node;
				curr_node++;
			}
			else if (curr == "element")
			{
				int i = 0;
				while (slowa[i] == "") i++;
				while (slowa[i + 1] == "") i++;
				element.ID[0] = stod(slowa[i + 1]) - 1;

				while (slowa[i + 2] == "") i++;
				element.ID[1] = stod(slowa[i + 2]) - 1;

				while (slowa[i + 3] == "") i++;
				element.ID[2] = stod(slowa[i + 3]) - 1;

				while (slowa[i + 4] == "") i++;
				element.ID[3] = stod(slowa[i + 4]) - 1;

				grid.element[curr_elem] = element;
				curr_elem++;
			}
			else if (curr == "bc") {
				int i = 0;
				int j = 0;

				int* BC = new int[slowa.size()];

				while (j < slowa.size())
				{
					while (slowa[i + j] == "") i++;
					BC[j] = stoi(slowa[i + j]);
					j++;
				}

				for (int m = 0; m < grid.nE; m++) {
					for (int n = 0; n < 4; n++) {
						grid.node[grid.element[m].ID[n]].BC = false;
					}
				}

				for (int l = 0; l < grid.nE; l++) {
					for (int m = 0; m < 4; m++) {
						for (int n = 0; n < slowa.size(); n++) {
							if (grid.element[l].ID[m] + 1 == BC[n]) grid.node[grid.element[l].ID[m]].BC = true;
						}
					}
				}

				delete[] BC;
			}
			else {
				data.close();
			}
		}
	}
}

double dN1dKsi(double eta) {
	return -0.25 * (1 - eta);
}

double dN2dKsi(double eta) {
	return 0.25 * (1 - eta);
}

double dN3dKsi(double eta) {
	return 0.25 * (1 + eta);
}

double dN4dKsi(double eta) {
	return -0.25 * (1 + eta);
}

double dN1dEta(double ksi) {
	return -0.25 * (1 - ksi);
}

double dN2dEta(double ksi) {
	return -0.25 * (1 + ksi);
}

double dN3dEta(double ksi) {
	return 0.25 * (1 + ksi);
}

double dN4dEta(double ksi) {
	return 0.25 * (1 - ksi);
}

double* f_xk(int npc) {						// The function returns nodes depending on the number of integration points
	double* xk = new double[sqrt(npc)];

	if (npc == 4) {
		xk[0] = -1 / sqrt(3.0);
		xk[1] = 1 / sqrt(3.0);
	}
	else if (npc == 9) {
		xk[0] = -sqrt(3.0 / 5.0);
		xk[1] = 0.0;
		xk[2] = sqrt(3.0 / 5.0);
	}
	else if (npc == 16) {
		xk[0] = -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));;
		xk[1] = -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
		xk[2] = sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
		xk[3] = sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));;
	}
	else exit(-1);

	return xk;
}

double** dNdKsi(int npc) {

	double** dNdKsi = new double* [npc];

	for (int i = 0; i < npc; i++) {
		dNdKsi[i] = new double[4];
	}

	for (int i = 0; i < npc; i++) {
		for (int j = 0; j < 4; j++) {
			dNdKsi[i][j] = 0;
		}
	}

	int k = 0;
	for (int i = 0; i < npc; i += sqrt(npc)) {
		for (int j = 0; j < sqrt(npc); j++) {
			dNdKsi[i + j][0] = dN1dKsi(f_xk(npc)[k]);
			dNdKsi[i + j][1] = dN2dKsi(f_xk(npc)[k]);
			dNdKsi[i + j][2] = dN3dKsi(f_xk(npc)[k]);
			dNdKsi[i + j][3] = dN4dKsi(f_xk(npc)[k]);
		}
		k++;
	}

	delete[] f_xk(npc);

	return dNdKsi;
}

double** dNdEta(int npc) {

	double** dNdEta = new double* [npc];

	for (int i = 0; i < npc; i++) {
		dNdEta[i] = new double[4];
	}

	for (int i = 0; i < npc; i++) {
		for (int j = 0; j < 4; j++) {
			dNdEta[i][j] = 0;
		}
	}

	int i = 0;
	while (i < sqrt(npc)) {
		for (int j = 0; j < npc; j += sqrt(npc)) {
			dNdEta[i + j][0] = dN1dEta(f_xk(npc)[i]);
			dNdEta[i + j][1] = dN2dEta(f_xk(npc)[i]);
			dNdEta[i + j][2] = dN3dEta(f_xk(npc)[i]);
			dNdEta[i + j][3] = dN4dEta(f_xk(npc)[i]);
		}
		i++;
	}

	delete[] f_xk(npc);

	return dNdEta;
}

Jakobian jakobian(Jakobian jak, ElemUniv elem_univ, double* x, double* y, int p, int npc) {

	double dXdKsi = 0.0;
	double dYdKsi = 0.0;
	double dXdEta = 0.0;
	double dYdEta = 0.0;

	for (int i = 0; i < 4; i++) {
		dXdKsi += elem_univ.dNdKsi[p][i] * x[i];
		dYdKsi += elem_univ.dNdKsi[p][i] * y[i];
		dXdEta += elem_univ.dNdEta[p][i] * x[i];
		dYdEta += elem_univ.dNdEta[p][i] * y[i];
	}

	jak.J[0][0] = dXdKsi;
	jak.J[0][1] = dYdKsi;
	jak.J[1][0] = dXdEta;
	jak.J[1][1] = dYdEta;

	jak.J1[0][0] = dYdEta;
	jak.J1[0][1] = -dYdKsi;
	jak.J1[1][0] = -dXdEta;
	jak.J1[1][1] = dXdKsi;

	jak.detJ = dXdKsi * dYdEta - dYdKsi * dXdEta;

	return jak;
}

double** dNdX(int npc, ElemUniv elem_univ, Element elem) {
	double** dNdX = new double* [npc];

	for (int i = 0; i < npc; i++) {
		dNdX[i] = new double[4];
	}

	for (int i = 0; i < npc; i++) {
		for (int j = 0; j < 4; j++) {
			dNdX[i][j] = 0;
		}
	}

	for (int i = 0; i < npc; i++) {
		for (int j = 0; j < 4; j++) {
			dNdX[i][j] += 1 / elem.Jakobian[i].detJ * elem.Jakobian[i].J1[0][0] * elem_univ.dNdKsi[i][j];
			dNdX[i][j] += 1 / elem.Jakobian[i].detJ * elem.Jakobian[i].J1[0][1] * elem_univ.dNdEta[i][j];
		}
	}

	return dNdX;
}

double** dNdY(int npc, ElemUniv elem_univ, Element elem) {
	double** dNdY = new double* [npc];

	for (int i = 0; i < npc; i++) {
		dNdY[i] = new double[4];
	}

	for (int i = 0; i < npc; i++) {
		for (int j = 0; j < 4; j++) {
			dNdY[i][j] = 0;
		}
	}

	for (int i = 0; i < npc; i++) {
		for (int j = 0; j < 4; j++) {
			dNdY[i][j] += 1 / elem.Jakobian[i].detJ * elem.Jakobian[i].J1[1][0] * elem_univ.dNdKsi[i][j];
			dNdY[i][j] += 1 / elem.Jakobian[i].detJ * elem.Jakobian[i].J1[1][1] * elem_univ.dNdEta[i][j];
		}
	}

	return dNdY;
}

double* f_w(int npc) {					// A function that returns weights depending on the number of integration points
	double* w = new double[sqrt(npc)];

	if (npc == 4) {
		w[0] = 1.0;
		w[1] = 1.0;
	}
	else if (npc == 9) {
		w[0] = 5.0 / 9.0;
		w[1] = 8.0 / 9.0;
		w[2] = 5.0 / 9.0;
	}
	else if (npc == 16) {
		w[0] = (18.0 - sqrt(30)) / 36.0;
		w[1] = (18.0 + sqrt(30)) / 36.0;
		w[2] = (18.0 + sqrt(30)) / 36.0;
		w[3] = (18.0 - sqrt(30)) / 36.0;
	}
	else exit(-1);

	return w;
}

double** H(int npc, double k, Element elem, double** dNdX, double** dNdY, Solver& H_global) {

	double* w = new double[npc];
	int indeks = 0;
	for (int i = 0; i < sqrt(npc); i++) {
		for (int j = 0; j < sqrt(npc); j++) {
			w[indeks] = f_w(npc)[j] * f_w(npc)[i];
			indeks++;
		}
	}

	double H_pc[4][4] = { 0 };

	for (int i = 0; i < npc; i++) {
		for (int j = 0; j < 4; j++) {
			for (int l = 0; l < 4; l++) {
				H_pc[j][l] = k * (dNdX[i][j] * dNdX[i][l] + dNdY[i][j] * dNdY[i][l]) * elem.Jakobian[i].detJ;
				H_pc[j][l] *= w[i];
				elem.H[j][l] += H_pc[j][l];
			}
		}
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			H_global.HG[elem.ID[i]][elem.ID[j]] += elem.H[i][j];
			H_global.HG[elem.ID[i]][elem.ID[j]] += elem.Hbc[i][j];
		}
	}

	delete[] f_w(npc);
	delete[] w;

	return elem.H;
}

double N1(double ksi, double eta) {
	return 0.25 * (1 - ksi) * (1 - eta);
}

double N2(double ksi, double eta) {
	return 0.25 * (1 + ksi) * (1 - eta);
}

double N3(double ksi, double eta) {
	return 0.25 * (1 + ksi) * (1 + eta);
}

double N4(double ksi, double eta) {
	return 0.25 * (1 - ksi) * (1 + eta);
}

void surface(ElemUniv& elem_univ, int npc2) {

	double ksi = 0.0;
	double eta = 0.0;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < npc2; j++) {
			if (i == 0) {
				ksi = f_xk(npc2 * npc2)[j];
				eta = -1.0;
			}
			if (i == 1) {
				ksi = 1.0;
				eta = f_xk(npc2 * npc2)[j];
			}
			if (i == 2) {
				ksi = f_xk(npc2 * npc2)[j];
				eta = 1.0;
			}
			if (i == 3) {
				ksi = -1.0;
				eta = f_xk(npc2 * npc2)[j];
			}
			elem_univ.surface[i].N[j][0] = N1(ksi, eta);
			elem_univ.surface[i].N[j][1] = N2(ksi, eta);
			elem_univ.surface[i].N[j][2] = N3(ksi, eta);
			elem_univ.surface[i].N[j][3] = N4(ksi, eta);
		}
	}

	delete[] f_xk(npc2 * npc2);
}

double length(double x1, double x2, double y1, double y2) {
	return sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2));
}

double** Hbc(double alpha, int npc2, Grid grid, Element elem, ElemUniv elem_univ) {

	double Hbc_pc[4][4] = { 0.0 };

	for (int i = 0, j = 1; i < 4; i++, j++) {
		if ((i + 1) >= 4) j = 0;
		if (grid.node[elem.ID[i]].BC == true && grid.node[elem.ID[j]].BC == true) {
			double len = length(grid.node[elem.ID[i]].x, grid.node[elem.ID[j]].x, grid.node[elem.ID[i]].y, grid.node[elem.ID[j]].y);
			double detJ = len / 2.0;

			for (int k = 0; k < npc2; k++) {
				for (int m = 0; m < 4; m++) {
					for (int n = 0; n < 4; n++) {
						Hbc_pc[m][n] = alpha * (elem_univ.surface[i].N[k][m] * elem_univ.surface[i].N[k][n]) * detJ;
						Hbc_pc[m][n] *= f_w(npc2 * npc2)[k];
						elem.Hbc[m][n] += Hbc_pc[m][n];
					}
				}
			}
		}
	}

	delete f_w(npc2 * npc2);

	return elem.Hbc;
}

double* P(double alpha, double tot, int npc2, Grid grid, Element elem, ElemUniv elem_univ, Solver& P_global) {

	double P_pc[4] = { 0.0 };

	for (int i = 0, j = 1; i < 4; i++, j++) {
		if ((i + 1) >= 4) j = 0;
		if (grid.node[elem.ID[i]].BC == true && grid.node[elem.ID[j]].BC == true) {
			double len = length(grid.node[elem.ID[i]].x, grid.node[elem.ID[j]].x, grid.node[elem.ID[i]].y, grid.node[elem.ID[j]].y);
			double detJ = len / 2.0;

			for (int k = 0; k < npc2; k++) {
				for (int m = 0; m < 4; m++) {
					P_pc[m] = alpha * (elem_univ.surface[i].N[k][m] * tot) * detJ;
					P_pc[m] *= f_w(npc2 * npc2)[k];
					elem.P[m] += P_pc[m];
				}
			}
		}
	}

	for (int i = 0; i < 4; i++) {
		P_global.PG[elem.ID[i]] += elem.P[i];
	}

	delete f_w(npc2 * npc2);

	return elem.P;
}

void initialize(Grid& grid, Solver& solver, ElemUniv& elem_univ, int npc, int npc2) {		// Initializing objects and filling them initially with zeros

	for (int i = 0; i < grid.nE; i++) {
		grid.element[i].Jakobian = new Jakobian[npc];

		for (int j = 0; j < npc; j++) {
			grid.element[i].Jakobian[j].detJ = 0.0;
			for (int k = 0; k < 2; k++) {
				for (int m = 0; m < 2; m++) {
					grid.element[i].Jakobian[j].J[k][m] = 0.0;
					grid.element[i].Jakobian[j].J1[k][m] = 0.0;
				}
			}
		}

		grid.element[i].H = new double* [4];
		grid.element[i].Hbc = new double* [4];
		grid.element[i].C = new double* [4];
		for (int j = 0; j < 4; j++) {
			grid.element[i].H[j] = new double[4];
			grid.element[i].Hbc[j] = new double[4];
			grid.element[i].C[j] = new double[4];
		}
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				grid.element[i].H[j][k] = 0.0;
				grid.element[i].Hbc[j][k] = 0.0;
				grid.element[i].C[j][k] = 0.0;
			}
		}

		grid.element[i].P = new double[4];
		for (int j = 0; j < 4; j++) {
			grid.element[i].P[j] = 0.0;
		}
	}

	elem_univ.dNdKsi = new double* [npc];
	elem_univ.dNdEta = new double* [npc];
	elem_univ.N = new double* [npc];
	for (int i = 0; i < npc; i++) {
		elem_univ.dNdKsi[i] = new double[4];
		elem_univ.dNdEta[i] = new double[4];
		elem_univ.N[i] = new double[4];
	}
	for (int i = 0; i < npc; i++) {
		for (int j = 0; j < 4; j++) {
			elem_univ.dNdKsi[i][j] = 0.0;
			elem_univ.dNdEta[i][j] = 0.0;
			elem_univ.N[i][j] = 0.0;
		}
	}

	for (int i = 0; i < 4; i++) {
		elem_univ.surface[i].N = new double* [npc];
		for (int j = 0; j < npc2; j++) {
			elem_univ.surface[i].N[j] = new double[4];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < npc2; j++) {
			for (int l = 0; l < 4; l++) {
				elem_univ.surface[i].N[j][l] = 0.0;
			}
		}
	}

	solver.HG = new double* [grid.nN];
	solver.CG = new double* [grid.nN];
	for (int i = 0; i < grid.nN; i++) {
		solver.HG[i] = new double[grid.nN];
		solver.CG[i] = new double[grid.nN];
	}
	for (int i = 0; i < grid.nN; i++) {
		for (int j = 0; j < grid.nN; j++) {
			solver.HG[i][j] = 0.0;
			solver.CG[i][j] = 0.0;
		}
	}

	solver.PG = new double[grid.nN];
	for (int i = 0; i < grid.nN; i++) {
		solver.PG[i] = 0.0;
	}
}

double* ElimGauss(double** A, double* b, int n) {	// A function that solves a system of equations using the Gaussian elimination method

	double* X = new double[n];
	for (int i = 0; i < n; i++) {
		X[i] = 0.0;
	}

	for (int i = 0; i < n - 1; i++) {
		int max = i;
		for (int j = i + 1; j < n; j++) {
			if (fabs(A[j][i]) > fabs(A[max][i])) max = j;
		}

		if (max != i) {
			for (int k = 0; k < n; k++) {
				swap(A[i][k], A[max][k]);
			}
			swap(b[i], b[max]);
		}

		if (fabs(A[i][i]) < 1e-12) {
			delete[] X;
			return nullptr;
		}

		for (int j = i + 1; j < n; j++) {
			double wsp = A[j][i] / A[i][i];
			for (int k = i; k < n; k++) {
				A[j][k] -= wsp * A[i][k];
			}
			b[j] -= wsp * b[i];
		}
	}

	for (int i = n - 1; i >= 0; i--) {
		X[i] = b[i];
		for (int j = n - 1; j > i; j--) {
			X[i] -= A[i][j] * X[j];
		}
		X[i] = X[i] / A[i][i];
	}

	return X;
}

double** N(int npc) {

	double** N = new double* [npc];

	for (int i = 0; i < npc; i++) {
		N[i] = new double[4];
	}

	for (int i = 0; i < npc; i++) {
		for (int j = 0; j < 4; j++) {
			N[i][j] = 0;
		}
	}

	int k = 0;
	for (int i = 0; i < sqrt(npc); i++) {
		for (int j = 0; j < sqrt(npc); j++) {
			N[k][0] = N1(f_xk(npc)[j], f_xk(npc)[i]);
			N[k][1] = N2(f_xk(npc)[j], f_xk(npc)[i]);
			N[k][2] = N3(f_xk(npc)[j], f_xk(npc)[i]);
			N[k][3] = N4(f_xk(npc)[j], f_xk(npc)[i]);

			k++;
		}
	}

	delete[] f_xk(npc);

	return N;
}

double** C(int npc, double rho, double cp, Element elem, ElemUniv elem_univ, Solver& C_global) {

	double* w = new double[npc];
	int indeks = 0;
	for (int i = 0; i < sqrt(npc); i++) {
		for (int j = 0; j < sqrt(npc); j++) {
			w[indeks] = f_w(npc)[j] * f_w(npc)[i];
			indeks++;
		}
	}

	double C_pc[4][4] = { 0 };

	for (int i = 0; i < npc; i++) {
		for (int j = 0; j < 4; j++) {
			for (int l = 0; l < 4; l++) {
				C_pc[j][l] = rho * cp * (elem_univ.N[i][j] * elem_univ.N[i][l]) * elem.Jakobian[i].detJ;
				C_pc[j][l] *= w[i];
				elem.C[j][l] += C_pc[j][l];
			}
		}
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			C_global.CG[elem.ID[i]][elem.ID[j]] += elem.C[i][j];
		}
	}

	delete[] f_w(npc);
	delete[] w;

	return elem.C;
}


int main()
{
	int npc = 4;	// NUMBER OF INTEGRATION POINTS IN AN ELEMENT
	int npc2 = 2;	// NUMBER OF INTEGRATION POINTS ON EDGES

	GlobalData globaldata;
	Grid grid;
	ElemUniv elem_univ;
	Solver solver;

	string filename1 = "Test1_4_4.txt";
	string filename2 = "Test2_4_4_MixGrid.txt";
	string filename3 = "Test3_31_31_kwadrat.txt";

	readData(filename1, globaldata, grid);

	cout << "Input data file: " << filename1 << endl;
	cout << "Number of integration points in element: " << npc << endl;
	cout << "Number of integration points on edges: " << npc2 << endl << endl;

	cout << "SimulationTime: " << globaldata.SimulationTime << endl;
	cout << "SimulationStepTime: " << globaldata.SimulationStepTime << endl;
	cout << "Conductivity: " << globaldata.Conductivity << endl;
	cout << "Alfa: " << globaldata.Alfa << endl;
	cout << "Tot: " << globaldata.Tot << endl;
	cout << "InitialTemp: " << globaldata.InitialTemp << endl;
	cout << "Density: " << globaldata.Density << endl;
	cout << "SpecificHeat: " << globaldata.SpecificHeat << endl;
	cout << "Nodes number: " << globaldata.nN << endl;
	cout << "Elements number: " << globaldata.nE << endl;
	cout << endl;

	cout << "Nodes: " << endl;
	for (int i = 0; i < grid.nN; i++)
	{
		cout << i + 1 << ".\t" << setprecision(10) << grid.node[i].x << ",\t ";
		if (grid.node[i].x == 0) cout << "\t";
		cout << setprecision(10) << grid.node[i].y << endl;
	}
	cout << endl;
	cout << "Elements: " << endl;
	for (int i = 0; i < grid.nE; i++)
	{
		cout << i + 1 << ". ";
		for (int j = 0; j < 4; j++)
		{
			cout << grid.element[i].ID[j] << " ";
		}
		cout << endl;
	}
	cout << "\n\n\n";

	double k = globaldata.Conductivity;
	double alpha = globaldata.Alfa;
	double tot = globaldata.Tot;
	double rho = globaldata.Density;
	double cp = globaldata.SpecificHeat;
	double tau = globaldata.SimulationTime;
	double dTau = globaldata.SimulationStepTime;

	double x[4] = { 0 };
	double y[4] = { 0 };

	initialize(grid, solver, elem_univ, npc, npc2);
	elem_univ.dNdKsi = dNdKsi(npc);
	elem_univ.dNdEta = dNdEta(npc);
	elem_univ.N = N(npc);
	surface(elem_univ, npc2);

	double* t0 = new double[grid.nN];		// Initial temperature
	for (int i = 0; i < grid.nN; i++) {
		t0[i] = globaldata.InitialTemp;
	}

	double* t1 = new double[grid.nN];		// Support tables for storing results
	double* min = new double[tau / dTau];
	double* max = new double[tau / dTau];
	for (int i = 0; i < grid.nN; i++) {
		t1[i] = 0.0;
	}
	for (int i = 0; i < tau / dTau; i++) {
		min[i] = 0.0;
		max[i] = 0.0;
	}

	double** H_C = new double* [grid.nN];	// [H] + ([C] / dTau)
	double* C_P = new double[grid.nN];		// {P} + ([C] / dTau) * {t0}
	for (int i = 0; i < grid.nN; i++) {
		H_C[i] = new double[grid.nN];
	}
	for (int i = 0; i < grid.nN; i++) {
		for (int j = 0; j < grid.nN; j++) {
			H_C[i][j] = 0.0;
		}
		C_P[i] = 0.0;
	}

	for (int t = 0, nr = 0; t < tau; t += dTau, nr++) {		// Beginning of loop over time

		for (int i = 0; i < grid.nE; i++) {		// Beginning of loop over elements

			for (int j = 0; j < 4; j++) {
				x[j] = grid.node[grid.element[i].ID[j]].x;
				y[j] = grid.node[grid.element[i].ID[j]].y;
			}

			for (int l = 0; l < npc; l++) {
				grid.element[i].Jakobian[l] = jakobian(grid.element[i].Jakobian[l], elem_univ, x, y, l, npc);
			}

			grid.element[i].Hbc = Hbc(alpha, npc2, grid, grid.element[i], elem_univ);
			grid.element[i].H = H(npc, k, grid.element[i], dNdX(npc, elem_univ, grid.element[i]), dNdY(npc, elem_univ, grid.element[i]), solver);
			grid.element[i].P = P(alpha, tot, npc2, grid, grid.element[i], elem_univ, solver);
			grid.element[i].C = C(npc, rho, cp, grid.element[i], elem_univ, solver);

			/*cout << endl << "Matrix [Hbc] for element " << i + 1 << ": " << endl;				// Auxiliary printing
			for (int m = 0; m < 4; m++) {														// Matrix [Hbc] local
				for (int n = 0; n < 4; n++) {													// Matrix [H] local (after adding Hbc)
					cout << grid.element[i].Hbc[m][n] << "\t";									// Matrix [C] local
				}																				// Vector [P] local
				cout << endl;
			}
			for (int m = 0; m < 4; m++) {
				for (int n = 0; n < 4; n++) {
					grid.element[i].H[m][n] += grid.element[i].Hbc[m][n];
				}
			}
			cout << endl << "Matrix [H] for element " << i + 1 << " (after adding Hbc): " << endl;
			for (int m = 0; m < 4; m++) {
				for (int n = 0; n < 4; n++) {
					cout << grid.element[i].H[m][n] << "\t";
				}
				cout << endl;
			}
			cout << endl << "Matrix [C] for element " << i + 1 << ": " << endl;
			for (int m = 0; m < 4; m++) {
				for (int n = 0; n < 4; n++) {
					cout << grid.element[i].C[m][n] << "\t";
				}
				cout << endl;
			}
			cout << endl << "Vector [P] for element " << i + 1 << ": " << endl;
			for (int m = 0; m < 4; m++) {
				cout << grid.element[i].P[m] << "\t";
			}
			cout << "\n\n\n";*/

			for (int m = 0; m < 4; m++) {				// Zeros local matrices and vectors
				for (int n = 0; n < 4; n++) {
					grid.element[i].H[m][n] = 0.0;
					grid.element[i].Hbc[m][n] = 0.0;
					grid.element[i].C[m][n] = 0.0;
				}
				grid.element[i].P[m] = 0.0;
			}
		}

		/*cout << endl << endl;							// Auxiliary printing
		cout << "Matrix [H] global:" << endl;			// Matrix [H] global
		for (int i = 0; i < grid.nN; i++) {				// Matrix [C] global
			for (int j = 0; j < grid.nN; j++) {			// Vector [P] global
				cout << solver.HG[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << "Matrix [C] global:" << endl;
		for (int i = 0; i < grid.nN; i++) {
			for (int j = 0; j < grid.nN; j++) {
				cout << solver.CG[i][j];
				if (j < grid.nN - 1) cout << " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << "Vector [P] global:" << endl;
		for (int i = 0; i < grid.nN; i++) {
			cout << solver.PG[i] << " ";
		}
		cout << "\n\n\n";*/

		for (int i = 0; i < grid.nN; i++) {						// [H] + ([C] / dTau) ("A" in the equation)
			for (int j = 0; j < grid.nN; j++) {
				H_C[i][j] = solver.HG[i][j] + (solver.CG[i][j] / dTau);
			}
		}

		for (int i = 0; i < grid.nN; i++) {						// P + C/dT * T0 ("B" in the equation)
			for (int j = 0; j < grid.nN; j++) {
				C_P[i] += (solver.CG[i][j] / dTau) * t0[j];
			}
			C_P[i] += solver.PG[i];
		}

		/*cout << "H_C:" << endl;
		for (int m = 0; m < grid.nN; m++) {
			for (int n = 0; n < grid.nN; n++) {
				cout << H_C[m][n] << " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << "C_P "<< t + dTau<<":" << endl;
		for (int m = 0; m < grid.nN; m++) {
			cout << C_P[m] << " ";
		}
		cout << "\n\n";*/

		t1 = ElimGauss(H_C, C_P, grid.nN);				// Solving a system of equations (Ax + B = 0, where x = t1)

		min[nr] = t1[0];
		max[nr] = t1[0];
		for (int i = 0; i < grid.nN; i++) {				// Looking for the minimum and maximum temperature after a given time
			if (t1[i] > max[nr]) max[nr] = t1[i];
			if (t1[i] < min[nr]) min[nr] = t1[i];
		}

		cout << "Temperature over time " << t + dTau << "s:" << endl;
		for (int i = 0; i < grid.nN; i++) {
			cout << t1[i] << " ";

		}
		cout << endl;
		cout << "Min: " << min[nr] << "\tMax: " << max[nr] << endl;
		cout << "\n\n";

		for (int i = 0; i < grid.nN; i++) {				// Zeros global matrices and vectors
			for (int j = 0; j < grid.nN; j++) {
				solver.HG[i][j] = 0.0;
				solver.CG[i][j] = 0.0;
				H_C[i][j] = 0.0;
			}
			solver.PG[i] = 0.0;
			C_P[i] = 0.0;
		}

		for (int i = 0; i < grid.nN; i++) {
			t0[i] = t1[i];								// In the next iteration, t0 is assumed to be t1 from the previous iteration
		}
	}

	cout << "Summary of results (number of integration points - " << sqrt(npc) << "):" << endl;
	cout << "Time[s]\t\tMinTemp\t\tMaxTemp" << endl;
	for (int i = 0, nr = 0; i < tau; i += dTau, nr++) {
		cout << i + dTau << "\t\t" << min[nr] << "\t" << max[nr] << endl;
	}

	delete[] t0;
	delete[] t1;										// Clearing memory
	delete[] min;
	delete[] max;
	delete[] grid.element;
	delete[] grid.node;
	delete[] dNdKsi(npc);
	delete[] dNdEta(npc);
	delete[] elem_univ.dNdKsi;
	delete[] elem_univ.dNdEta;
	delete[] elem_univ.N;
	delete[] solver.HG;
	delete[] solver.PG;
	delete[] solver.CG;
	for (int i = 0; i < 4; i++) {
		delete[] elem_univ.surface[i].N;
	}

	return 0;
}
// Storage for cell information
// Created by James McClure
// Copyright 2008-2011
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

//#define DEBUG

using namespace std;

double PI = 3.141592653589793;
//******************************************************
// STRUCTURES TO STORE SPHERES
//******************************************************

// ........................................................
//   CELL STORAGE STRUCTURE FOR PARTICLES
//		- cell (icx,icy,icz) :
//			Lower Boundary = icx*lencx,icy*lency,icz*lencz
//			Upper Boundary = (icx+1)*lencx,(icy+1)*lency,(icz+1)*lencz
//		- particle is inside of cell if cell contains centroid
//		- particle data stored externally (centroid, geometry)
// ........................................................
struct CellStorage {
	CellStorage(int &nx, int &ny, int &nz, double &lenx, double &leny, double &lenz);
	~CellStorage();
	// Set up cells to store all spheres
	
	int ncells,ncx,ncy,ncz;		// # cells
	int maxcell;				// max # particles per cell
	double lencx,lency,lencz;	// cell lengths
	
	int * DATA;					// store particles
	int * COUNT;				// actual # particles per cell

	int & CellCount(int icx, int icy, int icz) {
		return COUNT[icz*ncx*ncy+icy*ncx+icx]; 
	}
	
	int & CellEntry(int icx, int icy, int icz, int index) {
		i = icz*ncx*ncy+icy*ncx+icx;
		return DATA[maxcell*i+index];
	}
	
	void Reset();
	
private:
	int i;
};
//........ constructor .....................................
CellStorage::CellStorage(int &nx, int &ny, int &nz, double &lenx, double &leny, double &lenz)
{
	ncx = nx;					// # cells in x direction
	ncy = ny;					// # cells in y direction
	ncz = nz;					// # cells in z direction
	ncells = ncx*ncy*ncz;		// total # cells
	
	lencx = lenx / ncx;			// cell length in x direction
	lency = leny / ncy;			// cell length in y direction
	lencz = lenz / ncz;			// cell length in z direction
	
	maxcell = 10000;
	
	COUNT = new int [ncells];	
	DATA = new int [maxcell*ncells];	
}
//........ destructor .....................................
CellStorage::~CellStorage()
{
	delete COUNT;
	delete DATA;
}

// ...... reset the cell count ............................
void CellStorage::Reset()
{
	for (i=0;i<ncells;i++) COUNT[i] = 0;
}

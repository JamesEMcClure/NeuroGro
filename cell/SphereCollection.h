// store a collection of spheres
#include "cell/CellStorage.h"

// ........................................................
//	STORAGE FOR A COLLECTION OF SPHERES
//		- store particle centroid, radius	(cx,cy,cz,Radius)
//		- store the distance to move each particle (dx,dy,dz)
// ........................................................
struct SphereCollection {
	SphereCollection(int &n);
	~SphereCollection();
	int N;					// number of particles
	//.............................................
	// | centroid | radius | 
	double * data;		
	//.......... Access ..................
	// center of mass / centroid for particle i
	double & cx(int i) { return data[7*i];}
	double & cy(int i) { return data[7*i+1];}
	double & cz(int i){ return data[7*i+2];}
	// distance to move particle i
	double & dx(int i) {return data[7*i+3];}
	double & dy(int i){ return data[7*i+4];}
	double & dz(int i){ return data[7*i+5];}
	// cylinder length
	double & Radius(int i){ return data[7*i+6];}	
//	double & Length(int &i){ return data[8*i+7];}
	
	// ..... Initialize a system of spheres 
	void Initialize(double &Lx, double &Ly, double &Lz, double &porosity, 
					double &mu, double &sig);
	
//	int GetContacts(int & i);
//	int & Contact(int &i) {	return work[i];}
private:
	int i;
};

//........ constructor..........................
SphereCollection::SphereCollection(int &n)
{
	N = n;
	data = new double [7*N];
}

//........ destructor ...........................
SphereCollection::~SphereCollection()
{
	delete data;
}

// ...... initialization .........................
void SphereCollection::Initialize(double &Lx, double &Ly, double &Lz, double &porosity, 
				double &mu, double &sig)
{
	// double mu;
	double cxi,cyi,czi,r;
	srand((unsigned)time(0));
	// Expected volume for the 
	//  r = CUBE_ROOT( 3*VOLUME/4/PI) 
	mu = 0.333333333333*log(0.75*(1.0-porosity)*Lx*Ly*Lz/PI/N) - 1.5*sig;

	// r = exp( 0.333333333333333*log((1.0-porosity)*0.75*Lx*Ly*Lz/N/PI));
	for (i=0;i<N;i++){
//		cout << "INITIALIZE SPHERE " << 1.0*rand()/RAND_MAX << endl;
		// Uniformly distribute the centroids
		cxi = rand()*Lx/RAND_MAX;
		cyi = rand()*Ly/RAND_MAX;
		czi = rand()*Lz/RAND_MAX;
		cx(i) = cxi;
		cy(i) = cyi;
		cz(i) = czi;
		// Generate radii from lognormal distribution, truncate lower & upper  1.5*sig
		r = exp(sqrt(-2*log(1.0*rand()/RAND_MAX))*cos(2.0*PI*rand()/RAND_MAX)*sqrt(sig)+mu);
		while (r < exp(mu)/pow(exp(1.5*sqrt(sig)),2) || r > exp(mu)*pow(exp(1.5*sqrt(sig)),2) ){
			r = exp(sqrt(-2*log(1.0*rand()/RAND_MAX))*cos(2.0*PI*rand()/RAND_MAX)*sqrt(sig)+mu);
		}
		Radius(i) = r;
//		cout << "	Radius value: " << r << endl; 
	}
}

inline double CoordinationNumber(SphereCollection &Particles, double size,	
								 double Lx, double Ly, double Lz)
{
	double dist;
	int count;
	count = 0;
	for (int i=0;i<Particles.N;i++){
		// Compute the number of contacts for particle i
		for (int j=0;j<Particles.N;j++){
			// compute distance
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
					+ pow(Particles.cy(i)-Particles.cy(j),2)
					+ pow(Particles.cz(i)-Particles.cz(j),2))
					- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			// compute distance ( -x Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( +x Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			// compute distance ( -y Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			// compute distance ( +y Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			// compute distance ( -x Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( +x Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			// compute distance ( -y Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			// compute distance ( +y Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
		}
	}
	// Compute the averge coordination number 
	double cn = double (count) / double(Particles.N);
	return cn;
}
inline void Random_Displacement(SphereCollection &Particles, double size)
{
	int i;
	for (i=0;i<Particles.N;i++){
		Particles.cx(i) + size*(2*rand()/RAND_MAX - 1)*Particles.Radius(i);
		Particles.cy(i) + size*(2*rand()/RAND_MAX - 1)*Particles.Radius(i);
		Particles.cz(i) + size*(2*rand()/RAND_MAX - 1)*Particles.Radius(i);
	}
}

inline bool CHECK_OVERLAPS(SphereCollection &Particles, CellStorage &Storage, double tol)
{
	bool toReturn = false;
	int icx,icy,icz,jcx,jcy,jcz,number,count,i,j;
	int ii,jj;
	// Check for any overlaps
	for (icz=0; icz<Storage.ncz; icz++){
		for (icy=0; icy<Storage.ncy; icy++){
			for (icx=0; icx<Storage.ncx; icx++){
				// How many particles in this cell 
				number = Storage.CellCount(icx,icy,icz);
				for (i=0; i<number; i++){
					ii = Storage.CellEntry(icx,icy,icz,i);
					// Go over all the other cells 
					for (jcz=0; jcz<Storage.ncz; jcz++){
						for (jcy=0; jcy<Storage.ncy; jcy++){
							for (jcx=0; jcx<Storage.ncx; jcx++){
								count = Storage.CellCount(jcx,jcy,jcz);
								for (j=0; j<count; j++){
									jj = Storage.CellEntry(jcx,jcy,jcz,j);
									// Check for overlap between ii & jj
									if ( icx == jcx && icy == jcy && icz == jcz && ii==jj){
										// Same particle - don't count overlap with self
									}
									else {
										if ( sqrt( pow(Particles.cx(ii)-Particles.cx(jj),2)
												  +pow(Particles.cy(ii)-Particles.cy(jj),2)
												  +pow(Particles.cz(ii)-Particles.cz(jj),2) )
											- Particles.Radius(ii) - Particles.Radius(jj) > tol){
#ifdef DEBUG
											cout << "ERROR: Overlap between" << ii << " and " << jj << endl; 
#endif
											toReturn = true;
											
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	if ( toReturn == true){
		cout << "WARNING: Overlaps exist,  consider decreasing the number of cells!" << endl;
	}
	return toReturn;
}



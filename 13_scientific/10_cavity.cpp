#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <vector>
#include <cmath>
#include <mpi.h>

using namespace std;
typedef vector<vector<float>> matrix;

int main( int argc, char** argv ) {

  MPI_Init(&argc, &argv);
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nx = 41;
  int ny = 41;
  int nt = 500;
  int nit = 50;
  double dx = 2. / (nx - 1);
  double dy = 2. / (ny - 1);

  int sx = nx / size;
  if( nx % sx != 0 ) sx = sx + 1;
  int begin = rank * sx;
  int end = ( rank + 1 ) * sx;
  if( end > nx ) end = nx;

  double sdx = dx * dx;
  double sdy = dy * dy;

  double dt = .01;
  double rho = 1.;
  double nu = .02;

  float u[ny][nx];
  float v[ny][nx];
  float p[ny][nx];
  float b[ny][nx];
  float un[ny][nx];
  float vn[ny][nx];
  float pn[ny][nx];

  int length = end-begin;
  float up[length];
  float vp[length];
  float pp[length];
  float bp[length];
  float unp[length];
  float vnp[length];
  float pnp[length];

  for (int j=0; j<ny; j++) {
    for (int i=0; i<length; i++) {
      up[i] = 0;
      vp[i] = 0;
      pp[i] = 0;
      bp[i] = 0;
    }
  }

  for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
      u[j][i] = 0;
      v[j][i] = 0;
      p[j][i] = 0;
      b[j][i] = 0;
    }
  }
  ofstream ufile;
  ofstream vfile;
  ofstream pfile;
  if( rank == 0 ){
  	ufile.open("u.dat");
  	vfile.open("v.dat");
  	pfile.open("p.dat");
  }
  for (int n=0; n<nt; n++) {
    for (int j=1; j<ny-1; j++) {
      for (int i=begin; i<end; i++) {
        // Compute b[j][i]
	if( ( i == 0 ) | ( i == nx - 1 ) ){
		bp[i-begin] = b[j][i];
		continue;
	}
	bp[i-begin] = rho * (1 / dt *((u[j][i+1] - u[j][i-1]) / (2 * dx) + (v[j+1][i] - v[j-1][i]) / (2 * dy)) - pow( ((u[j][i+1] - u[j][i-1]) / (2 * dx) ) , 2 ) - 2 * ((u[j+1][i] - u[j-1][i]) / (2 * dy) * (v[j][i+1] - v[j][i-1]) / (2 * dx)) - pow( ((v[j+1][i] - v[j-1][i]) / (2 * dy)) , 2 ) );
      }
      MPI_Gather(&bp[0], length, MPI_FLOAT,&b[j][0], length, MPI_FLOAT,0,MPI_COMM_WORLD);
    }
    MPI_Bcast( &b[0][0], nx*ny, MPI_FLOAT, 0,MPI_COMM_WORLD );
    for (int it=0; it<nit; it++) {
      for (int j=0; j<ny; j++){
        for (int i=begin; i<end; i++)
	  pnp[i-begin] = p[j][i];
        MPI_Gather(&pnp[0], length, MPI_FLOAT,&pn[j][0], length, MPI_FLOAT,0,MPI_COMM_WORLD);
      }
      MPI_Bcast( &pn[0][0], nx*ny, MPI_FLOAT, 0,MPI_COMM_WORLD );
      for (int j=1; j<ny-1; j++) {
        for (int i=begin; i<end; i++) {
	  // Compute p[j][i]
	  if( ( i == 0 ) | ( i == nx - 1 ) ){
		  pp[i-begin] = p[j][i];
		  continue;
	  }
	  pp[i-begin] = (sdy * (pn[j][i+1] + pn[j][i-1]) + sdx * (pn[j+1][i] + pn[j-1][i]) - b[j][i] * sdx * sdy) / (2 * (sdx + sdy));
	}
	MPI_Gather(&pp[0], length, MPI_FLOAT,&p[j][0], length, MPI_FLOAT,0,MPI_COMM_WORLD);
      }
      MPI_Bcast( &p[0][0], nx*ny, MPI_FLOAT, 0,MPI_COMM_WORLD );
      for (int j=0; j<ny; j++) {
        // Compute p[j][0] and p[j][nx-1]
	p[j][nx-1] = p[j][nx-2];
	p[j][0] = p[j][1];
      }
      for (int i=0; i<nx; i++) {
	// Compute p[0][i] and p[ny-1][i]
	p[0][i] = p[1][i];
	p[ny-1][i] = 0;
      }
    }
    for (int j=0; j<ny; j++) {
      for (int i=begin; i<end; i++) {
        unp[i-begin] = u[j][i];
	vnp[i-begin] = v[j][i];
      }
      MPI_Gather(&unp[0], length, MPI_FLOAT,&un[j][0], length, MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Gather(&vnp[0], length, MPI_FLOAT,&vn[j][0], length, MPI_FLOAT,0,MPI_COMM_WORLD);
      // if( rank == 0 ) printf( "%i \n" , n , 4 , j );
    }
    MPI_Bcast( &un[0][0], nx*ny, MPI_FLOAT, 0,MPI_COMM_WORLD );
    MPI_Bcast( &vn[0][0], nx*ny, MPI_FLOAT, 0,MPI_COMM_WORLD );
    for (int j=1; j<ny-1; j++) {
      for (int i=begin; i<end; i++) {
	// Compute u[j][i] and v[j][i]
	if( ( i == 0 ) | ( i == nx - 1 ) ){
	       up[i-begin] = u[j][i];
	       vp[i-begin] = v[j][i];
	       continue;
	}
	up[i-begin] = un[j][i] - un[j][i] * dt / dx * (un[j][i] - un[j][i - 1])- vn[j][i] * dt / dy * (un[j][i] - un[j - 1][i])- dt / (2 * rho * dx) * (p[j][i+1] - p[j][i-1])+ nu * dt / sdx * (un[j][i+1] - 2 * un[j][i] + un[j][i-1])+ nu * dt / sdy * (un[j+1][i] - 2 * un[j][i] + un[j-1][i]);
        vp[i-begin] = vn[j][i] - un[j][i] * dt / dx * (vn[j][i] - vn[j][i - 1])- vn[j][i] * dt / dy * (vn[j][i] - vn[j - 1][i])- dt / (2 * rho * dx) * (p[j+1][i] - p[j-1][i])+ nu * dt / sdx * (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1])+ nu * dt / sdy * (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i]);
      }
      MPI_Gather(&up[0], length, MPI_FLOAT,&u[j][0], length, MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Gather(&vp[0], length, MPI_FLOAT,&v[j][0], length, MPI_FLOAT,0,MPI_COMM_WORLD);
      // if( rank == 0 ) printf( "%i \n" , n , 5 , j );
    }
    MPI_Bcast( &u[0][0], nx*ny, MPI_FLOAT, 0,MPI_COMM_WORLD );
    MPI_Bcast( &v[0][0], nx*ny, MPI_FLOAT, 0,MPI_COMM_WORLD );
    for (int j=0; j<ny; j++) {
      // Compute u[j][0], u[j][nx-1], v[j][0], v[j][nx-1]
      u[j][0]  = 0;
      u[j][nx-1] = 0;
      v[j][0]  = 0;
      v[j][nx-1] = 0;
    }
    for (int i=0; i<nx; i++) {
      // Compute u[0][i], u[ny-1][i], v[0][i], v[ny-1][i]
      u[0][i]  = 0;
      u[ny-1][i] = 1;
      v[0][i]  = 0;
      v[ny-1][i] = 0;
    }
    if ( (n % 10 == 0) & ( rank == 0 ) ){
      for (int j=0; j<ny; j++)
        for (int i=0; i<nx; i++)
          ufile << u[j][i] << " ";
      ufile << "\n";
      for (int j=0; j<ny; j++)
        for (int i=0; i<nx; i++)
          vfile << v[j][i] << " ";
      vfile << "\n";
      for (int j=0; j<ny; j++)
        for (int i=0; i<nx; i++)
          pfile << p[j][i] << " ";
      pfile << "\n";
    }
  }
  if( rank == 0 ){
  	ufile.close();
  	vfile.close();
  	pfile.close();
  }
  MPI_Finalize();
}

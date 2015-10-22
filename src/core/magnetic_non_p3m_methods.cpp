/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file magnetic_non_p3m_methods.cpp  All 3d non P3M methods to deal with the magnetic dipoles
 *   
 *  DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA
 *   Handling of a system of dipoles where no replicas 
 *   Assumes minimum image convention for those axis in which the 
 *   system is periodic as defined by setmd.
 *
 *   MDDS => Calculates dipole-dipole interaction of a perioidic system
 *   by explicitly summing the dipole-dipole interaction over several copies of the system
 *   Uses spherical summation order
 *
 */

#include "domain_decomposition.hpp"
#include "magnetic_non_p3m_methods.hpp"

#ifdef DIPOLES

double mu0 = 1.256637061e-6; // my0 = 4e-7*M_PI N/A^2

#ifdef EXCLUDED_VOLUME_FORCE
double evf_A = -1.0;
double evf_xi = -1.0;
double evf_cut = -1.0;
#endif

#ifdef MAGNETIC_CHARGE_SHEETS
// akoun 1984
// HARDCODED position, size and magnetization of magnetic charge sheets
// field and gradients from rectangular charge sheet with uniform magnetization simga = |J| along z
void fields_rectPole(double x, double y, double z, double a, double b, double* dhdx, double* dhdy, double* dhdz, double* h, double sigma, int div){
	double f,s,t,r,r2;
	double fpow, rmint, rmins, dr, st,rz,stdivrz,ff;
	
	//f = 1.0/(4.0*acos( -1.0 ));   acos(-1) == PI
	f = sigma/(4.0*M_PI);//*my0);

	//printf("relative distance %e %e %e\n", x,y,z);

	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			s = x - pow((-1),i) * a;
			t = y - pow((-1),j) * b;
			r2=s*s+t*t+z*z;
			r = sqrt(r2);
			fpow=f*pow((-1),(i+j));
			rmint=r-t;
			rmins=r-s;
			st=s*t;
			rz=r*z;
			stdivrz=(st)/(rz);
			ff=(1+(stdivrz*stdivrz))*rz*rz;
			
			if(div==0){
				// field values
				h[0] += fpow * log(rmint);      // log == ln
				h[1] += fpow * log(rmins);
				h[2] += fpow * atan(stdivrz);
			}
			else if (div==1){
				
				// nach x ableiten
				dr = s / r;
				dhdx[0] += fpow * dr / rmint;
				dhdx[1] += fpow * (dr - 1.0) / rmins;
				dhdx[2] += fpow * (t*rz - st*dr*z) / ff;
				
				// nach y ableiten
				dr = t / r;
				dhdy[0] += fpow * (dr - 1.0) / rmint;
				dhdy[1] += fpow * dr / rmins;
				dhdy[2] += fpow * (s*rz - st*dr*z) / ff;
				
				// nach z ableiten
				dr = z / r;
				dhdz[0] += fpow * dr / rmint;
				dhdz[1] += fpow * dr / rmins;
				dhdz[2] -= fpow * st * (dr*z+r) / ((stdivrz*stdivrz)+1) / (rz*rz);
			}
		}
	}
}

void fields_betweenRectPoles(double x, double y, double z, double magX, double magY, double magZ, double a, double b, double c, double* dhdx, double* dhdy, double* dhdz, double* h, double sigma, int div){
	
	double h1[3]={0,0,0},h2[3]={0,0,0};
	double dhdx1[3]={0,0,0}, dhdx2[3]={0,0,0}, dhdy1[3]={0,0,0}, dhdy2[3]={0,0,0}, dhdz1[3]={0,0,0}, dhdz2[3]={0,0,0};
	
	fields_rectPole(x-magX,y-magY, z-magZ-c, a,b,dhdx1, dhdy1, dhdz1, h1, sigma, div);
	fields_rectPole(x-magX,y-magY, z-magZ+c, a,b,dhdx2, dhdy2, dhdz2, h2, sigma, div);
	
	for(int i=0;i<3;i++){
		if(div==0){
			h[i]=h1[i]-h2[i];
		}
		else if(div==1){
			dhdx[i]=dhdx1[i]-dhdx2[i];
			dhdy[i]=dhdy1[i]-dhdy2[i];
			dhdz[i]=dhdz1[i]-dhdz2[i];
		}
		//dhdx[i]=dhdx1[i];
		//dhdy[i]=dhdy1[i];
		//dhdz[i]=dhdz1[i];
	}
}

void calculate_gradient_field(Particle* p1, double* pos, double* size, double sigma){
    double h[3]={0,0,0};
	double dhdx[3]={0,0,0}, dhdy[3]={0,0,0}, dhdz[3]={0,0,0};
	double p11[3];
	double dip[3]={0,0,0};
	int img[3];
	
	
	memcpy(p11, p1->r.p, 3*sizeof(double));
	memcpy(img, p1->l.i, 3*sizeof(int));
	fold_position(p11, img);
	
	// hardcoded magnet values
	//double magX = 50., magY = 10., magZ = 10.; // center of channel
	/*if(time_step < 1000000) magX = 20;
	if(time_step < 2000000) magX = 50;
	if(time_step < 3000000) magX = 70;*/
	//double a = 10, b = 10, c = 1e8; // charge sheets with dimensions 10x10, poles outside the channel
	//double sigma = 1;//0.05056;//100;//1.6;
	//double chi = 1.0;//0.226; //1.0;
	//double mSat= 1;//10632.1467;//0.124;//87500;  // saturated magnetic moment
	
	// Volume sphere = 4/3 pi r^3 --> r=0.5
	//double ch = 0.524 / (4*M_PI*1e-7);//chi * 0.524 / (4*M_PI*1e-7);
    //double chi = 0.0;
	double ch = 0.0;//0.524 / (4*M_PI*1e-7);//chi * 0.524 / (4*M_PI*1e-7);
    double chi = 0.0;
    double vol = 0.0;
	
	//printf("pos orig-fold: %e %e %e\n", p1->r.p[0]- p11[0], p1->r.p[1]-p11[1], p1->r.p[2]-p11[2]);
	
	//fields_betweenRectPoles(p11[0], p11[1], p11[2], magX, magY, magZ, a, b, c, dhdx, dhdy, dhdz, h, sigma,0);
	fields_betweenRectPoles(p11[0], p11[1], p11[2], pos[0], pos[1], pos[2], size[0], size[1], size[2], dhdx, dhdy, dhdz, h, sigma,0);

	//fields_betweenRectPoles(p1->r.p[0], p1->r.p[1], p1->r.p[2], magX, magY, magZ, a, b, c, dhdx, dhdy, dhdz, h, sigma, 0);
    
    //double ch = chi * 0.524 / (4*M_PI*1e-7); // Volume sphere = 4/3 pi r^3 --> r=0.5
    
    ch = 3*p1->p.susc/(3+p1->p.susc);
	vol = 4.188790205 * p1->p.radius * p1->p.radius * p1->p.radius; // vol = 4/3*pi*r^3
	chi = ch * vol / mu0; 
	
    //chi = ch*p1->p.susc;
    dip[0]= h[0]*chi;  // h equals the external B field
    dip[1]= h[1]*chi;
    dip[2]= h[2]*chi;
    
    // set dipole fields
    p1->r.dip[0] = dip[0]; 
    p1->r.dip[1] = dip[1];
    p1->r.dip[2] = dip[2];
    //if((int)time_step % 10000 == 0)printf("h: %e %e %e, dhdx:  %e %e %e, dhdy:  %e %e %e, \ndhdz:  %e %e %e, fg:  %e %e %e\n=====================", h[0], h[1],h[2], dhdx[0],dhdx[1],dhdx[2], dhdy[0],dhdy[1],dhdy[2], dhdz[0],dhdz[1],dhdz[2],fgx,fgy,fgz); 
}

// Calculate field and gradient force from permanent magnets
void calculate_moment_force(Particle* p1, double* pos, double* size, double sigma){
	double fgx, fgy, fgz;
    double h[3]={0,0,0};
	double dhdx[3]={0,0,0}, dhdy[3]={0,0,0}, dhdz[3]={0,0,0};
	double p11[3];
	double dip[3]={0,0,0},eDIP[3]={0,0,0},DIP=0,DIP2=0;
	int img[3];
	
	
	memcpy(p11, p1->r.p, 3*sizeof(double));
	memcpy(img, p1->l.i, 3*sizeof(int));
	fold_position(p11, img);
	
	// hardcoded magnet values
	//double magX = 12.5, magY = 12.5, magZ = 12.5; // center of channel
	/*if(time_step < 1000000) magX = 20;
	if(time_step >= 1000000 && time_step < 2000000) magX = 50;
	if(time_step >= 2000000 && time_step < 3000000) magX = 70;*/

	
	//double a = 12.5, b = 12.5, c = 1e8; // charge sheets with dimensions 10x10, poles outside the channel
	//double sigma = 1;//0.05056;//100;//1.6;
	//double chi = 1.0;//0.226; //1.0;
	//double mSat= 1;//10632.1467;//5567;//0.124;//87500;  // saturated magnetic moment
	
	//printf("pos orig-fold: %e %e %e\n", p1->r.p[0]- p11[0], p1->r.p[1]-p11[1], p1->r.p[2]-p11[2]);
	
	//fields_betweenRectPoles(p11[0], p11[1], p11[2], magX, magY, magZ, a, b, c, dhdx, dhdy, dhdz, h, sigma,1);
	fields_betweenRectPoles(p11[0], p11[1], p11[2], pos[0], pos[1], pos[2], size[0], size[1], size[2], dhdx, dhdy, dhdz, h, sigma,1);
	
	//fields_betweenRectPoles(p1->r.p[0], p1->r.p[1], p1->r.p[2], magX, magY, magZ, a, b, c, dhdx, dhdy, dhdz, h, sigma, 1);
    
    // calculate gradient force
    fgx = p1->r.dip[0]*dhdx[0] + p1->r.dip[1]*dhdx[1] + p1->r.dip[2]*dhdx[2];
	fgy = p1->r.dip[0]*dhdy[0] + p1->r.dip[1]*dhdy[1] + p1->r.dip[2]*dhdy[2];
	fgz = p1->r.dip[0]*dhdz[0] + p1->r.dip[1]*dhdz[1] + p1->r.dip[2]*dhdz[2];
    
    // set gradient force
    p1->f.f[0] +=fgx;
    p1->f.f[1] +=fgy;
    p1->f.f[2] +=fgz;
    
    //if((int)(sim_time*100)%300==0) printf("X %.2e fg:  %.2e %.2e %.2e, |fg| %.2e\n", p1->r.p[0],fgx,fgy,fgz,sqrt(fgx*fgx+fgy*fgy+fgz*fgz));
    
  //  if((int)time_step % 10000 == 0) printf("h: %e %e %e, dhdx:  %e %e %e, dhdy:  %e %e %e, \ndhdz:  %e %e %e, fg:  %e %e %e magX %e\n=====================", h[0], h[1],h[2], dhdx[0],dhdx[1],dhdx[2], dhdy[0],dhdy[1],dhdy[2], dhdz[0],dhdz[1],dhdz[2],fgx,fgy,fgz,magX); 
    //if((int)time_step % 100000 == 0) printf("h: %e %e %e, fg:  %e %e %e\n", h[0], h[1],h[2],fgx,fgy,fgz); 

}
#endif

#ifdef SOFTMAGNETIC
void calculate_dipolar_field(Particle* p1){
	double fgx, fgy, fgz;
    double h[3]={0,0,0};
	double p11[3], p22[3];
	double dip[3]={0,0,0},eDIP[3]={0,0,0},DIP=0,DIP2=0;
	int img[3];
	double dr[3]={0,0,0},r,r2,r3,r5,h2[3]={0,0,0};
    
    //double chi = 1.0;//0.226; //1.0;
	//double mSat= 1;//10632.1467;//5567;//0.124;//87500;  // saturated magnetic moment
    // M = chi*H, chi needs prefactor due to shape of sphere
    
    // Volume sphere = 4/3 pi r^3 --> r=0.5
    //double ch = 3*susc/(3+susc)
    double ch = 0.0;//0.524 / (4*M_PI*1e-7);//chi * 0.524 / (4*M_PI*1e-7);
    double chi = 0.0;
    double vol = 0.0;
    
    memcpy(p11, p1->r.p, 3*sizeof(double));
	memcpy(img, p1->l.i, 3*sizeof(int));
	fold_position(p11, img);
    
    //int flag=0;
    Particle* p;
    double dip2[3]={0,0,0}, ff=0.0, m1r=0.0,m2r=0.0;
    int c,i;
    // field from other beads
    for (c = 0; c < local_cells.n; c++) {
		// Iterate over all particles in this cell
		for(i=0;i<local_cells.cell[c]->n;i++) {
			// obtain field and gradients from magnetic charge sheets
			p = &local_cells.cell[c]->part[i];
			memcpy(p22, p->r.p, 3*sizeof(double));
			memcpy(img, p->l.i, 3*sizeof(int));
			fold_position(p22, img);
			// start with adding fields, when particle is the same in the loop
			//if(p->r.p[0]==p1->r.p[0] && p->r.p[1]==p1->r.p[1] && p->r.p[2]==p1->r.p[2]) {
			if(p11[0]==p22[0] && p11[1]==p22[1] && p11[2]==p22[2]) {
				//flag=1;
				//continue;
			}
			else{
			//if(flag==1){
				// Distance between particles
				get_mi_vector(dr,p1->r.p,p->r.p);
				r2=dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
				r=sqrt(r2);
				r3=r2*r;
				r5=r3*r2;
				ff=1.0/(4.0*M_PI*r5);
				m2r=p->r.dip[0]*dr[0]+p->r.dip[1]*dr[1]+p->r.dip[2]*dr[2];
				m1r=p1->r.dip[0]*dr[0]+p1->r.dip[1]*dr[1]+p1->r.dip[2]*dr[2];
				
			
				for (int k=0;k<3;k++) {
					h[k] += ff*(3*dr[k]*m2r-p->r.dip[k]*r2);				// Summe h von anderen beads
					h2[k] += ff*(3*dr[k]*m1r-dip[k]*r2);
				}
				
				ch = 3*p->p.susc/(3+p->p.susc);
				vol = 4.188790205 * p->p.radius * p->p.radius * p->p.radius; // vol = 4/3*pi*r^3
				chi = ch * vol / mu0;  // field values in B units, hence dividing by mu0
				p->r.dip[0]+= h2[0]*chi;  // add dipolar field to other particle; h equals the B field
				p->r.dip[1]+= h2[1]*chi;
				p->r.dip[2]+= h2[2]*chi;
			}
		}
	  }
	  
	ch = 3*p1->p.susc/(3+p1->p.susc);
	vol = 4.188790205 * p1->p.radius * p1->p.radius * p1->p.radius; // vol = 4/3*pi*r^3
	chi = ch * vol / mu0; 
	 
    //chi = ch*p1->p.susc;
    
    // m = vol*M = vol*ch*H --> ch = 3*susc/(3+susc) due to sphere
    dip[0]= h[0]*chi;  // h equals the B field
	dip[1]= h[1]*chi;
    dip[2]= h[2]*chi;
    
    p1->r.dip[0] += dip[0]; 
    p1->r.dip[1] += dip[1];
    p1->r.dip[2] += dip[2];
    
     DIP2=p1->r.dip[0]*p1->r.dip[0]+p1->r.dip[1]*p1->r.dip[1]+p1->r.dip[2]*p1->r.dip[2];
	 DIP=sqrt(DIP2);
	 
	//printf("|m| %e ", DIP);
	if(DIP > p1->p.sat){
		for(int j=0;j<3; j++){
			eDIP[j]=p1->r.dip[j]*1/DIP;  //Einheitsvektor
			p1->r.dip[j]=eDIP[j]*p1->p.sat;
		}
	}
    //if((int)time_step%10000==0)printf("m: %e %e %e\n", p1->r.dip[0],p1->r.dip[1],p1->r.dip[2]); 
}
#endif

// Calculates dipolar energy and/or force between two particles
double calc_dipole_dipole_ia(Particle* p1, Particle *p2, int force_flag)
{
  double u,r,pe1,pe2,pe3,pe4,r3,r5,r2,r7,a,b,cc,d,ab;
#ifdef ROTATION
  double bx,by,bz,ax,ay,az; 
#endif
  double  ffx,ffy,ffz;
  double dr[3];
  
#ifdef EXCLUDED_VOLUME_FORCE
  double ff, m,m1, m2,mm,aa,aa2,aa4,ex, fffx,fffy,fffz, ffff[3],ffff2,ffff3,fi[3],fe[3],fffffi,fffffe, multiplier;
#endif
	
  // Distance between particles
  get_mi_vector(dr,p1->r.p,p2->r.p);

  // Powers of distance
  r2=dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
  r=sqrt(r2);
  r3=r2*r;
  r5=r3*r2;
  r7=r5*r2;
 
  // Dot products
  pe1=p1->r.dip[0]*p2->r.dip[0]+p1->r.dip[1]*p2->r.dip[1]+p1->r.dip[2]*p2->r.dip[2];
  pe2=p1->r.dip[0]*dr[0]+p1->r.dip[1]*dr[1]+p1->r.dip[2]*dr[2];
  pe3=p2->r.dip[0]*dr[0]+p2->r.dip[1]*dr[1]+p2->r.dip[2]*dr[2];
#ifdef EXCLUDED_VOLUME_FORCE 
// dipole interaction force fom yung1998, additional prefactor to the original dawaanr force 
  pe4=3.0/r5*1e-7; // additional prefactor mu0/(4*pi)=1e-7 
#else
  pe4=3.0/r5
#endif

  // Energy, if requested
  u= coulomb.Dprefactor* ( pe1/r3 -   pe4*pe2*pe3);

  // Force, if requested
  if(force_flag) { 
    a=pe4*pe1;
    b=-15.0*pe2*pe3/r7;
    ab =a+b;
    cc=pe4*pe3;
    d=pe4*pe2;

    //  Result
    ffx=ab*dr[0]+cc*p1->r.dip[0]+d*p2->r.dip[0];
    ffy=ab*dr[1]+cc*p1->r.dip[1]+d*p2->r.dip[1];
    ffz=ab*dr[2]+cc*p1->r.dip[2]+d*p2->r.dip[2];
    
#ifdef EXCLUDED_VOLUME_FORCE
// mu0 is not fully implemented in calc_dipole_dipole_ia, take care in excluded volume force calculatoin

// Krishnamurthy2008 and Melle2003, same magnetization for both dipoles
#ifdef EXCLUDED_VOLUME_FORCE_V1
	aa = p1->p.radius + p2->p.radius;
	if(r < (evf_cut + aa)){
		m=p1->p.dipm;
		m2=m*m;

		aa2=aa*aa;
		aa4=aa2*aa2;

		ff = evf_A * m2 / aa4;
		ex = exp(- evf_xi * (r / r12 - 1));

		ffx += ff * dr[0] * ex;
		ffy += ff * dr[1] * ex;
		ffz += ff * dr[2] * ex;
	}
#endif

// adapted to yung1998 formulation
#ifdef EXCLUDED_VOLUME_FORCE_V2
	aa = p1->p.radius + p2->p.radius;
	if(r < (evf_cut + aa)){
		//m1=p1->p.dipm;
		//m2=p2->p.dipm;
		//mm=m1*m2;

		aa2=aa*aa;
		aa4=aa2*aa2;

		ff = evf_A * 3.0e-7 * pe1 / (aa4*aa);
		ex = exp(- evf_xi * (r / aa - 1.0));

		//printf("dawaanr dist = %e ; fi = %e, fe = %e, sum %e\n", r, ffy, ff * dr[1] * e, ffy + ff * dr[1] * e);
		//printf("dawaanr dist = %e, fi = %e %e %e, fe = %e %e %e, diff %e %e %e \n", r, ffx, ffy, ffz,ff * dr[0] * ex,ff * dr[1] * ex,ff * dr[2] * ex, ffx - ff * dr[0] * ex, ffy - ff * dr[1] * ex, ffz- ff * dr[2] * ex );
		//if(r < 1.2 && r > 0.8)
			printf("dawaanr dist = %e, fiy = %e, fey = %e, diffy %e \n", r, ffy,ff * dr[1] * ex, ffy - ff * dr[1] * ex );

		ffx += ff * dr[0] * ex;
		ffy += ff * dr[1] * ex;
		ffz += ff * dr[2] * ex;

		//fe[0]=fffx;
		//fe[1]=fffy;
		//fe[2]=fffz;

		//get_mi_vector(ffff,fe,fi);
		//ffff2=ffff[0]*ffff[0]+ffff[1]*ffff[1]+ffff[2]*ffff[2];
		//ffff3=sqrt(ffff2);

		//fffffe=sqrt(fffx*fffx+fffy*fffy+fffz*fffz);
		//if(fffffi>fffffe){
		//	if((int)(sim_time*10000)%300==0)printf("|fe| %2e |sum i+e| %.2e\n", fffffe,1.0*sqrt((ffx+fffx)*(ffx+fffx)+(ffy+fffy)*(ffy+fffy)+(ffz+fffz)*(ffz+fffz)));
		//}else{
		//	if((int)(sim_time*10000)%300==0)printf("|fe| %2e |sum i+e| %.2e\n", fffffe,-1.0*sqrt((ffx+fffx)*(ffx+fffx)+(ffy+fffy)*(ffy+fffy)+(ffz+fffz)*(ffz+fffz)));
		//}

		//ffx += fffx;
		//ffy += fffy;
		//ffz += fffz;
		//if((int)(sim_time*1000)%300==0)printf("fe = %.2e %.2e %.2e |fe| %2e sum %.2e %.2e %.2e |sum| %.2e\n", fffx, fffy, fffz, sqrt(fffx*fffx+fffy*fffy+fffz*fffz),ffx,ffy,ffz,sqrt(ffx*ffx+ffy*ffy+ffz*ffz));


		//printf("fe = %e\n", ff * dr[1] * e);
	}
#endif

// adaptive EVF, reverse magnetic interaction force, should work similar than exluded volume force
#ifdef EXCLUDED_VOLUME_FORCE_V3
	aa = p1->p.radius + p2->p.radius;
	if(r < (evf_cut + aa)){
		multiplier=1.0-((r-aa)/evf_cut);
		//if(multiplier>0.95 && multiplier < 1.05) multiplier = 1.0;

		//if((int)(sim_time*100)%300==0)printf("%d d = %.2e ; fi = %.2e %.2e %.2e |fi| %.2e; multi %.2e fe = %.2e %.2e %.2e |fe| %.2e |diff| %.2e\n", (int)(sim_time*10000), rr, ffx, ffy, ffz, sqrt(ffx*ffx+ffy*ffy+ffz*ffz), multiplier, ffx*multiplier, ffy*multiplier, ffz*multiplier, sqrt((ffx*multiplier)*(ffx*multiplier)+(ffy*multiplier)*(ffy*multiplier)+(ffz*multiplier)*(ffz*multiplier)), sqrt(ffx*ffx+ffy*ffy+ffz*ffz)-sqrt((ffx*multiplier)*(ffx*multiplier)+(ffy*multiplier)*(ffy*multiplier)+(ffz*multiplier)*(ffz*multiplier)));
		//if(multiplier>0.95 && multiplier < 1.05 && (int)(sim_time*10)%100==0)printf("%d d = %.2e ; fi = %.2e %.2e %.2e |fi| %.2e; multi %.2e fe = %.2e %.2e %.2e |fe| %.2e\n", (int)(sim_time*10000), rr, ffx, ffy, ffz, sqrt(ffx*ffx+ffy*ffy+ffz*ffz), multiplier, ffx*multiplier, ffy*multiplier, ffz*multiplier, sqrt((ffx*multiplier)*(ffx*multiplier)+(ffy*multiplier)*(ffy*multiplier)+(ffz*multiplier)*(ffz*multiplier)));

		//if(multiplier>0.95 && multiplier < 1.05){ //set force 0 in direct contact to avoid wobbling
		//	fe[0] = ffx;
		//	fe[1] = ffy;
		//	fe[2] = ffz;
		//}
		//else{
			ffx += (ffx*multiplier);
			ffy += (ffy*multiplier);
			ffz += (ffz*multiplier);
		//}
		/*ffx -= fe[0];
		ffy -= fe[1];
		ffz -= fe[2];*/
		//if((int)(sim_time*10000)%300==0)printf("%d d = %.2e ; fi = %.2e %.2e %.2e |fi| %.2e; multi %.2e fe = %.2e %.2e %.2e |fe| %.2e |diff| %.2e\n", (int)(sim_time*10000), rr, fi[0], fi[1],fi[2], sqrt(fi[0]*fi[0]+fi[1]*fi[1]+fi[2]*fi[2]), multiplier, fe[0], fe[1], fe[2], sqrt((fe[0])*(fe[0])+(fe[1])*(fe[1])+(fe[2])*(fe[2])),sqrt(fi[0]*fi[0]+fi[1]*fi[1]+fi[2]*fi[2])-sqrt((fe[0])*(fe[0])+(fe[1])*(fe[1])+(fe[2])*(fe[2])));

		//ffx -= (ffx*multiplier*dr[0]/rr);
		//ffy -= (ffy*multiplier*dr[1]/rr);
		//ffz -= (ffz*multiplier*dr[2]/rr);

		//if((int)(sim_time*1000)%100==0)printf("multi %.2e fe = %.2e %.2e %.2e |fe| %.2e\n", multiplier, ffx*multiplier, ffy*multiplier, ffz*multiplier, ffz, sqrt((ffx*multiplier)*(ffx*multiplier)+(ffy*multiplier)*(ffy*multiplier)+(ffz*multiplier)*(ffz*multiplier)));
	}
#endif
#endif
    
    // Add the force to the particles
    p1->f.f[0] +=coulomb.Dprefactor*ffx;
    p1->f.f[1] +=coulomb.Dprefactor*ffy;
    p1->f.f[2] +=coulomb.Dprefactor*ffz;
    p2->f.f[0] -=coulomb.Dprefactor*ffx;
    p2->f.f[1] -=coulomb.Dprefactor*ffy;
    p2->f.f[2] -=coulomb.Dprefactor*ffz;
//    if (p1->p.identity==248)
//    {
//      printf("xxx %g %g %g\n", dr[0],dr[1],dr[2]);
//      printf("%d %g %g %g - %g %g %g\n",p2->p.identity,ffx,ffy,ffz,p2->r.p[0],p2->r.p[1],p2->r.p[2]);
//     }
//    if (p2->p.identity==248)
 //   {
//      printf("xxx %g %g %g\n", dr[0],dr[1],dr[2]);
//      printf("%d %g %g %g - %g %g %g\n",p1->p.identity,-ffx,-ffy,-ffz,p1->r.p[0],p1->r.p[1],p1->r.p[2]);
//     }

      // Torques
#ifdef ROTATION
    ax=p1->r.dip[1]*p2->r.dip[2]-p2->r.dip[1]*p1->r.dip[2];
    ay=p2->r.dip[0]*p1->r.dip[2]-p1->r.dip[0]*p2->r.dip[2];
    az=p1->r.dip[0]*p2->r.dip[1]-p2->r.dip[0]*p1->r.dip[1];
    
    bx=p1->r.dip[1]*dr[2]-dr[1]*p1->r.dip[2];
    by=dr[0]*p1->r.dip[2]-p1->r.dip[0]*dr[2];
    bz=p1->r.dip[0]*dr[1]-dr[0]*p1->r.dip[1];
    
    p1->f.torque[0]+=coulomb.Dprefactor*(-ax/r3+bx*cc);
    p1->f.torque[1]+=coulomb.Dprefactor *(-ay/r3+by*cc);
    p1->f.torque[2]+=coulomb.Dprefactor*(-az/r3+bz*cc);
    
    // 2nd particle     
    bx=p2->r.dip[1]*dr[2]-dr[1]*p2->r.dip[2];
    by=dr[0]*p2->r.dip[2]-p2->r.dip[0]*dr[2];
    bz=p2->r.dip[0]*dr[1]-dr[0]*p2->r.dip[1];
	     
    p2->f.torque[0]+=coulomb.Dprefactor* (ax/r3+bx*d);
    p2->f.torque[1]+=coulomb.Dprefactor*(ay/r3+by*d);
    p2->f.torque[2]+=coulomb.Dprefactor* (az/r3+bz*d);
#endif
  }    
	
  // Return energy
  return u;
}

/* =============================================================================
                  DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA                
   =============================================================================
*/

double dawaanr_calculations(int force_flag, int energy_flag)
{
  double u; 
  int i,j,c,cc;
  
  // HARDCODED position, size and magnetization of magnetic charge sheets
  #ifdef MAGNETIC_CHARGE_SHEETS
  double pos[3]={50.,10.,10.}; // center of channel
  double size[3]={10,10,1000}; // charge sheets with dimensions 20x20, located outside the channel z=+-1000 Lm
  double sigma=1.0;
  #endif
  
  	//double magX = 50., magY = 10., magZ = 10.; // center of channel
	/*if(time_step < 1000000) magX = 20;
	if(time_step < 2000000) magX = 50;
	if(time_step < 3000000) magX = 70;*/
	//double a = 10, b = 10, c = 1e8; // charge sheets with dimensions 10x10, poles outside the channel
  
  
  if(n_nodes!=1) {fprintf(stderr,"error:  DAWAANR is just for one cpu .... \n"); errexit();}
  if(!(force_flag) && !(energy_flag) ) {fprintf(stderr," I don't know why you call dawaanr_caclulations with all flags zero \n"); return 0;}
  
  // calculate field from charge sheets and set dip in particles
  #ifdef MAGNETIC_CHARGE_SHEETS
  // Iterate over all cells
  for (c = 0; c < local_cells.n; c++) {
    // Iterate over all particles in this cell
    for(i=0;i<local_cells.cell[c]->n;i++) {
      // obtain field from magnetic charge sheets
		if( local_cells.cell[c]->part[i].p.sat == 0 ) 
			continue;
       calculate_gradient_field(&local_cells.cell[c]->part[i], pos, size, sigma);
	}
  }
  #endif
  
  // calculate field from surrounding dipoles and add it to dip in particles
  #ifdef SOFTMAGNETIC
  // Iterate over all cells
  for (c = 0; c < local_cells.n; c++) {
    // Iterate over all particles in this cell
    for(i=0;i<local_cells.cell[c]->n;i++) {
      // obtain field from surrounding dipoles
      	if( local_cells.cell[c]->part[i].p.sat == 0 ) 
			continue;
       calculate_dipolar_field(&local_cells.cell[c]->part[i]);
	}
  }
  #endif
  
  // calculate gradient force based on dip and gradient field of permanent magnet
  #ifdef MAGNETIC_CHARGE_SHEETS
   // Iterate over all cells
  for (c = 0; c < local_cells.n; c++) {
    // Iterate over all particles in this cell
    for(i=0;i<local_cells.cell[c]->n;i++) {
      // obtain gradient field from magnetic charge sheets
      	if( local_cells.cell[c]->part[i].p.sat == 0 ) 
			continue;
       calculate_moment_force(&local_cells.cell[c]->part[i], pos, size, sigma);
	}
  }
  #endif
  
  // Variable to sum up the energy
  u=0;

  // Iterate over all cells
  for (c = 0; c < local_cells.n; c++) {
    // Iterate over all particles in this cell
    for(i=0;i<local_cells.cell[c]->n;i++) {
      // If the particle has no dipole moment, ignore it
      if( local_cells.cell[c]->part[i].p.dipm == 0.0 ) 
        continue;
      
      // Consider interaction of this particle with the others in the same cell
      for (j=i+1;j<local_cells.cell[c]->n; j++)	{
        // If the particle has no dipole moment, ignore it
        if( local_cells.cell[c]->part[j].p.dipm == 0.0 ) 
          continue;
        // Calculate energy and/or force between the particles
	u+=calc_dipole_dipole_ia(&local_cells.cell[c]->part[i],&local_cells.cell[c]->part[j],force_flag);
      }

      // Calculate the ia between this particles and the particles in the 
      // other cells:
      // Iterate over all remaining cells
      for (cc = c+1; cc < local_cells.n; cc++) {
	// Iterate over the particles in this cell
	for (j=0;j<local_cells.cell[cc]->n;j++) {
	  // If it doesn't have dipole moment, ignore
	  if( local_cells.cell[cc]->part[j].p.dipm == 0.0 ) 
	    continue;
        
	  // Calculate energy and/or force between the particles
	  u+=calc_dipole_dipole_ia(&local_cells.cell[c]->part[i],&local_cells.cell[cc]->part[j],force_flag);
	}
      }
    }
  }
  
  // Return energy
  return u;
}


/************************************************************/

/*=================== */
/*=================== */
/*=================== */
/*=================== */
/*=================== */
/*=================== */

/* =============================================================================
                  DIRECT SUM FOR MAGNETIC SYSTEMS               
   =============================================================================
*/

int  Ncut_off_magnetic_dipolar_direct_sum=0;

/************************************************************/


int  magnetic_dipolar_direct_sum_sanity_checks()
{
  /* left for the future , at this moment nothing to do */
  
  return 0;
}

/************************************************************/


double  magnetic_dipolar_direct_sum_calculations(int force_flag, int energy_flag) {
  Cell *cell;
  Particle *part;
  int i,c,np;
  double *x=NULL,  *y=NULL, *z=NULL;
  double *mx=NULL,  *my=NULL, *mz=NULL;
  double *fx=NULL,  *fy=NULL, *fz=NULL;
#ifdef ROTATION
  double *tx=NULL,  *ty=NULL, *tz=NULL;
#endif
  int dip_particles,dip_particles2;
  double ppos[3];
  int img[3];
  double u;

  
  if(n_nodes!=1) {fprintf(stderr,"error: magnetic Direct Sum is just for one cpu .... \n"); errexit();}
  if(!(force_flag) && !(energy_flag) ) {fprintf(stderr," I don't know why you call dawaanr_caclulations with all flags zero \n"); return 0;}

  x = (double *) Utils::malloc(sizeof(double)*n_part);
  y = (double *) Utils::malloc(sizeof(double)*n_part);
  z = (double *) Utils::malloc(sizeof(double)*n_part);

  mx = (double *) Utils::malloc(sizeof(double)*n_part);
  my = (double *) Utils::malloc(sizeof(double)*n_part);
  mz = (double *) Utils::malloc(sizeof(double)*n_part);
 
  if(force_flag) {
    fx = (double *) Utils::malloc(sizeof(double)*n_part);
    fy = (double *) Utils::malloc(sizeof(double)*n_part);
    fz = (double *) Utils::malloc(sizeof(double)*n_part);
 
#ifdef ROTATION   
    tx = (double *) Utils::malloc(sizeof(double)*n_part);
    ty = (double *) Utils::malloc(sizeof(double)*n_part);
    tz = (double *) Utils::malloc(sizeof(double)*n_part);
#endif  
  }
 
  dip_particles=0;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np   = cell->n;
    for(i=0;i<np;i++) {
      if( part[i].p.dipm != 0.0 ) {
       
	mx[dip_particles]=part[i].r.dip[0];
	my[dip_particles]=part[i].r.dip[1];
	mz[dip_particles]=part[i].r.dip[2];

	/* here we wish the coordinates to be folded into the primary box */
                  
	ppos[0]=part[i].r.p[0];	  
	ppos[1]=part[i].r.p[1];	  
	ppos[2]=part[i].r.p[2];	  
	img[0]=part[i].l.i[0];	  
	img[1]=part[i].l.i[1];	  
        img[2]=part[i].l.i[2];	  		  
        fold_position(ppos, img);
	 
	x[dip_particles]=ppos[0];
	y[dip_particles]=ppos[1];
	z[dip_particles]=ppos[2];
	 

	if(force_flag) {
	  fx[dip_particles]=0;
	  fy[dip_particles]=0;
	  fz[dip_particles]=0;
	 
#ifdef ROTATION
	  tx[dip_particles]=0;
	  ty[dip_particles]=0;
	  tz[dip_particles]=0;
#endif
	} 
	 
	dip_particles++;
	  	 
      }
    }
  }
 
  /*now we do the calculations */
   
  { /* beginning of the area of calculation */
    int nx,ny,nz,i,j;
    double r,rnx,rny,rnz,pe1,pe2,pe3,r3,r5,r2,r7;
    double a,b,c,d;
#ifdef ROTATION
    double ax,ay,az,bx,by,bz;
#endif
    double rx,ry,rz;
    double rnx2,rny2;
    int NCUT[3],NCUT2;
	 
    for(i=0;i<3;i++){
      NCUT[i]=Ncut_off_magnetic_dipolar_direct_sum;
      if(PERIODIC(i) == 0)  {NCUT[i]=0;}  
    }
    NCUT2=Ncut_off_magnetic_dipolar_direct_sum*Ncut_off_magnetic_dipolar_direct_sum;
	     
	   
    u=0;
     
    fprintf(stderr,"Magnetic Direct sum takes long time. Done of %d: \n",dip_particles);

    for(i=0;i<dip_particles;i++){
      fprintf(stderr,"%d\r",i);
      for(j=0;j<dip_particles;j++){
	pe1=mx[i]*mx[j]+my[i]*my[j]+mz[i]*mz[j];
	rx=x[i]-x[j];
	ry=y[i]-y[j];
	rz=z[i]-z[j];
           
	for(nx=-NCUT[0];nx<=NCUT[0];nx++){
	  rnx=rx+nx*box_l[0]; 
	  rnx2=rnx*rnx;
	  for(ny=-NCUT[1];ny<=NCUT[1];ny++){
	    rny=ry+ny*box_l[1];
	    rny2=rny*rny;
	    for(nz=-NCUT[2];nz<=NCUT[2];nz++){
	      if( !(i==j && nx==0 && ny==0 && nz==0) ) {
		if(nx*nx+ny*ny +nz*nz<=  NCUT2){
		  rnz=rz+nz*box_l[2]; 
		  r2=rnx2+rny2+rnz*rnz;
		  r=sqrt(r2);
		  r3=r2*r;
		  r5=r3*r2;
		  r7=r5*r2;
			    
   
		  pe2=mx[i]*rnx+my[i]*rny+mz[i]*rnz;
		  pe3=mx[j]*rnx+my[j]*rny+mz[j]*rnz;
    
		  /*fprintf(stderr,"--------------------------------\n");
		    fprintf(stderr,"ij: %d %d\n",i,j);
		    fprintf(stderr,"xyz[i]: %lf %lf %lf\n",x[i],y[i],z[i]);
		    fprintf(stderr,"xyz[j]: %lf %lf %lf\n",x[j],y[j],z[j]);
		    fprintf(stderr,"mu xyz[i]: %lf %lf %lf\n",mx[i],my[i],mz[i]);
		    fprintf(stderr,"mu xyz[j]: %lf %lf %lf\n",mx[j],my[j],mz[j]);
		    fprintf(stderr,"rnxyz:  %lf %lf %lf\n",rnx,rny,rnz);
		    fprintf(stderr,"--------------------------------\n");*/
 
    
		  //Energy ............................   
    
		  u+= pe1/r3 - 3.0*pe2*pe3/r5;
	     
		  if(force_flag) {
		    //force ............................
		    a=mx[i]*mx[j]+my[i]*my[j]+mz[i]*mz[j];
		    a=3.0*a/r5;
		    b=-15.0*pe2*pe3/r7;
		    c=3.0*pe3/r5;
		    d=3.0*pe2/r5;
	     
		    fx[i]+=(a+b)*rnx+c*mx[i]+d*mx[j];
		    fy[i]+=(a+b)*rny+c*my[i]+d*my[j];
		    fz[i]+=(a+b)*rnz+c*mz[i]+d*mz[j];
	     
#ifdef ROTATION
		    //torque ............................
		    c=3.0/r5*pe3;
		    ax=my[i]*mz[j]-my[j]*mz[i];
		    ay=mx[j]*mz[i]-mx[i]*mz[j];
		    az=mx[i]*my[j]-mx[j]*my[i];
	     
		    bx=my[i]*rnz-rny*mz[i];
		    by=rnx*mz[i]-mx[i]*rnz;
		    bz=mx[i]*rny-rnx*my[i];
	     
		    tx[i]+=-ax/r3+bx*c;
		    ty[i]+=-ay/r3+by*c;
		    tz[i]+=-az/r3+bz*c;
#endif	  
		  } /* of force_flag  */
	     
		} }/* of nx*nx+ny*ny +nz*nz< NCUT*NCUT   and   !(i==j && nx==0 && ny==0 && nz==0) */
	    }/* of  for nz */
          }/* of  for ny  */
        }/* of  for nx  */
      }}   /* of  j and i  */ 


    fprintf(stderr,"done \n");
  }/* end of the area of calculation */
    
    
    
  /* set the forces, and torques of the particles within Espresso */
  if(force_flag) {
   
    dip_particles2=0;
    for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c];
      part = cell->part;
      np   = cell->n;
      for(i=0;i<np;i++) {
	if( part[i].p.dipm != 0.0 ) {
	 
	  part[i].f.f[0]+=coulomb.Dprefactor*fx[dip_particles2];
	  part[i].f.f[1]+=coulomb.Dprefactor*fy[dip_particles2];
	  part[i].f.f[2]+=coulomb.Dprefactor*fz[dip_particles2];
	
#ifdef ROTATION 
	  part[i].f.torque[0]+=coulomb.Dprefactor*tx[dip_particles2];
	  part[i].f.torque[1]+=coulomb.Dprefactor*ty[dip_particles2];
	  part[i].f.torque[2]+=coulomb.Dprefactor*tz[dip_particles2];
#endif
	  dip_particles2++;
	  	 
	}
      }
    }
   
    /* small checking */
    if(dip_particles != dip_particles2) { fprintf(stderr,"magnetic direct sum calculations: error mismatch of particles \n"); errexit();}
  } /*of if force_flag */
  
  /* free memory used */

  free(x);
  free(y);
  free(z);
  free(mx);
  free(my);
  free(mz);
 
  if(force_flag) {
    free(fx);
    free(fy);
    free(fz);
#ifdef ROTATION
    free(tx);
    free(ty);
    free(tz);
#endif
  }
 
  return 0.5*coulomb.Dprefactor*u;
} 
 
int dawaanr_set_params()
{
  if (n_nodes > 1) {
    return ES_ERROR;
  }
  if (coulomb.Dmethod != DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA ) {
    set_dipolar_method_local(DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA);
  } 
  // also necessary on 1 CPU, does more than just broadcasting
  mpi_bcast_coulomb_params();

  return ES_OK;
}

#ifdef EXCLUDED_VOLUME_FORCE
int dawaanr_set_params_evf(double A, double xi, double cut)
{
  if (n_nodes > 1) {
    return ES_ERROR;
  }

  evf_A = A;
  evf_xi = xi;
  evf_cut = cut;

  if (coulomb.Dmethod != DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA ) {
    coulomb.Dmethod = DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA;
  }
  // also necessary on 1 CPU, does more than just broadcasting
  mpi_bcast_coulomb_params();

  return ES_OK;
}
#endif

int mdds_set_params(int n_cut)
{
  if (n_nodes > 1) {
    return ES_ERROR;  
  }
  
  Ncut_off_magnetic_dipolar_direct_sum = n_cut;
  
  if (Ncut_off_magnetic_dipolar_direct_sum == 0) {
    fprintf(stderr,"Careful:  the number of extra replicas to take into account during the direct sum calculation is zero \n");
  }
  
  if (coulomb.Dmethod != DIPOLAR_DS  && coulomb.Dmethod != DIPOLAR_MDLC_DS) {
    set_dipolar_method_local(DIPOLAR_DS);
  }  
  
  // also necessary on 1 CPU, does more than just broadcasting
  mpi_bcast_coulomb_params();
  return ES_OK;
}

#endif

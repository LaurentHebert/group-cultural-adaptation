/**
 * @file   tevol_source_gca.hpp
 * @brief  ODE system for group SIS with group selection mechanism for transmission rates
 *
 * Source code. All parameters passed as arguments, but specify and compile to change precision or output format.
 * g++ -std=c++11 -O3 -o tevol_source_gca ./tevol_source_gca.cpp $(gsl-config --cflags) $(gsl-config --libs)
 *
 * @author  LHD
 * @since   2021-02-06
 */

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>

#include <boost/multi_array.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "dyn_gca.hpp"

double binomial(int n, int k, double p);

using namespace std;

int main(int argc, const char *argv[]) {
		 
	//Model parameters	
	double beta = atof(argv[1]); //basic diffusion rate of innovation 
	double gamma = atof(argv[2]); //recovery rate of innovation
	double rho = atof(argv[3]); //coupling between groups
	double b = atof(argv[4]); //collective benefit per adopter
	double c = atof(argv[5]); //cost per level of promotion
	double mu = atof(argv[6]); //mutation rate of adoption
	double epsilon = atof(argv[7]); //random initial condition
	int n = atoi(argv[8]); //group size
	int L = atoi(argv[9]); //number of levels
    bool timeseries = atoi(argv[10]); //1 for time series, 0 for final state

    if(argc<6) {cerr << "Requires bunch of parameters: basic diffusion rate (beta),\n recovery rate (gamma),"
                    << "\n group coupling (rho),\n adaptation benefit (b),\n promotion cost (c),"
                    << "\n initial adoption (epsilon),\n group size (n),\n number of levels (l_max), \n 1 for time series, 0 for final state.\n" << endl; return 0;}

    const int dim1 = L+1;
    const int dim2 = n+1;
    Sparam param = {beta, gamma, rho, b, c, mu, dim1, dim2};

    // Integrator parameters
    double t = 0;
    double dt = 1e-4;
    double t_step = 1.0;
    const double eps_abs = 1e-8;
    const double eps_rel = 1e-6;

    // Setting initial conditions
    typedef boost::multi_array<double,2> mat_type;
    typedef mat_type::index index;
    mat_type y(boost::extents[dim1][dim2]);
    fill(y.data(),y.data()+y.num_elements(),0.0);

    // Uniform initial conditions
    double lastI = 0.0;
    double lastIL = 0.0;
	for(int l=0; l<dim1; ++l) {
        for(int i=0; i<=n; ++i) {
            y[l][i] = (1.0/(1.0*(L+1)))*binomial(n,i,epsilon); //C_n,i
            lastI += 1.0*i*y[l][i];
        }
    }
    double newI = 0.0;
    double newIL = 0.0;
    double newS = 1.0;

    // Define GSL odeiv parameters
    const gsl_odeiv_step_type * step_type = gsl_odeiv_step_rkf45;
    gsl_odeiv_step * step = gsl_odeiv_step_alloc (step_type, dim1*dim2);
    gsl_odeiv_control * control = gsl_odeiv_control_y_new (eps_abs,eps_rel);
    gsl_odeiv_evolve * evolve = gsl_odeiv_evolve_alloc (dim1*dim2);
    gsl_odeiv_system sys = {dydt, NULL, dim1*dim2, &param};
	
	//Integration
    int status(GSL_SUCCESS);
    double diff = 1.0;
    double diffL = 1.0;
    for (double t_target = t+t_step; t_target < 10000 || diff > 1e-10; t_target += t_step ) { //stop by time and difference
        while (t < t_target) {
            status = gsl_odeiv_evolve_apply (evolve,control,step,&sys,&t,t_target,&dt,y.data());
            if (status != GSL_SUCCESS) {
				cout << "SNAFU" << endl;
                break;
			}
        } // end while

        //measure adoption
        vector<double> newIvec(dim1, 0.0); vector<double> popvec(dim1, 0.0); newI=0.0;
        for(int l=0; l<dim1; ++l) for(int i=0; i<=n; ++i) {
            newIvec[l] += 1.0*i*y[l][i]; popvec[l] += y[l][i]; newI += 1.0*i*y[l][i]; }
        
        //timeseries output
        if(timeseries) {   
          cout << t;
		  for(int l=0; l<dim1; ++l) cout << " " << newIvec[l]/popvec[l] << " " << popvec[l];
          cout << "\n";
        }
        diff = abs(newIvec[0]/popvec[0] - lastI);
        diffL = abs(newIvec[dim1-1]/popvec[dim1-1] - lastIL);
        lastI=newIvec[0]/popvec[0];
        lastIL=newIvec[dim1-1]/popvec[dim1-1];
	} //end while

    //calculate fitnesses
    vector<double> Zvec_(dim1, 0.0); vector<double> popvec_(dim1, 0.0);
    for(int l=0; l<dim1; ++l) {
        for(int i=0; i<dim2; ++i) {
            Zvec_[l] += exp(b*i-c*l)*y[l][i]; //should this really be proportional to C_i,l?
            popvec_[l] += y[l][i];
        }
        Zvec_[l] /= popvec_[l];
    }
    double Z_ = accumulate(Zvec_.begin(), Zvec_.end(), 0.0);

    //final state output
    
    if(!timeseries) {
      vector<double> Ivec(dim1, 0.0); vector<double> pop(dim1, 0.0); newI=0.0;
      for(int l=0; l<dim1; ++l) for(int i=0; i<=n; ++i) {
        Ivec[l] += 1.0*i*y[l][i]; pop[l] += y[l][i]; newI += 1.0*i*y[l][i];
      }
      cout << beta << " " << b << " " << c << " " << rho;
      for(int l=0; l<dim1; ++l) cout << " " << Ivec[l]/pop[l] << " " << pop[l];
      for(int l=0; l<dim1; ++l) cout << " " << Zvec_[l]/Z_;
      cout << "\n";
    }
    

    //output fitness distribution
    /*for(int l=0; l<dim1; ++l) {
        cout << beta << " " << b << " " << c << " " << rho << " " << l;
        for(int i=0; i<dim2; ++i) {
            int g = 0;
            for(; g<int(1000*y[l][i]); ++g) cout << " " << exp(b*i-c*l);
            for(; g<1000; ++g) cout << " " << "nan";
        }
        cout << "\n";
    }*/
    

    // localization in cliques
    /*double xi = 0.0;
    double It = 0.0;
	for(int s=2; s<dim; ++s) for(int i=0; i<=s; ++i) It += 1.0*i*y[s][i];
    for(int s=2; s<dim; ++s) {
        double In = 0.0;
        for(int i=0; i<=s; ++i) In += 1.0*i*y[s][i];
        xi += (In/It)*(In/It);
    }*/
    
    /*cout << beta << " " << lastI;
    for(int n=dim-1; n>1; n-=10) {
        lastI=0.0;
        for(int i=0; i<=n; ++i) lastI += 1.0*norms*i*y[n][i]/(1.0*n*pow(1.0*n,-gammas));
	    cout << " " << lastI;
    }
    cout << "\n";
	*/
    cout.flush();

    // Free memory
    gsl_odeiv_evolve_free(evolve);
    gsl_odeiv_control_free(control);
    gsl_odeiv_step_free(step);
    
    return 0;
}

double binomial(int n, int k, double p)
	{
	//check
	if(k>n || k<0 || n<0) return 0.0;
		
	double result=0.L;
	double postfix = pow(p,k)*pow(1.L-p,n-k);
	if(n<2*k) k=n-k;
	
	if(k<=100&&n<=200)
		{
		double numerator=1;
		for(int i=0; i<k; i++) numerator = numerator*((double)(n-i));
		double denominator=1;
		for(int i=2; i<=k; i++) denominator = denominator*((double)i);
		result = (double)(numerator/denominator)*postfix;
		}
	else
		{
		if(k<100||n<300)
			{
			result=1.L;
			for(int i=1; i<=k; i++) result = result*((double)(n-k+i))/((double)i);
			result = result*postfix;
			}
		else if(n<1000)
			{
			double nd=(double)n, kd=(double)k, nkd=(double)(n-k);
			result = nd*log(nd)-nd+0.5L*log(2.L*M_PI*nd)+1.L/(12.L*nd)-1.L/(360.L*pow(nd,3))+1.L/(1260.L*pow(nd,5));
			result -= kd*log(kd)-kd + 0.5L*log(2.L*M_PI*kd) + 1.L/(12.L*kd) - 1.L/(360.L*pow(kd,3)) + 1.L/(1260.L*pow(kd,5));
			result -= nkd*log(nkd)-nkd + 0.5L*log(2.L*M_PI*nkd) + 1.L/(12.L*nkd) - 1.L/(360.L*pow(nkd,3)) + 1.L/(1260.L*pow(nkd,5));
			result = exp(result);
			result = result*postfix;
			}
			
		else result = (1.L/sqrt(2.L*M_PI*(double)n*p*(1.L-p)))*exp(-pow(k-n*p,2)/(2.L*(double)n*p*(1.L-p)));
		}
	if(result!=result) result=0.L;
	if(isinf(result)) result=0.L;
	return result;
	}

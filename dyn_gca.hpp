#ifndef DYN_GCA_INCLUDED
#define DYN_GCA_INCLUDED

#include <boost/multi_array.hpp>

using namespace std;

/**
 * @file   dyn_gca.hpp
 * @brief  ODE system for group SIS with group selection mechanism for transmission rates
 *
 * @author  LHD
 * @since   2021-02-06
 */

struct Sparam {
    const double beta;
    const double gamma;
    const double rho;
    const double b;
    const double c;
    const double mu;
    const int dim1;
    const int dim2;
}; // parameter structure

//********** function dydt definition **************************************************************
int dydt(double t, const double y[], double f[], void * param) {
// ODE system for interacting contagions

    // Cast parameters
    Sparam& p = *static_cast<Sparam* >(param);

    // Create multi_array reference to y and f
    typedef boost::multi_array_ref<const double,2> CSTmatref_type;
    typedef boost::multi_array_ref<double,2> matref_type;
    typedef CSTmatref_type::index indexref;
    CSTmatref_type yref(y,boost::extents[p.dim1][p.dim2]);
    matref_type fref(f,boost::extents[p.dim1][p.dim2]);

    // Calculate mean-field coupling and observed fitness landscape
    vector<double> Zvec(p.dim1, 0.0); vector<double> popvec(p.dim1, 0.0); double R=0.0;
    for(int l=0; l<p.dim1; ++l) {
        for(int i=0; i<p.dim2; ++i) {
            Zvec[l] += exp(p.b*i-p.c*l)*yref[l][i]; //should this really be proportional to C_i,l?
            popvec[l] += yref[l][i];
            R += p.rho*i*yref[l][i];
        }
        if(popvec[l]>0.0) Zvec[l] /= popvec[l];
    }
    double Z = accumulate(Zvec.begin(), Zvec.end(), 0.0);
    int n=p.dim2-1;
    int L=p.dim1-1;

    // Compute derivatives
    for(int l=0; l<p.dim1; ++l) {
        for(int i=0; i<p.dim2; ++i) {
            fref[l][i] = -1.0*p.gamma*i*yref[l][i] - p.beta*l*(i+R)*(n-i)*yref[l][i];
            if(i>0) fref[l][i] += p.beta*l*((i-1)+R)*(n-i+1)*yref[l][i-1];
            if(i<n) fref[l][i] += p.gamma*(i+1)*yref[l][i+1];
            if(l>0) fref[l][i] += p.rho*yref[l-1][i]*(Zvec[l]/Zvec[l-1]+p.mu) - p.rho*yref[l][i]*(Zvec[l-1]/Zvec[l]+p.mu);
            if(l<L) fref[l][i] += p.rho*yref[l+1][i]*(Zvec[l]/Zvec[l+1]+p.mu) - p.rho*yref[l][i]*(Zvec[l+1]/Zvec[l]+p.mu);
        }
    }

    return GSL_SUCCESS;

} //********** end function dydt definition ********************************************************

#endif // DYN_GCA_HPP_INCLUDED

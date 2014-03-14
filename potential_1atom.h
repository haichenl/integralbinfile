/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef _psi_src_lib_libmints_potential_1Atom_h_
#define _psi_src_lib_libmints_potential_1Atom_h_

#include <vector>
#include <libmints/typedefs.h>

#include <libciomr/libciomr.h>

#include <libmints/mints.h>
#include <libmints/cdsalclist.h>
#include <libmints/potential.h>

#include <physconst.h>

namespace psi {

    class Matrix;
    class BasisSet;
    class GaussianShell;
    class ObaraSaikaTwoCenterVIRecursion;
    class ObaraSaikaTwoCenterVIDerivRecursion;
    class ObaraSaikaTwoCenterVIDeriv2Recursion;
    class OneBodyAOInt;
    class IntegralFactory;
    class SphericalTransform;
    class OneBodySOInt;
    class CdSalcList;

/*! \ingroup MINTS
 *  \class PotentialInt
 *  \brief Computes potential integrals.
 * Use an IntegralFactory to create this object.
 */
class PotentialInt_1Atom : public PotentialInt
{
    /// Computes integrals between two shell objects.
    void compute_pair(const GaussianShell& s1,
                                const GaussianShell& s2)
    {
        int ao12;
        int am1 = s1.am();
        int am2 = s2.am();
        int nprim1 = s1.nprimitive();
        int nprim2 = s2.nprimitive();
        double A[3], B[3];
        A[0] = s1.center()[0];
        A[1] = s1.center()[1];
        A[2] = s1.center()[2];
        B[0] = s2.center()[0];
        B[1] = s2.center()[1];
        B[2] = s2.center()[2];
    
        int izm = 1;
        int iym = am1 + 1;
        int ixm = iym * iym;
        int jzm = 1;
        int jym = am2 + 1;
        int jxm = jym * jym;
    
        // compute intermediates
        double AB2 = 0.0;
        AB2 += (A[0] - B[0]) * (A[0] - B[0]);
        AB2 += (A[1] - B[1]) * (A[1] - B[1]);
        AB2 += (A[2] - B[2]) * (A[2] - B[2]);
    
        memset(buffer_, 0, s1.ncartesian() * s2.ncartesian() * sizeof(double));
    
        double ***vi = potential_recur_->vi();
    
        double** Zxyzp = Zxyz_->pointer();
        int ncharge = Zxyz_->rowspi()[0];
    
        for (int p1=0; p1<nprim1; ++p1) {
            double a1 = s1.exp(p1);
            double c1 = s1.coef(p1);
            for (int p2=0; p2<nprim2; ++p2) {
                double a2 = s2.exp(p2);
                double c2 = s2.coef(p2);
                double gamma = a1 + a2;
                double oog = 1.0/gamma;
    
                double PA[3], PB[3];
                double P[3];
    
                P[0] = (a1*A[0] + a2*B[0])*oog;
                P[1] = (a1*A[1] + a2*B[1])*oog;
                P[2] = (a1*A[2] + a2*B[2])*oog;
                PA[0] = P[0] - A[0];
                PA[1] = P[1] - A[1];
                PA[2] = P[2] - A[2];
                PB[0] = P[0] - B[0];
                PB[1] = P[1] - B[1];
                PB[2] = P[2] - B[2];
    
                double over_pf = exp(-a1*a2*AB2*oog) * sqrt(M_PI*oog) * M_PI * oog * c1 * c2;
    
                // Loop over atoms of basis set 1 (only works if bs1_ and bs2_ are on the same
                // molecule)
                
                // spring:
                int atom=curr_atom;
                //
                
                double PC[3];

                double Z = Zxyzp[atom][0];

                PC[0] = P[0] - Zxyzp[atom][1];
                PC[1] = P[1] - Zxyzp[atom][2];
                PC[2] = P[2] - Zxyzp[atom][3];

                // Do recursion
                potential_recur_->compute(PA, PB, PC, gamma, am1, am2);

                ao12 = 0;
                for(int ii = 0; ii <= am1; ii++) {
                    int l1 = am1 - ii;
                    for(int jj = 0; jj <= ii; jj++) {
                        int m1 = ii - jj;
                        int n1 = jj;
                        /*--- create all am components of sj ---*/
                        for(int kk = 0; kk <= am2; kk++) {
                            int l2 = am2 - kk;
                            for(int ll = 0; ll <= kk; ll++) {
                                int m2 = kk - ll;
                                int n2 = ll;

                                // Compute location in the recursion
                                int iind = l1 * ixm + m1 * iym + n1 * izm;
                                int jind = l2 * jxm + m2 * jym + n2 * jzm;

                                buffer_[ao12++] += -vi[iind][jind][0] * over_pf * Z;

//                                fprintf(outfile, "ao12=%d, vi[%d][%d][0] = %20.14f, over_pf = %20.14f, Z = %f\n", ao12-1, iind, jind, vi[iind][jind][0], over_pf, Z);
                            }
                        }
                    }
                }
                
            }
        }
    }
    /// Computes integrals between two shell objects.
    //void compute_pair_deriv1(const GaussianShell&, const GaussianShell& );
    //void compute_pair_deriv2(const GaussianShell&, const GaussianShell& );

protected:
	  int curr_atom;
    /// Recursion object that does the heavy lifting.
    //ObaraSaikaTwoCenterVIRecursion* potential_recur_;

    /// Matrix of coordinates/charges of partial charges
    //SharedMatrix Zxyz_;

public:
    /// Constructor. Assumes nuclear centers/charges as the potential
    PotentialInt_1Atom(int input_curr_atom, std::vector<SphericalTransform>& st, boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2, int deriv) :
    PotentialInt(st, bs1, bs2, deriv)
    {
    	  curr_atom = input_curr_atom;
    }
    
    virtual ~PotentialInt_1Atom()
    {
    }

    /// Computes the first derivatives and stores them in result
    //virtual void compute_deriv1(std::vector<SharedMatrix > &result);

    /// Computes the second derivatives and store them in result
    //virtual void compute_deriv2(std::vector<SharedMatrix>& result);

    /// Set the field of charges
    //void set_charge_field(SharedMatrix Zxyz) { Zxyz_ = Zxyz; }

    /// Get the field of charges
    //SharedMatrix charge_field() const { return Zxyz_; }

    /// Does the method provide first derivatives?
    //bool has_deriv1() { return true; }
};

}

#endif

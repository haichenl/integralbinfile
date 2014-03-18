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

#ifndef _psi_src_lib_libmints_integral_binfile_h_
#define _psi_src_lib_libmints_integral_binfile_h_

#include <boost/shared_ptr.hpp>
#include <vector>

#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/mints.h>
#include <libint/libint.h>
#include "potential_1atom.h"
#include "potential_ptq.h"

/*! \def INT_NCART(am)
    Gives the number of cartesian functions for an angular momentum.
*/
#define INT_NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)
/*! \def INT_PURE(am)
    Gives the number of spherical functions for an angular momentum.
*/
#define INT_NPURE(am) (2*(am)+1)
/*! \def INT_NFUNC(pu,am)
    Gives the number of functions for an angular momentum based on pu.
*/
#define INT_NFUNC(pu,am) ((pu)?INT_NPURE(am):INT_NCART(am))
/*! \def INT_CARTINDEX(am,i,j)
    Computes offset index for cartesian function.
*/
#define INT_CARTINDEX(am,i,j) (((i) == (am))? 0 : (((((am) - (i) + 1)*((am) - (i)))>>1) + (am) - (i) - (j)))
/*! \def INT_ICART(a, b, c)
    Given a, b, and c compute a cartesian offset.
*/
#define INT_ICART(a, b, c) (((((((a)+(b)+(c)+1)<<1)-(a))*((a)+1))>>1)-(b)-1)
/*! \def INT_IPURE(l, m)
    Given l and m compute a pure function offset.
*/
#define INT_IPURE(l, m) ((l)+(m))

namespace psi {

class BasisSet;
class GaussianShell;
class OneBodyAOInt;
class OneBodySOInt;
class TwoBodyAOInt;
class ThreeCenterOverlapInt;
class Symmetry;
class CartesianIter;
class RedundantCartesianIter;
class RedundantCartesianSubIter;
class ShellRotation;
class SymmetryOperation;
class SOTransform;
class SOBasisSet;
class CorrelationFactor;

/*! \ingroup MINTS */
class IntegralFactory_binfile : public IntegralFactory
{
protected:

public:
    /** Initialize IntegralFactory object given a BasisSet for each center. */
    IntegralFactory_binfile(boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2,
                    boost::shared_ptr<BasisSet> bs3, boost::shared_ptr<BasisSet> bs4) : IntegralFactory(bs1, bs2, bs3, bs4)
    {
    }
    /** Initialize IntegralFactory object given a BasisSet for two centers. Becomes (bs1 bs2 | bs1 bs2). */
    IntegralFactory_binfile(boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2) : IntegralFactory(bs1, bs2)
    {
    }
    /** Initialize IntegralFactory object given a BasisSet for two centers. Becomes (bs1 bs1 | bs1 bs1). */
    IntegralFactory_binfile(boost::shared_ptr<BasisSet> bs1) : IntegralFactory(bs1)
    {
    }

    virtual ~IntegralFactory_binfile()
    {
    }
    
    virtual OneBodyAOInt* ao_potential_1atom(int curr_atom, int deriv=0)
    {
        return new PotentialInt_1Atom(curr_atom, spherical_transforms_, bs1_, bs2_, deriv);
    }
    
    virtual OneBodyAOInt* ao_potential_ptq(SharedVector ptq_zxyz, int deriv=0)
    {
        return new PotentialInt_ptq(ptq_zxyz, spherical_transforms_, bs1_, bs2_, deriv);
    }
    
};

}

#endif

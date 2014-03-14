#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.h>

// spring
#include "inout.h"
#include "integral_1atom.h"
//

INIT_PLUGIN

namespace psi{ namespace integralbinfile {

extern "C"
int read_options(std::string name, Options &options)
{
    if (name == "INTEGRALBINFILE"|| options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /*- Whether to compute two-electron integrals -*/
        options.add_bool("DO_TEI", true);
    }

    return true;
}

extern "C"
PsiReturnType integralbinfile(Options &options)
{
	  // spring:
	  Binary_ofstream integral_bin_file("integral.bin");
	  //
	  
    int print = options.get_int("PRINT");
    int doTei = options.get_bool("DO_TEI");

    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();

    // Form basis object:
    // Create a basis set parser object.
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    // Construct a new basis set.
    boost::shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, molecule, "BASIS");

    // The integral factory oversees the creation of integral objects
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory
            (aoBasis, aoBasis, aoBasis, aoBasis));
    
    // N.B. This should be called after the basis has been built, because the geometry has not been
    // fully initialized until this time.
    molecule->print();
    int nbf[] = { aoBasis->nbf() };
    double nucrep = molecule->nuclear_repulsion_energy();
    fprintf(outfile, "\n    Nuclear repulsion energy: %16.8f\n\n", nucrep);

    // The matrix factory can create matrices of the correct dimensions...
    boost::shared_ptr<MatrixFactory> factory(new MatrixFactory);
    factory->init_with(1, nbf, nbf);

    // Form the one-electron integral objects from the integral factory
    boost::shared_ptr<OneBodyAOInt> sOBI(integral->ao_overlap());
    boost::shared_ptr<OneBodyAOInt> tOBI(integral->ao_kinetic());
    boost::shared_ptr<OneBodyAOInt> vOBI(integral->ao_potential());
    // Form the one-electron integral matrices from the matrix factory
    SharedMatrix sMat(factory->create_matrix("Overlap"));
    SharedMatrix tMat(factory->create_matrix("Kinetic"));
    SharedMatrix vMat(factory->create_matrix("Potential"));
    SharedMatrix hMat(factory->create_matrix("One Electron Ints"));
    // Compute the one electron integrals, telling each object where to store the result
    sOBI->compute(sMat);
    tOBI->compute(tMat);
    vOBI->compute(vMat);

    sMat->print();
    tMat->print();
    vMat->print();
    
    // spring: output
    
    // number of atoms
    int natom = molecule->natom();
    integral_bin_file << natom;
    
    // overlap 
    integral_bin_file << sMat;
    
    // kinetic 
    integral_bin_file << tMat;
    
    // potential
    for(int iatom = 0; iatom < natom; iatom++)
    {
    	  boost::shared_ptr<IntegralFactory_1Atom> integral_1atom(new IntegralFactory_1Atom(iatom, aoBasis, aoBasis, aoBasis, aoBasis));
    	  SharedMatrix vMat_temp(factory->create_matrix("Potential_1Atom"));
        boost::shared_ptr<OneBodyAOInt> v1AtomOBI(integral_1atom->ao_potential());
        v1AtomOBI->compute(vMat_temp);
        integral_bin_file << vMat_temp;
    }

    // Form h = T + V by first cloning T and then adding V
    hMat->copy(tMat);
    hMat->add(vMat);
    hMat->print();

    if(doTei){
         fprintf(outfile, "\n  Two-electron Integrals\n\n");

        // Now, the two-electron integrals
        boost::shared_ptr<TwoBodyAOInt> eri(integral->eri());
        // The buffer will hold the integrals for each shell, as they're computed
        const double *buffer = eri->buffer();
        // The iterator conveniently lets us iterate over functions within shells
        AOShellCombinationsIterator shellIter = integral->shells_iterator();
        int count=0;
        for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
            // Compute quartet
            eri->compute_shell(shellIter);
            // From the quartet get all the integrals
            AOIntegralsIterator intIter = shellIter.integrals_iterator();
            for (intIter.first(); intIter.is_done() == false; intIter.next()) {
                int p = intIter.i();
                int q = intIter.j();
                int r = intIter.k();
                int s = intIter.l();
                fprintf(outfile, "\t(%2d %2d | %2d %2d) = %20.15f\n",
                    p, q, r, s, buffer[intIter.index()]);
                
                // spring:
                integral_bin_file << p << q << r << s;
                integral_bin_file << buffer[intIter.index()];
                //
                
                ++count;
            }
        }
        fprintf(outfile, "\n\tThere are %d unique integrals\n\n", count);
    }

    return Success;
}

}} // End Namespaces

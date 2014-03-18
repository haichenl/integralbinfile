#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.h>

// spring start
#include "integral_binfile.h"
#include "binaryio.h"
// spring end

INIT_PLUGIN

namespace psi{ namespace integralbinfile {

extern "C"
int read_options(std::string name, Options &options)
{
    if (name == "INTEGRALBINFILE"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        /*- Whether to compute two-electron integrals -*/
        options.add_bool("DO_TEI", true);
        /*- Name of the output binary file -*/
        options.add_str("BINFILENAME", "integral.bin");
        /*- Whether to compute internal environment potential integrals -*/
        options.add_bool("DO_ENV_IN", false);
        /*- The number of internal environmental point charges -*/
        options.add_int("NUM_PTQ_IN", 0);
        /*- zxyz array of internal environmental point charges -*/
        options.add("ZXYZ_PTQ_IN", new ArrayType());
        /*- Whether to compute external environment potential integrals -*/
        options.add_bool("DO_ENV_EX", false);
        /*- The number of external environmental point charges -*/
        options.add_int("NUM_PTQ_EX", 0);
        /*- zxyz array of external environmental point charges -*/
        options.add("ZXYZ_PTQ_EX", new ArrayType());
    }
    
    return true;
}

void set_ZxyzMat(SharedMatrix ZxyzMat, boost::shared_ptr<double[]> Zxyz_array, int num_ptq) 
{
    for(int i = 0; i < num_ptq; i++) {
        for(int j = 0; j < 4; j++) {
            ZxyzMat->set(i, j, Zxyz_array[i*4+j]);
        }
    }
}

void write_env(std::string InOrEx, boost::shared_ptr<MatrixFactory> factory, boost::shared_ptr<IntegralFactory_binfile> integral, Binary_ofstream &integral_bin_file, Options &options)
{
    std::string do_env_string = "DO_ENV_";
    do_env_string += InOrEx;
    int doEnv = options.get_bool(do_env_string.c_str());
    integral_bin_file << doEnv;
    if(doEnv) {
        std::string num_ptq_string = "NUM_PTQ_";
        num_ptq_string += InOrEx;
        int num_ptq = options.get_int(num_ptq_string.c_str());
        integral_bin_file << num_ptq;
        std::string zxyz_ptq_string = "ZXYZ_PTQ_";
        zxyz_ptq_string += InOrEx;
        boost::shared_ptr<double[]> Zxyz_array(options.get_double_array(zxyz_ptq_string.c_str()));
        boost::shared_ptr<MatrixFactory> Zxyzfactory(new MatrixFactory);
        SharedMatrix ZxyzMat = Zxyzfactory->create_shared_matrix("Zxyz", num_ptq, 4);
        set_ZxyzMat(ZxyzMat, Zxyz_array, num_ptq);
        for(int i = 0; i < num_ptq; i++) {
            SharedVector curr_ptqvec(ZxyzMat->get_row(0, i));
            boost::shared_ptr<OneBodyAOInt> ptqOBI(integral->ao_potential_ptq(curr_ptqvec));
            SharedMatrix vMat_ptq(factory->create_matrix("Potential_ptq"));
            ptqOBI->compute(vMat_ptq);
            integral_bin_file<<vMat_ptq;
        }
    }
}

extern "C"
PsiReturnType integralbinfile(Options &options)
{
    // spring start 
    std::string binfilename = options.get_str("BINFILENAME");
    Binary_ofstream integral_bin_file(binfilename.c_str());
    // spring end
    
    int print = options.get_int("PRINT");
    int doTei = options.get_bool("DO_TEI");

    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();

    // Form basis object:
    // Create a basis set parser object.
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    // Construct a new basis set.
    boost::shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, molecule, "BASIS");

    // The integral factory oversees the creation of integral objects
    boost::shared_ptr<IntegralFactory_binfile> integral(new IntegralFactory_binfile
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
    
    // spring start: output integrals into the binary file
    
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
        SharedMatrix vMat_temp(factory->create_matrix("Potential_1Atom"));
        boost::shared_ptr<OneBodyAOInt> v1AtomOBI(integral->ao_potential_1atom(iatom));
        v1AtomOBI->compute(vMat_temp);
        integral_bin_file << vMat_temp;
    }
    // spring end

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
                
                // spring start 
                integral_bin_file << p << q << r << s;
                integral_bin_file << buffer[intIter.index()];
                // spring end
                
                ++count;
            }
        }
        fprintf(outfile, "\n\tThere are %d unique integrals\n\n", count);
    }
    
    // spring start 
    // nuclear repulsion energy
    double e_nuc = molecule->nuclear_repulsion_energy();
    integral_bin_file << e_nuc;
    // spring end
    
    // spring start
    // whether we do internal environment calculation 
    write_env("IN", factory, integral, integral_bin_file, options);
    // whether we do external environment calculation 
    write_env("EX", factory, integral, integral_bin_file, options);
    // spring end

    return Success;
}

}} // End Namespaces

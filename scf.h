#include <psi4-dec.h>
#include <libmints/typedefs.h>

namespace psi{ 
    class Options;
    namespace integralbinfile{

class SCF{
  public:
    /// The constuctor
    SCF(Options &options);
    /// The destuctor
    ~SCF();
    /// Computes the SCF energy, and returns it
    double compute_energy();
    SharedMatrix get_orb();
  protected:
    /// The options object, to interact with the input file
    Options &options_;
    /// The amount of information to print to the output file
    int print_;
    /// The number of doubly occupied orbitals
    int ndocc_;
    /// The number of symmetrized spin orbitals
    int nso_;
    /// The maximum number of iterations
    int maxiter_;
    /// The nuclear repulsion energy
    double e_nuc_;
    /// The convergence criterion for the density
    double d_convergence_;
    /// The convergence criterion for the energy
    double e_convergence_;
    /// The one electron integrals
    SharedMatrix H_;
    /// The overlap matrix
    SharedMatrix S_;
    /// The inverse square root of the overlap matrix
    SharedMatrix X_;
    /// The Fock Matrix
    SharedMatrix F_;
    /// The transformed Fock matrix
    SharedMatrix Ft_;
    /// The MO coefficients
    SharedMatrix C_;
    /// The density matrix
    SharedMatrix D_;
    /// A 4d array containing all two electron integrals
    double ****tei_;
    /// Computes the electronic part of the SCF energy, and returns it
    double compute_electronic_energy();
    /// Sets up the integrals object
    void init_integrals();
    /// Forms the density matrix from the MO coefficients
    void form_density();
    /// Initializes a 4 dimensional array, setting the elements to zero
    void init_matrix(double****& matrix, int dim1, int dim2, int dim3, int dim4);
    /// Frees the memory used for a 4D array
    void free_matrix(double****& matrix, int dim1, int dim2, int dim3, int dim4);
};

}} //End namespaces

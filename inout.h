#include <fstream>

using namespace std;

class
Binary_ifstream: private ifstream
{
public:
  ///
#ifdef WIN32
  Binary_ifstream(const char* filename): ifstream(filename, std::ios::binary)
  {
  }
#else
  Binary_ifstream(const char* filename): ifstream(filename)
  {
  }
#endif
  ///
  bool operator!()
  {
    return ifstream::operator!();
  }
  /*
  template<class X>
  Binary_ifstream& operator>>(const X& x)
  {
  read( (char*) &x, sizeof(X));
  return *this;
  }
  */
  ///
  Binary_ifstream& operator>>(const double& x)
  {
    read( (char*) &x, sizeof(double));
    return *this;
  }
  ///
  Binary_ifstream& operator>>(const int& x)
  {
    read( (char*) &x, sizeof(int));
    return *this;
  }
  ///
  Binary_ifstream& operator>>(const unsigned& x)
  {
    read( (char*) &x, sizeof(unsigned));
    return *this;
  }
  ///
  Binary_ifstream& operator>>(const std::complex<double>& x)
  {
    read( (char*) &x, sizeof(std::complex<double>));
    return *this;
  }
  ///
  Binary_ifstream& operator>>(bool& x)
  {
    int iboolean=0;
    read( (char*) &iboolean, sizeof(int));
    x = (iboolean==1)? true : false; 
    return *this;
  }
  /// std::ifstream::read()
  Binary_ifstream& read(char* pch, int nCount)
  {
    ifstream::read(pch, nCount);
    return *this;
  }
  /// set file pointer to the beginning of the file
  /// added Oct. 30, 2001 by Angela
  void reset()
  {
    // FIXME EEM: does not work with g++-3.0
    //    ifstream::rdbuf()->seekpos(0);
  }
  ///
  Binary_ifstream& operator>>(std::string& st)
  {
    int length;
    *this >> length;
    char* ctemp = new char[length+1];
    read(ctemp, length);
    std::string result(ctemp, length);
    st = result;
    return *this;
  }
  ///
  bool good()
  {
    return ifstream::good();
  }
  ///
  void close()
  {
    ifstream::close();
  }
  ///
  void open(const char* szName)
  {
    //void open( const char* szName, int nMode = ios::out, int nProt = filebuf::openprot );
    ifstream::open( szName, std::ios::binary);
  }
};


//@Man:
//@Memo: Binary output file (wrapper of ofstream)
//@Doc:
class
Binary_ofstream: private ofstream
{
public:
  ///
#ifdef WIN32
  Binary_ofstream(const char* filename): ofstream(filename, std::ios::binary)
  {
  }
#else
  Binary_ofstream(const char* filename): ofstream(filename)
  {
  }
#endif
  ///
  bool operator!()
  {
    return ofstream::operator!();
  }
  /* may be too dangerous, since it can output classes that shouldn't be
  stored in this way.
  template<class X>
  Binary_ofstream& operator<<(const X& x)
  {
  write( (char*) &x, sizeof(X));
  return *this;
  }
  */
  ///
  Binary_ofstream& operator<<(const float& x)
  {
    write( (char*) &x, sizeof(float));
    return *this;
  }
  ///
  Binary_ofstream& operator<<(const double& x)
  {
    write( (char*) &x, sizeof(double));
    return *this;
  }
  ///
  Binary_ofstream& operator<<(const int& x)
  {
    write( (char*) &x, sizeof(int));
    return *this;
  }
  ///
  // CML 3-7-14
  Binary_ofstream& operator<<(const long unsigned int& x)
  {
    write( (char*) &x, sizeof(long unsigned int));
    return *this;
  }
  ///
  Binary_ofstream& operator<<(const unsigned& x)
  {
    write( (char*) &x, sizeof(unsigned));
    return *this;
  }
  ///
  Binary_ofstream& operator<<(const std::complex<double>& x)
  {
    write( (char*) &x, sizeof(std::complex<double>));
    return *this;
  }
  ///
  Binary_ofstream& operator<<(const bool& x)
  {
    int ibool = x ? 1 : 0;
    write( (char*) &ibool, sizeof(int));
    return *this;
  }
  /// std::ofstream::write()
  Binary_ofstream& write(const char* pch, int nCount)
  {
    ofstream::write(pch, nCount);
    return *this;
  }
  ///
  Binary_ofstream& operator<<(const std::string& st)
  {
    *this << (int) st.size();
    write(st.c_str(), (int) st.size() );
    return *this;
  }
  Binary_ofstream& operator<<(const psi::SharedMatrix m)
  {
    *this << m->nirrep();
    for (int ih = 0; ih < m->nirrep(); ih++)
      {
      	*this << m->rowdim(ih) << m->coldim(ih);
      	write((char*) m->get_const_pointer(ih),sizeof(double)*m->size(ih));
      }
  }
  /// ofstream::close()
  void close()
  {
    ofstream::close();
  }
  /// ofstream::open( szName, std::ios::binary)
  void open(const char* szName)
  {
    //void open( const char* szName, int nMode = ios::out, int nProt = filebuf::openprot );
    ofstream::open( szName, std::ios::binary);
  }
};
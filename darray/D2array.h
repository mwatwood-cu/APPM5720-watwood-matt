// A simple one dimensional array class that can use arbitrary index range

class D2array
{
  public:
    D2array( int nc, int ibeg, int iend, int jbeg, int jend);
    D2array( int ibeg, int iend, int jbeg, int jend);
    D2array();
    ~D2array() {if( m_data != 0 ) delete[] m_data;}
    void define( int nc, int ibeg, int iend, int jbeg, int jend);
    void define( int ibeg, int iend, int jbeg, int jend);
    inline bool in_range( int c, int i, int j)
    {return 1 <= c && c <= m_nc && m_ib <= i && i <= m_ie && m_jb <= j && j <= m_je;}
// overload parenthesis operator 
    inline double& operator()( int c, int i, int j){
            // Turn on array bounds check 
#ifdef BZ_DEBUG  
        try
        {  
            if (!in_range(c,i,j)) throw 10;
        }
        catch(int e) 
        {
            std::cout <<
                "Error Index (c,i,j) = (" << c << "," <<
                i << "," << j << ") not in range 1<= c <= " <<
                m_nc << " " << m_ib << " <= i <= " << m_ie <<
		" " << m_jb << " <= j <= " << m_je <<std::endl;
        }
#endif
       return m_data[m_base+m_offc*c+m_offi*i+m_offj*j];
    }
    inline double& operator()( int i, int j)
    {
#ifdef BZ_DEBUG
        try{
            if (!in_range(1,i,j)) throw 10;
        }
        catch(int e)
        {
            std::cout << "Error Index (c,i,j) = (" << 1
                      << "," << i << "," << j << ") not in range 1<= c <= "
                      << m_nc << " "
                      << m_ib << " <= i <= " << m_ie <<" "
		      << m_jb << " <= j <= " << m_je << std::endl;
        }
#endif
    return m_data[m_base+m_offi*i+m_offc+m_offj*j];
    }
    inline bool is_defined(){return m_data != NULL;}
    inline double* c_ptr(){return m_data;}
    static bool m_corder;
    int m_ib, m_ie, m_jb, m_je;
    ssize_t m_base;
    size_t m_offi, m_offj, m_offc, m_npts;
    void define_offsets();
    void set_value( double scalar );
    void copy( const D2array& u );
    int m_nc, m_ni, m_nj;
  private:
   double* m_data;
};

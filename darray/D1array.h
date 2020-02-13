// A simple one dimensional array class that can use arbitrary index range

class D1array
{
  public:
    D1array( int nc, int ibeg, int iend);
    D1array( int ibeg, int iend);
    D1array();
    ~D1array() {if( m_data != 0 ) delete[] m_data;}
    void define( int nc, int ibeg, int iend);
    void define( int ibeg, int iend);
    inline bool in_range( int c, int i)
    {return 1 <= c && c <= m_nc && m_ib <= i && i <= m_ie;}
// overload parenthesis operator 
    inline double& operator()( int c, int i){
            // Turn on array bounds check 
#ifdef BZ_DEBUG  
        try
        {  
            if (!in_range(c,i)) throw 10;
        }
        catch(int e) 
        {
            std::cout <<
                "Error Index (c,i) = (" << c << "," <<
                i << ") not in range 1<= c <= " <<
                m_nc << " " << m_ib << " <= i <= " << m_ie;
        }
#endif
       return m_data[m_base+m_offc*c+m_offi*i];
    }
    inline double& operator()( int i)
    {
#ifdef BZ_DEBUG
        try{
            if (!in_range(1,i)) throw 10;
        }
        catch(int e)
        {
            std::cout << "Error Index (c,i,j,k) = (" << 1
                      << "," << i << ") not in range 1<= c <= "
                      << m_nc << " "
                      << m_ib << " <= i <= " << m_ie;
        }
#endif
    return m_data[m_base+m_offi*i+m_offc];
    }
    inline bool is_defined(){return m_data != NULL;}
    inline double* c_ptr(){return m_data;}
    static bool m_corder;
    int m_ib, m_ie;
    ssize_t m_base;
    size_t m_offi, m_offc, m_npts;
    void define_offsets();
    void set_value( double scalar );
    void copy( const D1array& u );
    int m_nc, m_ni;
  private:
   double* m_data;
};

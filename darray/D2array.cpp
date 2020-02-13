#include <iostream>

#include "D2array.h"

using namespace std;

// Default value of array ordering
bool D2array::m_corder = false;

// Constructors
//-----------------------------------------------------------------------
D2array::D2array( int nc, int ibeg, int iend, int jbeg, int jend)
{
    m_nc = nc;
    m_ib = ibeg;
    m_ie = iend;
    m_jb = jbeg;
    m_je = jend;
    m_ni = m_ie-m_ib+1;
    m_nj = m_je-m_jb+1;
    if( m_nc*m_ni*m_nj > 0 )
        m_data = new double[m_nc*m_ni*m_nj];
    else
        m_data = NULL;
    define_offsets();
}
// Single variable array
D2array::D2array(int ibeg, int iend, int jbeg, int jend)
{
    m_nc = 1;
    m_ib = ibeg;
    m_ie = iend;
    m_jb = jbeg;
    m_je = jend;
    m_ni = m_ie-m_ib+1;
    m_nj = m_je-m_jb+1;
    if( m_nc*m_ni*m_nj > 0 )
        m_data = new double[m_nc*m_ni*m_nj];
    else
        m_data = NULL;
    define_offsets();
}

//-----------------------------------------------------------------------

//Default constructor
D2array::D2array()
{
    m_nc = m_ib = m_ie = m_jb = m_je = 0;
    m_data = NULL;
}

//
void D2array::define( int nc, int ibeg, int iend, int jbeg, int jend)
{
    if( m_data != NULL )
        delete[] m_data;

    m_nc = nc;
    m_ib = ibeg;
    m_ie = iend;
    m_jb = jbeg;
    m_je = jend;
    m_ni = m_ie-m_ib+1;
    m_nj = m_je-m_jb+1;
    if( m_nc*m_ni*m_nj > 0 )
        m_data = new double[m_nc*m_ni*m_nj];
    else
        m_data = NULL;
    define_offsets();
}

void D2array::define(int ibeg, int iend, int jbeg, int jend)
{
    if( m_data != NULL )
        delete[] m_data;

    m_nc = 1;
    m_ib = ibeg;
    m_ie = iend;
    m_jb = jbeg;
    m_je = jend;
    m_ni = m_ie-m_ib+1;
    m_nj = m_je-m_jb+1;
    if( m_nc*m_ni*m_nj > 0 )
        m_data = new double[m_nc*m_ni*m_nj];
    else
        m_data = NULL;
    define_offsets();
}



//-----------------------------------------------------------------------
void D2array::set_value(double scalar)
{
    for( size_t i = 0 ; i < m_npts ; i++ )
        m_data[i] = scalar;
}

//-----------------------------------------------------------------------
void D2array::define_offsets()
{
    m_npts = static_cast<size_t>(m_ni)*static_cast<size_t>(m_nj)*m_nc;
    if( m_corder )
    {
        m_offc = m_ni*m_nj;
        m_offi = m_nj;
        m_offj = 1;
    }
    else
    {
        m_offc = 1;
        m_offi = m_nc;
        m_offj = m_nc*m_ni;
    }
    m_base = -m_offc -(m_offi*m_ib)-(m_offj*m_jb);
}

//-----------------------------------------------------------------------
void D2array::copy( const D2array& u )
{
    if( m_data != NULL )
        delete[] m_data;
    
    m_nc = u.m_nc;
    m_ib = u.m_ib;
    m_ie = u.m_ie;
    m_jb = u.m_jb;
    m_je = u.m_je;
    m_ni = m_ie-m_ib+1;
    m_nj = m_je-m_jb+1;
    if( m_nc*m_ni*m_nj > 0 )
    {
        m_data = new double[m_nc*m_ni*m_nj];
        for( int i = 0 ; i < m_nc*m_ni*m_nj ; i++ )
            m_data[i] = u.m_data[i];
    }
    else
        m_data = NULL;
    define_offsets();
}



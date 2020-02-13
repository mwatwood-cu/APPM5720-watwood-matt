#include <iostream>

#include "D1array.h"

using namespace std;

// Default value of array ordering
bool D1array::m_corder = false;

// Constructors
//-----------------------------------------------------------------------
D1array::D1array( int nc, int ibeg, int iend)
{
    m_nc = nc;
    m_ib = ibeg;
    m_ie = iend;
    m_ni = m_ie-m_ib+1;
    if( m_nc*m_ni > 0 )
        m_data = new double[m_nc*m_ni];
    else
        m_data = NULL;
    define_offsets();
}
// Single variable array
D1array::D1array(int ibeg, int iend)
{
    m_nc = 1;
    m_ib = ibeg;
    m_ie = iend;
    m_ni = m_ie-m_ib+1;
    if( m_nc*m_ni > 0 )
        m_data = new double[m_nc*m_ni];
    else
        m_data = NULL;
    define_offsets();
}

//-----------------------------------------------------------------------

//Default constructor
D1array::D1array()
{
    m_nc = m_ib = m_ie = 0;
    m_data = NULL;
}

//
void D1array::define( int nc, int ibeg, int iend)
{
    if( m_data != NULL )
        delete[] m_data;

    m_nc = nc;
    m_ib = ibeg;
    m_ie = iend;
    m_ni = m_ie-m_ib+1;
    if( m_nc*m_ni > 0 )
        m_data = new double[m_nc*m_ni];
    else
        m_data = NULL;
    define_offsets();
}

void D1array::define(int ibeg, int iend)
{
    if( m_data != NULL )
        delete[] m_data;

    m_nc = 1;
    m_ib = ibeg;
    m_ie = iend;
    m_ni = m_ie-m_ib+1;
    if( m_nc*m_ni > 0 )
        m_data = new double[m_nc*m_ni];
    else
        m_data = NULL;
    define_offsets();
}



//-----------------------------------------------------------------------
void D1array::set_value(double scalar)
{
    for( size_t i = 0 ; i < m_npts ; i++ )
        m_data[i] = scalar;
}

//-----------------------------------------------------------------------
void D1array::define_offsets()
{
    m_npts = static_cast<size_t>(m_ni)*m_nc;
    if( m_corder )
    {
            // (i,c)=i-ib+ni*(c-1)
        m_offc = m_ni;
        m_offi = 1;
    }
    else
    {
            // (c,i)=c-1+nc*(i-ib)
        m_offc = 1;
        m_offi = m_nc;
    }
    m_base = -(m_offc)-(m_offi*m_ib);
}

//-----------------------------------------------------------------------
void D1array::copy( const D1array& u )
{
    if( m_data != NULL )
        delete[] m_data;
    
    m_nc = u.m_nc;
    m_ib = u.m_ib;
    m_ie = u.m_ie;
    m_ni = m_ie-m_ib+1;
    if( m_nc*m_ni > 0 )
    {
        m_data = new double[m_nc*m_ni];
        for( int i = 0 ; i < m_nc*m_ni ; i++ )
            m_data[i] = u.m_data[i];
    }
    else
        m_data = NULL;
    define_offsets();
}



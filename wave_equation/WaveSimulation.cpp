
#include <iostream>
#include "WaveSimulation.h"
#include <cmath>

using namespace std;

//Constructors where any that are not the most general just call the most general constructor.

WaveSimulation::WaveSimulation(int grid_points, double x_length, double y_length)
{
    Initialize(grid_points, grid_points, 0, x_length, 0, y_length);
}

WaveSimulation::WaveSimulation(int grid_points, double x_start, double x_end, double y_start, double y_end)
{
    Initialize(grid_points, grid_points, x_start, x_end, y_start, y_end);
}

WaveSimulation::WaveSimulation(int grid_points_x, int grid_points_y, double x_length, double y_length)
{
    Initialize(grid_points_x, grid_points_y, 0, x_length, 0, y_length);
}

WaveSimulation::WaveSimulation(int grid_points_x, int grid_points_y, double x_start, double x_end, double y_start, double y_end)
{
    Initialize(grid_points_x, grid_points_y, x_start, x_end, y_start, y_end);
}

void WaveSimulation::Initialize(int grid_points_x, int grid_points_y, double x_start, double x_end, double y_start, double y_end)
{
    //Saving Input ranges.
    m_x_grid_points = grid_points_x;
    m_y_grid_points = grid_points_y;
    
    m_x_start = x_start;
    m_y_start = y_start;
    m_x_end = x_end;
    m_y_end = y_end;

    //Creating arrays
    m_solution_array.define(1, -1, m_y_grid_points, -1, m_x_grid_points);
    m_solution_array.set_value(0);
    for(int i=-1; i<=m_y_grid_points;i++) {
        for(int j=-1; j<=m_x_grid_points; j++) {
            double y_value =(double)i/m_y_grid_points*(m_y_end-m_y_start);
            m_solution_array(i,j) = y_value*y_value;
        }
    }
    
    m_laplacian.define(1, 0, m_y_grid_points-1, 0, m_x_grid_points-1);  
    m_laplacian.set_value(2);
    
    m_speed_sound.define(1, -1, m_y_grid_points, -1, m_x_grid_points);
    m_speed_sound.set_value(1);

    m_time_steps = 1;
}

//Private function implementations

void WaveSimulation::Calculate_Laplacian()
{
    Calculate_Laplacian(0, m_x_grid_points, 0, m_y_grid_points);
}

void WaveSimulation::Calculate_Laplacian(int x_start, int x_finish, int y_start, int y_finish)
{
    double h_x = (double)(m_x_end-m_x_start)/m_x_grid_points;
    double h_y = (double)(m_y_end-m_y_start)/m_y_grid_points;
    for(int i=y_start; i<y_finish; i++) {
        for(int j=x_start; j<x_finish; j++) {
            //First do the x direction then add the y direction
            m_laplacian(i,j) = 1/(2*h_x*h_x)*Calculate_Derivative(1, i, j)
                              +1/(2*h_y*h_y)*Calculate_Derivative(0, i, j);
        }
     }

}

double WaveSimulation::Calculate_Derivative(int x_dir, int i, int j)
{
    double c_2_p1, c_2, c_2_m1;
    double first_piece, second_piece, third_piece;
    if(x_dir == 1)
    {
        c_2_p1 = m_speed_sound(i+1,j)*m_speed_sound(i+1,j);
        c_2 = m_speed_sound(i,j)*m_speed_sound(i,j);
        c_2_m1 = m_speed_sound(i-1,j)*m_speed_sound(i-1,j);
        first_piece = (c_2_p1+c_2)*m_solution_array(i+1,j);
        second_piece = (c_2_p1+2*c_2+c_2_m1)*m_solution_array(i,j);
        third_piece = (c_2+c_2_m1)*m_solution_array(i-1,j); 
    }
    else 
    {
        c_2_p1 = m_speed_sound(i,j+1)*m_speed_sound(i,j+1);
        c_2 = m_speed_sound(i,j)*m_speed_sound(i,j);
        c_2_m1 = m_speed_sound(i,j-1)*m_speed_sound(i,j-1);
        first_piece = (c_2_p1+c_2)*m_solution_array(i,j+1);
        second_piece = (c_2_p1+2*c_2+c_2_m1)*m_solution_array(i,j);
        third_piece = (c_2+c_2_m1)*m_solution_array(i,j-1); 
    }
    return first_piece-second_piece+third_piece;
}

void WaveSimulation::print_all()
{
    cout<<"Speed of Sound"<<endl;
    for(int i=-1; i<10; i++) {
        for(int j=-1; j<10; j++) {
            cout<<m_speed_sound(i,j)<<", ";
        }
        cout<<endl;
    }

    cout<<"Laplacian"<<endl;
    for(int i=0; i<10; i++) {
        for(int j=0; j<10; j++) {
            cout<<m_laplacian(i,j)<<", ";
        }
        cout<<endl;
    }
   
    cout<<"Solution Data"<<endl;
    for(int i=-1; i<10; i++) {
        for(int j=-1; j<10; j++) {
            cout<<m_solution_array(i,j)<<", ";
        }
        cout<<endl;
    }
}

void WaveSimulation::Verbose_Test()
{
    cout<<"Testing u = sin(4x+3)cos(3y+5) and c = 1 for all (x,y)"<<endl;

   double y_value;
   double x_value;
    for(int i=-1; i<=m_y_grid_points;i++) {
        for(int j=-1; j<=m_x_grid_points; j++) {
            x_value =(double)j/m_x_grid_points*(m_x_end-m_x_start);
            y_value =(double)i/m_y_grid_points*(m_y_end-m_y_start);
            m_solution_array(i,j) = sin(4*x_value+3)*cos(3*y_value+5);
        }
    } 
   D2array exact;
   exact.define(1, 0, m_y_grid_points, 0, m_x_grid_points);
   exact.set_value(0);
   for(int i=0; i<m_y_grid_points;i++) {
        for(int j=0; j<m_x_grid_points; j++) {
            y_value =(double)i/m_y_grid_points*(m_y_end-m_y_start);
            x_value =(double)j/m_x_grid_points*(m_x_end-m_x_start);
            exact(i,j) = -16*sin(4*x_value+3)*cos(3*y_value+5)
                        -9*cos(3*y_value+5)*sin(4*x_value+3);
        }
    } 
   
   Calculate_Laplacian();

   double l1_error = Calculate_Error(exact, m_laplacian, 1);
   double linf_error = Calculate_Error(exact, m_laplacian, 0);
   double l2_error = Calculate_Error(exact, m_laplacian, 2);

   cout<<"L1 Norm Error: "<<l1_error<<endl;
   cout<<"L2 Norm Error: "<<l2_error<<endl;
   cout<<"L-Inf Norm Error: "<<linf_error<<endl;
   
}

double WaveSimulation::Calculate_Error(D2array& solution, D2array& approx, int type)
{
    double error = 0;
    double h = (m_x_end-m_x_start)/m_x_grid_points;
    if(type ==0)
    {
        double tmp_error = 0;
        for(int i=0; i<m_y_grid_points;i++) {
            for(int j=0; j<m_x_grid_points; j++) {
                tmp_error += abs(solution(i,j)-approx(i,j))*h;
            }
            error = tmp_error>error?tmp_error:error;
            tmp_error = 0;
        }
    }
    else if(type ==1)
    {
        double tmp_error = 0;
        for(int j=0; j<m_x_grid_points;j++) {
            for(int i=0; i<m_y_grid_points; i++) {
               tmp_error += abs(solution(i,j)-approx(i,j))*h;
            }
            error = tmp_error>error?tmp_error:error;
            tmp_error = 0;
        }
    }
    else
    {
        for(int i=0; i<m_y_grid_points;i++) {
            for(int j=0; j<m_x_grid_points; j++) {
                double diff = solution(i,j)-approx(i,j);
                error += pow(diff,2);
            }
        }
        error = pow(error*h, 0.5);
    }
    return error;
}

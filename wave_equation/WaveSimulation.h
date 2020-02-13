//A class to hold the needed values for a wave equation model
//
#include "../darray/D2array.h"
#include "../darray/D1array.h"
class WaveSimulation 
{
    public:
      //Class Constructors and Destructors
      WaveSimulation(int grid_points, double x_length, double y_length);
      WaveSimulation(int grid_points, double x_start, double x_end, double y_start, double y_end);
      WaveSimulation(int grid_points_x, int grid_points_y, double x_length, double y_length); 
      WaveSimulation(int grid_points_x, int grid_points_y, double x_start, double x_end, double y_start, double y_end);
      //No pointers so  this is simple destructor
      ~WaveSimulation(){}
      
      void print_all();

      //Will calculate for all locations in domain
      void Calculate_Laplacian();
      //Will only calculate a given area. 
      void Calculate_Laplacian(int x_start, int x_finish, int y_start, int y_finish);
      
      void Verbose_Test();

      double Calculate_Error(D2array& solution, D2array& approx, int order);

    private:
      //Helper function to simplify constructors
      void Initialize(int grid_points_x, int grid_points_y, double x_start, double x_end, double y_start, double y_end);

      //Helper function
      double Calculate_Derivative(int x_dir, int i, int j);

      //Member variables for domain size and shape
      double m_x_start;
      double m_x_end;
      double m_y_start;
      double m_y_end;
      int m_x_grid_points;
      int m_y_grid_points;
      
      //Member variables for solution
      D2array m_solution_array;

      //Member variables for discretization parameters
      D2array m_laplacian;

      //Member variables for parameters
      D2array m_speed_sound;
     
      //Member variables for simulation time
      int m_time_steps;


};

/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/boundary_velocity/function_dike.h>
#include <aspect/utilities.h>
#include <aspect/global.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace BoundaryVelocity
  {
    template <int dim>
    FunctionDike<dim>::FunctionDike ()
      :
      boundary_velocity_function (dim)
    {}




    template <int dim>
	double
	FunctionDike<dim>::
	read_from_txt(const std::string path,
	              const unsigned int mpi_communicator)const
	{
		double temp_input = 0;
		if(this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == mpi_communicator)
		{
			std::ifstream data_from_txt(path.c_str());
			std::string line;
			getline(data_from_txt, line);
			while(line == "")
			{
				data_from_txt.close();
				sleep(1.);
				std::ifstream data_from_txt(path.c_str());
				getline(data_from_txt, line);
			}
			std::stringstream ss(line);
			ss >> temp_input;
			data_from_txt.close();
			return temp_input;
		} 
		return temp_input;
	}	

    template <int dim>
    Tensor<1,dim>
    FunctionDike<dim>::
    boundary_velocity (const types::boundary_id boundary_indicator,
                       const Point<dim> &position) const
    {
      Tensor<1,dim> velocity;

      Utilities::NaturalCoordinate<dim> point =
        this->get_geometry_model().cartesian_to_other_coordinates(position, coordinate_system);

      //20251201 判定岩墙启动情况
	  int communicator_id = 0;
	  unsigned int out_comm_id = 0;
	  while (out_comm_id < 5000)
	  {
		  if (Utilities::MPI::this_mpi_process (this->get_mpi_communicator()) == out_comm_id)
		  {
			  communicator_id = out_comm_id;
			  break;
		  } 
		  out_comm_id ++;
	  } 
	  //20260109 关闭岩墙功能时不读取各项数据，按照普通的function运行
	  const double total_melt_extraction = dike_mode ? read_from_txt(melt_info_path + "mantle_total_melt_extraction.txt", communicator_id) : 0.;
	  //Unit:km
	  //岩墙宽度不固定时需要使用
	  const double dike_ymax_km = dike_mode ? constant_dike_length ? 0. : read_from_txt(melt_info_path + "dike_emplace_location_ymax.txt", communicator_id) : 0.;
      const double dike_ymin_km = dike_mode ? constant_dike_length ? 0. : read_from_txt(melt_info_path + "dike_emplace_location_ymin.txt", communicator_id) : 0.;
	  //Unit : m
	  const double dike_length_value = constant_dike_length ? dike_length : (dike_ymax_km - dike_ymin_km) * 1e3;
	  //岩墙长度不固定时，可能出现为0的情况（如初始步或无明显软流圈上涌）
      const bool dike_start = dike_mode && dike_length_value > 0 && total_melt_extraction * 1e6 >= dike_width * dike_length_value;
      for (unsigned int d=0; d<dim; ++d)
	  {
		  //暂时只支持2D模型和笛卡尔坐标系。如有3D或其他坐标系需求，请按需修改程序。
		  if(dike_start && coordinate_system_string == "cartesian" && dim == 2 && d == 1 && boundary_indicator == this->get_geometry_model().translate_symbolic_boundary_name_to_id ("bottom"))
		  {
			  const double temp_velocity = boundary_velocity_function.value(Utilities::convert_array_to_point<dim>(point.get_coordinates()), d);
			  //修正规则：根据岩墙长度，时间步长和M值确定当前步额外扩张面积，并在底边界进行修正。
			  //固定岩墙长度模式下，可能会出现岩墙长度不满足设置值的情况，导致一定的误差。
			  
			  //Unit:km^2
			  //M值不固定时需要使用
			  const double current_melt_extraction = constant_M ? 0. : read_from_txt(melt_info_path + "mantle_current_melt_extraction.txt", communicator_id);
			  //Unit : m/s
			  const double full_spreading_rate = spreading_rate / year_in_seconds;
			  //Unit : none
			  const double M_value = constant_M ? M_for_dike : current_melt_extraction * 1e6 / (dike_length_value * full_spreading_rate * this->get_timestep());
			  //M * full_spreading_rate * this->get_timestep()：额外扩张区域宽度。dike_length：额外扩张区域长度
			  //求解扩张面积后需要换算为m/yr，因此要除 this->get_timestep() / year_in_seconds
			  //化简后变为M_value * dike_length_value * full_spreading_rate * year_in_seconds
			  //修正的速度平均分摊在底边界上
			  const double modify_bottom_velocity = M_value * dike_length_value * full_spreading_rate * year_in_seconds / max_x;
			  velocity[d] = temp_velocity - modify_bottom_velocity;
		  }
		  else
			  velocity[d] = boundary_velocity_function.value(Utilities::convert_array_to_point<dim>(point.get_coordinates()), d);
	  }

      if (use_spherical_unit_vectors)
        velocity = Utilities::Coordinates::spherical_to_cartesian_vector(velocity, position);

      // ASPECT always wants things in MKS system. however, as described
      // in the documentation of this class, we interpret the formulas
      // given to this plugin as meters per year if the global flag
      // for using years instead of seconds is given. so if someone
      // write "5" in their parameter file and sets the flag, then this
      // means "5 meters/year" and we need to convert it to the ASPECT
      // time system by dividing by the number of seconds per year
      // to get MKS units
      if (this->convert_output_to_years())
        return velocity / year_in_seconds;
      else
        return velocity;
    }


    template <int dim>
    void
    FunctionDike<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        boundary_velocity_function.set_time (this->get_time() / year_in_seconds);
      else
        boundary_velocity_function.set_time (this->get_time());
    }


    template <int dim>
    void
    FunctionDike<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary velocity model");
      {
        prm.enter_subsection("Function");
        {
          prm.declare_entry ("Coordinate system", "cartesian",
                             Patterns::Selection ("cartesian|spherical|depth"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
                             "are interpreted as r,phi or r,phi,theta in 2D/3D "
                             "respectively with theta being the polar angle. `depth' "
                             "will create a function, in which only the first "
                             "parameter is non-zero, which is interpreted to "
                             "be the depth of the point.");
          prm.declare_entry ("Use spherical unit vectors", "false",
                             Patterns::Bool (),
                             "Specify velocity as $r$, $\\phi$, and $\\theta$ components "
                             "instead of $x$, $y$, and $z$. Positive velocities point up, east, "
                             "and north (in 3D) or out and clockwise (in 2D). "
                             "This setting only makes sense for spherical geometries.");

          Functions::ParsedFunction<dim>::declare_parameters (prm, dim);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    FunctionDike<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary velocity model");
      {
        prm.enter_subsection("Function");
        {
          coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
		  coordinate_system_string = prm.get ("Coordinate system");
          use_spherical_unit_vectors = prm.get_bool("Use spherical unit vectors");
          if (use_spherical_unit_vectors)
            AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical,
                         ExcMessage ("Spherical unit vectors should not be used "
                                     "when geometry model is not spherical."));
        }
        try
          {
            boundary_velocity_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Boundary velocity model.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
          }
        prm.leave_subsection();
      }
      prm.leave_subsection();
	  
	  //20251201
	  //根据新增的速度修正功能所需，读取输入参数
	  prm.enter_subsection("Geometry model");
	  {
		  prm.enter_subsection("Box");
		  {
			  max_x = prm.get_double ("X extent");
		  }
		  prm.leave_subsection();
	  }
	  prm.leave_subsection();
	
	  prm.enter_subsection("Postprocess");
	  {
		  prm.enter_subsection("Melt extraction");
		  {
			  //20251104
			  melt_info_path = this->get_output_directory() + "melt_info/";
			  dike_mode = prm.get_bool ("Use dike mode");
			  constant_M = prm.get_bool ("Constant M");
			  M_for_dike = prm.get_double ("M for dike");
			  dike_width = prm.get_double ("Width of dike(m)");
			  dike_length = prm.get_double ("Length of dike(m)");
			  spreading_rate = prm.get_double ("Full spreading rate of model(m/yr)");
			  //20251201
			  constant_dike_length = prm.get_bool ("Constant dike length");
		  }
		  prm.leave_subsection();
	  }
	  prm.leave_subsection();
	}
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryVelocity
  {
    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(FunctionDike,
                                            "function dike",
                                            "Implementation of a model in which the boundary "
                                            "velocity is given in terms of an explicit formula "
                                            "that is elaborated in the parameters in section "
                                            "``Boundary velocity model|Function''. The format of these "
                                            "functions follows the syntax understood by the "
                                            "muparser library, see Section~\\ref{sec:muparser-format}."
                                            "\n\n"
                                            "The formula you describe in the mentioned "
                                            "section is a semicolon separated list of velocities "
                                            "for each of the $d$ components of the velocity vector. "
                                            "These $d$ formulas are interpreted as having units "
                                            "m/s, unless the global input parameter ``Use "
                                            "years in output instead of seconds'' is set, in "
                                            "which case we interpret the formula expressions "
                                            "as having units m/year."
                                            "\n\n"
                                            "Likewise, since the symbol $t$ indicating time "
                                            "may appear in the formulas for the prescribed "
                                            "velocities, it is interpreted as having units "
                                            "seconds unless the global parameter above has "
                                            "been set.")
  }
}

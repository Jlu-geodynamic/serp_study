/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/rift/include/comp/lithosphere_rift.h>
#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/rift/include/temp/lithosphere_rift.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/box.h>
#include <aspect/boundary_temperature/interface.h>
#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/visco_plastic_grainsize.h>
#include <aspect/heating_model/interface.h>

#include <cmath>

namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    LithosphereRift<dim>::LithosphereRift ()
    {}


    template <int dim>
    void
    LithosphereRift<dim>::
    initialize ()
    {

      // Check that the required radioactive heating model ("compositional heating") is used
      const std::vector<std::string> &heating_models = this->get_heating_model_manager().get_active_heating_model_names();
      AssertThrow(std::find(heating_models.begin(), heating_models.end(), "compositional heating") != heating_models.end(),
                  ExcMessage("The lithosphere with rift initial temperature plugin requires the compositional heating plugin."));

      // Check that the required material model ("visco plastic") is used
      //AssertThrow((dynamic_cast<MaterialModel::ViscoPlasticStrain<dim> *> (const_cast<MaterialModel::Interface<dim> *>(&this->get_material_model()))) != 0,
      AssertThrow((dynamic_cast<MaterialModel::ViscoPlasticGrainSize<dim> *> (const_cast<MaterialModel::Interface<dim> *>(&this->get_material_model()))) != 0,
                  ExcMessage("The lithosphere with rift initial temperature plugin requires the viscoplastic material model plugin."));
    }


	//20240701 计算初始温度时，将所有地幔放在一起，并使用mantle_upper的参数
	//20240702 克拉通区域的初始温度计算按照克拉通厚度调整。
	//卡拉通整体温度结构会更冷。
    template <int dim>
    double
    LithosphereRift<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // Determine coordinate system
      const bool cartesian_geometry = dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model()) != NULL ? true : false;

      // Get the distance to the line segments along a path parallel to the surface
      //InitialComposition::LithosphereRift<dim> *ic = dynamic_cast<InitialComposition::LithosphereRift<dim> *> (const_cast<InitialComposition::Interface<dim> *>(&this->get_initial_composition()));

      double distance_to_rift_axis = 1e23;
      Point<dim-1> surface_position;
      std::pair<double, unsigned int> distance_to_L_polygon;
      const std::list<std::unique_ptr<InitialComposition::Interface<dim> > > & initial_composition_objects = this->get_initial_composition_manager().get_active_initial_composition_conditions();
      for (typename std::list<std::unique_ptr<InitialComposition::Interface<dim> > >::const_iterator it = initial_composition_objects.begin(); it != initial_composition_objects.end(); ++it)
        if ( InitialComposition::LithosphereRift<dim> *ic = dynamic_cast<InitialComposition::LithosphereRift<dim> *> ((*it).get()))
          {
            surface_position = ic->surface_position(position, cartesian_geometry);
            distance_to_rift_axis = ic->distance_to_rift(surface_position);
            distance_to_L_polygon = ic->distance_to_polygon(surface_position);
          }

      // Compute the local thickness of the upper crust, lower crust and mantle part of the lithosphere
      // based on the distance from the rift axis.
	  std::vector<double> local_thicknesses(4);
      
	  //20240704
	  //新增选项：边界内为克拉通，或边界外为克拉通
	  //20240705
	  //克拉通和普通岩石圈的初始温度单独计算
	  //20260127
	  //膏盐模式添加一层亏损地幔
	  const bool center_is_cratontic = add_cratontic && center_cratontic && surface_position[0] >= cratontic_boundary[0] && surface_position[0] < cratontic_boundary[1];
	  const bool side_is_cratontic = add_cratontic && !center_cratontic &&(surface_position[0] < cratontic_boundary[0] || surface_position[0] >= cratontic_boundary[1]);
	  bool in_cratontic = false;
	  if (salt_mode)
	  {
		  local_thicknesses[0] = ((0.5+0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*polygon_thicknesses[distance_to_L_polygon.second][0]
                                  +(0.5-0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*thicknesses[0])
                                 *(!blend_rift_and_polygon && distance_to_L_polygon.first > 0.-2.*sigma_polygon ? 1. :
                                   (1.0 - A[0] * std::exp((-std::pow(distance_to_rift_axis,2)/(2.0*std::pow(sigma_rift,2))))));
          local_thicknesses[1] = ((0.5+0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*polygon_thicknesses[distance_to_L_polygon.second][1]
                                  +(0.5-0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*thicknesses[1])
                                 *(!blend_rift_and_polygon && distance_to_L_polygon.first > 0.-2.*sigma_polygon ? 1. :
                                   (1.0 - A[1] * std::exp((-std::pow(distance_to_rift_axis,2)/(2.0*std::pow(sigma_rift,2))))));
		  local_thicknesses[2] = ((0.5+0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*polygon_thicknesses[distance_to_L_polygon.second][2]
                                  +(0.5-0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*thicknesses[2])
                                 *(!blend_rift_and_polygon && distance_to_L_polygon.first > 0.-2.*sigma_polygon ? 1. :
                                   (1.0 - A[2] * std::exp((-std::pow(distance_to_rift_axis,2)/(2.0*std::pow(sigma_rift,2))))));
	      const double mantle_l_thickness = ((0.5 + 0.5 * std::tanh(distance_to_L_polygon.first / sigma_polygon)) * polygon_thicknesses[distance_to_L_polygon.second][3]
                                 + (0.5 - 0.5 * std::tanh(distance_to_L_polygon.first / sigma_polygon)) * thicknesses[3])
                                  * (!blend_rift_and_polygon && distance_to_L_polygon.first > 0. - 2. * sigma_polygon ? 1. :
                                   (1.0 - A[3] * std::exp((-std::pow(distance_to_rift_axis, 2) / (2.0 * std::pow(sigma_rift, 2))))));
		  const double mantle_d_thickness = ((0.5 + 0.5 * std::tanh(distance_to_L_polygon.first / sigma_polygon)) * polygon_thicknesses[distance_to_L_polygon.second][4]
                                 + (0.5 - 0.5 * std::tanh(distance_to_L_polygon.first / sigma_polygon)) * thicknesses[4])
                                  * (!blend_rift_and_polygon && distance_to_L_polygon.first > 0. - 2. * sigma_polygon ? 1. :
                                   (1.0 - A[4] * std::exp((-std::pow(distance_to_rift_axis, 2) / (2.0 * std::pow(sigma_rift, 2))))));
		  local_thicknesses[3] = mantle_l_thickness + mantle_d_thickness;
	  }
	  else if (!salt_mode && center_is_cratontic || side_is_cratontic)
	  {
		  in_cratontic = true;
		  local_thicknesses[0] = cratontic_layer_thicknesses[0];
		  local_thicknesses[1] = cratontic_layer_thicknesses[1];
		  //20240716
		  local_thicknesses[2] = exclude_mantle_lower ? cratontic_layer_thicknesses[2] + cratontic_layer_thicknesses[3]
		                                              : cratontic_layer_thicknesses[2] + cratontic_layer_thicknesses[3] + cratontic_layer_thicknesses[4];
		  local_thicknesses[3] = 0.;
	  }
	  else
	  {
		  local_thicknesses[0] = ((0.5+0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*polygon_thicknesses[distance_to_L_polygon.second][0]
                                  +(0.5-0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*thicknesses[0])
                                 *(!blend_rift_and_polygon && distance_to_L_polygon.first > 0.-2.*sigma_polygon ? 1. :
                                   (1.0 - A[0] * std::exp((-std::pow(distance_to_rift_axis,2)/(2.0*std::pow(sigma_rift,2))))));
          local_thicknesses[1] = ((0.5+0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*polygon_thicknesses[distance_to_L_polygon.second][1]
                                  +(0.5-0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*thicknesses[1])
                                 *(!blend_rift_and_polygon && distance_to_L_polygon.first > 0.-2.*sigma_polygon ? 1. :
                                   (1.0 - A[1] * std::exp((-std::pow(distance_to_rift_axis,2)/(2.0*std::pow(sigma_rift,2))))));
	      //20240701
	      const double mantle_upper = ((0.5+0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*polygon_thicknesses[distance_to_L_polygon.second][2]
                                  +(0.5-0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*thicknesses[2])
                                 *(!blend_rift_and_polygon && distance_to_L_polygon.first > 0.-2.*sigma_polygon ? 1. :
                                   (1.0 - A[2] * std::exp((-std::pow(distance_to_rift_axis,2)/(2.0*std::pow(sigma_rift,2))))));
	      const double mantle_middle = ((0.5+0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*polygon_thicknesses[distance_to_L_polygon.second][3]
                                  +(0.5-0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*thicknesses[3])
                                 *(!blend_rift_and_polygon && distance_to_L_polygon.first > 0.-2.*sigma_polygon ? 1. :
                                   (1.0 - A[3] * std::exp((-std::pow(distance_to_rift_axis,2)/(2.0*std::pow(sigma_rift,2))))));
	      const double mantle_lower = ((0.5+0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*polygon_thicknesses[distance_to_L_polygon.second][4]
                                  +(0.5-0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*thicknesses[4])
                                 *(!blend_rift_and_polygon && distance_to_L_polygon.first > 0.-2.*sigma_polygon ? 1. :
                                   (1.0 - A[4] * std::exp((-std::pow(distance_to_rift_axis,2)/(2.0*std::pow(sigma_rift,2))))));
			
		  local_thicknesses[2] = exclude_mantle_lower ? mantle_upper + mantle_middle
			                                          : mantle_upper + mantle_middle + mantle_lower;
		  local_thicknesses[3] = 0.;
		  
          
	  }
      //20241105
	  //计算温度结构时，如果有等于0的区域，给他们赋予一个非常小的值
	  for (unsigned int a = 0; a < local_thicknesses.size(); ++a)
		  if (local_thicknesses[a] == 0)
			  local_thicknesses[a] = local_thicknesses[a] + 1;

      const double depth = this->get_geometry_model().depth(position);

      return temperature(depth, in_cratontic, local_thicknesses);
    }


    template <int dim>
    double
    LithosphereRift<dim>::
    temperature (const double depth,
	             const bool in_cratontic,
                 const std::vector<double> layer_thicknesses) const
    {
	  if(salt_mode)
	  {
		  //盐岩模式下不考虑克拉通相关的问题
		  const double LAB_temp = LAB_isotherm[0];
		  
		  const double m1 = heat_productivities[0] * pow(layer_thicknesses[0], 2);
		  const double m2 = heat_productivities[1] * pow(layer_thicknesses[1], 2);
		  const double m3 = heat_productivities[2] * pow(layer_thicknesses[2], 2);
		  const double m4 = conductivities[0] * conductivities[1] * conductivities[2];
		  const double m5 = conductivities[1] * conductivities[2] * densities[0];
		  const double m6 = conductivities[0] * conductivities[1] * densities[2];
		  const double m7 = conductivities[0] * conductivities[2] * densities[1];
		  const double m8 = heat_productivities[0] * layer_thicknesses[0];
		  const double m9 = heat_productivities[1] * layer_thicknesses[1];
		  const double m10 = heat_productivities[2] * layer_thicknesses[2];
		  const double m11 = heat_productivities[1] * layer_thicknesses[0] * layer_thicknesses[1];
		  const double m12 = heat_productivities[2] * layer_thicknesses[0] * layer_thicknesses[2];
		  const double m13 = heat_productivities[2] * layer_thicknesses[1] * layer_thicknesses[2];
		  const double m14 = heat_productivities[1] * layer_thicknesses[0] * pow(layer_thicknesses[1], 2);
		  const double m15 = heat_productivities[2] * layer_thicknesses[0] * pow(layer_thicknesses[2], 2);
		  
		  const double n1 = m1 * m5 + 2 * m8 * layer_thicknesses[1] * m5 + 2 * m8 * layer_thicknesses[2] * m5 + 2 * m8 * layer_thicknesses[3] * m5 + m2 * m7 + 2 * m9 * layer_thicknesses[2] * m7 + 2 * m9 * layer_thicknesses[3] * m7 + m3 * m6 + 2 * m10 * layer_thicknesses[3] * m6 - 2 * T0 * m4 + 2 * LAB_temp * m4;
		  const double n2 = 2 * layer_thicknesses[0] * m4 + 2 * layer_thicknesses[1] * m4 + 2 * layer_thicknesses[2] * m4 + 2 * layer_thicknesses[3] * m4;
		  const double n3 = -m1 * m5 + m2 * m7 + 2 * m9 * layer_thicknesses[2] * m7 + 2 * m9 * layer_thicknesses[3] * m7 + m3 * m6 + 2 * m10 * layer_thicknesses[3] * m6 - 2 * T0 * m4 + 2 * LAB_temp * m4;
		  const double n4 = -m1 * m5 - 2 * m11 * m7 - m2 * m7 + m3 * m6 + 2 * m10 * layer_thicknesses[3] * m6 - 2 * T0 * m4 + 2 * LAB_temp * m4;
		  const double n5 = -m1 * m5 - 2 * m11 * m7 - m2 * m7 - 2 * m12 * m6 - 2 * m13 * m6 - m3 * m6 - 2 * T0 * m4 + 2 * LAB_temp * m4;
		  const double n6 = m1 * layer_thicknesses[1] * m5 + m1 * layer_thicknesses[2] * m5 + m1 * layer_thicknesses[3] * m5 + m14 * m7 + 2 * m11 * layer_thicknesses[2] * m7 + 2 * m11 * layer_thicknesses[3] * m7 + m15 * m6 + 2 * m12 * layer_thicknesses[3] * m6 + 2 * T0 * layer_thicknesses[1] * m4 + 2 * T0 * layer_thicknesses[2] * m4 + 2 * T0 * layer_thicknesses[3] * m4 + 2 * LAB_temp * layer_thicknesses[0] * m4;
		  const double n7 = m1 * layer_thicknesses[2] * m5 + m1 * layer_thicknesses[3] * m5 + 2 * m11 * layer_thicknesses[2] * m7 + 2 * m11 * layer_thicknesses[3] * m7 + m2 * layer_thicknesses[2] * m7 + m2 * layer_thicknesses[3] * m7 + m15 * m6 + 2 * m12 * layer_thicknesses[3] * m6 + heat_productivities[2] * layer_thicknesses[1] * pow(layer_thicknesses[2], 2) * m6 + 2 * m13 * layer_thicknesses[3] * m6 + 2 * T0 * layer_thicknesses[2] * m4 + 2 * T0 * layer_thicknesses[3] * m4 + 2 * LAB_temp * layer_thicknesses[0] * m4 + 2 * LAB_temp * layer_thicknesses[1] * m4;
		  const double n8 = m1 * layer_thicknesses[3] * m5 + 2 * m11 * layer_thicknesses[3] * m7 + m2 * layer_thicknesses[3] * m7 + 2 * m12 * layer_thicknesses[3] * m6 + 2 * m13 * layer_thicknesses[3] * m6 + m3 * layer_thicknesses[3] * m6 + 2 * T0 * layer_thicknesses[3] * m4 + 2 * LAB_temp * layer_thicknesses[0] * m4 + 2 * LAB_temp * layer_thicknesses[1] * m4 + 2 * LAB_temp * layer_thicknesses[2] * m4;
		  
		  const double C1 = n1 / n2;
          const double C3 = n3 / n2;
          const double C5 = n4 / n2;
          const double C7 = n5 / n2;
          const double T1 = n6 / n2;
          const double T2 = n7 / n2;
          const double T3 = n8 / n2;
		  
		  // Temperature in layer 1
          if (depth < layer_thicknesses[0])
              return -(densities[0] * heat_productivities[0]) / (2 * conductivities[0]) * std::pow(depth, 2) + C1 * depth + T0;

          // Temperature in layer 2
          else if (depth < layer_thicknesses[0] + layer_thicknesses[1])
              return -(densities[1] * heat_productivities[1]) / (2 * conductivities[1]) * std::pow(depth - layer_thicknesses[0], 2.) + C3 * (depth - layer_thicknesses[0]) + T1;

          // Temperature in layer 3
          else if (depth < layer_thicknesses[0] + layer_thicknesses[1] + layer_thicknesses[2])
              return -(densities[2] * heat_productivities[2]) / (2 * conductivities[2]) * std::pow(depth - layer_thicknesses[0] - layer_thicknesses[1], 2.) + C5 * (depth - layer_thicknesses[0] - layer_thicknesses[1]) + T2;

          // Temperature in layer 4
          else if (depth < layer_thicknesses[0] + layer_thicknesses[1] + layer_thicknesses[2] + layer_thicknesses[3])
              return -(densities[3] * heat_productivities[3]) / (2 * conductivities[3]) * std::pow(depth - layer_thicknesses[0] - layer_thicknesses[1] - layer_thicknesses[2], 2.) + C7 * (depth - layer_thicknesses[0] - layer_thicknesses[1] - layer_thicknesses[2]) + T3;

          else if (use_compensation_depth && depth < compensation_depth)
              return LAB_temp;

          else
              return 10. * LAB_temp;

	  }
	  else
	  {
		  // Compute some constants
		  const double a = 0.5*densities[0]*heat_productivities[0]*layer_thicknesses[0] + 0.5*densities[1]*heat_productivities[1]*layer_thicknesses[1] + conductivities[0]/layer_thicknesses[0]*T0;
		  const double b = 1./(conductivities[0]/layer_thicknesses[0]+conductivities[1]/layer_thicknesses[1]);
		  const double c = in_cratontic? 0.5*densities[1]*heat_productivities[1]*layer_thicknesses[1] + conductivities[2]/layer_thicknesses[2]*LAB_isotherm[1] : 0.5*densities[1]*heat_productivities[1]*layer_thicknesses[1] + conductivities[2]/layer_thicknesses[2]*LAB_isotherm[0];
		  const double d = 1./(conductivities[1]/layer_thicknesses[1]+conductivities[2]/layer_thicknesses[2]);

		  //Temperature at boundary between layer 1 and 2
		  const double T1 = (a*b + conductivities[1]/layer_thicknesses[1]*c*d*b) / (1.-(conductivities[1]*conductivities[1])/(layer_thicknesses[1]*layer_thicknesses[1])*d*b);
		  //Temperature at boundary between layer 2 and 3
		  const double T2 = (c + conductivities[1]/layer_thicknesses[1]*T1) * d;

		  // Temperature in layer 1
		  if (depth < layer_thicknesses[0])
			return -0.5*densities[0]*heat_productivities[0]/conductivities[0]*std::pow(depth,2) + (0.5*densities[0]*heat_productivities[0]*layer_thicknesses[0]/conductivities[0] + (T1-T0)/layer_thicknesses[0])*depth + T0;
		  // Temperature in layer 2
		  else if (depth < layer_thicknesses[0]+layer_thicknesses[1])
			return -0.5*densities[1]*heat_productivities[1]/conductivities[1]*std::pow(depth-layer_thicknesses[0],2.) + (0.5*densities[1]*heat_productivities[1]*layer_thicknesses[1]/conductivities[1] + (T2-T1)/layer_thicknesses[1])*(depth-layer_thicknesses[0]) + T1;
		  // Temperature in layer 3
		  else if (depth < layer_thicknesses[0]+layer_thicknesses[1]+layer_thicknesses[2])
			return in_cratontic ? (LAB_isotherm[1]-T2)/layer_thicknesses[2] *(depth-layer_thicknesses[0]-layer_thicknesses[1]) + T2 : (LAB_isotherm[0]-T2)/layer_thicknesses[2] *(depth-layer_thicknesses[0]-layer_thicknesses[1]) + T2;
		  // Return a constant sublithospheric temperature of 10*LAB_isotherm.
		  // This way we can combine the continental geotherm with an adiabatic profile from the input file
		  // using the "minimum" operator on the "List of initial temperature models"
		  else if (use_compensation_depth && depth < compensation_depth)
			return in_cratontic ? 10 * LAB_isotherm[1] : LAB_isotherm[0];
		  else
			return 10.*LAB_isotherm[0];
	  }
      
    }


    template <int dim>
    void
    LithosphereRift<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Lithosphere with rift");
        {

          prm.declare_entry ("Surface temperature", "273.15",
                             Patterns::Double (0),
                             "The value of the surface temperature. Units: Kelvin.");
          prm.declare_entry ("LAB isotherm temperature", "1673.15",
                             Patterns::List(Patterns::Double(0)),
                             "The value of the isothermal boundary temperature assumed at the LAB "
                             "and up to the reference depth . Units: Kelvin.");
          prm.declare_entry ("Use temperature compensation depth", "false",
                             Patterns::Bool (),
                             "Whether or not to use a compensation depth to which the LAB isotherm temperature "
                             "is prescribed when no lithosphere is present. Below this depth, 10 times the isotherm "
                             "is prescribed . Units: -.");
          prm.declare_entry ("Temperature compensation depth", "100e3",
                             Patterns::Double (0),
                             "The depth of the temperature compensation depth, "
                             "i.e. the depth up to which the LAB isotherm is prescribed when not in the lithosphere. Units: Kelvin.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    LithosphereRift<dim>::parse_parameters (ParameterHandler &prm)
    {
      unsigned int n_fields = 0;
      prm.enter_subsection ("Compositional fields");
      {
        n_fields = prm.get_integer ("Number of fields");
      }
      prm.leave_subsection();

      prm.enter_subsection ("Initial composition model");
      {
        prm.enter_subsection("Lithosphere with rift");
        {
          sigma_rift           = prm.get_double ("Standard deviation of Gaussian rift geometry");
          sigma_polygon        = prm.get_double ("Half width of polygon smoothing");
          blend_rift_and_polygon = prm.get_bool ("Blend polygons and rifts");
          A                    = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Amplitude of Gaussian rift geometry"))),
                                                                         5,
                                                                         "Amplitude of Gaussian rift geometry");
          thicknesses = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Layer thicknesses"))),
                                                                5,
                                                                "Layer thicknesses");
          // Split the string into the separate polygons
          const std::vector<std::string> temp_thicknesses = Utilities::split_string_list(prm.get("Lithospheric polygon layer thicknesses"),';');
          const unsigned int n_polygons = temp_thicknesses.size();
          polygon_thicknesses.resize(n_polygons);
          for (unsigned int i_polygons = 0; i_polygons < n_polygons; ++i_polygons)
            {
              polygon_thicknesses[i_polygons] = Utilities::string_to_double(Utilities::split_string_list(temp_thicknesses[i_polygons],','));
              AssertThrow(polygon_thicknesses[i_polygons].size()==5, ExcMessage ("The number of layer thicknesses should be equal to 4."));
            }
			
		  //20240702
		  add_cratontic = prm.get_bool ("Add cratontic lithosphere");
	      cratontic_boundary = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Segment of cratontic lithosphere boundary"))),
                                                                2,
                                                                "Segment of cratontic lithosphere boundary");
		  cratontic_layer_thicknesses = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cratontic lithosphere layer thicknesses"))),
                                                                5,
                                                                "Cratontic lithosphere layer thicknesses");
		  center_cratontic = prm.get_bool ("Inside the boundaries is cratontic lithosphere");
		  //20240716
		  exclude_mantle_lower = prm.get_bool ("Exclude phase 'mantle_lower' from mantle material");
		  //20260114
		  //新增盐岩模式，新模式初始温度计算四层，其他模式计算三层
		  salt_mode = prm.get_bool ("Salt mode");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Lithosphere with rift");
        {
          LAB_isotherm = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("LAB isotherm temperature"))),
                                                                2,
                                                                "LAB isotherm temperature");
          T0 = prm.get_double ("Surface temperature");
          use_compensation_depth = prm.get_bool ("Use temperature compensation depth");
          compensation_depth = prm.get_double ("Temperature compensation depth");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();


      // Retrieve the indices of the fields that represent the lithospheric layers.
	  //20260112 屏蔽报错
	  //20260114 考虑多模式切换，屏蔽报错，后期根据不同模式设置报错
      //AssertThrow(this->introspection().compositional_name_exists("upper"),ExcMessage("We need a compositional field called 'upper' representing the upper crust."));
      //AssertThrow(this->introspection().compositional_name_exists("lower"),ExcMessage("We need a compositional field called 'lower' representing the lower crust."));
      //AssertThrow(this->introspection().compositional_name_exists("mantle_upper"),ExcMessage("We need a compositional field called 'mantle_upper' representing the lithospheric part of the mantle."));
	  //AssertThrow(this->introspection().compositional_name_exists("mantle_middle"),ExcMessage("We need a compositional field called 'mantle_middle' representing the heterogeneous part of the mantle."));
	  //AssertThrow(this->introspection().compositional_name_exists("mantle_lower"),ExcMessage("We need a compositional field called 'mantle_lower' representing the heterogeneous part of the mantle."));

      // For now, we assume a 3-layer system with an upper crust, lower crust and lithospheric mantle
	  //20260112 现在可删除mantle_middle或mantle_lower
	  //20260114 根据不同模式切换
	  //20260127 膏盐模式添加一层亏损地幔
      const unsigned int id_layer_1 = salt_mode ? this->introspection().compositional_index_for_name("sediment_init") : this->introspection().compositional_index_for_name("upper");
      const unsigned int id_layer_2 = salt_mode ? this->introspection().compositional_index_for_name("salt") : this->introspection().compositional_index_for_name("lower");
      const unsigned int id_layer_3 = salt_mode ? this->introspection().compositional_index_for_name("crust") : this->introspection().compositional_index_for_name("mantle_upper");
	  const unsigned int id_layer_4 = salt_mode ? this->introspection().compositional_index_for_name("mantle_L") : this->introspection().compositional_name_exists("mantle_middle") ? this->introspection().compositional_index_for_name("mantle_middle") : 9999;
	  const unsigned int id_layer_5 = salt_mode ? this->introspection().compositional_name_exists("mantle_D") ? this->introspection().compositional_index_for_name("mantle_D") : 9999 : this->introspection().compositional_name_exists("mantle_lower") ? this->introspection().compositional_index_for_name("mantle_lower") : 9999;
      // Retrieve other material properties set in different sections such that there
      // is no need to set them twice.

      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Compositional heating");
        {
          // The heating model compositional heating prefixes an entry for the background material
          const std::vector<double> temp_heat_productivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Compositional heating values"))),
                                                               n_fields+1,
                                                               "Compositional heating values");
          // This sets the heat productivity in W/m3 units
          heat_productivities.push_back(temp_heat_productivities[id_layer_1+1]);
          heat_productivities.push_back(temp_heat_productivities[id_layer_2+1]);
          heat_productivities.push_back(temp_heat_productivities[id_layer_3+1]);
		  if(salt_mode || this->introspection().compositional_name_exists("mantle_middle"))
			  heat_productivities.push_back(temp_heat_productivities[id_layer_4+1]);
	      if((salt_mode && this->introspection().compositional_name_exists("mantle_D"))||this->introspection().compositional_name_exists("mantle_lower"))
			  heat_productivities.push_back(temp_heat_productivities[id_layer_5+1]);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Visco Plastic Grain Size");
        {
          // The material model viscoplastic prefixes an entry for the background material
          const std::vector<double> temp_densities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Densities"))),
                                                     n_fields+1,
                                                     "Densities");
          const std::vector<double> temp_thermal_diffusivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal diffusivities"))),
                                                                 n_fields+1,
                                                                 "Thermal diffusivities");
          const std::vector<double> temp_heat_capacities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Heat capacities"))),
                                                           n_fields+1,
                                                           "Heat capacities");

          densities.push_back(temp_densities[id_layer_1+1]);
          densities.push_back(temp_densities[id_layer_2+1]);
          densities.push_back(temp_densities[id_layer_3+1]);
		  if(salt_mode || this->introspection().compositional_name_exists("mantle_middle"))
		      densities.push_back(temp_densities[id_layer_4+1]);
	      if((salt_mode && this->introspection().compositional_name_exists("mantle_D"))||this->introspection().compositional_name_exists("mantle_lower"))
		      densities.push_back(temp_densities[id_layer_5+1]);

          // Thermal diffusivity kappa = k/(rho*cp), so thermal conducitivity k = kappa*rho*cp
          conductivities.push_back(temp_thermal_diffusivities[id_layer_1+1] * densities[0] * temp_heat_capacities[id_layer_1+1]);
          conductivities.push_back(temp_thermal_diffusivities[id_layer_2+1] * densities[1] * temp_heat_capacities[id_layer_2+1]);
          conductivities.push_back(temp_thermal_diffusivities[id_layer_3+1] * densities[2] * temp_heat_capacities[id_layer_3+1]);
		  if(salt_mode || this->introspection().compositional_name_exists("mantle_middle"))
		      conductivities.push_back(temp_thermal_diffusivities[id_layer_4+1] * densities[3] * temp_heat_capacities[id_layer_4+1]);
	      if((salt_mode && this->introspection().compositional_name_exists("mantle_D"))||this->introspection().compositional_name_exists("mantle_lower"))
		      conductivities.push_back(temp_thermal_diffusivities[id_layer_5+1] * densities[4] * temp_heat_capacities[id_layer_5+1]);

          //20260112 屏蔽报错，现在相态不一定为5个
		  // To obtain the radioactive heating rate in W/kg, we divide the volumetric heating rate by density
          //AssertThrow(heat_productivities.size() == 5 && densities.size() == 5 && conductivities.size() == 5,
                      //ExcMessage("The entries for density, conductivity and heat production do not match with the expected number of layers (5)."))

          //for (unsigned int i = 0; i<5; ++i)
		  for (unsigned int i = 0; i< heat_productivities.size(); ++i)
            heat_productivities[i] /= densities[i];
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
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(LithosphereRift,
                                              "lithosphere with rift",
                                              "This is a temperature initial condition that "
                                              "computes a continental geotherm based on the "
                                              "conductive equations of Turcotte and Schubert Ch. 4.6. "
                                              "The geotherm can be computed for any number of crustal "
                                              "layers, for each of which a density, heat production and thermal "
                                              "conductivity should be supplied. "
                                              "Make sure the top and bottom temperatures of the lithosphere "
                                              "agree with temperatures set in for example the temperature "
                                              "boundary conditions.")
  }
}

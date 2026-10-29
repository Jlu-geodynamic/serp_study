/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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


#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/rift/include/geom/lithosphere_rift.h>
#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/rift/include/comp/lithosphere_rift.h>
#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/rift/include/temp/lithosphere_rift.h>
#include <aspect/geometry_model/box.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/utilities.h>
#include <boost/lexical_cast.hpp>
#include <aspect/geometry_model/initial_topography_model/ascii_data.h>
#include <aspect/geometry_model/initial_topography_model/function.h>
#include <aspect/geometry_model/initial_topography_model/prm_polygon.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>
#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <aspect/plugins.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/function_lib.h>
namespace aspect
{
  namespace InitialTopographyModel
  {
    template <int dim>
    void
    LithosphereRift<dim>::
    initialize ()
    {
      // Compute the maximum topography amplitude based on isostasy.
      // Assume the reference density is representative for each layer (despite temperature dependence)
	  
      //20241031
	  //使用新规则计算高程
	  //base面通过layer thickness的输入值计算得出
	  //原理：艾利重力均衡，将软流圈视为液体，通过计算V排确定不同板片间的高程差
	  //计算温度影响
	  
	  // For now, we assume a 3-layer system with an upper crust, lower crust and lithospheric mantle
	  //20260112 现在可删除mantle_middle和mantle_lower
	  //20260127 不适用于膏盐模式，未来会考虑进行适配
      const unsigned int id_upper = this->introspection().compositional_index_for_name("upper");
      const unsigned int id_lower = this->introspection().compositional_index_for_name("lower");
      const unsigned int id_mantle_upper = this->introspection().compositional_index_for_name("mantle_upper");
	  const unsigned int id_mantle_middle = this->introspection().compositional_name_exists("mantle_middle") ? this->introspection().compositional_index_for_name("mantle_middle") : 9999;
	  const unsigned int id_mantle_lower = this->introspection().compositional_name_exists("mantle_lower") ? this->introspection().compositional_index_for_name("mantle_lower") : 9999;
	  
	  
	  densities.push_back(temp_densities[0]);
	  densities.push_back(temp_densities[id_upper+1]);
	  densities.push_back(temp_densities[id_lower+1]);
	  densities.push_back(temp_densities[id_mantle_upper+1]);
	  if(this->introspection().compositional_name_exists("mantle_middle"))
	      densities.push_back(temp_densities[id_mantle_middle+1]);
      if(this->introspection().compositional_name_exists("mantle_lower"))
	      densities.push_back(temp_densities[id_mantle_lower+1]);
  
	  thermal_expansivities.push_back(temp_thermal_expansivities[0]);
	  thermal_expansivities.push_back(temp_thermal_expansivities[id_upper+1]);
	  thermal_expansivities.push_back(temp_thermal_expansivities[id_lower+1]);
	  thermal_expansivities.push_back(temp_thermal_expansivities[id_mantle_upper+1]);
	  if(this->introspection().compositional_name_exists("mantle_middle"))
	      thermal_expansivities.push_back(temp_thermal_expansivities[id_mantle_middle+1]);
      if(this->introspection().compositional_name_exists("mantle_lower"))
	      thermal_expansivities.push_back(temp_thermal_expansivities[id_mantle_lower+1]);
	  
	  // This sets the heat productivity in W/m3 units
      heat_productivities.push_back(temp_heat_productivities[id_upper+1]);
      heat_productivities.push_back(temp_heat_productivities[id_lower+1]);
      heat_productivities.push_back(temp_heat_productivities[id_mantle_upper+1]);
	  if(this->introspection().compositional_name_exists("mantle_middle"))
          heat_productivities.push_back(temp_heat_productivities[id_mantle_middle+1]);
      if(this->introspection().compositional_name_exists("mantle_lower"))
	      heat_productivities.push_back(temp_heat_productivities[id_mantle_lower+1]);
	  for(unsigned int a = 0; a < heat_productivities.size(); ++a)
	  {
		  heat_productivities[a] /= densities[a];
	  }		  
													   
	  conductivities.push_back(temp_thermal_diffusivities[id_upper+1] * densities[1] * temp_heat_capacities[id_upper+1]);
	  conductivities.push_back(temp_thermal_diffusivities[id_lower+1] * densities[2] * temp_heat_capacities[id_lower+1]);
	  conductivities.push_back(temp_thermal_diffusivities[id_mantle_upper+1] * densities[3] * temp_heat_capacities[id_mantle_upper+1]);
	  if(this->introspection().compositional_name_exists("mantle_middle"))
	      conductivities.push_back(temp_thermal_diffusivities[id_mantle_middle+1] * densities[4] * temp_heat_capacities[id_mantle_middle+1]);
	  if(this->introspection().compositional_name_exists("mantle_lower"))
	      conductivities.push_back(temp_thermal_diffusivities[id_mantle_lower+1] * densities[5] * temp_heat_capacities[id_mantle_lower+1]);
	  
	  std::vector<double> base_thickness_for_temp(3);
	  std::vector<double> subduction_thickness_for_temp(3);
	  for (unsigned int a = 0; a < 3; a++) 
	  {
		  if (a == 2)
		  {
			  base_thickness_for_temp[a] = thicknesses[2] + thicknesses[3] + thicknesses[4] == 0 ? 
			                               1. : thicknesses[2] + thicknesses[3] + thicknesses[4];
			  subduction_thickness_for_temp[a] = cratontic_layer_thicknesses[2] + cratontic_layer_thicknesses[3] + cratontic_layer_thicknesses[4] == 0 ?
                                               	 1. : cratontic_layer_thicknesses[2] + cratontic_layer_thicknesses[3] + cratontic_layer_thicknesses[4];
		  }
		  else
		  {
			  base_thickness_for_temp[a] = thicknesses[a] == 0 ? 1. : thicknesses[a];
		      subduction_thickness_for_temp[a] = cratontic_layer_thicknesses[a] == 0 ? 1. : cratontic_layer_thicknesses[a];
		  }
	  }

	  const double a_base = 0.5*densities[1]*heat_productivities[0]*base_thickness_for_temp[0] + 0.5*densities[2]*heat_productivities[1]*base_thickness_for_temp[1] + conductivities[0]/base_thickness_for_temp[0]*T0;
      const double b_base = 1./(conductivities[0]/base_thickness_for_temp[0]+conductivities[1]/base_thickness_for_temp[1]);
      const double c_base = 0.5*densities[2]*heat_productivities[1]*base_thickness_for_temp[1] + conductivities[2]/base_thickness_for_temp[2]*LAB_isotherm[0];
      const double d_base = 1./(conductivities[1]/base_thickness_for_temp[1]+conductivities[2]/base_thickness_for_temp[2]);

      //Temperature at boundary between layer 1 and 2
      const double T1_base = (a_base*b_base + conductivities[1]/base_thickness_for_temp[1]*c_base*d_base*b_base) / (1.-(conductivities[1]*conductivities[1])/(base_thickness_for_temp[1]*base_thickness_for_temp[1])*d_base*b_base);
      //Temperature at boundary between layer 2 and 3
      const double T2_base = (c_base + conductivities[1]/base_thickness_for_temp[1]*T1_base) * d_base;
	  
	  //岩石圈分界面温度，按照线性增长计算
	  const double T3_base = thicknesses[2] + thicknesses[3] + thicknesses[4] != 0 ? T2_base + (LAB_isotherm[0] - T2_base) * (thicknesses[2] / (thicknesses[2] + thicknesses[3] + thicknesses[4])) : LAB_isotherm[0];
	  const double T4_base = thicknesses[3] + thicknesses[4] != 0 ? T3_base + (LAB_isotherm[0] - T3_base) * (thicknesses[3] / (thicknesses[3] + thicknesses[4])) : LAB_isotherm[0];
	  
	  
      const double a_subduction = 0.5*densities[1]*heat_productivities[0]*subduction_thickness_for_temp[0] + 0.5*densities[2]*heat_productivities[1]*subduction_thickness_for_temp[1] + conductivities[0]/subduction_thickness_for_temp[0]*T0;
      const double b_subduction = 1./(conductivities[0]/subduction_thickness_for_temp[0]+conductivities[1]/subduction_thickness_for_temp[1]);
      const double c_subduction = 0.5*densities[2]*heat_productivities[1]*subduction_thickness_for_temp[1] + conductivities[2]/subduction_thickness_for_temp[2]*LAB_isotherm[1];
      const double d_subduction = 1./(conductivities[1]/subduction_thickness_for_temp[1]+conductivities[2]/subduction_thickness_for_temp[2]);

      //Temperature at boundary between layer 1 and 2
      const double T1_subduction = (a_subduction*b_subduction + conductivities[1]/subduction_thickness_for_temp[1]*c_subduction*d_subduction*b_subduction) / (1.-(conductivities[1]*conductivities[1])/(subduction_thickness_for_temp[1]*subduction_thickness_for_temp[1])*d_subduction*b_subduction);
      //Temperature at boundary between layer 2 and 3
      const double T2_subduction = (c_subduction + conductivities[1]/subduction_thickness_for_temp[1]*T1_subduction) * d_subduction;
	  
	  const double T3_subduction = cratontic_layer_thicknesses[2] + cratontic_layer_thicknesses[3] + cratontic_layer_thicknesses[4] != 0 ? T2_subduction + (LAB_isotherm[1] - T2_subduction) * (cratontic_layer_thicknesses[2] / (cratontic_layer_thicknesses[2] + cratontic_layer_thicknesses[3] + cratontic_layer_thicknesses[4])) : LAB_isotherm[1];
	  const double T4_subduction = cratontic_layer_thicknesses[3] + cratontic_layer_thicknesses[4] != 0 ? T3_subduction + (LAB_isotherm[1] - T3_subduction) * (cratontic_layer_thicknesses[3] / (cratontic_layer_thicknesses[3] + cratontic_layer_thicknesses[4])) : LAB_isotherm[1];
	  
	  //20241122
	  //计算交界面密度
	  //地表密度：参考密度 *（1+Tp*热膨胀系数）
	  //后续密度：参考密度 *（1+（Tp-交界面温度）*热膨胀系数）
	  //温度单位为摄氏度
	  //温度超过Tp时即为模型输入的参考密度
	  double test_G_base = 0;
	  double test_G_2 = 0;
	  for(unsigned int a = 1; a <= base_thickness_for_temp[0] + base_thickness_for_temp[1] + base_thickness_for_temp[2]; a++)
	  {
		  //按照模型的顺序计算浮力
		  //1.计算温度；2.计算密度；3.计算浮力
		  //每100m计算一次
		  double temp_t = 0;
		  double temp_den = 0;
		  if(a % 100 == 0)
		  {
			  if (a > 0 && a <= base_thickness_for_temp[0])
			  {
				  temp_t = -0.5*densities[1]*heat_productivities[0]/conductivities[0]*std::pow(a,2) + (0.5*densities[1]*heat_productivities[0]*base_thickness_for_temp[0]/conductivities[0] + (T1_base-T0)/base_thickness_for_temp[0])*a + T0;
			      temp_den = densities[1] * (1 + std::max(0., Tp - temp_t) * thermal_expansivities[1]);
			      test_G_base += temp_den * 100;
				  temp_t = 0;
				  temp_den = 0;
			  }
			  else if (a > base_thickness_for_temp[0] && a <= base_thickness_for_temp[0] + base_thickness_for_temp[1])
			  {
				  temp_t = -0.5*densities[2]*heat_productivities[1]/conductivities[1]*std::pow(a-base_thickness_for_temp[0],2.) + (0.5*densities[2]*heat_productivities[1]*base_thickness_for_temp[1]/conductivities[1] + (T2_base-T1_base)/base_thickness_for_temp[1])*(a-base_thickness_for_temp[0]) + T1_base;
			      temp_den = densities[2] * (1 + std::max(0., Tp - temp_t) * thermal_expansivities[2]);
			      test_G_base += temp_den * 100;
				  temp_t = 0;
				  temp_den = 0;
			  }
			  else
			  {
				  temp_t = (LAB_isotherm[0]-T2_base)/base_thickness_for_temp[2] *(a-base_thickness_for_temp[0]-base_thickness_for_temp[1]) + T2_base;
			      temp_den = densities[3] * (1 + std::max(0., Tp - temp_t) * thermal_expansivities[3]);
			      test_G_base += temp_den * 100;
				  temp_t = 0;
				  temp_den = 0;
			  }
		  }
	  }
	  
	  for(unsigned int a = 1; a <= subduction_thickness_for_temp[0] + subduction_thickness_for_temp[1] + subduction_thickness_for_temp[2]; a++)
	  {
		  //按照模型的顺序计算浮力
		  //1.计算温度；2.计算密度；3.计算浮力
		  //每100m计算一次
		  double temp_t = 0;
		  double temp_den = 0;
		  if(a % 100 == 0)
		  {
			  if (a > 0 && a <= subduction_thickness_for_temp[0])
			  {
				  temp_t = -0.5*densities[1]*heat_productivities[0]/conductivities[0]*std::pow(a,2) + (0.5*densities[1]*heat_productivities[0]*subduction_thickness_for_temp[0]/conductivities[0] + (T1_subduction-T0)/subduction_thickness_for_temp[0])*a + T0;
			      temp_den = densities[1] * (1 + std::max(0., Tp - temp_t) * thermal_expansivities[1]);
			      test_G_2 += temp_den * 100;
				  temp_t = 0;
				  temp_den = 0;
			  }
			  else if (a > subduction_thickness_for_temp[0] && a <= subduction_thickness_for_temp[0] + subduction_thickness_for_temp[1])
			  {
				  temp_t = -0.5*densities[2]*heat_productivities[1]/conductivities[1]*std::pow(a-subduction_thickness_for_temp[0],2.) + (0.5*densities[2]*heat_productivities[1]*subduction_thickness_for_temp[1]/conductivities[1] + (T2_subduction-T1_subduction)/subduction_thickness_for_temp[1])*(a-subduction_thickness_for_temp[0]) + T1_subduction;
			      temp_den = densities[2] * (1 + std::max(0., Tp - temp_t) * thermal_expansivities[2]);
			      test_G_2 += temp_den * 100;
				  temp_t = 0;
				  temp_den = 0;
			  }
			  else
			  {
				  temp_t = (LAB_isotherm[1]-T2_subduction)/subduction_thickness_for_temp[2] *(a-subduction_thickness_for_temp[0]-subduction_thickness_for_temp[1]) + T2_subduction;
			      temp_den = densities[3] * (1 + std::max(0., Tp - temp_t) * thermal_expansivities[3]);
			      test_G_2 += temp_den * 100;
				  temp_t = 0;
				  temp_den = 0;
			  }
		  }
	  }
	  
	  
	  const double G_layer_base = 0.5 * 
	                            ( densities[1] * (2 + (std::max(0., Tp - T0) + std::max(0., Tp - T1_base)) * thermal_expansivities[1]) * thicknesses[0] 
	                            + densities[2] * (2 + (std::max(0., Tp - T1_base) + std::max(0., Tp - T2_base)) * thermal_expansivities[2]) * thicknesses[1] 
								+ densities[3] * (2 + (std::max(0., Tp - T2_base) + std::max(0., Tp - T3_base)) * thermal_expansivities[3]) * thicknesses[2] 
								+ densities[4] * (2 + (std::max(0., Tp - T3_base) + std::max(0., Tp - T4_base)) * thermal_expansivities[4]) * thicknesses[3] 
								+ densities[5] * (2 + std::max(0., Tp - T4_base) * thermal_expansivities[5]) * thicknesses[4]);
	  
	  const double G_layer_2 = 0.5 * 
	                            ( densities[1] * (2 + (std::max(0., Tp - T0) + std::max(0., Tp - T1_subduction)) * thermal_expansivities[1]) * cratontic_layer_thicknesses[0] 
	                            + densities[2] * (2 + (std::max(0., Tp - T1_subduction) + std::max(0., Tp - T2_subduction)) * thermal_expansivities[2]) * cratontic_layer_thicknesses[1] 
								+ densities[3] * (2 + (std::max(0., Tp - T2_subduction) + std::max(0., Tp - T3_subduction)) * thermal_expansivities[3]) * cratontic_layer_thicknesses[2] 
								+ densities[4] * (2 + (std::max(0., Tp - T3_subduction) + std::max(0., Tp - T4_subduction)) * thermal_expansivities[4]) * cratontic_layer_thicknesses[3] 
								+ densities[5] * (2 + std::max(0., Tp - T4_subduction) * thermal_expansivities[5]) * cratontic_layer_thicknesses[4]);
	  //h_base含义：Layer thicknesses结构的岩石圈如果泡在软流圈里，露出软流圈的高度
	  h_base = 0;
      //h_base = (thicknesses[0] + thicknesses[1] + thicknesses[2] + thicknesses[3] + thicknesses[4]) - G_layer_base / densities[0];
	  h_base = (thicknesses[0] + thicknesses[1] + thicknesses[2] + thicknesses[3] + thicknesses[4]) - test_G_base / densities[0];
	  //const double h_layer_2 = (cratontic_layer_thicknesses[0] + cratontic_layer_thicknesses[1] + cratontic_layer_thicknesses[2] + cratontic_layer_thicknesses[3] + cratontic_layer_thicknesses[4]) - G_layer_2 / densities[0];
	  const double h_layer_2 = (cratontic_layer_thicknesses[0] + cratontic_layer_thicknesses[1] + cratontic_layer_thicknesses[2] + cratontic_layer_thicknesses[3] + cratontic_layer_thicknesses[4]) - test_G_2 / densities[0];
	  //泡在软流圈里的高度相减，便是不同厚度岩石圈的高程差，此处假定Layer thicknesses结构的岩石圈高程为0
	  const double delta_h = h_layer_2 - h_base;
      
	  this->get_pcout() << "   Base height: " << h_base << " m" << std::endl;
      this->get_pcout() << "   Subduction plate height: " << h_layer_2 << " m" << std::endl;
	  this->get_pcout() << "   Subduction plate initial topography of polygon: " << delta_h << " m" << std::endl;
    }


    template <int dim>
    double
    LithosphereRift<dim>::
    value (const Point<dim-1> &position) const
    {
      // When cartesian, position contains x(,y); when spherical, position contains lon(,lat) (in degrees);
      // Turn into a Point<dim-1>
      Point<dim-1> surface_position;
      for (unsigned int d=0; d<dim-1; ++d)
        surface_position[d] = position[d];
	
	  

      // Get the distance to the line segments along a path parallel to the surface
      double distance_to_rift_axis = 1e23;
      std::pair<double,unsigned int> distance_to_L_polygon;
     // const std::list<std::unique_ptr<InitialComposition::Interface<dim> > > initial_composition_objects = this->get_initial_composition_manager().get_active_initial_composition_conditions();
      const std::list<std::unique_ptr<InitialComposition::Interface<dim> > > & initial_composition_objects = this->get_initial_composition_manager().get_active_initial_composition_conditions();
	  
	  //20240702
	  //克拉通地势变化与克拉通左侧边界区域保持一致。
	  //20240704
	  //新增选项：边界内为克拉通，或边界外为克拉通
	  //20241031 俯冲模式下这里没有用处，暂留
	  const bool center_is_cratontic = add_cratontic && center_cratontic && surface_position[0] >= cratontic_boundary[0] && surface_position[0] < cratontic_boundary[1];
	  const bool side_is_cratontic = add_cratontic && !center_cratontic &&(surface_position[0] < cratontic_boundary[0] || surface_position[0] >= cratontic_boundary[1]);
	  for (typename std::list<std::unique_ptr<InitialComposition::Interface<dim> > >::const_iterator it = initial_composition_objects.begin(); it != initial_composition_objects.end(); ++it)
	  if ( InitialComposition::LithosphereRift<dim> *ic = dynamic_cast<InitialComposition::LithosphereRift<dim> *> ((*it).get()))
		{
		  distance_to_rift_axis = ic->distance_to_rift(surface_position);
		  distance_to_L_polygon = ic->distance_to_polygon(surface_position);
		}
	  // Compute the topography based on distance to the rift and distance to the polygon
      std::vector<double> local_thicknesses(5);
	  if (center_is_cratontic || side_is_cratontic)
	  {
		  local_thicknesses = cratontic_layer_thicknesses;
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
          local_thicknesses[2] = ((0.5+0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*polygon_thicknesses[distance_to_L_polygon.second][2]
                              +(0.5-0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*thicknesses[2])
                             *(!blend_rift_and_polygon && distance_to_L_polygon.first > 0.-2.*sigma_polygon ? 1. :
                               (1.0 - A[2] * std::exp((-std::pow(distance_to_rift_axis,2)/(2.0*std::pow(sigma_rift,2))))));
	      local_thicknesses[3] = ((0.5+0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*polygon_thicknesses[distance_to_L_polygon.second][3]
                              +(0.5-0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*thicknesses[3])
                             *(!blend_rift_and_polygon && distance_to_L_polygon.first > 0.-2.*sigma_polygon ? 1. :
                               (1.0 - A[3] * std::exp((-std::pow(distance_to_rift_axis,2)/(2.0*std::pow(sigma_rift,2))))));
	      local_thicknesses[4] = ((0.5+0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*polygon_thicknesses[distance_to_L_polygon.second][4]
                              +(0.5-0.5*std::tanh(distance_to_L_polygon.first/sigma_polygon))*thicknesses[4])
                             *(!blend_rift_and_polygon && distance_to_L_polygon.first > 0.-2.*sigma_polygon ? 1. :
                               (1.0 - A[4] * std::exp((-std::pow(distance_to_rift_axis,2)/(2.0*std::pow(sigma_rift,2))))));
	  }
	  
	  //20241107
	  //评估高程时计算温度影响
	  //计算温度时地幔物质归为一类
      std::vector<double> local_thickness_for_temp(3);
	  for (unsigned int a = 0; a < 3; a++) 
	  {
		  if (a == 2)
			  local_thickness_for_temp[a] = local_thicknesses[2] + local_thicknesses[3] + local_thicknesses[4] == 0 ? 
			                               1. : local_thicknesses[2] + local_thicknesses[3] + local_thicknesses[4];
		  else
			  local_thickness_for_temp[a] = local_thicknesses[a] == 0 ? 1. : local_thicknesses[a];
	  }
	  const double a_local = 0.5*densities[0]*heat_productivities[0]*local_thickness_for_temp[0] + 0.5*densities[1]*heat_productivities[1]*local_thickness_for_temp[1] + conductivities[0]/local_thickness_for_temp[0]*T0;
      const double b_local = 1./(conductivities[0]/local_thickness_for_temp[0]+conductivities[1]/local_thickness_for_temp[1]);
      const double c_local = 0.5*densities[1]*heat_productivities[1]*local_thickness_for_temp[1] + conductivities[2]/local_thickness_for_temp[2]*LAB_isotherm[0];
      const double d_local = 1./(conductivities[1]/local_thickness_for_temp[1]+conductivities[2]/local_thickness_for_temp[2]);
	  const double T1_local = (a_local*b_local + conductivities[1]/local_thickness_for_temp[1]*c_local*d_local*b_local) / (1.-(conductivities[1]*conductivities[1])/(local_thickness_for_temp[1]*local_thickness_for_temp[1])*d_local*b_local);
      const double T2_local = (c_local + conductivities[1]/local_thickness_for_temp[1]*T1_local) * d_local;
	  //岩石圈分界面温度，按照线性增长计算
	  const double T3_local = local_thicknesses[2] + local_thicknesses[3] + local_thicknesses[4] != 0 ? T2_local + (LAB_isotherm[0] - T2_local) * (local_thicknesses[2] / (local_thicknesses[2] + local_thicknesses[3] + local_thicknesses[4])) : LAB_isotherm[0];
	  const double T4_local = local_thicknesses[3] + local_thicknesses[4] != 0 ? T3_local + (LAB_isotherm[0] - T3_local) * (local_thicknesses[3] / (local_thicknesses[3] + local_thicknesses[4])) : LAB_isotherm[0];

      //确定LAB温度
	  const double local_T_LAB = center_is_cratontic || side_is_cratontic ? LAB_isotherm[1] : LAB_isotherm[0];
	  // The local lithospheric column
	  //岩石圈分界面温度，按照线性增长计算
      const double G_local = 0.5 * 
	                            ( densities[1] * (2 + (2 * local_T_LAB - T0 - T1_local) * thermal_expansivities[1]) * local_thicknesses[0] 
	                            + densities[2] * (2 + (2 * local_T_LAB - T1_local - T2_local) * thermal_expansivities[2]) * local_thicknesses[1] 
								+ densities[3] * (2 + (2 * local_T_LAB - T2_local - T3_local) * thermal_expansivities[3]) * local_thicknesses[2] 
								+ densities[4] * (2 + (2 * local_T_LAB - T3_local - T4_local) * thermal_expansivities[4]) * local_thicknesses[3] 
								+ densities[5] * (2 + (local_T_LAB - T4_local) * thermal_expansivities[5]) * local_thicknesses[4]);
	  // The total local lithosphere thickness
      const double sum_local_thicknesses = std::accumulate(local_thicknesses.begin(), local_thicknesses.end(),0);
	  //测试：所有区域使用新规则计算高程，基准面设定为Layer thicknesses的岩石圈
	  //const double h_local = center_is_cratontic || side_is_cratontic ? sum_local_thicknesses - G_local / densities[0] - 195 : sum_local_thicknesses - G_local / densities[0];
	  const double h_local = sum_local_thicknesses - G_local / densities[0];

	  return h_local - h_base;
    }


    template <int dim>
    double
    LithosphereRift<dim>::
    max_topography () const
    {
      return maximum_topography;
    }


    template <int dim>
    void
    LithosphereRift<dim>::
    declare_parameters (ParameterHandler &)
    {
    }



    template <int dim>
    void
    LithosphereRift<dim>::parse_parameters (ParameterHandler &prm)
    {
      unsigned int n_fields;
	  Tp = prm.get_double ("Adiabatic surface temperature");
	  
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
		  cratontic_layer_thicknesses = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cratontic lithosphere layer thicknesses"))),
                                                                5,
                                                                "Cratontic lithosphere layer thicknesses");
          // Split the string into the separate polygons
          const std::vector<std::string> temp_thicknesses = Utilities::split_string_list(prm.get("Lithospheric polygon layer thicknesses"),';');
          const unsigned int n_polygons = temp_thicknesses.size();
          polygon_thicknesses.resize(n_polygons);
          for (unsigned int i_polygons = 0; i_polygons < n_polygons; ++i_polygons)
            {
              polygon_thicknesses[i_polygons] = Utilities::string_to_double(Utilities::split_string_list(temp_thicknesses[i_polygons],','));
              AssertThrow(polygon_thicknesses[i_polygons].size()==5, ExcMessage ("The number of layer thicknesses should be equal to 5."));
            }
		  //20240702
		  add_cratontic = prm.get_bool ("Add cratontic lithosphere");
	      cratontic_boundary = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Segment of cratontic lithosphere boundary"))),
                                                                2,
                                                                "Segment of cratontic lithosphere boundary");
		  center_cratontic = prm.get_bool ("Inside the boundaries is cratontic lithosphere");
		  
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
	  
	  
	  prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Compositional heating");
        {
          // The heating model compositional heating prefixes an entry for the background material
          temp_heat_productivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Compositional heating values"))),
                                                               n_fields+1,
                                                               "Compositional heating values");
          
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Visco Plastic Grain Size");
        {
          // The material model viscoplastic prefixes an entry for the background material
          temp_densities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Densities"))),
                                                                   n_fields+1,
                                                                   "Densities");
		  temp_thermal_expansivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal expansivities"))),
                                                                   n_fields+1,
                                                                   "Thermal expansivities");
          temp_thermal_diffusivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal diffusivities"))),
                                                                 n_fields+1,
                                                                 "Thermal diffusivities");
          temp_heat_capacities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Heat capacities"))),
                                                           n_fields+1,
                                                           "Heat capacities");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
	  prm.enter_subsection("Geometry model");
		{
			prm.enter_subsection("Box");
			{
				max_x = prm.get_double ("X extent");
			}
			prm.leave_subsection ();
		}
		prm.leave_subsection ();
		
		prm.enter_subsection ("Initial temperature model");
        {
          prm.enter_subsection("Lithosphere with rift");
          {
            LAB_isotherm = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("LAB isotherm temperature"))),
                                                                2,
                                                                "LAB isotherm temperature");
            T0 = prm.get_double ("Surface temperature");
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
  namespace InitialTopographyModel
  {
    ASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(LithosphereRift,
                                             "lithosphere with rift",
                                             "An initial topography model that defines the initial topography "
                                             "as constant inside each of a set of polylineal parts of the "
                                             "surface. The polylines, and their associated surface elevation, "
                                             "are defined in the `Geometry model/Initial topography/Prm polyline' "
                                             "section.")
  }
}

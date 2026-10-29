/*
  Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

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


#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/heat/latent_heat_viscoplastic_melt.h>
#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/visco_plastic_grainsize.h>
#include <aspect/material_model/interface.h>


namespace aspect
{
  namespace HeatingModel
  {
	  
	  
	
	  template <int dim>
	  double
	  LatentHeatViscoPlasticMelt<dim>::
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
	  bool
	  LatentHeatViscoPlasticMelt<dim>::
	  dike_depend(const double posx,
	              const double posy,
	              const unsigned int mpi_communicator)const
	  {

		  //确定岩墙位置
		  //20251201 新增岩墙长度不固定的判定情况
		  const unsigned int emplace_xmax = read_from_txt(melt_info_path + "dike_emplace_location_xmax.txt", mpi_communicator);
		  const unsigned int emplace_ymin = read_from_txt(melt_info_path + "dike_emplace_location_ymin.txt", mpi_communicator);
		  const double total_melt_extract = read_from_txt(melt_info_path + "mantle_total_melt_extraction.txt", mpi_communicator);
		  const double emplace_ymax = constant_dike_length ? 0. : read_from_txt(melt_info_path + "dike_emplace_location_ymax.txt", mpi_communicator);
		  const double dike_length_value = constant_dike_length ? dike_length : (emplace_ymax - emplace_ymin) * 1e3;
		  //判定条件：在岩墙内，且累计提取熔融超过设置的岩墙面积
		  //20260108 岩墙x位置判定修改：由最左侧端点改为中心点
		  return dike_length_value > 0 
		         //&& posx < emplace_xmax * 1e3 && posx >= emplace_xmax * 1e3 - dike_width
				 && posx < emplace_xmax * 1e3 + 0.5 * dike_width && posx >= emplace_xmax * 1e3 - 0.5 * dike_width
		         && posy >= emplace_ymin * 1e3 && posy < emplace_ymin * 1e3 + dike_length_value
		         && total_melt_extract * 1e6 >= dike_width * dike_length_value;
	  }	
	
    template <int dim>
    void
    LatentHeatViscoPlasticMelt<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
		  heating_model_outputs.heating_source_terms[q] = 0.0;
		  heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
		  heating_model_outputs.rates_of_temperature_change[q] = 0.0;
		  if (this->introspection().compositional_name_exists("melt_fraction") &&  this->get_timestep_number() > 0)
		  {

			  //当时间步长过小时，熔融造成的降温会非常严重
			  //因此将最小的时间步长限制在2500yr
			  const double time = std::max(this->get_timestep(), 2500 * 3600 * 24 * 365.25);
			  
			  //1.熔融吸热导致的温度变化
			  //在程序中与放射热的形式相同
			  //等式右侧单位 
			  //mantle_melt_temp_reduction/crustal_melt_temp_reduction: J*kg-1*K-1
			  //melt_fractions: none
			  //temperature: K
			  //densities: kg*m-3
			  //time: s
				  
		      //等式左侧单位
		      //heating_source_terms：J*s-1*m-3 (Watt*m-3)，
		      //从物理意义的角度出发，heating_source_terms 反映辐射热对温度的影响（radioactive heat production，Hr）
			  //在模型当中，adiabatic heat、latent heat、radioactive heat production 三项为并列关系
			  //因此 adiabatic heat、latent heat 在计算过程中被转换为和 radioactive heat production 相同的单位。
			  
			  //2.洋壳放热导致的温度变化
			  //根据不同熔融机制形成的产物进行区分

			  //20260112
			  const double last_melt_fractions = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("melt_fraction")];
			  const double melt_fractions_change = material_model_outputs.reaction_terms[q][this->introspection().compositional_index_for_name("melt_fraction")];
			  //混合熔融时，根据成分判定当前位置的熔融类型
			  const double target_mantle_fraction_change = include_mantle_melting && this->introspection().compositional_name_exists("new_crust") ? 
														    material_model_outputs.reaction_terms[q][this->introspection().compositional_index_for_name("new_crust")]
														  : 0.;
			  const double target_crust_fraction_change = include_crustal_melting && this->introspection().compositional_name_exists("core_complex") ? 
															material_model_outputs.reaction_terms[q][this->introspection().compositional_index_for_name("core_complex")]
														  : 0.;
			  
			  //判定熔融形式，混合型与独立型的计算规则有区别
			  if (include_mantle_melting && !include_crustal_melting)
				  heating_model_outputs.heating_source_terms[q] = last_melt_fractions > 1e-3 && melt_fractions_change > 1e-5 ?
			                                                      - mantle_melt_temp_reduction
			                                                      * melt_fractions_change
		                                                          * material_model_inputs.temperature[q] 
														          * material_model_outputs.densities[q]
														          / time
																: 0.;
				  
			  else if (!include_mantle_melting && include_crustal_melting)
				  heating_model_outputs.heating_source_terms[q] = last_melt_fractions > 1e-3 && melt_fractions_change > 1e-5 ?
				                                                  - crustal_melt_temp_reduction
			                                                      * melt_fractions_change
		                                                          * material_model_inputs.temperature[q] 
														          * material_model_outputs.densities[q]
														          / time
																: 0.;
			  //20260114 地幔相与地壳相占比计算移动至此处
              //盐岩模式不启用，判定方式：查询"salt"相是否存在	  
			  else if (!this->introspection().compositional_name_exists("salt") && include_mantle_melting && include_crustal_melting)
			  {
				  
				  const double sediment_1_value = this->introspection().compositional_name_exists("sediment_1") ? material_model_inputs.composition[q][this->introspection().compositional_index_for_name("sediment_1")] : 0.;
				  const double sediment_2_value = this->introspection().compositional_name_exists("sediment_2") ? material_model_inputs.composition[q][this->introspection().compositional_index_for_name("sediment_2")] : 0.;
				  const double crust_fraction = include_crustal_melting ? 
												material_model_inputs.composition[q][this->introspection().compositional_index_for_name("upper")]
											  + material_model_inputs.composition[q][this->introspection().compositional_index_for_name("lower")]
											  + material_model_inputs.composition[q][this->introspection().compositional_index_for_name("new_crust")]
											  : 0.;
				  const double mantle_fraction = include_mantle_melting ? 
											   1. - sediment_1_value
											   - sediment_2_value
											   - material_model_inputs.composition[q][this->introspection().compositional_index_for_name("upper")]
											   - material_model_inputs.composition[q][this->introspection().compositional_index_for_name("lower")]
											   - material_model_inputs.composition[q][this->introspection().compositional_index_for_name("mantle_upper")]
											   - target_mantle_fraction_change
											   - target_crust_fraction_change
											   : 0.;
				  if (mantle_fraction >= crust_fraction)
					  heating_model_outputs.heating_source_terms[q] = last_melt_fractions > 1e-3 && melt_fractions_change > 1e-5 ?
			                                                      - mantle_melt_temp_reduction
			                                                      * melt_fractions_change
		                                                          * material_model_inputs.temperature[q] 
														          * material_model_outputs.densities[q]
														          / time
																: 0.;
				  else
					  heating_model_outputs.heating_source_terms[q] = last_melt_fractions > 1e-3 && melt_fractions_change > 1e-5 ?
				                                                  - crustal_melt_temp_reduction
			                                                      * melt_fractions_change
		                                                          * material_model_inputs.temperature[q] 
														          * material_model_outputs.densities[q]
														          / time
																: 0.;
			  }
				  
			  else
				  heating_model_outputs.heating_source_terms[q] = 0.;
			  
              //2.洋壳放热导致的温度变化
			  //根据不同熔融机制形成的产物进行区分
			  //为了不与熔融吸热过程冲突，放热过程不使用ifelse，均为独立判定，仅在新生物质的成分满足条件时赋值

			  //仅有地幔熔融
			  if (heating_process && !dike_mode && !include_crustal_melting && target_mantle_fraction_change > 0.1)
				  heating_model_outputs.heating_source_terms[q] = mantle_melt_temp_reduction 
				                                                  * target_mantle_fraction_change
		                                                          * material_model_inputs.temperature[q] 
														          * material_model_outputs.densities[q]
														          / time;
			  //仅有地壳熔融
			  if (heating_process && !include_mantle_melting && target_crust_fraction_change > 0.1)
				  heating_model_outputs.heating_source_terms[q] = crustal_melt_temp_reduction 
				                                                  * target_crust_fraction_change
		                                                          * material_model_inputs.temperature[q] 
														          * material_model_outputs.densities[q]
														          / time;
			  //混合熔融
			  //暂定取成分高的熔融信息
			  if (heating_process && include_crustal_melting && include_mantle_melting && (target_mantle_fraction_change > 0.1 || target_crust_fraction_change > 0.1))
			  {
				  const double target_mantle_change_for_composite = target_mantle_fraction_change > 0.1 ? target_mantle_fraction_change : 0.;
				  const double target_crust_change_for_composite = target_crust_fraction_change > 0.1 ? target_crust_fraction_change : 0.;
				  heating_model_outputs.heating_source_terms[q] = target_mantle_change_for_composite >= target_crust_change_for_composite ?
				                                                  mantle_melt_temp_reduction 
				                                                  * target_mantle_fraction_change
		                                                          * material_model_inputs.temperature[q] 
														          * material_model_outputs.densities[q]
														          / time
																: crustal_melt_temp_reduction 
				                                                  * target_crust_fraction_change
		                                                          * material_model_inputs.temperature[q] 
														          * material_model_outputs.densities[q]
														          / time;

			  }

		      /* const double posx = material_model_inputs.position[q][0];
              const double posy = material_model_inputs.position[q][1];
		      //20251104 判定是否在岩墙内
		      bool in_dike = false;
		      if(dike_mode)
		      {
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
			  in_dike = dike_depend(posx, posy, communicator_id);
		      } */
			  //20251114 新增判定条件，确定洋壳变化是否非0
			  //20260112 删除M值判定
			  //20260112 测试：仅使用相变判定是否加热
			  //if (heating_process && this->get_parameters().enable_additional_stokes_rhs && in_dike && target_mantle_fraction_change > 0.1)
		      if (heating_process && this->get_parameters().enable_additional_stokes_rhs && dike_mode && target_mantle_fraction_change > 0.1)
			  {
				  heating_model_outputs.heating_source_terms[q] = mantle_melt_temp_reduction 
				                                                  * target_mantle_fraction_change
		                                                          * material_model_inputs.temperature[q] 
														          * material_model_outputs.densities[q]
														          / time;
			  }  
			  
		}
    }
  }

    template <int dim>
    void
    LatentHeatViscoPlasticMelt<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Latent heat viscoplastic melt");
        {
          prm.declare_entry ("Change of entropy on mantle melting", "400.",
                             Patterns::Double (),
                             "Change of entropy on mantle melting "
                             "Units: J*kg-1*K-1.");
		  prm.declare_entry ("Change of entropy on crustal melting", "250.",
                             Patterns::Double (),
                             "Change of entropy on crustal melting "
                             "Units: J*kg-1*K-1.");
		  //20251201 新增变量，控制是否进行放热过程
		  prm.declare_entry ("Use heating process of magma intrusion", "true",
                             Patterns::Bool (),
                             "Use heating process of magma intrusion"
                             "Units: none.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    LatentHeatViscoPlasticMelt<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Latent heat viscoplastic melt");
        {
          mantle_melt_temp_reduction = prm.get_double ("Change of entropy on mantle melting");
		  crustal_melt_temp_reduction = prm.get_double ("Change of entropy on crustal melting");
		  heating_process = prm.get_bool ("Use heating process of magma intrusion");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
	  
	  prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Visco Plastic Grain Size");
        {
			include_mantle_melting = prm.get_bool ("Include mantle melting");
			include_crustal_melting = prm.get_bool("Include crustal melting");
			//input_total_force = prm.get_double ("Total force in one side");
		}
		prm.leave_subsection();
	  }
	  prm.leave_subsection();
	  
	  //20251104
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
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(LatentHeatViscoPlasticMelt,
                                  "latent heat viscoplastic melt",
                                  "Implementation of a standard model for latent heat "
                                  "of melting. This assumes that there is a compositional field "
                                  "called porosity, and it uses the reaction term of this field "
                                  "(the fraction of material that melted in the current time step) "
                                  "multiplied by a constant entropy change for melting all "
                                  "of the material as source term of the heating model.\n"
                                  "If there is no field called porosity, the heating terms are 0.")
  }
}

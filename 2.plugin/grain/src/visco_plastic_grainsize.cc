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


#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/visco_plastic_grainsize.h>
#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/postprocess/melt_extraction.h>
#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/rift/include/comp/lithosphere_rift.h>



#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/signaling_nan.h>
#include <aspect/newton.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/geometry_model/box.h>

#include <iostream>
#include <dirent.h>
#include <sys/stat.h>

#include <mpi.h>
// 此模型以material model中的visco plastic模型为基础修改而来，grain size相关计算参考aterial model中的grain size model模型
// 原文件位置以及名称：material_model/visco_plastic.cc 
//                     material_model/rheology/visco_plastic.cc 
//					   material_model/rheology/diffusion_creep.cc 
//					   material_model/grain_size.cc
namespace aspect
{
  namespace MaterialModel
  {
	//20250417
	template <int dim>
    void
    ViscoPlasticGrainSize<dim>::initialize()
    {
	 //20260310 添加地幔相和玄武岩相开关	
	 double delta_temp = 0;
     double n_temperature = 0;
	 double delta_press = 0;
	 double n_pressure = 0;
     if(use_perplex_table && perplex_include_mantle_phase)
	 {
		  delta_temp = 0;
		  n_temperature = 0;
		  delta_press = 0;
		  n_pressure = 0;
		  std::string temp;
		  //不含碳
		  perplex_table_normal = std::make_unique<MaterialUtilities::Lookup::PerplexReader>(data_directory + normal_material_file_name, false, this->get_mpi_communicator());
		  //获取perplex范围，当温压条件位于范围之外时，使用参考密度
		  std::istringstream in_normal(Utilities::read_and_distribute_file_content(data_directory + normal_material_file_name, this->get_mpi_communicator()));
		  std::getline(in_normal, temp); // get next line, table file name
		  std::getline(in_normal, temp); // get next line, dimension of table
		  std::getline(in_normal, temp); // get next line, either T(K) or P(bar)
		  for (unsigned int i=0; i<2; i++)
			{
			  std::string natural_variable;
			  in_normal >> natural_variable;

			  if (natural_variable == "T(K)")
				{
				  std::getline(in_normal, temp);
				  in_normal >> normal_table_min_T;
				  std::getline(in_normal, temp);
				  in_normal >> delta_temp;
				   std::getline(in_normal, temp);
				  in_normal >> n_temperature;
				  std::getline(in_normal, temp); // get next line, either T(K), P(bar) or number of columns
				}
			  else if (natural_variable == "P(bar)")
				{
				  std::getline(in_normal, temp);
				  in_normal >> normal_table_min_P;
				  normal_table_min_P *= 1e5;  // conversion from [bar] to [Pa]
				  std::getline(in_normal, temp);
				  in_normal >> delta_press;
				  delta_press *= 1e5; // conversion from [bar] to [Pa]
				  std::getline(in_normal, temp);
				  in_normal >> n_pressure;
				  std::getline(in_normal, temp); // get next line, either T(K), P(bar) or number of columns
				}
			}
		  normal_table_max_T = normal_table_min_T + (n_temperature-1) * delta_temp;
		  normal_table_max_P = normal_table_min_P + (n_pressure-1) * delta_press;
		  
		  //含碳
		  perplex_table_carbon = std::make_unique<MaterialUtilities::Lookup::PerplexReader>(data_directory + carbon_material_file_name, false, this->get_mpi_communicator());
		  //获取perplex范围，当温压条件位于范围之外时，使用参考密度
		  delta_temp = 0;
		  n_temperature = 0;
		  delta_press = 0;
		  n_pressure = 0;
		  std::istringstream in_carbon(Utilities::read_and_distribute_file_content(data_directory + carbon_material_file_name, this->get_mpi_communicator()));
		  std::getline(in_carbon, temp); // get next line, table file name
		  std::getline(in_carbon, temp); // get next line, dimension of table
		  std::getline(in_carbon, temp); // get next line, either T(K) or P(bar)
		  for (unsigned int i=0; i<2; i++)
			{
			  std::string natural_variable;
			  in_carbon >> natural_variable;

			  if (natural_variable == "T(K)")
				{
				  std::getline(in_carbon, temp);
				  in_carbon >> carbon_table_min_T;
				  std::getline(in_carbon, temp);
				  in_carbon >> delta_temp;
				   std::getline(in_carbon, temp);
				  in_carbon >> n_temperature;
				  std::getline(in_carbon, temp); // get next line, either T(K), P(bar) or number of columns
				}
			  else if (natural_variable == "P(bar)")
				{
				  std::getline(in_carbon, temp);
				  in_carbon >> carbon_table_min_P;
				  carbon_table_min_P *= 1e5;  // conversion from [bar] to [Pa]
				  std::getline(in_carbon, temp);
				  in_carbon >> delta_press;
				  delta_press *= 1e5; // conversion from [bar] to [Pa]
				  std::getline(in_carbon, temp);
				  in_carbon >> n_pressure;
				  std::getline(in_carbon, temp); // get next line, either T(K), P(bar) or number of columns
				}
			}
		  carbon_table_max_T = carbon_table_min_T + (n_temperature-1) * delta_temp;
		  carbon_table_max_P = carbon_table_min_P + (n_pressure-1) * delta_press;
		  
		  //含水
		  perplex_table_wet = std::make_unique<MaterialUtilities::Lookup::PerplexReader>(data_directory + wet_material_file_name, false, this->get_mpi_communicator());
		  //获取perplex范围，当温压条件位于范围之外时，使用参考密度
		  delta_temp = 0;
		  n_temperature = 0;
		  delta_press = 0;
		  n_pressure = 0;
		  std::istringstream in_wet(Utilities::read_and_distribute_file_content(data_directory + wet_material_file_name, this->get_mpi_communicator()));
		  std::getline(in_wet, temp); // get next line, table file name
		  std::getline(in_wet, temp); // get next line, dimension of table
		  std::getline(in_wet, temp); // get next line, either T(K) or P(bar)
		  for (unsigned int i=0; i<2; i++)
			{
			  std::string natural_variable;
			  in_wet >> natural_variable;

			  if (natural_variable == "T(K)")
				{
				  std::getline(in_wet, temp);
				  in_wet >> wet_table_min_T;
				  std::getline(in_wet, temp);
				  in_wet >> delta_temp;
				   std::getline(in_wet, temp);
				  in_wet >> n_temperature;
				  std::getline(in_wet, temp); // get next line, either T(K), P(bar) or number of columns
				}
			  else if (natural_variable == "P(bar)")
				{
				  std::getline(in_wet, temp);
				  in_wet >> wet_table_min_P;
				  wet_table_min_P *= 1e5;  // conversion from [bar] to [Pa]
				  std::getline(in_wet, temp);
				  in_wet >> delta_press;
				  delta_press *= 1e5; // conversion from [bar] to [Pa]
				  std::getline(in_wet, temp);
				  in_wet >> n_pressure;
				  std::getline(in_wet, temp); // get next line, either T(K), P(bar) or number of columns
				}
			}
		  wet_table_max_T = wet_table_min_T + (n_temperature-1) * delta_temp;
		  wet_table_max_P = wet_table_min_P + (n_pressure-1) * delta_press;
	 }
	 
	 //20260310 玄武岩相
	 if(use_perplex_table && perplex_include_basalt_phase)
	 {
		  delta_temp = 0;
		  n_temperature = 0;
		  delta_press = 0;
		  n_pressure = 0;
		  std::string temp;
		  //玄武岩
		  perplex_table_basalt = std::make_unique<MaterialUtilities::Lookup::PerplexReader>(data_directory + basalt_material_file_name, false, this->get_mpi_communicator());
		  //获取perplex范围，当温压条件位于范围之外时，使用参考密度
		  std::istringstream in_basalt(Utilities::read_and_distribute_file_content(data_directory + basalt_material_file_name, this->get_mpi_communicator()));
		  std::getline(in_basalt, temp); // get next line, table file name
		  std::getline(in_basalt, temp); // get next line, dimension of table
		  std::getline(in_basalt, temp); // get next line, either T(K) or P(bar)
		  for (unsigned int i=0; i<2; i++)
			{
			  std::string natural_variable;
			  in_basalt >> natural_variable;

			  if (natural_variable == "T(K)")
				{
				  std::getline(in_basalt, temp);
				  in_basalt >> basalt_table_min_T;
				  std::getline(in_basalt, temp);
				  in_basalt >> delta_temp;
				   std::getline(in_basalt, temp);
				  in_basalt >> n_temperature;
				  std::getline(in_basalt, temp); // get next line, either T(K), P(bar) or number of columns
				}
			  else if (natural_variable == "P(bar)")
				{
				  std::getline(in_basalt, temp);
				  in_basalt >> basalt_table_min_P;
				  normal_table_min_P *= 1e5;  // conversion from [bar] to [Pa]
				  std::getline(in_basalt, temp);
				  in_basalt >> delta_press;
				  delta_press *= 1e5; // conversion from [bar] to [Pa]
				  std::getline(in_basalt, temp);
				  in_basalt >> n_pressure;
				  std::getline(in_basalt, temp); // get next line, either T(K), P(bar) or number of columns
				}
			}
		  basalt_table_max_T = basalt_table_min_T + (n_temperature-1) * delta_temp;
		  basalt_table_max_P = basalt_table_min_P + (n_pressure-1) * delta_press;
	 }
    }

    template <int dim>
    bool
    ViscoPlasticGrainSize<dim>::
    is_yielding (const double pressure,
                 const double temperature,
                 const std::vector<double> &composition,
                 const SymmetricTensor<2,dim> &strain_rate) const
    {
      /* The following returns whether or not the material is plastically yielding
       * as documented in evaluate.
       */
      bool plastic_yielding = false;

      MaterialModel::MaterialModelInputs <dim> in (/*n_evaluation_points=*/1,
                                                                           this->n_compositional_fields());
      unsigned int i = 0;

      in.pressure[i] = pressure;
      in.temperature[i] = temperature;
      in.composition[i] = composition;
      in.strain_rate[i] = strain_rate;

      const std::vector<double> temp_volume_fractions
        = MaterialUtilities::compute_composition_fractions(composition,
                                                           rheology->get_volumetric_composition_mask());
      //修正成分判断
	  const std::vector<double> volume_fractions = rheology->fractions_modify(temp_volume_fractions);

      const IsostrainViscosities isostrain_viscosities
        = rheology->calculate_isostrain_viscosities(in, i, volume_fractions);

      std::vector<double>::const_iterator max_composition
        = std::max_element(volume_fractions.begin(),volume_fractions.end());

      plastic_yielding = isostrain_viscosities.composition_yielding[std::distance(volume_fractions.begin(),
                                                                                  max_composition)];

      return plastic_yielding;
    }



    template <int dim>
    bool
    ViscoPlasticGrainSize<dim>::
    is_yielding(const MaterialModelInputs<dim> &in) const
    {
      Assert(in.n_evaluation_points() == 1, ExcInternalError());
      const std::vector<double> temp_volume_fractions = MaterialUtilities::compute_composition_fractions(in.composition[0], rheology->get_volumetric_composition_mask());
	  //修正成分判断
	  const std::vector<double> volume_fractions = rheology->fractions_modify(temp_volume_fractions);

      /* The following handles phases in a similar way as in the 'evaluate' function.
       * Results then enter the calculation of plastic yielding.
       */
      std::vector<double> phase_function_values(phase_function.n_phase_transitions(), 0.0);

      if (phase_function.n_phase_transitions() > 0)
        {
          const double gravity_norm = this->get_gravity_model().gravity_vector(in.position[0]).norm();

          double reference_density;
          if (this->get_adiabatic_conditions().is_initialized())
            {
              reference_density = this->get_adiabatic_conditions().density(in.position[0]);
            }
          else
            {
              EquationOfStateOutputs<dim> eos_outputs_all_phases (this->n_compositional_fields()+1+phase_function.n_phase_transitions());
              equation_of_state.evaluate(in, 0, eos_outputs_all_phases);
              reference_density = eos_outputs_all_phases.densities[0];
            }

          MaterialUtilities::PhaseFunctionInputs<dim> phase_inputs(in.temperature[0],
                                                                   in.pressure[0],
                                                                   this->get_geometry_model().depth(in.position[0]),
                                                                   gravity_norm*reference_density,
                                                                   numbers::invalid_unsigned_int);

          for (unsigned int j=0; j < phase_function.n_phase_transitions(); j++)
            {
              phase_inputs.phase_index = j;
              phase_function_values[j] = phase_function.compute_value(phase_inputs);
            }
        }

      /* The following returns whether or not the material is plastically yielding
       * as documented in evaluate.
       */
      const IsostrainViscosities isostrain_viscosities = rheology->calculate_isostrain_viscosities(in, 0, volume_fractions, phase_function_values, phase_function.n_phase_transitions_for_each_composition());

      std::vector<double>::const_iterator max_composition = std::max_element(volume_fractions.begin(), volume_fractions.end());
      const bool plastic_yielding = isostrain_viscosities.composition_yielding[std::distance(volume_fractions.begin(), max_composition)];

      return plastic_yielding;
    }
	
	template <int dim>
	bool 
	ViscoPlasticGrainSize<dim>::
	get_convey_flag(const std::string path,
	                const bool melt_convey,
	                const unsigned int crood_x,
			        const unsigned int crood_y_in_flag_container,
	                const unsigned int mpi_communicator)const
	{
		
		unsigned int flag_out = 0;
		if(melt_convey && this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == mpi_communicator)
		{
			std::stringstream get_name;
		    get_name << crood_x;
		    std::string temp_txt_name = get_name.str();
			std::ifstream data_from_txt((path + temp_txt_name + ".txt").c_str());
			if (!data_from_txt)
			{
				data_from_txt.close();
				sleep(0.01);
				std::ifstream data_from_txt((path + temp_txt_name + ".txt").c_str());
			}
			std::string flag;
			getline(data_from_txt, flag);
				
            double flag_single; 		
            std::vector<double> flag_temp;			
            std::stringstream ss(flag);
            while(ss >> flag_single)
			{
				flag_temp.push_back(flag_single);
			}
			data_from_txt.close();
			flag_out = flag_temp[crood_y_in_flag_container];
		} 
		if (flag_out == 1)
			return true;
		else 
			return false;
	}
	
	template <int dim>
	std::vector<double>
    ViscoPlasticGrainSize<dim>::
	melt_convey_function(const bool convey_start_new_crust,
	                     const bool convey_start_core_complex,
						 const std::vector<double> &raw_volume_fractions)const
	{
		std::vector<double> convey_temp_volume_fractions = raw_volume_fractions;
		unsigned int convey_target_mantle_melting_id = 0;
		unsigned int convey_target_crustal_melting_id = 0;
		//只有地幔熔融
		if (!include_crustal_melting && convey_start_new_crust && this->introspection().compositional_name_exists("new_crust"))
		{
			convey_target_mantle_melting_id = this->introspection().compositional_index_for_name("new_crust") + 1;
			//输送到目标区域 
			convey_temp_volume_fractions[convey_target_mantle_melting_id] = 1.;
			for (unsigned int m = 0; m < raw_volume_fractions.size(); m++)
			{
				if (m != convey_target_mantle_melting_id)
					convey_temp_volume_fractions[m] = 0;
			}
		}
		//只有地壳熔融
		else if (!include_mantle_melting && convey_start_core_complex && this->introspection().compositional_name_exists("core_complex"))
		{
			convey_target_crustal_melting_id = this->introspection().compositional_index_for_name("core_complex") + 1;
			//输送到目标区域 
			convey_temp_volume_fractions[convey_target_crustal_melting_id] = 1.;
			for (unsigned int m = 0; m < raw_volume_fractions.size(); m++)
			{
				if (m != convey_target_crustal_melting_id)
					convey_temp_volume_fractions[m] = 0;
			}
		}
		//混合熔融
		//两种产物不互相代替
		else if (include_crustal_melting && include_mantle_melting && this->introspection().compositional_name_exists("new_crust") && this->introspection().compositional_name_exists("core_complex"))
		{
			convey_target_mantle_melting_id = this->introspection().compositional_index_for_name("new_crust") + 1;
			convey_target_crustal_melting_id = this->introspection().compositional_index_for_name("core_complex") + 1;
			//新生地壳
			if (convey_start_new_crust)
			{
				convey_temp_volume_fractions[convey_target_mantle_melting_id] = 1. - convey_temp_volume_fractions[convey_target_crustal_melting_id];
				for (unsigned int m = 0; m < raw_volume_fractions.size(); m++)
				{
					if (m != convey_target_mantle_melting_id && m != convey_target_crustal_melting_id)
					convey_temp_volume_fractions[m] = 0;
				}
			}
			//核杂岩
			if (convey_start_core_complex)
			{
				convey_temp_volume_fractions[convey_target_crustal_melting_id] = 1. - convey_temp_volume_fractions[convey_target_mantle_melting_id];
				for (unsigned int m = 0; m < raw_volume_fractions.size(); m++)
				{
					if (m != convey_target_mantle_melting_id && m != convey_target_crustal_melting_id)
					convey_temp_volume_fractions[m] = 0;
				}
			}
			//暂定规定不能同时在一个网格内产生两种产物
			//如果有需求的话，可以根据情况添加
		}
		else
			convey_temp_volume_fractions = raw_volume_fractions;
		
		 //熔体传输过程结束后，按照当前模型的设置进行组分修正
		const std::vector<double>out_volume_fractions = rheology->fractions_modify(convey_temp_volume_fractions);
		//输出
		return out_volume_fractions;
	} 
	
    template <int dim>
    std::array<double, 6>
    ViscoPlasticGrainSize<dim>::
	melted_judgement(const double temperature,
					 const double pressure,
					 const std::vector<double> &composition,
					 const std::vector<double> &premelt_volume_fractions)const
	{
	   //20250701
	   //调整成分场调用规则，便于在prm文件中删除不需要的场
	   //20240706
	   const std::vector<double> temp_composition = composition;
	   const double judge_temp_max_melt = this->introspection().compositional_name_exists("max_melt") ? std::abs(composition[this->introspection().compositional_index_for_name("max_melt")]) : 0.;
	   //20240701
	   //地幔分为上、中、下三层，默认下地幔和软流圈可发生熔融。
	   //后续可根据应用情况调整程序。
	   //20240905
	   //地壳熔融，判定是否在地壳内
	   const double judge_upper_value = this->introspection().compositional_name_exists("upper") ? premelt_volume_fractions[this->introspection().compositional_index_for_name("upper") + 1] : 0.;
	   const double judge_lower_value = this->introspection().compositional_name_exists("lower") ? premelt_volume_fractions[this->introspection().compositional_index_for_name("lower") + 1] : 0.;
	   //const double mantle_upper_value = this->introspection().compositional_name_exists("mantle_upper") ? premelt_volume_fractions[this->introspection().compositional_index_for_name("mantle_upper") + 1] : 0.;
	   const double judge_mantle_middle_value = this->introspection().compositional_name_exists("mantle_middle") ? premelt_volume_fractions[this->introspection().compositional_index_for_name("mantle_middle") + 1] : 0.;
	   //20260127 新增膏盐模式的亏损地幔
	   const double judge_mantle_lower_value = this->introspection().compositional_name_exists("mantle_lower") ? premelt_volume_fractions[this->introspection().compositional_index_for_name("mantle_lower") + 1] : this->introspection().compositional_name_exists("mantle_D") ? premelt_volume_fractions[this->introspection().compositional_index_for_name("mantle_D") + 1] : 0.;
	   const double judge_new_crust_value = this->introspection().compositional_name_exists("new_crust") ? premelt_volume_fractions[this->introspection().compositional_index_for_name("new_crust") + 1] : 0.;
	   const double judge_last_melt_fraction = this->introspection().compositional_name_exists("melt_fraction") ? std::abs(composition[this->introspection().compositional_index_for_name("melt_fraction")]) : 0.;
	   //默认软流圈地幔含水量最多
	   //含水量最多降低到30ppm
	   const double temp_water_content = this->introspection().compositional_name_exists("water_content") ? std::max(30. , std::abs(composition[this->introspection().compositional_index_for_name("water_content")])) : 0.;
	   //20240904 地壳成分占比
	   const double crust_fractions = judge_upper_value + judge_lower_value + judge_new_crust_value;
	   //20240904 地壳熔融程度
	   double crust_melt_fraction = 0.;
	   //20240701 地幔成分占比
	   const double mantle_fractions = judge_mantle_lower_value + judge_mantle_middle_value + premelt_volume_fractions[0];
	   //地幔熔融计算，迭代相关参数
	   //输出熔体分数，当温压条件无法发生熔融时，输出为0
	   double iterate_output_melt_fraction = 0;
	   //迭代初始步使用的初始熔体分数
	   double iterate_melt_fraction = judge_last_melt_fraction;
	   //迭代最大差值
	   double iterate_accuracy = 1e-5;
	   //迭代初始步使用的含水量
	   double iterate_water_content = temp_water_content;
	   //迭代过程中判定值为负的次数
	   double iterate_negative_time = 0;
	   //迭代过程中熔融过多的次数
	   double iterate_excess_time = 0;
	   //迭代次数记录
	   double iterate_time = 0;
	   //判定是否需要进行迭代
	   double iterate_depend = iterate_accuracy;
	   
	   //20250113 为最大熔融设置变量
	   const double max_carbon_melt_fraction = 0.02;
	   double T_carbon = 0;								  
	   //仅判断地幔物质的熔融程度
	   //20240703
	   //在压力小于10GPa（~330km）时计算熔融程度
	   //地幔熔融
	  if (include_mantle_melting && mantle_fractions > 0 && this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && pressure <= 10e9)
	  {
		    const double melt_CO2_wt = include_carbonate_melting ? X_CO2 : 0.;
		    const double T_carbon_2G = 19.21 * melt_CO2_wt + 1491.37 * std::log((100 - 0.86 * melt_CO2_wt) / 100);
		    const double T_carbon_3G = 27.04 * melt_CO2_wt + 1490.75 * std::log((100 - 1.18 * melt_CO2_wt) / 100);
		    const double T_carbon_4G = 31.90 * melt_CO2_wt + 1469.92 * std::log((100 - 1.31 * melt_CO2_wt) / 100);
		    const double T_carbon_5G = -5.01 * melt_CO2_wt + 1514.84 * std::log((100 + 1.23 * melt_CO2_wt) / 100);
		    if (pressure <= 2e9)
			 //  小于2GPa的区域，取2GPa的值
			    T_carbon = T_carbon_2G;
		    else if (pressure > 2e9 && pressure <= 3e9)
			    T_carbon = T_carbon_2G + (T_carbon_3G - T_carbon_2G) * (pressure - 2e9) / 1e9;
		    else if (pressure > 3e9 && pressure <= 4e9)
			    T_carbon = T_carbon_3G + (T_carbon_4G - T_carbon_3G) * (pressure - 3e9) / 1e9;
		    else if (pressure > 4e9 && pressure <= 5e9)
			    T_carbon = T_carbon_4G + (T_carbon_5G - T_carbon_4G) * (pressure - 4e9) / 1e9;
		   //对于大于5GPa（~165km以下）的区域，使用4~5GPa间的斜率
		    else
			    T_carbon = T_carbon_5G + (T_carbon_5G - T_carbon_4G) * (pressure - 5e9) / 1e9;
		   const double mantle_lower_fraction_in_melt = judge_mantle_lower_value + judge_mantle_middle_value > 0.02 ? 
		                                                judge_mantle_lower_value + judge_mantle_middle_value : 0;
		   //20250122 调整为30%含量
		   const double have_lower_mantle = mantle_lower_fraction_in_melt > 0.3 ? 1 : 0;
		   //20250109
		   //规定最大含碳熔融不超过2%
		   const double carbon_content = judge_temp_max_melt < max_carbon_melt_fraction ? 1 : 0;								  
		   //20240712
		   //使用熔融程度限制
	       //采用差值准则进行迭代计算
	       //首先判断是否需要进入迭代循环
		   const double X_water_0 = include_hydrous_melting ? temp_water_content * (61.2 / 1000) / D_H2O : 0.;
		   const double T_water_0 = K * std::pow(X_water_0 , gamma);
		   const double T_solidus_0  = A1 + 273.15
									+ A2 * pressure
									+ A3 * pressure * pressure
									- T_water_0
									- T_carbon * carbon_content * have_lower_mantle
									- mantle_lower_fraction_in_melt * hetero_temperature;
		   
		   while (temperature > T_solidus_0 && std::abs(iterate_depend) >= iterate_accuracy)
		   {
			   //迭代计数
			   iterate_time += 1;
			   //每次循环前恢复含水量
			   iterate_water_content = temp_water_content;
			   //使用更新过的迭代熔体分数计算含水量变化，以对比差值
			   //含水量最多降低到30ppm
			   
			   iterate_water_content = std::max(30. , iterate_water_content * std::pow(1 - std::max(0., iterate_melt_fraction - judge_last_melt_fraction), (1 - D_H2O) / D_H2O));
			   double X_water = include_hydrous_melting ? iterate_water_content * (61.2 / 1000) / D_H2O : 0.;
		       double T_water = K * std::pow(X_water , gamma);
			   // anhydrous melting of peridotite after Katz, 2003
		       double T_solidus  = A1 + 273.15
									  + A2 * pressure
								  	  + A3 * pressure * pressure
									  - T_water
									  - T_carbon * carbon_content * have_lower_mantle
									  - mantle_lower_fraction_in_melt * hetero_temperature;
			   //仅在温度大于T_solidus时继续计算
				if (temperature >= T_solidus)
				  {
					  double T_lherz_liquidus = B1 + 273.15
												  + B2 * pressure
												  + B3 * pressure * pressure
												  - T_water
												  - T_carbon * carbon_content * have_lower_mantle
												  - mantle_lower_fraction_in_melt * hetero_temperature;
					  double T_liquidus = C1 + 273.15
												+ C2 * pressure
												+ C3 * pressure * pressure
												- T_water
												- T_carbon * carbon_content * have_lower_mantle
												- mantle_lower_fraction_in_melt * hetero_temperature;
					  if (temperature > T_lherz_liquidus)
					  {
						  iterate_output_melt_fraction = 1.0;
						  break;
					  }
						
					  else
						iterate_output_melt_fraction = std::pow((temperature - T_solidus) / (T_lherz_liquidus - T_solidus),beta);

					  // melt fraction after melting of all clinopyroxene
					  double R_cpx = r1 + r2 * std::max(0.0, pressure);
					  double F_max = M_cpx / R_cpx;

					  if ( (!include_carbonate_melting 
					      || include_carbonate_melting && ((have_lower_mantle == 1 && T_carbon == 0 && carbon_content == 0) || have_lower_mantle == 0))
						  && iterate_output_melt_fraction > F_max && temperature < T_liquidus)
						{
						  double T_max = std::pow(F_max,1/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
						  iterate_output_melt_fraction = F_max + (1 - F_max) * pow((temperature - T_max) / (T_liquidus - T_max),beta);
						  break;
						}
						//判定是否需要继续迭代
						iterate_depend = iterate_output_melt_fraction - iterate_melt_fraction;
				  }

				  if(include_carbonate_melting && T_carbon > 0 && carbon_content > 0 && have_lower_mantle == 1 && iterate_output_melt_fraction >= max_carbon_melt_fraction)
				  {
					  iterate_output_melt_fraction = max_carbon_melt_fraction;
					  break;
				  }
				  
				  //迭代超过30次时,强制退出循环
				  //不包含含水熔融时，不进行循环计算
				  if (iterate_time > 30 || X_water == 0. || !include_hydrous_melting || std::abs(iterate_depend) < iterate_accuracy)
					  break;
					   
				  //差值小于0时，代表等式右侧熔融过多，减小熔融值
				  else if (iterate_depend < 0)
				   {
					   iterate_negative_time += 1;
					   iterate_melt_fraction =std::max( 0.5 * iterate_melt_fraction, iterate_melt_fraction + iterate_depend / (iterate_negative_time / 20 + 1));
				   }
				   //当温度小于固相线温度时，代表等式右侧熔融过多，减小熔融值
				   else if (temperature < T_solidus)
				  {
					   iterate_excess_time += 1;
					   iterate_melt_fraction =std::max(0.5 * iterate_melt_fraction, iterate_melt_fraction - iterate_depend / (iterate_excess_time / 20 + 1));
				  }
				  //其他情况代表等式右侧熔融过少，增加熔融值
				   else
					   iterate_melt_fraction += iterate_depend / (iterate_time / 20 + 1);

		   }
	   }
	 //地壳熔融
	 if (include_crustal_melting && crust_fractions > 0 && this->simulator_is_past_initialization() && this->get_timestep_number() > 0)
	 {
		 const double T_solidus_crust = A1_crust + 273.15
									  + A2_crust * pressure
								  	  + A3_crust * pressure * pressure;
		 if (temperature >= T_solidus_crust)
		 {
			 const double T_liquidus_crust = B1_crust + 273.15
									       + B2_crust * pressure
								  	       + B3_crust * pressure * pressure;
			 const double T_ss = (temperature - 0.5 * (T_solidus_crust + T_liquidus_crust)) / (T_liquidus_crust - T_solidus_crust);
			 crust_melt_fraction = 0.5 + T_ss + (std::pow(T_ss, 2) - 0.25) * (0.4256 + 2.988 * T_ss);
		 }
	 }
	 //循环内使用右侧熔融值计算含水量，在走出循环后，需要使用左侧熔融值进行含水量计算
	 //注：两侧差值小于1e-5，因此含水量差距不会很大
	 double out_water_content = temp_water_content;
	 double out_max_melt_mantle = include_mantle_melting && mantle_fractions > 0 ? judge_temp_max_melt : 0.;
	 double out_max_melt_crust = include_crustal_melting && crust_fractions > 0 ? judge_temp_max_melt : 0.;
	 if (iterate_output_melt_fraction > 0)
	 {
		 //含水量最多降低到30ppm
		 out_water_content = include_hydrous_melting ? std::max(30. , out_water_content * std::pow(1 - std::max(0., iterate_output_melt_fraction - judge_last_melt_fraction), (1 - D_H2O) / D_H2O)) : temp_water_content;
		 out_max_melt_mantle = std::max(iterate_output_melt_fraction, out_max_melt_mantle);
	 }
	 if (crust_melt_fraction > 0)
		 out_max_melt_crust = std::max(crust_melt_fraction, out_max_melt_crust);
	 //计算密度变化，地壳与地幔的计算方案不同，需分开处理
	 double out_density_change = 0.;
	 if (include_mantle_melting && ! include_crustal_melting)
		 out_density_change = out_max_melt_mantle * 72.6;
	 else if (!include_mantle_melting && include_crustal_melting)
		 out_density_change = out_max_melt_crust * gamma_crust * crustal_melt_temp_reduction;
	 //考虑到二者混合的概率较低，暂时取成分含量高的区域的熔融信息。
	 else
		 out_density_change = mantle_fractions >= crust_fractions ? 
	                          out_max_melt_mantle * 72.6
							: out_max_melt_crust * gamma_crust * crustal_melt_temp_reduction;
	 return {iterate_output_melt_fraction, out_max_melt_mantle, out_water_content, crust_melt_fraction, out_max_melt_crust, out_density_change};
	}
	
    template <int dim>
    void
    ViscoPlasticGrainSize<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      //20250902
	  //添加RHS
	  MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>
      *force = out.template get_additional_output<AdditionalMaterialOutputsStokesRHS<dim> >();
	  //create_additional_material_model_outputs(out);

	  // Store which components do not represent volumetric compositions (e.g. strain components).
      const ComponentMask volumetric_compositions = rheology->get_volumetric_composition_mask();

      EquationOfStateOutputs<dim> eos_outputs (this->n_compositional_fields()+1);
      EquationOfStateOutputs<dim> eos_outputs_all_phases (this->n_compositional_fields()+1+phase_function.n_phase_transitions());

      std::vector<double> average_elastic_shear_moduli (in.n_evaluation_points());

      // Store value of phase function for each phase and composition
      // While the number of phases is fixed, the value of the phase function is updated for every point
      std::vector<double> phase_function_values(phase_function.n_phase_transitions(), 0.0);

	  //熔体运输初始化
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
	  const double posx = in.position[0][0];
	  const double posy = in.position[0][1];
	  const unsigned int crood_x = std::floor( posx / ((max_x / repetitions_x) / std::pow(2, melt_fraction_refinement)));
	  const unsigned int crood_y = std::floor( posy / ((max_y / repetitions_y) / std::pow(2, melt_fraction_refinement)));
	  const int crood_y_in_flag_container = crood_y - min_crood_y;
	  
      // Loop through all requested points
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
      {
		  //20251104
	      //判定是否在岩墙内
	      bool in_dike = false;
		  //20251201 M值和phi值计算放于此处
		  double M_value = 0.;
		  double phi_dike = 0.; 
		  //20251114
		  //岩墙周围区域也不累积应变
		  bool strain_accumulate_banning = false;
		  bool convey_start_new_crust = false;
		  const double posx_dike = in.position[i][0];
		  const double posy_dike = in.position[i][1];
		  if(dike_mode)
		  {
			  //确定岩墙位置
			  const double emplace_xmax = rheology->read_from_txt(melt_info_path + "dike_emplace_location_xmax.txt", communicator_id);
			  const double emplace_ymin = rheology->read_from_txt(melt_info_path + "dike_emplace_location_ymin.txt", communicator_id);
			  const double total_melt_extract = rheology->read_from_txt(melt_info_path + "mantle_total_melt_extraction.txt", communicator_id);
			  
			  //20251201
			  const double emplace_ymax = constant_dike_length ? 0. : rheology->read_from_txt(melt_info_path + "dike_emplace_location_ymax.txt", communicator_id);
			  const double current_melt_extract = constant_M ? 0. : rheology->read_from_txt(melt_info_path + "mantle_current_melt_extraction.txt", communicator_id);
			  //判定条件：在岩墙内，且累计提取熔融超过设置的岩墙面积
			  //20251115 临时加入温度限制,累积应变调整不作温度限制
			  //20251118 取消温度限制
			  //20251201 岩墙长度判定
			  const double dike_length_value = constant_dike_length ? dike_length : (emplace_ymax - emplace_ymin) * 1e3;
			  //20251201 M值和phi值计算
			  M_value = constant_M ? M_for_dike : current_melt_extract * 1e6 / (dike_length_value * spreading_rate * this->get_timestep() / year_in_seconds);
			  //测试year_in_seconds能否在此处使用
			  phi_dike = M_value * spreading_rate / (dike_width * year_in_seconds);
			  //20260108 岩墙x位置判定修改：由最左侧端点改为中心点
			  in_dike = dike_length_value > 0 
			              //&& posx_dike < emplace_xmax * 1e3 && posx_dike >= emplace_xmax * 1e3 - dike_width
						  && posx_dike < emplace_xmax * 1e3 + 0.5 * dike_width && posx_dike >= emplace_xmax * 1e3 - 0.5 * dike_width
						  && posy_dike >= emplace_ymin * 1e3 && posy_dike < emplace_ymin * 1e3 + dike_length_value
						  && total_melt_extract * 1e6 >= dike_width * dike_length_value;
			  //考虑到岩墙长度较长，Y（Z）方向增加的不累积应变区域按照岩墙宽度设定
			  //20260108 不累积应变去扩张范围改为定值1.5e3
			  strain_accumulate_banning = dike_length_value > 0 
			              //&& posx_dike < emplace_xmax * 1e3 + dike_width && posx_dike >= emplace_xmax * 1e3 - 2 * dike_width
						  && posx_dike < emplace_xmax * 1e3 + 0.5 * dike_width + 1.5e3 && posx_dike >= emplace_xmax * 1e3 - 0.5 * dike_width - 1.5e3
						  && posy_dike >= emplace_ymin * 1e3 - 1.5e3 && posy_dike < emplace_ymin * 1e3 + dike_length_value + 1.5e3
						  && total_melt_extract * 1e6 >= dike_width * dike_length_value;
			  convey_start_new_crust = in_dike && melt_convey;
		  }
		  else
			  convey_start_new_crust = include_mantle_melting && crood_y_in_flag_container >= 0
									   ? get_convey_flag(path_mantle_melting, melt_convey, crood_x, crood_y_in_flag_container, communicator_id)
									   : false;

		  const bool convey_start_core_complex = !dike_mode && include_crustal_melting && crood_y_in_flag_container >= 0
											   ? get_convey_flag(path_crustal_melting, melt_convey, crood_x, crood_y_in_flag_container, communicator_id)
											   : false;

          // First compute the equation of state variables and thermodynamic properties
          equation_of_state.evaluate(in, i, eos_outputs_all_phases);

          const double gravity_norm = this->get_gravity_model().gravity_vector(in.position[i]).norm();
          const double reference_density = (this->get_adiabatic_conditions().is_initialized())
                                           ?
                                           this->get_adiabatic_conditions().density(in.position[i])
                                           :
                                           eos_outputs_all_phases.densities[0];

          // The phase index is set to invalid_unsigned_int, because it is only used internally
          // in phase_average_equation_of_state_outputs to loop over all existing phases
          MaterialUtilities::PhaseFunctionInputs<dim> phase_inputs(in.temperature[i],
                                                                   in.pressure[i],
                                                                   this->get_geometry_model().depth(in.position[i]),
                                                                   gravity_norm*reference_density,
                                                                   numbers::invalid_unsigned_int);

          // Compute value of phase functions
          for (unsigned int j=0; j < phase_function.n_phase_transitions(); j++)
            {
              phase_inputs.phase_index = j;
              phase_function_values[j] = phase_function.compute_value(phase_inputs);
            }

          // Average by value of gamma function to get value of compositions
          phase_average_equation_of_state_outputs(eos_outputs_all_phases,
                                                  phase_function_values,
                                                  phase_function.n_phase_transitions_for_each_composition(),
                                                  eos_outputs);

          const std::vector<double> temp_volume_fractions = MaterialUtilities::compute_composition_fractions(in.composition[i], volumetric_compositions);
		  //修正成分判断
		  //依据当前的设定，熔融程度超过2%的物质将被认为发生熔融作用
		  //后续将根据具体情况修正判定逻辑
		   const std::vector<double> volume_fractions = rheology->fractions_modify(temp_volume_fractions);
		   const std::vector<double> convey_volume_fractions = use_melt_model && melt_convey ? 
		                                              melt_convey_function(convey_start_new_crust, convey_start_core_complex, volume_fractions)
		                                            : volume_fractions;
	      //20250701
		  //调整成分场调用规则
		  const double mantle_upper_value = this->introspection().compositional_name_exists("mantle_upper") ? volume_fractions[this->introspection().compositional_index_for_name("mantle_upper") + 1] : 0.;
	      const double mantle_middle_value = this->introspection().compositional_name_exists("mantle_middle") ? volume_fractions[this->introspection().compositional_index_for_name("mantle_middle") + 1] : 0.;
	      const double mantle_lower_value = this->introspection().compositional_name_exists("mantle_lower") ? volume_fractions[this->introspection().compositional_index_for_name("mantle_lower") + 1] : 0.;
		  //判定熔融程度
		  //本模型守恒公式遵循visco-plastic模式,即不考虑流体的相关问题
		  //熔融程度降低时，判断条件待定
		  //20240905                             地幔熔融程度  地幔最大熔融程度  当前含水量  地壳熔融程度  地壳最大熔融程度  密度变化
		  std::array<double, 6> melt_factors = {       0,             0,             0,           0,              0,             0    };
		  
		  if (use_melt_model)
			  melt_factors = melted_judgement(in.temperature[i], std::abs(in.pressure[i]), in.composition[i], volume_fractions);
		  
		  // not strictly correct if thermal expansivities are different, since we are interpreting
          // these compositions as volume fractions, but the error introduced should not be too bad.
		  //熔融后当前位置的密度、热力学参数在此处更新
		  const double temp_density = MaterialUtilities::average_value (volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic);
		  //20241227 修正温度比Tp高的区域的软流圈密度
		  //20250901 添加参数：控制是否启用修正高Tp区域的的软流圈密度
		  const double density_Asthen_above_tp = density_modify_above_tp ? densities[0] * (1 - std::max(0., (in.temperature[i] - Tp)) * thermal_expansivities[0]) : temp_density;
		  const double density_for_Asten = in.temperature[i] >=  Tp ? density_Asthen_above_tp : temp_density;
		  //20250701
	      //调整成分场调用规则
		  //20260807 密度调整不包含mantle_lower
		  //const double Asten_phase = exclude_mantle_lower ? volume_fractions[0] + mantle_lower_value : volume_fractions[0];
		  const double Asten_phase = volume_fractions[0];
		  //20250417
		  //perplex修正密度
		  //在表格范围外的密度按照原规则计算，即由参考密度、温度和热膨胀系数控制
		  //20250512
		  //perplex新增湿地幔相
		  //20260310 perplex添加地幔相和玄武岩相开关
		  double out_density_perplex = 0;
		  //20260310 将洋壳相挪至perplex前
		  const double new_crust_phase = this->introspection().compositional_name_exists("new_crust") ? convey_volume_fractions[this->introspection().compositional_index_for_name("new_crust") + 1] : 0.;
		  if(use_perplex_table)
		  {
			  //20250701
			  //调整成分场调用规则
			  const double density_perplex_normal_mantle = perplex_include_mantle_phase ? perplex_table_normal->density(in.temperature[i],in.pressure[i]) : 0.;
		      const double density_perplex_carbon_mantle = perplex_include_mantle_phase ? perplex_table_carbon->density(in.temperature[i],in.pressure[i]) : 0.;
			  const double density_perplex_wet_mantle = perplex_include_mantle_phase ? perplex_table_wet->density(in.temperature[i],in.pressure[i]) : 0.;
			  //20260310 玄武岩相调用规则：将P,T限制在表格内
			  //后续考虑将地幔相调整为相同规则
			  //注：此为不得已的方式，尽量让模型的P,T条件在输入表格的范围内
			  const double perplex_basalt_T = std::max(basalt_table_min_T, std::min(basalt_table_max_T, in.temperature[i]));
			  const double perplex_basalt_P = std::max(basalt_table_min_P, std::min(basalt_table_max_P, in.pressure[i]));
			  const double density_perplex_basalt = perplex_include_basalt_phase ? perplex_table_basalt->density(perplex_basalt_T, perplex_basalt_P) : 0.;
			  if (perplex_include_mantle_phase && in.temperature[i] >= normal_table_min_T && in.temperature[i] <= normal_table_max_T
				&& in.pressure[i] >= normal_table_min_P && in.pressure[i] <= normal_table_max_P
				&& (mantle_upper_value > 0.5 || mantle_middle_value > 0.5))
				  out_density_perplex = density_perplex_normal_mantle;
			  else if(perplex_include_mantle_phase && in.temperature[i] >= carbon_table_min_T && in.temperature[i] <= carbon_table_max_T
					  && in.pressure[i] >= carbon_table_min_P && in.pressure[i] <= carbon_table_max_P
					  && mantle_lower_value > 0.5)
				  out_density_perplex = density_perplex_carbon_mantle;
		      else if(perplex_include_mantle_phase && in.temperature[i] >= wet_table_min_T && in.temperature[i] <= wet_table_max_T
				&& in.pressure[i] >= wet_table_min_P && in.pressure[i] <= wet_table_max_P
				&& volume_fractions[0] > 0.5)
				  out_density_perplex = density_perplex_wet_mantle;
			  else if(perplex_include_basalt_phase && new_crust_phase >= 0.2)
				  out_density_perplex = density_perplex_basalt;
			  else
				  out_density_perplex = temp_density;
		  }
		  
		  //20250121
		  const double density_for_output = Asten_phase > 0.5 ? density_for_Asten : temp_density;
		  //20260127 新增蛇纹岩化相关的密度更新，暂定蛇纹岩化最多可减少400kg/m3的密度，后续可作为输入值
		  //20260310 将洋壳相挪至perplex前
		  //const double new_crust_phase = this->introspection().compositional_name_exists("new_crust") ? convey_volume_fractions[this->introspection().compositional_index_for_name("new_crust") + 1] : 0.;
		  const double serp_value = this->introspection().compositional_name_exists("serp") ? in.composition[i][this->introspection().compositional_index_for_name("serp")] : 0.;
		  const double density_change_serp = serp_value > 0.1 && new_crust_phase < 0.2 ? serp_density_change * std::min(1., serp_value) : 0.;
		  //const double density_for_output = temp_density;
		  //20240614 用最大熔融程度计算密度变化
		  //20250428
		  //20260310 在不启用perplex地幔相调整时，需对洋壳相含量进行判定
		  if(use_perplex_table && (perplex_include_mantle_phase || (!perplex_include_mantle_phase && perplex_include_basalt_phase && new_crust_phase >= 0.2)))
			  out.densities[i] = out_density_perplex;
		  else
			  out.densities[i] = use_melt_model? density_for_output - melt_factors[5] - density_change_serp : density_for_output - density_change_serp;
          out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.thermal_expansion_coefficients, MaterialUtilities::arithmetic);
          out.specific_heat[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.specific_heat_capacities, MaterialUtilities::arithmetic);

          //20250701
	      //调整成分场调用规则
		  //20260121 由于蛇纹岩化统计需要累计塑性应变信息，将此部分移动至主函数内
		  const double plastic_strain = this->introspection().compositional_name_exists("plastic_strain") ? in.composition[i][this->introspection().compositional_index_for_name("plastic_strain")] : 0.;
          if (define_conductivities == false)
            {
              double thermal_diffusivity = 0.0;

              for (unsigned int j=0; j < volume_fractions.size(); ++j)
                thermal_diffusivity += volume_fractions[j] * thermal_diffusivities[j];
			//20250701
			//调整成分场调用规则
			//20260121 由于蛇纹岩化统计需要累计塑性应变信息，将此部分移动至主函数内
			//const double plastic_strain = this->introspection().compositional_name_exists("plastic_strain") ? in.composition[i][this->introspection().compositional_index_for_name("plastic_strain")] : 0.;
			//判断是否需要进行导热系数k的修正
			//20251204 测试：岩墙内启用k_modify_active
			const bool k_modify_active = k_modify_order && (plastic_strain > 0.1 || in_dike);
			//20241230 所有地幔相物质热导率受温压控制
			//热导率取α*ρ*C与k(T,P)间较小的值
			//数据来源：Clauser and Huenges，1995；Liao and Gerya，2013&2014
			//可在prm文件中设置是否启用
			const double k_profile_mantle_T_P = (0.73 + 1293 / (in.temperature[i] + 77)) * std::exp(4e-5 * in.pressure[i] / 1e6);
			//20250701
			//调整成分场调用规则
			//20251115 热液循环不缩小范围
			//20251118 测试：岩墙内热液循环不作限制
			const double mantle_phase_for_k = volume_fractions[0] + mantle_upper_value + mantle_middle_value + mantle_lower_value;
			  //if (!k_modify_active || (k_modify_active && in_dike))
			  if (!k_modify_active)
			  {
				  // Thermal conductivity at the given positions. If the temperature equation uses
                  // the reference density profile formulation, use the reference density to
                  // calculate thermal conductivity. Otherwise, use the real density. If the adiabatic
                  // conditions are not yet initialized, the real density will still be used.
				  //20250122
				  //温压条件启动时，对所有地幔相使用
                  if (this->get_parameters().formulation_temperature_equation ==
                      Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile &&
                      this->get_adiabatic_conditions().is_initialized())
                    out.thermal_conductivities[i] = use_k_modified_by_T_P && mantle_phase_for_k > 0.2 ? 
					                                //std::min(thermal_diffusivity * out.specific_heat[i] * this->get_adiabatic_conditions().density(in.position[i]), k_profile_mantle_T_P)
													k_profile_mantle_T_P
												  : thermal_diffusivity * out.specific_heat[i] * this->get_adiabatic_conditions().density(in.position[i]);
                  else
                    out.thermal_conductivities[i] = use_k_modified_by_T_P && mantle_phase_for_k > 0.2 ? 
				                                    //std::min(thermal_diffusivity * out.specific_heat[i] * out.densities[i], k_profile_mantle_T_P)
													k_profile_mantle_T_P
											      : thermal_diffusivity * out.specific_heat[i] * out.densities[i];
			  }
			  //20251104
			  //热液循环在dike内不使用
			  //20251114，增大热液循环取消范围
			  //20251118 测试：岩墙内热液循环不作限制
			  //else if (k_modify_active && !in_dike)
			  else if (k_modify_active)
			  {
				  const bool use_reference_strainrate = (this->get_timestep_number() == 0) && (in.strain_rate[i].norm() <= std::numeric_limits<double>::min());
				  double edot_ii;
			      if (use_reference_strainrate)
			        edot_ii = ref_strain_rate;
		          else
			        edot_ii = std::max(std::sqrt(std::fabs(second_invariant(in.strain_rate[i]))), min_strain_rate);
				  const double active_factor = 1 + (k_factor - 1) * (1 - std::exp(- edot_ii / critial_edot_ii));
				  if (this->get_parameters().formulation_temperature_equation ==
                      Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile &&
                      this->get_adiabatic_conditions().is_initialized())
                    out.thermal_conductivities[i] = use_k_modified_by_T_P && mantle_phase_for_k > 0.2 ? 
					                                //std::min(active_factor * thermal_diffusivity * out.specific_heat[i] * this->get_adiabatic_conditions().density(in.position[i]), active_factor * k_profile_mantle_T_P)
													active_factor * k_profile_mantle_T_P
												  : active_factor * thermal_diffusivity * out.specific_heat[i] * this->get_adiabatic_conditions().density(in.position[i]);
                  else
                    out.thermal_conductivities[i] = use_k_modified_by_T_P && mantle_phase_for_k > 0.2 ? 
				                                    //std::min(active_factor * thermal_diffusivity * out.specific_heat[i] * out.densities[i], active_factor * k_profile_mantle_T_P)
													active_factor * k_profile_mantle_T_P
				                                  : active_factor * thermal_diffusivity * out.specific_heat[i] * out.densities[i];
			  }

            }
          else
            {
              // Use thermal conductivity values specified in the parameter file, if this
              // option was selected.
              out.thermal_conductivities[i] = MaterialUtilities::average_value (volume_fractions, thermal_conductivities, MaterialUtilities::arithmetic);
            }

          out.compressibilities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.compressibilities, MaterialUtilities::arithmetic);
          out.entropy_derivative_pressure[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.entropy_derivative_pressure, MaterialUtilities::arithmetic);
          out.entropy_derivative_temperature[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.entropy_derivative_temperature, MaterialUtilities::arithmetic);

          // Compute the effective viscosity if requested and retrieve whether the material is plastically yielding



          bool plastic_yielding = false;
		  
          if (in.requests_property(MaterialProperties::viscosity))
            {
              // Currently, the viscosities for each of the compositional fields are calculated assuming
              // isostrain amongst all compositions, allowing calculation of the viscosity ratio.
              // TODO: This is only consistent with viscosity averaging if the arithmetic averaging
              // scheme is chosen. It would be useful to have a function to calculate isostress viscosities.
              const IsostrainViscosities isostrain_viscosities =
                rheology->calculate_isostrain_viscosities(in, i, volume_fractions, phase_function_values, phase_function.n_phase_transitions_for_each_composition());
              // The isostrain condition implies that the viscosity averaging should be arithmetic (see above).
              // We have given the user freedom to apply alternative bounds, because in diffusion-dominated
              // creep (where n_diff=1) viscosities are stress and strain-rate independent, so the calculation
              // of compositional field viscosities is consistent with any averaging scheme.
			  
			  //20260112 调整岩墙内粘度，由流变学计算部分移至此处，目的为减少岩墙判定次数
              //out.viscosities[i] = MaterialUtilities::average_value(volume_fractions, isostrain_viscosities.composition_viscosities, rheology->viscosity_averaging);
			  const double temp_viscosities = MaterialUtilities::average_value(volume_fractions, isostrain_viscosities.composition_viscosities, rheology->viscosity_averaging);
			  out.viscosities[i] = in_dike ? std::min(temp_viscosities , 1e21) : temp_viscosities;
			  
              // Decide based on the maximum composition if material is yielding.
              // This avoids for example division by zero for harmonic averaging (as plastic_yielding
              // holds values that are either 0 or 1), but might not be consistent with the viscosity
              // averaging chosen.
              std::vector<double>::const_iterator max_composition = std::max_element(volume_fractions.begin(),volume_fractions.end());
              plastic_yielding = isostrain_viscosities.composition_yielding[std::distance(volume_fractions.begin(),max_composition)];

              // Compute viscosity derivatives if they are requested
              if (MaterialModel::MaterialModelDerivatives<dim> *derivatives =
                    out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >())
                rheology->compute_viscosity_derivatives(i, volume_fractions, isostrain_viscosities.composition_viscosities, in, out, phase_function_values, phase_function.n_phase_transitions_for_each_composition());
            }

          // Now compute changes in the compositional fields (i.e. the accumulated strain).
          for (unsigned int c = 0; c < in.composition[i].size(); ++c)
              out.reaction_terms[i][c] = 0.0;

          // Calculate changes in strain invariants and update the reaction terms
		  //20251110
		  //新增岩墙相关
		  //20251114 扩大不累积应变的区域
		  if(dike_mode)
			  //rheology->dikemode_fill_reaction_outputs(in, i, rheology->min_strain_rate, plastic_yielding, in_dike, strain_rate_dev, out);
			  rheology->dikemode_fill_reaction_outputs(in, i, rheology->min_strain_rate, plastic_yielding, strain_accumulate_banning, in.strain_rate[i], out);
		  else
			  rheology->strain_rheology.fill_reaction_outputs(in, i, rheology->min_strain_rate, plastic_yielding, out);
		  
		  //20260403 新增塑性应变线性恢复
		  //规则：在应变速率小于1e-15的区域，根据dt线性恢复
		  if(plastic_linear_healing && this->introspection().compositional_name_exists("plastic_strain") && this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && in.requests_property(MaterialProperties::reaction_terms))
		  {
			  const double edot_ii = std::max(std::sqrt(std::fabs(second_invariant(in.strain_rate[i]))), min_strain_rate);
			  const double plastic_healing = edot_ii >= 1e-15 ? 0. : plastic_healing_rate * this->get_timestep();
			  const double temp_plastic_strain = plastic_strain + out.reaction_terms[i][this->introspection().compositional_index_for_name("plastic_strain")];
			  //20260415 临时添加
			  if(plastic_healing_proportion == 0.)
				  out.reaction_terms[i][this->introspection().compositional_index_for_name("plastic_strain")] += std::max(-temp_plastic_strain, -plastic_healing);
			  else
				  out.reaction_terms[i][this->introspection().compositional_index_for_name("plastic_strain")] += -(temp_plastic_strain * this->get_timestep() / (this->get_timestep() + plastic_healing_proportion));
			  
		  }
          
          if (this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && in.requests_property(MaterialProperties::reaction_terms))
		  {
			  //20240617
			  //20260114 根据模式判定
			  const double posx = in.position[i][0];
			  const unsigned int layer_1_id = salt_mode ? this->introspection().compositional_index_for_name("sediment_init") : this->introspection().compositional_index_for_name("upper");
			  const unsigned int layer_2_id = salt_mode ? this->introspection().compositional_index_for_name("salt") : this->introspection().compositional_index_for_name("lower");
			  //20240701
			  //20260112 可删除mantle_middle或mantle_lower
			  const unsigned int layer_3_id = salt_mode ? this->introspection().compositional_index_for_name("crust") : this->introspection().compositional_index_for_name("mantle_upper");
			  const unsigned int layer_4_id = salt_mode ? this->introspection().compositional_index_for_name("mantle_L") : this->introspection().compositional_name_exists("mantle_middle") ? this->introspection().compositional_index_for_name("mantle_middle") : 9999;
			  const unsigned int layer_5_id = salt_mode && this->introspection().compositional_name_exists("mantle_D") ? this->introspection().compositional_index_for_name("mantle_D") : !salt_mode && this->introspection().compositional_name_exists("mantle_lower")? this->introspection().compositional_index_for_name("mantle_lower") : 9999;
			  const unsigned int water_content_id = this->introspection().compositional_index_for_name("water_content");
			  //20260203 可删除rift_grain_size
			  const unsigned int grain_size_id = this->introspection().compositional_name_exists("rift_grain_size") ? this->introspection().compositional_index_for_name("rift_grain_size") : 9999;
			  //20241017 填充规则修改
			  //20241029 暂时默认俯冲输入方向为右侧，输入物质的结构按照cratontic_thicknesses的结构设定
			  //如需调整，可新增输入判定条件
			  //20250109 添加判定：是否含有克拉通成分
			  //20250109 现程序默认将两侧物质重置，在可视化时可以将重置区域截掉
			  //如有单侧的需求，可按情况添加判定条件
			  const std::vector<double> input_layer_thickness = add_cratontic ? cratontic_thicknesses : layer_thickness;
			  //if (input_from_flank && max_y - posy <= input_depth && (posx <= 10e3 || max_x - posx <= 10e3))
			  const double depth = this->get_geometry_model().depth(in.position[i]);
		      if (input_from_flank && depth <= input_depth && (posx <= 10e3 || max_x - posx <= 10e3))
			  {
				  
				  const double layer_1 = input_layer_thickness[0];
	              const double layer_2 = input_layer_thickness[0] + input_layer_thickness[1];
	              const double layer_3 = input_layer_thickness[0] + input_layer_thickness[1] + input_layer_thickness[2];
	              const double layer_4 = input_layer_thickness[0] + input_layer_thickness[1] + input_layer_thickness[2] + input_layer_thickness[3];
	              const double layer_5 = input_layer_thickness[0] + input_layer_thickness[1] + input_layer_thickness[2] + input_layer_thickness[3] + input_layer_thickness[4];
				  for (unsigned int a = 0; a < in.composition[i].size(); ++a) 
				  {
					  if (a == water_content_id)
						  //填入含水量，除背景值外默认为50
						  {
							  if(depth <= layer_5)
								  out.reaction_terms[i][a] = 50. - in.composition[i][a];
							  else
								  out.reaction_terms[i][a] = init_water_content[0] - in.composition[i][a];
						  }
					  else if(a == grain_size_id)
					      //20241104 粒度修正在粒度计算中完成
						  out.reaction_terms[i][a] = 0.;
					  else if(
					           //上地壳
							   //20250109 深度等于0的区域也填充上地壳组分
					          //(a == upper_id && depth > 0 && depth <= layer_1)
							  (a == layer_1_id && depth <= layer_1)
					         ||
						       //下地壳
						      (a == layer_2_id && depth > 0 && depth > layer_1 && depth <= layer_2)
                             ||
							   //上地幔
				              (a == layer_3_id && depth > 0 && depth > layer_2 && depth <= layer_3)
							 ||
							   //中地幔
							  (a == layer_4_id && depth > 0 && depth > layer_3 && depth <= layer_4)
							 ||
							   //下地幔
							  (a == layer_5_id && depth > 0 && depth > layer_4 && depth <= layer_5)
							 )
						  out.reaction_terms[i][a] = 1. - in.composition[i][a];
					  else
						  //两侧输入时不填入沉积物质
						  out.reaction_terms[i][a] = 0. - in.composition[i][a];
				  }
			  }
			  //模型在compositional fields中新增“rift_grain_size”作为存储粒径以及实现可视化的媒介，初始化结束后，在此处更新每一步粒径的变更情况。
		      //初始化部分在rift/source/comp/lithosphere_rift.cc中。
			  //20260203 可删除rift_grain_size
			  if(this->introspection().compositional_name_exists("rift_grain_size"))
				  rheology->fill_grainsize_outputs(in, out, i, plastic_yielding, volume_fractions, phase_function_values, phase_function.n_phase_transitions_for_each_composition());
			  //模型的熔融状态在此处更新
			  //20241230新增
			  //当模型底部为物质输入窗口时，对模型底边界含水量修正
			  //20260506 修改：当熔融和含水熔融同时关闭时修正含水量
			  if ((!use_melt_model && !include_hydrous_melting) || (use_melt_model && !input_from_flank && posy <= banned_depth_meter))
			  {
				  double out_water_no_melt = 0;
				  for (unsigned int a = 0; a < volume_fractions.size(); a++)
				  {
					  out_water_no_melt += volume_fractions[a] * init_water_content[a];
				  }
					  
				  out.reaction_terms[i][water_content_id] = out_water_no_melt - in.composition[i][water_content_id];
			  }
			  
              //在熔融状态发生变化时，更新相关参数
			  else if (use_melt_model)
			  {
			      const unsigned int melt_id = this->introspection().compositional_index_for_name("melt_fraction");
				  const unsigned int melt_record_id = this->introspection().compositional_index_for_name("melt_record");
				  const unsigned int max_melt_id = this->introspection().compositional_index_for_name("max_melt");
				  
				  //更新熔融状态
				  double out_melt_fraction = 0.;
				  double out_max_melt = 0.;
				  //只有地幔熔融
				  if (include_mantle_melting && !include_crustal_melting)
				  {
					  out_melt_fraction = melt_factors[0];
					  out_max_melt = melt_factors[1];
				  }
				  //只有地壳熔融
				  //20260114 地壳熔融暂时与岩盐模式不适配
				  else if (!salt_mode && !include_mantle_melting && include_crustal_melting) 
				  {
					  out_melt_fraction = melt_factors[3];
					  out_max_melt = melt_factors[4];
				  }
				  //混合熔融
				  //需要进一步测试，主要的问题在moho面区域
				  //考虑到二者混合的概率较低，暂时取成分含量高的区域的熔融信息。
				  //考虑到相变的问题，最大熔融取值情况可能需更改
				  else if (!salt_mode && include_mantle_melting && include_crustal_melting)
				  {
					  const double crust_fraction = volume_fractions[layer_1_id + 1] + volume_fractions[layer_2_id + 1];
					  //20260112
					  const double mantle_middle_fraction = this->introspection().compositional_name_exists("mantle_middle") ? volume_fractions[layer_4_id + 1] : 0.;
					  const double mantle_lower_fraction = this->introspection().compositional_name_exists("mantle_lower") ? volume_fractions[layer_5_id + 1] : 0.;
					  //const double mantle_fraction = volume_fractions[layer_4_id + 1] + volume_fractions[layer_5_id + 1] + volume_fractions[0];
					  const double mantle_fraction = mantle_middle_fraction + mantle_lower_fraction + volume_fractions[0];
					  out_melt_fraction = mantle_fraction >= crust_fraction ? melt_factors[0] : melt_factors[3];
					  out_max_melt = mantle_fraction >= crust_fraction ? melt_factors[1] : melt_factors[4];
				  }
				  else
				  {
					  out_melt_fraction = 0.;
				      out_max_melt = 0.;
				  }
				  //20241105
				  //物质输入区域不启用熔融功能
				  out.reaction_terms[i][melt_id] = input_from_flank && max_x - posx <= 10e3 ? 
				                                   0. - in.composition[i][melt_id]
				                                 : out_melt_fraction - in.composition[i][melt_id];
				  out.reaction_terms[i][max_melt_id] = input_from_flank && max_x - posx <= 10e3 ? 
				                                       0. - in.composition[i][max_melt_id]
				                                     : out_max_melt - in.composition[i][max_melt_id];
				  out.reaction_terms[i][melt_record_id] = input_from_flank && max_x - posx <= 10e3 ? 
				                                          0. - in.composition[i][melt_record_id]
														: out_melt_fraction - in.composition[i][melt_id] - in.composition[i][melt_record_id];
				  //含水量更新
				  //20260117 修改逻辑
				  //if (include_hydrous_melting && posy > banned_depth_meter && (!input_from_flank || (input_from_flank && max_x - posx >= 10e3)))
				  if (include_hydrous_melting && posy > banned_depth_meter && (!input_from_flank || (input_from_flank && posx > 10e3 && max_x - posx >= 10e3)))
				  {
					  out.reaction_terms[i][water_content_id] = melt_factors[2] - in.composition[i][water_content_id];
					  //20240614
				      //out.reaction_terms[i][CO2_content_id] = melt_factors[3] - in.composition[i][CO2_content_id];
				  }
				  
				  //熔体运输功能
				  if(melt_convey)
				  {
					  for (unsigned int d = 0; d < in.composition[i].size(); d++)
					  {
						  if (convey_volume_fractions[d+1] != volume_fractions[d+1])
							  out.reaction_terms[i][d] = convey_volume_fractions[d+1] - in.composition[i][d];
					  }
				  }
			  }
			  
			  //20260121 新增蛇纹岩化统计，仅发生在地幔相，不参与流变学计算
			  if(include_serp)
			  {
				  //默认开启蛇纹岩化统计时存在serp与d_serp
				  const unsigned int serp_id = this->introspection().compositional_index_for_name("serp");
				  const unsigned int d_serp_id = this->introspection().compositional_index_for_name("d_serp");
				  //计算当前区域地幔相占比，蛇纹岩化仅允许在地幔物质内进行
				  //使用洋壳形成判定后的相态占比
				  //软流圈
				  const double serp_mantle_0 = convey_volume_fractions[0];
				  //第一层地幔，盐岩模式与其他模式进行区分
				  const double serp_mantle_1 = salt_mode ? convey_volume_fractions[layer_4_id + 1] : convey_volume_fractions[layer_3_id + 1];
				  //其他两层地幔仅出现在非盐岩模式，后续可根据实际情况修改
				  const double serp_mantle_2 = !salt_mode && layer_4_id <= convey_volume_fractions.size() ? convey_volume_fractions[layer_4_id + 1] : 0.;
				  //20260127 膏盐模式新增一层亏损地幔，现在只需判断layer_5_id的值
				  const double serp_mantle_3 = layer_5_id <= convey_volume_fractions.size() ? convey_volume_fractions[layer_5_id + 1] : 0.;
				  
				  const double serp_mantle_phase = serp_mantle_0 + serp_mantle_1 + serp_mantle_2 + serp_mantle_3;
				  //获取当前蛇纹岩化程度
				  const double temp_serp = std::min(1., in.composition[i][serp_id]);
				  //计算当前蛇纹岩化速度，累计塑性应变之前已经提取过，单位1/s
				  const double serp_edot_ii = std::max(std::sqrt(std::fabs(second_invariant(in.strain_rate[i]))), min_strain_rate);
				  //20260814 新增Christiansen et al., 2026, NC蛇纹岩化计算
				  double current_d_serp = 0.;
				  switch(serp_pattern)
				  {
					  case 1:
					  {
						  //20260814 Liu et al.,2023,Geology
						  //20260127 添加洋壳相关判定,new_crust_phase赋值在密度输出部分
						  current_d_serp = in.temperature[i] <= 673 && serp_mantle_phase >= 0.8 && new_crust_phase < 0.2 && serp_edot_ii > serp_strain_rate_crit && plastic_strain > serp_plastic_strain_crit 
										  ? serp_A * std::exp(-serp_B * std::pow(in.temperature[i] - 273 - serp_C, 2)) : 0.;
						  break;
					  }
					  case 2:
					  {
						  //20260814 Christiansen et al., 2026, NC，速率相关计算参考文章中加利福尼亚区域的方辉橄榄岩结果
						  //仅在蛇纹岩化窗口内计算，其他情况直接置0
						  if(in.temperature[i] >= 125 + 273 && in.temperature[i] <= 325 + 273 && plastic_strain > serp_plastic_strain_crit && serp_mantle_phase >= 0.8 && new_crust_phase < 0.2)
						  {
							  //应变速率相关项，反映每日水流量的影响
							  //文章中2e6 L/d的结果为基准值
							  //应变速率min_strain_rate_christiansen对应文章中2e6 L/d的结果
							  //应变速率serp_strain_rate_crit对应文章中5e7 L/d的结果
							  double fatcor_waterflow = 0.;
							  if(serp_edot_ii > serp_strain_rate_crit)
								  //当应变速率大于等于1e-14时，认为蛇纹岩化速率到达顶峰
								  //不考虑水流量过高导致的蛇纹岩化速率下降
								  fatcor_waterflow = 6.;
							  else if(serp_edot_ii > min_strain_rate_christiansen && serp_edot_ii <= serp_strain_rate_crit)
								  //应变速率范围内，根据应变速率的指数线性增长
								  fatcor_waterflow = 1. + 5 * (std::log10(std::min(serp_edot_ii, serp_strain_rate_crit)) - std::log10(min_strain_rate_christiansen)) / (std::log10(serp_strain_rate_crit) - std::log10(min_strain_rate_christiansen));
							  else
								  //其他情况认为不发生蛇纹岩化反应
								  fatcor_waterflow = 0.;
							  //速率计算，与T相关的部分根据Christiansen et al., 2026, NC的结果进行曲线拟合
							  current_d_serp = fatcor_waterflow * std::pow(10, serp_A_christiansen * std::pow(in.temperature[i] - 273, 2) + serp_B_christiansen * (in.temperature[i] - 273) + serp_C_christiansen);
						  }
						  break;
					  }
					  default:
					  {
						  AssertThrow(false, ExcNotImplemented());
						  break;
					  }
				  }
				  
				  //计算当前步累计蛇纹岩化，蛇纹岩化速率可在paraview中根据此参数计算得出
				  //蛇纹岩化程度不超过1
				  const double temp_d_serp_timestep = current_d_serp * this->get_timestep();
				  const double d_serp_timestep = temp_d_serp_timestep + temp_serp > 1 ? 1 - temp_serp : temp_d_serp_timestep;
				  //每年蛇纹岩化速率
				  const double d_serp_yr = d_serp_timestep * year_in_seconds / this->get_timestep();
				  //输出
				  //20260126 挤压背景下边界区域置0
				  out.reaction_terms[i][d_serp_id] = input_from_flank && depth <= input_depth && (posx <= 10e3 || max_x - posx <= 10e3) ? 0. -  in.composition[i][d_serp_id] : d_serp_yr - in.composition[i][d_serp_id];
				  out.reaction_terms[i][serp_id] = input_from_flank && depth <= input_depth && (posx <= 10e3 || max_x - posx <= 10e3) ?  0. -  in.composition[i][serp_id] : d_serp_timestep;
			  }
		  }
          // Fill plastic outputs if they exist.
          rheology->fill_plastic_outputs(i,volume_fractions,plastic_yielding,in,out, phase_function_values, phase_function.n_phase_transitions_for_each_composition());
		  
		  

          if (rheology->use_elasticity)
            {
              // Compute average elastic shear modulus
              average_elastic_shear_moduli[i] = MaterialUtilities::average_value(volume_fractions,
                                                                                 rheology->elastic_rheology.get_elastic_shear_moduli(),
                                                                                 rheology->viscosity_averaging);

              // Fill the material properties that are part of the elastic additional outputs
              if (ElasticAdditionalOutputs<dim> *elastic_out = out.template get_additional_output<ElasticAdditionalOutputs<dim> >())
                {
                  elastic_out->elastic_shear_moduli[i] = average_elastic_shear_moduli[i];
                }
            }
		  if(this->get_parameters().enable_additional_stokes_rhs && force)
		  {
			  force->rhs_u[i][1] = 0.;
			  force->rhs_p[i] = 0.;
			  if(in_dike)
				  //20251201 负值代表膨胀
				  force->rhs_p[i] = -phi_dike;
			  else
				  force->rhs_p[i] = 0.;
		  }
		  
        }
      // If we use the full strain tensor, compute the change in the individual tensor components.
      rheology->strain_rheology.compute_finite_strain_reaction_terms(in, out);

      if (rheology->use_elasticity)
        {
          rheology->elastic_rheology.fill_elastic_force_outputs(in, average_elastic_shear_moduli, out);
          rheology->elastic_rheology.fill_reaction_outputs(in, average_elastic_shear_moduli, out);
        }
    }



    template <int dim>
    double
    ViscoPlasticGrainSize<dim>::
    reference_viscosity () const
    {
      return rheology->ref_visc;
    }



    template <int dim>
    bool
    ViscoPlasticGrainSize<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible();
    }



    template <int dim>
    double ViscoPlasticGrainSize<dim>::
    get_min_strain_rate () const
    {
      return rheology->min_strain_rate;
    }

    template <int dim>
    void
    ViscoPlasticGrainSize<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Visco Plastic Grain Size");
        {
          MaterialUtilities::PhaseFunction<dim>::declare_parameters(prm);

          EquationOfState::MulticomponentIncompressible<dim>::declare_parameters (prm);

          Rheology::ViscoPlasticGrainSize<dim>::declare_parameters(prm);

          // Equation of state parameters
          prm.declare_entry ("Thermal diffusivities", "0.8e-6",
                             Patterns::List(Patterns::Double (0.)),
                             "List of thermal diffusivities, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  "
                             "Units: \\si{\\meter\\squared\\per\\second}.");
          prm.declare_entry ("Define thermal conductivities","false",
                             Patterns::Bool (),
                             "Whether to directly define thermal conductivities for each compositional field "
                             "instead of calculating the values through the specified thermal diffusivities, "
                             "densities, and heat capacities. ");
          prm.declare_entry ("Thermal conductivities", "3.0",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal conductivities, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
		  prm.declare_entry ("Thermal conductivities modification","true",
                             Patterns::Bool (),
                             "Whether to modify thermal conductivities of under-sea-water place");
          prm.declare_entry ("Factor of thermal conductivities modification", "14",
                             Patterns::Double(0),
                             "If thermal conductivities of current location needs to be modefied, multiply the factor."
                             "Units: none.");
		  prm.declare_entry ("Critial strain rate", "3e-14",
                             Patterns::Double(0),
							 "Critial strain rate."
                             "Units: \\si{\\none\\per\\second}.");
		//熔融新增
		  prm.declare_entry ("Use melt model", "false",
                            Patterns::Bool (),
                             "Use melt model "
                             "Units: none.");		
		  prm.declare_entry	("Include mantle melting", "false",
                            Patterns::Bool (),
                             "Include mantle melting."
                             "Units: none.");	
		  prm.declare_entry ("A1", "1085.7",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the solidus "
                             "of peridotite. "
                             "Units: \\si{\\degreeCelsius}.");
          prm.declare_entry ("A2", "1.329e-7",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal}.");
          prm.declare_entry ("A3", "-5.1e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
          prm.declare_entry ("B1", "1475.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the lherzolite "
                             "liquidus used for calculating the fraction "
                             "of peridotite-derived melt. "
                             "Units: \\si{\\degreeCelsius}.");
          prm.declare_entry ("B2", "8.0e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-derived melt. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal}.");
          prm.declare_entry ("B3", "-3.2e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
          prm.declare_entry ("C1", "1780.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the liquidus "
                             "of peridotite. "
                             "Units: \\si{\\degreeCelsius}.");
          prm.declare_entry ("C2", "4.50e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal}.");
          prm.declare_entry ("C3", "-2.0e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
		  prm.declare_entry ("beta", "1.5",
                             Patterns::Double (),
                             "Exponent of the melting temperature in "
                             "the melt fraction calculation. "
                             "Units: non-dimensional.");
		  prm.declare_entry ("r1", "0.5",
                             Patterns::Double (),
                             "Constant in the linear function that "
                             "approximates the clinopyroxene reaction "
                             "coefficient. "
                             "Units: non-dimensional.");
          prm.declare_entry ("r2", "8e-11",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the linear function that approximates "
                             "the clinopyroxene reaction coefficient. "
                             "Units: \\si{\\per\\pascal}.");
		  prm.declare_entry ("Mass fraction cpx", "0.15",
                             Patterns::Double (),
                             "Mass fraction of clinopyroxene in the "
                             "peridotite to be molten. "
                             "Units: non-dimensional.");
		  prm.declare_entry("Water distribution coefficient DH2O", "0.01",
		                    Patterns::Double (),
							"We use a constant distribution coefficient to"
							"control the loss rate of water."
							"Units:none");
		  prm.declare_entry ("Solidus change of heterogeneous mantle", "100",
                            Patterns::Double (),
                             "Solidus change of heterogeneous mantle "
                             "Units: °C");				 
          prm.declare_entry ("Initial water content", "50",
                           Patterns::List(Patterns::Double (0.)),
                           "List of initial water content, "
                           "for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "Units: \\si{\\part\\per\\million}.");
		  prm.declare_entry	("Use melt convey function", "false",
                            Patterns::Bool (),
                             "Use melt convey function."
                             "Units: none.");	 
		  prm.declare_entry	("Include hydrous melting", "false",
                            Patterns::Bool (),
                             "Include hydrous melting."
                             "Units: none.");	
		  prm.declare_entry("K", "4.3e-2",
		                    Patterns::Double (),
							"Practor of water content in melt"
							"Units:none");
		  prm.declare_entry("gamma", "0.75",
		                    Patterns::Double (),
							"Expont of water content in melt"
							"Units:none");
		  prm.declare_entry	("Include carbonate melting", "false",
                            Patterns::Bool (),
                             "Include carbonate melting."
                             "Units: none.");	
		  prm.declare_entry ("CO2 content in the melt", "10",
                            Patterns::Double (),
                             "CO2 content in the melt "
                             "Units: wt%");
          prm.declare_entry	("Include crustal melting", "false",
                            Patterns::Bool (),
                             "Include crustal melting."
                             "Units: none.");	
		  prm.declare_entry ("A1 for crustal melting", "720.0",
                             Patterns::Double (),
                             "A1 for crustal melting "
                             "Units: \\si{\\degreeCelsius}.");
          prm.declare_entry ("A2 for crustal melting", "-1.2e-7",
                             Patterns::Double (),
                             "A2 for crustal melting "
                             "Units: \\si{\\degreeCelsius\\per\\pascal}.");
          prm.declare_entry ("A3 for crustal melting", "1.2e-16",
                             Patterns::Double (),
                             "A3 for crustal melting "
                             "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
          prm.declare_entry ("B1 for crustal melting", "1220.0",
                             Patterns::Double (),
                             "B1 for crustal melting "
                             "Units: \\si{\\degreeCelsius}.");
          prm.declare_entry ("B2 for crustal melting", "-1.2e-7",
                             Patterns::Double (),
                             "B2 for crustal melting "
                             "Units: \\si{\\degreeCelsius\\per\\pascal}.");
          prm.declare_entry ("B3 for crustal melting", "1.6e-16",
                             Patterns::Double (),
                             "B3 for crustal melting "
                             "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
		  prm.declare_entry ("Crustal melting expansion factor", "0.13",
		                     Patterns::Double (),
							 "Crustal melting expansion factor"
							 "Units: none.");
							 
		  //20240617
		  prm.declare_entry	("Put material input channel to flank area", "false",
                            Patterns::Bool (),
                             "Put material input channel to flank area."
                             "Units: none.");	 
		  prm.declare_entry ("Input channel depth", "120e3",
                            Patterns::Double (),
                             "Input channel depth."
                             "Units: m");				
          //20241230
          prm.declare_entry	("Conductivities modification of mantle phase via T and P", "false",
                            Patterns::Bool (),
                             "Conductivities modification of mantle phase via T and P."
                             "Units: none.");
          //20250417
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/entropy-table/opxtable/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the `data/' subdirectory of ASPECT.");
          prm.declare_entry ("Normal material file name", "normal_material_table.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file name of the normal material data.");  
		  prm.declare_entry ("Carbon material file name", "carbon_material_table.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file name of the carbon material data.");  
		  //20250512
		  prm.declare_entry ("Wet material file name", "wet_material_table.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file name of the wet material data.");  
		  //20250428
		  prm.declare_entry ("Use perplex table", "false",
                            Patterns::Bool (),
                             "Use perplex table "
                             "Units: none.");		
		  //20260310 perplex添加玄武岩相
		  prm.declare_entry ("Perplex include mantle phase", "false",
                            Patterns::Bool (),
                             "Perplex include mantle phase"
                             "Units: none.");	
		  prm.declare_entry ("Perplex include basalt phase", "false",
                            Patterns::Bool (),
                             "Perplex include basalt phase"
                             "Units: none.");	
		  prm.declare_entry ("Basalt material file name", "basalt_material_table.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file name of the basalt material data.");  
		  //20250901
		  prm.declare_entry	("Density modification when temperature above Tp", "false",
                            Patterns::Bool (),
                             "Density modification when temperature above Tp."
                             "Units: none.");
		  prm.declare_entry ("Total force in one side", "1e12",
                             Patterns::Double (),
                             "Total force in one side "
                             "Units: \\Newton\\}.");
							 
		  //20260121 新增蛇纹岩化统计，仅发生在地幔相，不参与流变学计算
		  prm.declare_entry ("Use serpentinization statistics", "false",
                            Patterns::Bool (),
                             "Use serpentinization statistics "
                             "Units: none.");
		  //20260814 新增蛇纹岩相关：蛇纹岩化模式选择
		  prm.declare_entry ("Serp patterns", "liu2023",
                           Patterns::Selection("liu2023|christiansen2026"),
                           "Select what type of viscosity law to use between liu2023 "
                           "and christiansen2026 options.");
		  prm.declare_entry ("A for serpentinization", "1e-10",
                             Patterns::Double (),
                             "A for serpentinization "
                             "Units: \\none\\per\\second\\}.");
		  prm.declare_entry ("B for serpentinization", "2.5e-4",
                             Patterns::Double (),
                             "B for serpentinization "
                             "Units: \\none\\per\\°C\\square\\}.");
		  prm.declare_entry ("C for serpentinization", "270",
                             Patterns::Double (),
                             "C for serpentinization "
                             "Units: \\°C\\}.");
		  //20260814 新增蛇纹岩相关：Christiansen2026模式下氢气产率系数A，B，C
		  prm.declare_entry ("A for Christiansen2026", "2.2725e-5",
                             Patterns::Double (),
                             "A for Christiansen2026 "
                             "Units: \\none\\per\\°C\\square\\}.");
		  prm.declare_entry ("B for Christiansen2026", "-1.37e-2",
                             Patterns::Double (),
                             "B for Christiansen2026 "
                             "Units: \\none\\per\\°C\\}.");
		  prm.declare_entry ("C for Christiansen2026", "-12.5017",
                             Patterns::Double (),
                             "C for Christiansen2026 "
                             "Units: \\none\\}.");
		  prm.declare_entry ("Strain rate criteria for serpentinization", "1e-14",
                             Patterns::Double (),
                             "Strain rate criteria for serpentinization "
                             "Units: \\none\\per\\second\\}.");
		  prm.declare_entry ("Plastic strain criteria for serpentinization", "0.1",
                             Patterns::Double (),
                             "Strain rate criteria for serpentinization "
                             "Units: \\none\\per\\second\\}.");
		  //20260814 新增蛇纹岩相关：Christiansen2026模式下最小应变速率限制
		  //其他应变相关参数通用
		  prm.declare_entry ("Min strain rate criteria for Christiansen2026", "1e-16",
                             Patterns::Double (),
                             "Min strain rate criteria for Christiansen2026 "
                             "Units: \\none\\per\\second\\}.");
		  //20260320 新增参数，控制蛇纹岩化对内摩擦角的弱化程度
		  prm.declare_entry ("Plastic wekening factor for serpentinization", "0.1",
                             Patterns::Double (),
                             "Plastic wekening factor for serpentinization "
                             "Units: \\none\\}.");
		  //20260327 新增参数，控制蛇纹岩化导致的密度变化
		  prm.declare_entry ("Density decrease for serpentinization", "400",
                             Patterns::Double (),
                             "Density decrease for serpentinization "
                             "Units: \\kilogram\\per\\meter\\square\\}.");
		  //20260403 新增参数，控制是否开启塑性应变线性弱化
		  prm.declare_entry ("Use linear plastic strain healing", "false",
                            Patterns::Bool (),
                             "Use linear plastic strain healing "
                             "Units: none.");
		 //单位：1/s，默认参数：150Ma恢复5累计塑性应变，可根据需要调整
		  prm.declare_entry ("Plastic strain healing rate", "1.06e-15",
                             Patterns::Double (),
                             "Plastic strain healing rate "
                             "Units: \\per\\second\\}.");
		 //20260415 临时添加 按比例弱化，输入值等于0时不启用
		 prm.declare_entry ("Plastic strain healing rate (proportion)", "0.",
                             Patterns::Double (),
                             "Plastic strain healing rate "
                             "Units: none.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
	  
	  
    }



    template <int dim>
    void
    ViscoPlasticGrainSize<dim>::parse_parameters (ParameterHandler &prm)
    {
      // increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;
	  Tp = prm.get_double ("Adiabatic surface temperature");
	  prm.enter_subsection("Geometry model");
	  {
		  prm.enter_subsection("Box");
		  {
			  max_x = prm.get_double ("X extent");
			  repetitions_x = prm.get_integer ("X repetitions");
			  max_y = prm.get_double ("Y extent");
			  repetitions_y = prm.get_integer ("Y repetitions");
		  }
		  prm.leave_subsection();
	  }
	  prm.leave_subsection();
	  
	  prm.enter_subsection ("Initial composition model");
      {
        prm.enter_subsection("Lithosphere with rift");
        {
			layer_thickness = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Layer thicknesses"))),
																5,
																"Layer thicknesses");
			cratontic_thicknesses = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cratontic lithosphere layer thicknesses"))),
																5,
																"Cratontic lithosphere layer thicknesses");
			add_cratontic = prm.get_bool ("Add cratontic lithosphere");
			cratontic_boundary = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Segment of cratontic lithosphere boundary"))),
																2,
																"Segment of cratontic lithosphere boundary");
			exclude_mantle_lower = prm.get_bool ("Exclude phase 'mantle_lower' from mantle material");
		    center_cratontic = prm.get_bool ("Inside the boundaries is cratontic lithosphere");
			//20260114
			//新增盐岩模式，新模式初始温度计算四层，其他模式计算三层
			salt_mode = prm.get_bool ("Salt mode");
		}
		prm.leave_subsection();
	  }
	  prm.leave_subsection();
	  
	  prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Latent heat viscoplastic melt");
        {
			crustal_melt_temp_reduction = prm.get_double ("Change of entropy on crustal melting");
		}
		prm.leave_subsection();
	  }
	  prm.leave_subsection();
	  
	  //20240614
		prm.enter_subsection("Mesh refinement");
		{
			melt_fraction_refinement = prm.get_integer("Initial adaptive refinement") + prm.get_integer("Initial global refinement");
		}
		prm.leave_subsection ();
		
	  prm.enter_subsection("Postprocess");
	  {
		  prm.enter_subsection("Melt extraction");
		  {
			  depth_max = prm.get_double ("Max depth of melt convey function");
			  //20240614
			  //melt_fraction_refinement = prm.get_integer ("Mesh refinement of melt fraction");
			  
			  min_crood_y = std::ceil((max_y - depth_max - 5e3) / ((max_y / repetitions_y) / std::pow(2, melt_fraction_refinement)));
			  path_mantle_melting = this->get_output_directory() + "convey_flag_mantle_melting/";
			  path_crustal_melting = this->get_output_directory() + "convey_flag_crustal_melting/";
			  //20251104
			  melt_info_path = this->get_output_directory() + "melt_info/";
			  dike_mode = prm.get_bool ("Use dike mode");
			  constant_M = prm.get_bool ("Constant M");
			  M_for_dike = prm.get_double ("M for dike");
			  dike_width = prm.get_double ("Width of dike(m)");
			  dike_length = prm.get_double ("Length of dike(m)");
			  spreading_rate = prm.get_double ("Full spreading rate of model(m/yr)");
			  //20251201
			  constant_dike_length = prm.get_bool("Constant dike length");
		  }
		  prm.leave_subsection();
	  }
	  prm.leave_subsection();

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Visco Plastic Grain Size");
        {
          // Phase transition parameters
          phase_function.initialize_simulator (this->get_simulator());
          phase_function.parse_parameters (prm);

          std::vector<unsigned int> n_phase_transitions_for_each_composition
          (phase_function.n_phase_transitions_for_each_composition());

          // We require one more entry for density, etc as there are phase transitions
          // (for the low-pressure phase before any transition).
          for (unsigned int &n : n_phase_transitions_for_each_composition)
            n += 1;

          // Equation of state parameters
          equation_of_state.initialize_simulator (this->get_simulator());
          equation_of_state.parse_parameters (prm,
                                              std::make_shared<std::vector<unsigned int>>(n_phase_transitions_for_each_composition));


          thermal_diffusivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal diffusivities"))),
                                                                          n_fields,
                                                                          "Thermal diffusivities");

          define_conductivities = prm.get_bool ("Define thermal conductivities");

          thermal_conductivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivities"))),
                                                                           n_fields,
                                                                           "Thermal conductivities");
		  //20241227
		  densities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Densities"))),
                                                                           n_fields,
                                                                           "Densities");
		  
		  thermal_expansivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal expansivities"))),
                                                                           n_fields,
                                                                           "Thermal expansivities");
																		   
		  //k:Thermal conductivities
		  k_modify_order = prm.get_bool ("Thermal conductivities modification");
		  k_factor = prm.get_double ("Factor of thermal conductivities modification");
		  critial_edot_ii = prm.get_double ("Critial strain rate");
		  //defiened in rheology/visco_plastic_grainsize.cc
		  ref_strain_rate = prm.get_double("Reference strain rate");
		  min_strain_rate = prm.get_double("Minimum strain rate");
		  
		  //熔融新增
		  use_melt_model  = prm.get_bool ("Use melt model");
		  include_mantle_melting = prm.get_bool ("Include mantle melting");
		  A1              = prm.get_double ("A1");
          A2              = prm.get_double ("A2");
          A3              = prm.get_double ("A3");
          B1              = prm.get_double ("B1");
          B2              = prm.get_double ("B2");
          B3              = prm.get_double ("B3");
          C1              = prm.get_double ("C1");
          C2              = prm.get_double ("C2");
          C3              = prm.get_double ("C3");
          r1              = prm.get_double ("r1");
          r2              = prm.get_double ("r2");
          beta            = prm.get_double ("beta");
		  M_cpx           = prm.get_double ("Mass fraction cpx");
		  
		  //异质地幔对固相线的影响
		  hetero_temperature = prm.get_double ("Solidus change of heterogeneous mantle");
		  //变量"Depth of banning grain calculation"在rheology/diffusion_creep_grainsize.cc中声明
		  banned_depth_meter = prm.get_double("Depth of banning grain calculation");
		  
		  init_water_content = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Initial water content"))),
                                                                            n_fields,
                                                                            "Initial water content");
		  //含水熔融（地幔）
		  include_hydrous_melting = prm.get_bool("Include hydrous melting");
		  D_H2O           = prm.get_double ("Water distribution coefficient DH2O");
		  K = prm.get_double ("K");
		  gamma = prm.get_double ("gamma");
		  
		  //含碳熔融（地幔）
		  include_carbonate_melting = prm.get_bool("Include carbonate melting");
		  X_CO2           = prm.get_double ("CO2 content in the melt");
		  
		  
		  //地壳熔融
		  include_crustal_melting = prm.get_bool("Include crustal melting");
		  A1_crust                = prm.get_double ("A1 for crustal melting");
		  A2_crust                = prm.get_double ("A2 for crustal melting");
		  A3_crust                = prm.get_double ("A3 for crustal melting");
		  B1_crust                = prm.get_double ("B1 for crustal melting");
		  B2_crust                = prm.get_double ("B2 for crustal melting");
		  B3_crust                = prm.get_double ("B3 for crustal melting");
		  gamma_crust             = prm.get_double ("Crustal melting expansion factor");
		  
	      //熔融搬运相关参数																
		  melt_convey = prm.get_bool("Use melt convey function");
		  
          rheology = std_cxx14::make_unique<Rheology::ViscoPlasticGrainSize<dim>>();
          rheology->initialize_simulator (this->get_simulator());
		  //20251201 新增 constant_dike_length
          rheology->parse_parameters(prm, 
		                             max_x,
									 max_y,
									 add_cratontic,
									 center_cratontic,
									 exclude_mantle_lower,
									 dike_mode,
									 constant_M,
									 constant_dike_length,
									 M_for_dike,
									 dike_width,
									 dike_length,
									 spreading_rate,
									 layer_thickness,
									 cratontic_thicknesses,
									 cratontic_boundary,
		                             std::make_shared<std::vector<unsigned int>>(n_phase_transitions_for_each_composition));
		  
		  //20240617
		  input_from_flank = prm.get_bool ("Put material input channel to flank area");
		  input_depth = prm.get_double ("Input channel depth");
		  init_grain_size = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Initial Grain size"))),
                                                                            n_fields,
                                                                            "Initial Grain size");
		  //20241230
		  use_k_modified_by_T_P = prm.get_bool ("Conductivities modification of mantle phase via T and P");
		  
		  //20250417
		  data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
		  normal_material_file_name = prm.get ("Normal material file name");
		  carbon_material_file_name = prm.get ("Carbon material file name");
		  use_perplex_table = prm.get_bool ("Use perplex table");
		  //20250512
		  wet_material_file_name = prm.get ("Wet material file name");
		  //20260310 perplex添加玄武岩相
		  perplex_include_mantle_phase = prm.get_bool ("Perplex include mantle phase");
		  perplex_include_basalt_phase = prm.get_bool ("Perplex include basalt phase");
		  basalt_material_file_name = prm.get ("Basalt material file name");
		  //20250901
		  density_modify_above_tp = prm.get_bool ("Density modification when temperature above Tp");
		  
		  //20251201 弃用，暂时保留
		  input_total_force = prm.get_double ("Total force in one side");
		  
		  //20260121 新增蛇纹岩化统计，仅发生在地幔相，不参与流变学计算
		  include_serp = prm.get_bool ("Use serpentinization statistics");
		  serp_A = prm.get_double ("A for serpentinization");
		  serp_B = prm.get_double ("B for serpentinization");
		  serp_C = prm.get_double ("C for serpentinization");
		  serp_strain_rate_crit = prm.get_double ("Strain rate criteria for serpentinization");
		  serp_plastic_strain_crit = prm.get_double ("Plastic strain criteria for serpentinization");
		  //20260327 新增参数，控制蛇纹岩化导致的密度变化
		  serp_density_change = prm.get_double ("Density decrease for serpentinization");
		  //20260403
		  plastic_linear_healing = prm.get_bool ("Use linear plastic strain healing");
		  plastic_healing_rate = prm.get_double ("Plastic strain healing rate");
		  //20260415
		  plastic_healing_proportion = prm.get_double ("Plastic strain healing rate (proportion)");
		  
		  //20260814 新增蛇纹岩化新模式相关
		  if(prm.get ("Serp patterns") == "liu2023")
			  serp_pattern = 1;
		  else if(prm.get ("Serp patterns") == "christiansen2026")
			  serp_pattern = 2;
		  else
			  AssertThrow(false, ExcMessage("Not a valid Serp patterns"));
		  serp_A_christiansen = prm.get_double ("A for Christiansen2026");
		  serp_B_christiansen = prm.get_double ("B for Christiansen2026");
		  serp_C_christiansen = prm.get_double ("C for Christiansen2026");
		  min_strain_rate_christiansen = prm.get_double ("Min strain rate criteria for Christiansen2026");
		  AssertThrow(min_strain_rate_christiansen < serp_strain_rate_crit, ExcMessage("Serp parameter error, check 'crit strain rate' and 'min strain rtate for christiansen'."));

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::strain_rate | NonlinearDependence::compositional_fields;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
    }

    template <int dim>
    void
    ViscoPlasticGrainSize<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      rheology->create_plastic_outputs(out);

      if (rheology->use_elasticity)
        rheology->elastic_rheology.create_elastic_outputs(out);
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ViscoPlasticGrainSize,
                                   "visco plastic grain size",
                                   "An implementation of an incompressible visco(elastic)-plastic rheology "
                                   "with options for selecting dislocation creep, diffusion creep or "
                                   "composite viscous flow laws. Prior to yielding, one may select to "
                                   "modify the viscosity to account for viscoelastic effects. Plasticity "
                                   "limits viscous stresses through a Drucker Prager yield criterion. "
                                   "The implementation of this material model is based heavily on the "
                                   "`DiffusionDislocation' (Bob Myhill), `DruckerPrager' "
                                   "(Anne Glerum), and `Viscoelastic' (John Naliboff) material models. "
                                   "\n\n "
                                   "The viscosity for dislocation or diffusion creep is defined as "
                                   "\\[v = \\frac 12 A^{-\\frac{1}{n}} d^{\\frac{m}{n}} "
                                   "\\dot{\\varepsilon}_{ii}^{\\frac{1-n}{n}} "
                                   "\\exp\\left(\\frac{E + PV}{nRT}\\right)\\] "
                                   "where $A$ is the prefactor, $n$ is the stress exponent, "
                                   "$\\dot{\\varepsilon}_{ii}$ is the square root of the deviatoric "
                                   "strain rate tensor second invariant, $d$ is grain size, "
                                   "$m$ is the grain size exponent, $E$ is activation energy, "
                                   "$V$ is activation volume, $P$ is pressure, $R$ is the gas "
                                   "exponent and $T$ is temperature. "
                                   "This form of the viscosity equation is commonly used in "
                                   "geodynamic simulations. See, for example, Billen and Hirth "
                                   "(2007), G3, 8, Q08012. Significantly, other studies may use "
                                   "slightly different forms of the viscosity equation leading to "
                                   "variations in how specific terms are defined or combined. For "
                                   "example, the grain size exponent should always be positive in "
                                   "the diffusion viscosity equation used here, while other studies "
                                   "place the grain size term in the denominator and invert the sign "
                                   "of the grain size exponent. When examining previous work, one "
                                   "should carefully check how the viscous prefactor and grain size "
                                   "terms are defined. "
                                   "\n\n "
                                   "One may select to use the diffusion ($v_{\\text{diff}}$; $n=1$, $m\neq 0$), "
                                   "dislocation ($v_{\\text{disl}}$, $n>1$, $m=0$) or composite "
                                   "$\\frac{v_{\\text{diff}} v_{\\text{disl}}}{v_{\\text{diff}}+v_{\\text{disl}}}$ equation form. "
                                   "\n\n "
                                   "The diffusion and dislocation prefactors can be weakened with a factor "
                                   "between 0 and 1 according to the total or the viscous strain only. "
                                   "\n\n "
                                   "Viscosity is limited through one of two different `yielding' mechanisms. "
                                   "\n\n"
                                   "The first plasticity mechanism limits viscous stress through a "
                                   "Drucker Prager yield criterion, where the yield stress in 3D is  "
                                   "$\\sigma_y = \\frac{6C\\cos(\\phi) + 2P\\sin(\\phi)} "
                                   "{\\sqrt{3}(3+\\sin(\\phi))}$ "
                                   "and "
                                   "$\\sigma_y = C\\cos(\\phi) + P\\sin(\\phi)$ "
                                   "in 2D. Above, $C$ is cohesion and $\\phi$  is the angle of "
                                   "internal friction.  Note that the 2D form is equivalent to the "
                                   "Mohr Coulomb yield surface.  If $\\phi$ is 0, the yield stress "
                                   "is fixed and equal to the cohesion (Von Mises yield criterion). "
                                   "When the viscous stress ($2v{\\varepsilon}_{ii}$) exceeds "
                                   "the yield stress, the viscosity is rescaled back to the yield "
                                   "surface: $v_{y}=\\sigma_{y}/(2{\\varepsilon}_{ii})$. "
                                   "This form of plasticity is commonly used in geodynamic models. "
                                   "See, for example, Thieulot, C. (2011), PEPI 188, pp. 47-68. "
                                   "\n\n"
                                   "The user has the option to linearly reduce the cohesion and "
                                   "internal friction angle as a function of the finite strain magnitude. "
                                   "The finite strain invariant or full strain tensor is calculated through "
                                   "compositional fields within the material model. This implementation is "
                                   "identical to the compositional field finite strain plugin and cookbook "
                                   "described in the manual (author: Gassmoeller, Dannberg). If the user selects to track "
                                   "the finite strain invariant ($e_{ii}$), a single compositional field tracks "
                                   "the value derived from $e_{ii}^t = (e_{ii})^{(t-1)} + \\dot{e}_{ii}\\; dt$, where $t$ and $t-1$ "
                                   "are the current and prior time steps, $\\dot{e}_{ii}$ is the second invariant of the "
                                   "strain rate tensor and $dt$ is the time step size. In the case of the "
                                   "full strain tensor $F$, the finite strain magnitude is derived from the "
                                   "second invariant of the symmetric stretching tensor $L$, where "
                                   "$L = F [F]^T$. The user must specify a single compositional "
                                   "field for the finite strain invariant or multiple fields (4 in 2D, 9 in 3D) "
                                   "for the finite strain tensor. These field(s) must be the first listed "
                                   "compositional fields in the parameter file. Note that one or more of the finite strain "
                                   "tensor components must be assigned a non-zero value initially. This value can be "
                                   "be quite small (e.g., 1.e-8), but still non-zero. While the option to track and use "
                                   "the full finite strain tensor exists, tracking the associated compositional fields "
                                   "is computationally expensive in 3D. Similarly, the finite strain magnitudes "
                                   "may in fact decrease if the orientation of the deformation field switches "
                                   "through time. Consequently, the ideal solution is track the finite strain "
                                   "invariant (single compositional) field within the material and track "
                                   "the full finite strain tensor through particles."
                                   "When only the second invariant of the strain is tracked, one has the option to "
                                   "track the full strain or only the plastic strain. In the latter case, strain is only tracked "
                                   "in case the material is plastically yielding, i.e. the viscous stress > yield stress. "
                                   "\n\n"
                                   "Viscous stress may also be limited by a non-linear stress limiter "
                                   "that has a form similar to the Peierls creep mechanism. "
                                   "This stress limiter assigns an effective viscosity "
                                   "$\\sigma_{\\text{eff}} = \\frac{\\tau_y}{2\\varepsilon_y} "
                                   "{\\frac{\\varepsilon_{ii}}{\\varepsilon_y}}^{\\frac{1}{n_y}-1}$ "
                                   "Above $\\tau_y$ is a yield stress, $\\varepsilon_y$ is the "
                                   "reference strain rate, $\\varepsilon_{ii}$ is the strain rate "
                                   "and $n_y$ is the stress limiter exponent.  The yield stress, "
                                   "$\\tau_y$, is defined through the Drucker Prager yield criterion "
                                   "formulation. This method of limiting viscous stress has been used "
                                   "in various forms within the geodynamic literature \\cite{chri92,vavv02,cibi13,cibi15}."
                                   "When $n_y$ is 1, it essentially becomes a linear viscosity model, "
                                   "and in the limit $n_y\\rightarrow \\infty$ it converges to the "
                                   "standard viscosity rescaling method (concretely, values $n_y>20$ "
                                   "are large enough)."
                                   "\n\n "
                                   "The visco-plastic rheology described above may also be modified to include "
                                   "viscoelastic deformation, thus producing a viscoelastic plastic constitutive "
                                   "relationship. "
                                   "\n\n "
                                   "The viscoelastic rheology behavior takes into account the elastic shear "
                                   "strength (e.g., shear modulus), while the tensile and volumetric "
                                   "strength (e.g., Young's and bulk modulus) are not considered. The "
                                   "model is incompressible and allows specifying an arbitrary number "
                                   "of compositional fields, where each field represents a different "
                                   "rock type or component of the viscoelastic stress tensor. The stress "
                                   "tensor in 2D and 3D, respectively, contains 3 or 6 components. The "
                                   "compositional fields representing these components must be named "
                                   "and listed in a very specific format, which is designed to minimize "
                                   "mislabeling stress tensor components as distinct 'compositional "
                                   "rock types' (or vice versa). For 2D models, the first three "
                                   "compositional fields must be labeled 'stress\\_xx', 'stress\\_yy' and 'stress\\_xy'. "
                                   "In 3D, the first six compositional fields must be labeled 'stress\\_xx', "
                                   "'stress\\_yy', 'stress\\_zz', 'stress\\_xy', 'stress\\_xz', 'stress\\_yz'. "
                                   "\n\n "
                                   "Combining this viscoelasticity implementation with non-linear viscous flow "
                                   "and plasticity produces a constitutive relationship commonly referred to "
                                   "as partial elastoviscoplastic (e.g., pEVP) in the geodynamics community. "
                                   "While extensively discussed and applied within the geodynamics "
                                   "literature, notable references include: "
                                   "Moresi et al. (2003), J. Comp. Phys., v. 184, p. 476-497. "
                                   "Gerya and Yuen (2007), Phys. Earth. Planet. Inter., v. 163, p. 83-105. "
                                   "Gerya (2010), Introduction to Numerical Geodynamic Modeling. "
                                   "Kaus (2010), Tectonophysics, v. 484, p. 36-47. "
                                   "Choi et al. (2013), J. Geophys. Res., v. 118, p. 2429-2444. "
                                   "Keller et al. (2013), Geophys. J. Int., v. 195, p. 1406-1442. "
                                   "\n\n "
                                   "The overview below directly follows Moresi et al. (2003) eqns. 23-38. "
                                   "However, an important distinction between this material model and "
                                   "the studies above is the use of compositional fields, rather than "
                                   "particles, to track individual components of the viscoelastic stress "
                                   "tensor. The material model will be updated when an option to track "
                                   "and calculate viscoelastic stresses with particles is implemented. "
                                   "\n\n "
                                   "Moresi et al. (2003) begins (eqn. 23) by writing the deviatoric "
                                   "rate of deformation ($\\hat{D}$) as the sum of elastic "
                                   "($\\hat{D_{e}}$) and viscous ($\\hat{D_{v}}$) components: "
                                   "$\\hat{D} = \\hat{D_{e}} + \\hat{D_{v}}$.  "
                                   "These terms further decompose into "
                                   "$\\hat{D_{v}} = \\frac{\\tau}{2\\eta}$ and "
                                   "$\\hat{D_{e}} = \\frac{\\overset{\\nabla}{\\tau}}{2\\mu}$, where "
                                   "$\\tau$ is the viscous deviatoric stress, $\\eta$ is the shear viscosity, "
                                   "$\\mu$ is the shear modulus and $\\overset{\\nabla}{\\tau}$ is the "
                                   "Jaumann corotational stress rate. This later term (eqn. 24) contains the "
                                   "time derivative of the deviatoric stress ($\\dot{\\tau}$) and terms that "
                                   "account for material spin (e.g., rotation) due to advection: "
                                   "$\\overset{\\nabla}{\\tau} = \\dot{\\tau} + {\\tau}W -W\\tau$. "
                                   "Above, $W$ is the material spin tensor (eqn. 25): "
                                   "$W_{ij} = \\frac{1}{2} \\left (\\frac{\\partial V_{i}}{\\partial x_{j}} - "
                                   "\\frac{\\partial V_{j}}{\\partial x_{i}} \\right )$. "
                                   "\n\n "
                                   "If plasticity is included, the deviatoric rate of deformation may be written as: "
                                   "$\\hat{D} = \\hat{D_{e}} + \\hat{D_{v}} + \\hat{D_{p}}$, where $\\hat{D_{p}}$ "
                                   "is the plastic component. $\\hat{D_{p}}$ decomposes to $\\frac{\\tau_{y}}{2\\eta_{y}}$, "
                                   "where $\\tau_{y}$ is the yield stress and $\\eta_{y}$ is the viscosity rescaled "
                                   "to the yield surface. "
                                   "The Jaumann stress-rate can also be approximated using terms from the "
                                   "previous time step ($t$) and current time step ($t + \\Delta t^{e}$): "
                                   "$\\smash[t]{\\overset{\\nabla}{\\tau}}^{t + \\Delta t^{e}} \\approx "
                                   "\\frac{\\tau^{t + \\Delta t^{e} - \\tau^{t}}}{\\Delta t^{e}} - "
                                   "W^{t}\\tau^{t} + \\tau^{t}W^{t}$. "
                                   "In this material model, the size of the time step above ($\\Delta t^{e}$) "
                                   "can be specified as the numerical time step size or an independent fixed time "
                                   "step. If the latter case is selected, the user has an option to apply a "
                                   "stress averaging scheme to account for the differences between the numerical "
                                   "and fixed elastic time step (eqn. 32). If one selects to use a fixed elastic time "
                                   "step throughout the model run, this can still be achieved by using CFL and "
                                   "maximum time step values that restrict the numerical time step to a specific time."
                                   "\n\n "
                                   "The formulation above allows rewriting the total rate of deformation (eqn. 29) as\n "
                                   "$\\tau^{t + \\Delta t^{e}} = \\eta_{eff} \\left ( "
                                   "2\\hat{D}^{t + \\triangle t^{e}} + \\frac{\\tau^{t}}{\\mu \\Delta t^{e}} + "
                                   "\\frac{W^{t}\\tau^{t} - \\tau^{t}W^{t}}{\\mu}  \\right )$. "
                                   "\n\n "
                                   "The effective viscosity (eqn. 28) is a function of the viscosity ($\\eta$), "
                                   "elastic time step size ($\\Delta t^{e}$) and shear relaxation time "
                                   "($ \\alpha = \\frac{\\eta}{\\mu} $): "
                                   "$\\eta_{eff} = \\eta \\frac{\\Delta t^{e}}{\\Delta t^{e} + \\alpha}$ "
                                   "The magnitude of the shear modulus thus controls how much the effective "
                                   "viscosity is reduced relative to the initial viscosity. "
                                   "\n\n "
                                   "Elastic effects are introduced into the governing Stokes equations through "
                                   "an elastic force term (eqn. 30) using stresses from the previous time step: "
                                   "$F^{e,t} = -\\frac{\\eta_{eff}}{\\mu \\Delta t^{e}} \\tau^{t}$. "
                                   "This force term is added onto the right-hand side force vector in the "
                                   "system of equations. "
                                   "\n\n "
                                   "When plastic yielding occurs, the effective viscosity in equation 29 and 30 is the "
                                   "plastic viscosity (equation 36). If the current stress is below the plastic "
                                   "yield stress, the effective viscosity is still as defined in equation 28. "
                                   "During non-linear iterations, we define the current stress prior to yielding "
                                   "(e.g., value compared to yield stress) as "
                                   "$\\tau^{t + \\Delta t^{e}} = \\eta_{eff} \\left ( 2\\hat{D}^{t + \\triangle t^{e}} + "
                                   "\\frac{\\tau^{t}}{\\mu \\Delta t^{e}} \\right ) $"
                                   "\n\n "
                                   "Compositional fields can each be assigned individual values of "
                                   "thermal diffusivity, heat capacity, density, thermal "
                                   "expansivity and rheological parameters. "
                                   "\n\n "
                                   "If more than one compositional field is present at a given "
                                   "point, viscosities are averaged with an arithmetic, geometric "
                                   "harmonic (default) or maximum composition scheme. "
                                   "\n\n "
                                   "The value for the components of this formula and additional "
                                   "parameters are read from the parameter file in subsection "
                                   " 'Material model/Visco Plastic'.")
  }
}

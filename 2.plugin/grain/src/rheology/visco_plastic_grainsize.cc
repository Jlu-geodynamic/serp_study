/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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

#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/rheology/visco_plastic_grainsize.h>
#include <aspect/material_model/utilities.h>
#include <aspect/utilities.h>
#include <aspect/newton.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/fe/fe_values.h>

// 此模型以material model中的visco plastic模型为基础修改而来，grain size相关计算参考aterial model中的grain size model模型
// 原文件位置以及名称：material_model/visco_plastic.cc 
//                     material_model/rheology/visco_plastic.cc 
//					   material_model/rheology/diffusion_creep.cc 
//					   material_model/grain_size.cc
namespace aspect
{
  namespace MaterialModel
  {
    namespace
    {
      std::vector<std::string> make_plastic_additional_outputs_names_grain_size()
      {
        std::vector<std::string> names;
        names.emplace_back("current_cohesions");
        names.emplace_back("current_friction_angles");
        names.emplace_back("plastic_yielding");
        names.emplace_back("(1/df)/(1/df+1/ds)");
		//Asthenosphere
        names.emplace_back("volume_0");
		//viscosity in dislcation creep
		names.emplace_back("eta_dislocation");
		//viscosity in diffusion creep
		names.emplace_back("eta_diffusion");
		//effective viscosity in viscous
		names.emplace_back("eta_for_viscous");
		//plastic stress
		names.emplace_back("plastic_stress");
		//viscous stress
		names.emplace_back("viscous_stress");
		//yield_strength
		names.emplace_back("yield_strength");
        return names;
      }
    }

    template <int dim>
    PlasticAdditionalOutputsGrainSize<dim>::PlasticAdditionalOutputsGrainSize(const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_plastic_additional_outputs_names_grain_size()),
      cohesions(n_points, numbers::signaling_nan<double>()),
      friction_angles(n_points, numbers::signaling_nan<double>()),
      yielding(n_points, numbers::signaling_nan<double>()),
      creep_visco_ratio(n_points, numbers::signaling_nan<double>()),
      volume_0(n_points, numbers::signaling_nan<double>()),
	  eta_dislocation(n_points, numbers::signaling_nan<double>()),
	  eta_diffusion(n_points, numbers::signaling_nan<double>()),
	  eta_for_viscous(n_points, numbers::signaling_nan<double>()),
	  plastic_stress(n_points, numbers::signaling_nan<double>()),
	  viscous_stress(n_points, numbers::signaling_nan<double>()),
	  yield_strength(n_points, numbers::signaling_nan<double>())
    {}



    template <int dim>
    std::vector<double>
     PlasticAdditionalOutputsGrainSize<dim>::get_nth_output(const unsigned int idx) const
    {
      AssertIndexRange (idx, 10);
      switch (idx)
        {
          case 0:
            return cohesions;
          case 1:
            return friction_angles;
          case 2:
            return yielding;
          case 3:
              return creep_visco_ratio;
          case 4:
              return volume_0;
		  case 5:
		      return eta_dislocation;
		  case 6:
		      return eta_diffusion;
		  case 7:
		      return eta_for_viscous;
		  case 8:
		      return plastic_stress;
		  case 9:
		      return viscous_stress;
		  case 10:
		      return yield_strength;
          default:
            AssertThrow(false, ExcInternalError());
        }
      // We will never get here, so just return something
      return cohesions;
    }

    namespace Rheology
    {
      template <int dim>
      ViscoPlasticGrainSize<dim>::ViscoPlasticGrainSize ()
      {} 
	  
	  template <int dim>
	  double
	  ViscoPlasticGrainSize<dim>::
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
				  sleep(0.01);
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
	  ViscoPlasticGrainSize<dim>::
	  dike_depend(const double posx,
	              const double posy,
	              const unsigned int mpi_communicator)const
	  {
		  //确定岩墙位置
		  //20251201 新增岩墙长度不固定的判定情况
		  const double emplace_xmax = read_from_txt(melt_info_path + "dike_emplace_location_xmax.txt", mpi_communicator);
		  const double emplace_ymin = read_from_txt(melt_info_path + "dike_emplace_location_ymin.txt", mpi_communicator);
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
	  
	  //20251110
	  //岩墙相关的累积应变调整
	  //暂时未考虑healing的情况
	  template <int dim>
      std::array<double, 3>
      ViscoPlasticGrainSize<dim>::
      dikemode_compute_strain_weakening_factors(const unsigned int j,
                                                const std::vector<double> &composition) const
      {
        double viscous_weakening = 1.0;
        std::pair<double, double> brittle_weakening (1.0, 1.0);

        switch (dike_weakening_mechanism)
          {
            case dike_none:
            {
              break;
            }
            case dike_plastic_weakening_with_plastic_strain_only:
            {
              const unsigned int plastic_strain_index = this->introspection().compositional_index_for_name("plastic_strain");
              brittle_weakening = strain_rheology.calculate_plastic_weakening(composition[plastic_strain_index], j);
              break;
            }
            case dike_plastic_weakening_with_plastic_strain_and_viscous_weakening_with_viscous_strain:
            {
              const unsigned int plastic_strain_index = this->introspection().compositional_index_for_name("plastic_strain");
              brittle_weakening = strain_rheology.calculate_plastic_weakening(composition[plastic_strain_index], j);
              const unsigned int viscous_strain_index = this->introspection().compositional_index_for_name("viscous_strain");
              viscous_weakening = strain_rheology.calculate_viscous_weakening(composition[viscous_strain_index], j);
              break;
            }
            case dike_plastic_weakening_with_plastic_strain_and_viscous_weakening_with_grain_size:
            {
              const unsigned int plastic_strain_index = this->introspection().compositional_index_for_name("plastic_strain");
              brittle_weakening = strain_rheology.calculate_plastic_weakening(composition[plastic_strain_index], j);
              const unsigned int viscous_strain_index = this->introspection().compositional_index_for_name("viscous_strain_grainsize");
              viscous_weakening = strain_rheology.calculate_viscous_weakening(composition[viscous_strain_index], j);
              break;
            }
            default:
            {
              AssertThrow(false, ExcNotImplemented());
              break;
            }
          }

        const std::array<double, 3> weakening_factors = {{brittle_weakening.first,brittle_weakening.second,viscous_weakening}};

        return weakening_factors;

      }
	  
	  //20251110
	  //岩墙相关的累积应变调整
	  //暂时未考虑healing的情况
	  //20251124
	  //将岩墙内累积应变置0
	  template <int dim>
      void
      ViscoPlasticGrainSize<dim>::
      dikemode_fill_reaction_outputs(const MaterialModel::MaterialModelInputs<dim> &in,
                                     const int i,
                                     const double min_strain_rate,
                                     const bool plastic_yielding,
									 const bool in_dike,
									 const SymmetricTensor<2,dim> strain_rate_dev,
                                     MaterialModel::MaterialModelOutputs<dim> &out) const
	  {
		  if  (this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && in.requests_property(MaterialProperties::reaction_terms))
		  {
			  const double edot_ii = std::max(sqrt(std::fabs(second_invariant(deviator(strain_rate_dev)))),min_strain_rate);
			  const double dike_delta_e_ii = edot_ii * this->get_timestep();
			  if(dike_weakening_mechanism == dike_plastic_weakening_with_plastic_strain_only && plastic_yielding)
				  out.reaction_terms[i][this->introspection().compositional_index_for_name("plastic_strain")] =
					in_dike ?  0. - in.composition[i][this->introspection().compositional_index_for_name("plastic_strain")] : std::max(dike_delta_e_ii, -in.composition[i][this->introspection().compositional_index_for_name("plastic_strain")]);
			  if(dike_weakening_mechanism == dike_plastic_weakening_with_plastic_strain_and_viscous_weakening_with_viscous_strain)
			  {
				  if(plastic_yielding)
					  out.reaction_terms[i][this->introspection().compositional_index_for_name("plastic_strain")] =
						in_dike ?  0. - in.composition[i][this->introspection().compositional_index_for_name("plastic_strain")] : std::max(dike_delta_e_ii, -in.composition[i][this->introspection().compositional_index_for_name("plastic_strain")]);
				  else
					  out.reaction_terms[i][this->introspection().compositional_index_for_name("viscous_strain")] =
						in_dike ?  0. - in.composition[i][this->introspection().compositional_index_for_name("viscous_strain")] : std::max(dike_delta_e_ii, -in.composition[i][this->introspection().compositional_index_for_name("viscous_strain")]);
			  }
			  if(dike_weakening_mechanism == dike_plastic_weakening_with_plastic_strain_and_viscous_weakening_with_grain_size)
			  {
				  if(plastic_yielding)
					  out.reaction_terms[i][this->introspection().compositional_index_for_name("plastic_strain")] =
						in_dike ?  0. - in.composition[i][this->introspection().compositional_index_for_name("plastic_strain")] : std::max(dike_delta_e_ii, -in.composition[i][this->introspection().compositional_index_for_name("plastic_strain")]);
				  else
					  out.reaction_terms[i][this->introspection().compositional_index_for_name("viscous_strain_grainsize")] =
						in_dike ?  0. - in.composition[i][this->introspection().compositional_index_for_name("viscous_strain_grainsize")] : std::max(dike_delta_e_ii, -in.composition[i][this->introspection().compositional_index_for_name("viscous_strain_grainsize")]);
			  }
			  //'noninitial_plastic_strain'记录未修正岩墙内累积塑性应变的情况
			  //如模型添加了初始噪声，则此项不包含初始噪声
			  if (this->introspection().compositional_name_exists("noninitial_plastic_strain") && plastic_yielding == true)
				  out.reaction_terms[i][this->introspection().compositional_index_for_name("noninitial_plastic_strain")] =
					std::max(dike_delta_e_ii, -in.composition[i][this->introspection().compositional_index_for_name("noninitial_plastic_strain")]);
		  }
		  
	  }
	  
	  template <int dim>
	  std::vector<double>
	  ViscoPlasticGrainSize<dim>::
	  fractions_modify(const std::vector<double> &temp_volume_fractions) const
	  {
		  const std::vector<double>temp_volume = temp_volume_fractions;
		  //消除成分比例过小的部分
		  std::vector<double>out_volume(temp_volume_fractions.size(), 0);
          double temp_fractions = 0;
          double fractions_store = 0;
          double full_fractions = 0;
          unsigned int temp_max_index = 0;
          for (unsigned int k = 0; k < temp_volume_fractions.size(); ++k)
          {
              //1.填充fill_volume,记录各成分比例之和
              out_volume[k] = temp_volume[k];
              full_fractions += temp_volume[k];

              //2.将 < fraction_ratio的部分消除，并累加被消去的比重
              if (out_volume[k] < fraction_ratio)
              {
                  out_volume[k] = 0;
                  fractions_store += temp_volume[k];
              }

              //3.记录成分比例最大的部分
              //当比例相等时，认为序号较小的为最大值（background > sediment1 > sediment2 > upper > lower > mantle）
              if (temp_fractions < temp_volume[k])
              {
                  temp_fractions = temp_volume[k];
                  temp_max_index = k;
              }
          }
          //循环后，若比例之和 > 1，在fractions_store中消除 > 1的部分
          if (full_fractions > 1)
              fractions_store += 1 - full_fractions;
          //将消去的成分比例加到占比最大的区域（控制最大值为1）
          const unsigned int max_index = temp_max_index;
          out_volume[max_index] = std::min((out_volume[temp_max_index] + fractions_store) , 1e-0);
		  return out_volume;
	  }
	  
      template <int dim>
      IsostrainViscosities
      ViscoPlasticGrainSize<dim>::
      calculate_isostrain_viscosities (const MaterialModel::MaterialModelInputs<dim> &in,
                                       const unsigned int i,
                                       const std::vector<double> &volume_fractions,
                                       const std::vector<double> &phase_function_values,
                                       const std::vector<unsigned int> &n_phases_per_composition) const
      {
	    const double posx = in.position[i][0];
        const double posy = in.position[i][1];
		//20260127 蛇纹岩化相关：根据蛇纹岩化程度修改内摩擦角，蛇纹岩化程度大于0.1时，内摩擦角调整为输入值的0.1倍
		const double new_crust_phase = this->introspection().compositional_name_exists("new_crust") ? volume_fractions[this->introspection().compositional_index_for_name("new_crust") + 1] : 0.;
		const double serp_value = this->introspection().compositional_name_exists("serp") ? in.composition[i][this->introspection().compositional_index_for_name("serp")] : 0.;
		//判断模型初始化是否完成
        const bool initialization_past = this->simulator_is_past_initialization() && this->get_timestep_number() > 0;
		const double melt_fraction = this->introspection().compositional_name_exists("melt_fraction") ? std::abs(in.composition[i][this->introspection().compositional_index_for_name("melt_fraction")]) : 0.;
		
        IsostrainViscosities output_parameters;
        // Initialize or fill variables used to calculate viscosities
        output_parameters.composition_yielding.resize(volume_fractions.size(), false);
        output_parameters.composition_viscosities.resize(volume_fractions.size(), numbers::signaling_nan<double>());
        // Assemble stress tensor if elastic behavior is enabled
        SymmetricTensor<2,dim> stress_old = numbers::signaling_nan<SymmetricTensor<2,dim>>();
        if (use_elasticity == true)
          {
            for (unsigned int j=0; j < SymmetricTensor<2,dim>::n_independent_components; ++j)
				stress_old[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)] = in.composition[i][j];
          }
        // The first time this function is called (first iteration of first time step)
        // a specified "reference" strain rate is used as the returned value would
        // otherwise be zero.
        const bool use_reference_strainrate = (this->get_timestep_number() == 0) &&
                                              (in.strain_rate[i].norm() <= std::numeric_limits<double>::min());
        double edot_ii;
        if (use_reference_strainrate)
          edot_ii = ref_strain_rate;
        else
         // Calculate the square root of the second moment invariant for the deviatoric strain rate tensor.
          edot_ii = std::max(std::sqrt(std::fabs(second_invariant(in.strain_rate[i]))),
                             min_strain_rate);
        // Calculate viscosities for each of the individual compositional phases
        for (unsigned int j=0; j < volume_fractions.size(); ++j)
          {
            // Step 1: viscous behavior
            // Choice of activation volume depends on whether there is an adiabatic temperature
            // gradient used when calculating the viscosity. This allows the same activation volume
            // to be used in incompressible and compressible models.
            const double temperature_for_viscosity = in.temperature[i] + adiabatic_temperature_gradient_for_viscosity*in.pressure[i];
            AssertThrow(temperature_for_viscosity != 0, ExcMessage(
                          "The temperature used in the calculation of the visco-plastic rheology is zero. "
                          "This is not allowed, because this value is used to divide through. It is probably "
                          "being caused by the temperature being zero somewhere in the model. The relevant "
                          "values for debugging are: temperature (" + Utilities::to_string(in.temperature[i]) +
                          "), adiabatic_temperature_gradient_for_viscosity ("
                          + Utilities::to_string(adiabatic_temperature_gradient_for_viscosity) + ") and pressure ("
                          + Utilities::to_string(in.pressure[i]) + ")."));
            // Step 1a: compute viscosity from diffusion creep law
            //模型选择在最开始的阶段进行粒径变化计算，其中使用到的粘度是通过dislocation creep以及diffusion creep计算得到的。
			//使用的应变速率在输入中直接取得，不涉及后续的优化措施。
			//含水量最多降低到30 ppm
			//20260203 可删除rift_grain_size
			//20260225 蛇纹岩化区域参与塑性应变弱化
			double yield_stress_pre_grainsize_calculation = 0.;
			if(this->introspection().compositional_name_exists("rift_grain_size"))
			{
				const std::array<double, 3> weakening_factors_pre_grainsize_calculation = dike_mode ? dikemode_compute_strain_weakening_factors(j, in.composition[i]) : strain_rheology.compute_strain_weakening_factors(j, in.composition[i]);
			    const double current_cohesion_pre_grainsize_calculation = drucker_prager_parameters.cohesions[j] * weakening_factors_pre_grainsize_calculation[0];
                const double current_friction_pre_grainsize_calculation = serp_value > 0.1 && new_crust_phase < 0.2 ? drucker_prager_parameters.angles_internal_friction[j] * serp_weakening * weakening_factors_pre_grainsize_calculation[1]: drucker_prager_parameters.angles_internal_friction[j] * weakening_factors_pre_grainsize_calculation[1];
			    const double pressure_for_plasticity_pre_grainsize_calculation = std::max(in.pressure[i],0.0);
			    yield_stress_pre_grainsize_calculation = drucker_prager_plasticity.compute_yield_stress(current_cohesion_pre_grainsize_calculation,
                                                                                       current_friction_pre_grainsize_calculation,
                                                                                       pressure_for_plasticity_pre_grainsize_calculation,
                                                                                       drucker_prager_parameters.max_yield_stress);
			}
			
			const double water_content_for_viscosity = initialization_past ?
			                            std::max(30. , std::abs(in.composition[i][this->introspection().compositional_index_for_name("water_content")]))
			                           : init_water_content[j];
			//20260203 可删除rift_grain_size
			const double last_grain_size = this->introspection().compositional_name_exists("rift_grain_size") ? initialization_past ? 
			                               std::max(std::abs(in.composition[i][this->introspection().compositional_index_for_name("rift_grain_size")])
			                               , diffusion_creep_grainsize.minimum_grain_size) 
			                               : diffusion_creep_grainsize.initial_grain_size[j] : diffusion_creep_grainsize.initial_grain_size[j];
			//20240702
            const double current_grain_size = diffusion_creep_grainsize.compute_grainsize(edot_ii, in.pressure[i], temperature_for_viscosity, j,
                                                                                          last_grain_size, posx, posy, water_content_for_viscosity, 
																						  init_water_content[j], yield_stress_pre_grainsize_calculation, melt_fraction,
																						  phase_function_values, n_phases_per_composition);
            const double viscosity_diffusion = diffusion_creep_grainsize.compute_viscosity(in.pressure[i], temperature_for_viscosity, j, current_grain_size,
                                                                                 water_content_for_viscosity, init_water_content[j], phase_function_values,
                                                                                 n_phases_per_composition);
            // Step 1b: compute viscosity from dislocation creep law
            const double viscosity_dislocation = dislocation_creep_grainsize.compute_viscosity(edot_ii, in.pressure[i], temperature_for_viscosity, water_content_for_viscosity,
			                                                                         init_water_content[j], j,
                                                                                     phase_function_values,
                                                                                     n_phases_per_composition);
            // Step 1c: select what form of viscosity to use (diffusion, dislocation, fk, or composite)
            double viscosity_pre_yield = 0.0;
            switch (viscous_flow_law)
              {
                case diffusion:
                {
                  viscosity_pre_yield = viscosity_diffusion;
                  break;
                }
                case dislocation:
                {
                  viscosity_pre_yield = viscosity_dislocation;
                  break;
                }
                case frank_kamenetskii:
                {
                  viscosity_pre_yield = frank_kamenetskii_rheology->compute_viscosity(in.temperature[i], j);
                  break;
                }
                case composite:
                {
                  viscosity_pre_yield = (viscosity_diffusion * viscosity_dislocation)/
                                        (viscosity_diffusion + viscosity_dislocation);
                  break;
                }
                default:
                {
                  AssertThrow(false, ExcNotImplemented());
                  break;
                }
              }

            // Step 1d: compute viscosity from Peierls creep law and harmonically average with current viscosities
            if (use_peierls_creep)
              {
                const double viscosity_peierls = peierls_creep->compute_viscosity(edot_ii, in.pressure[i], temperature_for_viscosity, j,
                                                                                  phase_function_values,
                                                                                  n_phases_per_composition);
                viscosity_pre_yield = (viscosity_pre_yield * viscosity_peierls) / (viscosity_pre_yield + viscosity_peierls);
              }

            // Step 1e: multiply the viscosity by a constant (default value is 1)
            viscosity_pre_yield = constant_viscosity_prefactors.compute_viscosity(viscosity_pre_yield, j);
            // Step 2: calculate the viscous stress magnitude
            // and strain rate. If requested compute visco-elastic contributions.
            double current_edot_ii = edot_ii;
            if (use_elasticity)
              {
                const std::vector<double> &elastic_shear_moduli = elastic_rheology.get_elastic_shear_moduli();

                if (use_reference_strainrate == true)
                  current_edot_ii = ref_strain_rate;
                else
                  {
                    const double viscoelastic_strain_rate_invariant = elastic_rheology.calculate_viscoelastic_strain_rate(in.strain_rate[i],
                                                                      stress_old,
                                                                      elastic_shear_moduli[j]);

                    current_edot_ii = std::max(viscoelastic_strain_rate_invariant,
                                               min_strain_rate);
                  }
                // Step 2a: calculate viscoelastic (effective) viscosity
                viscosity_pre_yield = elastic_rheology.calculate_viscoelastic_viscosity(viscosity_pre_yield,
                                                                                        elastic_shear_moduli[j]);
              }
            // Step 2b: calculate current (viscous or viscous + elastic) stress magnitude
            double current_stress = 2. * viscosity_pre_yield * current_edot_ii;
            // Step 3: strain weakening
            // Step 3a: calculate strain weakening factors for the cohesion, friction, and pre-yield viscosity
            // If no brittle and/or viscous strain weakening is applied, the factors are 1.
            const std::array<double, 3> weakening_factors = dike_mode ? dikemode_compute_strain_weakening_factors(j, in.composition[i]) : strain_rheology.compute_strain_weakening_factors(j, in.composition[i]);
            // Step 3b: calculate weakened friction, cohesion, and pre-yield viscosity and adjust the current_stress accordingly
            const double current_cohesion = drucker_prager_parameters.cohesions[j] * weakening_factors[0];
			//20260127 蛇纹岩化相关：根据蛇纹岩化程度修改内摩擦角，蛇纹岩化程度大于0.1时，内摩擦角调整为输入值的0.1倍
			//注：未测试与岩墙模式的兼容性
			//20260225 蛇纹岩化区域参与塑性应变弱化
			const double current_friction = serp_value > 0.1 && new_crust_phase < 0.2 ? drucker_prager_parameters.angles_internal_friction[j] * serp_weakening * weakening_factors[1]: drucker_prager_parameters.angles_internal_friction[j] * weakening_factors[1];
            viscosity_pre_yield *= weakening_factors[2];
            current_stress *= weakening_factors[2];
            // Step 4: plastic yielding
            // Determine if the pressure used in Drucker Prager plasticity will be capped at 0 (default).
            // This may be necessary in models without gravity and the dynamic stresses are much higher
            // than the lithostatic pressure.
            double pressure_for_plasticity = in.pressure[i];
            if (allow_negative_pressures_in_plasticity == false)
              pressure_for_plasticity = std::max(in.pressure[i],0.0);
            // Step 4a: calculate Drucker-Prager yield stress
            const double yield_stress = drucker_prager_plasticity.compute_yield_stress(current_cohesion,
                                                                                       current_friction,
                                                                                       pressure_for_plasticity,
                                                                                       drucker_prager_parameters.max_yield_stress);
            // Step 4b: select if yield viscosity is based on Drucker Prager or stress limiter rheology
            double viscosity_yield = viscosity_pre_yield;
            switch (yield_mechanism)
              {
                case stress_limiter:
                {
                  //Step 4b-1: always rescale the viscosity back to the yield surface
                  const double viscosity_limiter = yield_stress / (2.0 * ref_strain_rate)
                                                   * std::pow((edot_ii/ref_strain_rate),
                                                              1./exponents_stress_limiter[j] - 1.0);
                  viscosity_yield = 1. / ( 1./viscosity_limiter + 1./viscosity_pre_yield);
                  break;
                }
                case drucker_prager:
                {
                  // Step 4b-2: if the current stress is greater than the yield stress,
                  // rescale the viscosity back to yield surface
                  if (current_stress >= yield_stress)
                    {
                      viscosity_yield = drucker_prager_plasticity.compute_viscosity(current_cohesion,
                                                                                    current_friction,
                                                                                    pressure_for_plasticity,
                                                                                    current_edot_ii,
                                                                                    drucker_prager_parameters.max_yield_stress,
                                                                                    viscosity_pre_yield);
                      output_parameters.composition_yielding[j] = true;
                    }
                  break;
                }
                default:
                {
                  AssertThrow(false, ExcNotImplemented());
                  break;
                }
              }
            // Step 5: limit the viscosity with specified minimum and maximum bounds
			//20250122 新增最小粘度切换
			//20250909 测试：使用应力边界时将，将边界区域设置为最大粘滞系数
			//20251117 测试：调整岩墙内的最大粘滞系数
			//20260112 岩墙内粘度调整放到上一级程序内
			const double min_vis_depend = this->get_time() >= switch_time_in_year * 3600 * 24 * 365.25 ? switch_min_visc : min_visc;
            //output_parameters.composition_viscosities[j] = in_dike ? std::min(std::max(viscosity_yield, min_vis_depend), 1e21) : std::min(std::max(viscosity_yield, min_vis_depend), max_visc);
			output_parameters.composition_viscosities[j] = std::min(std::max(viscosity_yield, min_vis_depend), max_visc);
          }
        return output_parameters;
      }

      template <int dim>
      void
      ViscoPlasticGrainSize<dim>::
      compute_viscosity_derivatives(const unsigned int i,
                                    const std::vector<double> &volume_fractions,
                                    const std::vector<double> &composition_viscosities,
                                    const MaterialModel::MaterialModelInputs<dim> &in,
                                    MaterialModel::MaterialModelOutputs<dim> &out,
                                    const std::vector<double> &phase_function_values,
                                    const std::vector<unsigned int> &n_phases_per_composition)const
      {
        MaterialModel::MaterialModelDerivatives<dim> *derivatives =
          out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >();
        if (derivatives != nullptr)
          {
            // compute derivatives if necessary
            std::vector<SymmetricTensor<2,dim> > composition_viscosities_derivatives(volume_fractions.size());
            std::vector<double> composition_dviscosities_dpressure(volume_fractions.size());
            const double finite_difference_accuracy = 1e-7;
            // A new material model inputs variable that uses the strain rate and pressure difference.
            MaterialModel::MaterialModelInputs<dim> in_derivatives = in;
            // For each independent component, compute the derivative.
            for (unsigned int component = 0; component < SymmetricTensor<2,dim>::n_independent_components; ++component)
              {
                const TableIndices<2> strain_rate_indices = SymmetricTensor<2,dim>::unrolled_to_component_indices (component);
                // components that are not on the diagonal are multiplied by 0.5, because the symmetric tensor
                // is modified by 0.5 in both symmetric directions (xy/yx) simultaneously and we compute the combined
                // derivative
                const SymmetricTensor<2,dim> strain_rate_difference = in.strain_rate[i]
                                                                      + std::max(std::fabs(in.strain_rate[i][strain_rate_indices]), min_strain_rate)
                                                                      * (component > dim-1 ? 0.5 : 1 )
                                                                      * finite_difference_accuracy
                                                                      * Utilities::nth_basis_for_symmetric_tensors<dim>(component);
                in_derivatives.strain_rate[i] = strain_rate_difference;
                std::vector<double> eta_component =
                  calculate_isostrain_viscosities(in_derivatives, i, volume_fractions,
                                                  phase_function_values, n_phases_per_composition).composition_viscosities;

                // For each composition of the independent component, compute the derivative.
                for (unsigned int composition_index = 0; composition_index < eta_component.size(); ++composition_index)
                  {
                    // compute the difference between the viscosity with and without the strain-rate difference.
                    double viscosity_derivative = eta_component[composition_index] - composition_viscosities[composition_index];
                    if (viscosity_derivative != 0)
                      {
                        // when the difference is non-zero, divide by the difference.
                        viscosity_derivative /= std::max(std::fabs(strain_rate_difference[strain_rate_indices]), min_strain_rate)
                                                * finite_difference_accuracy;
                      }
                    composition_viscosities_derivatives[composition_index][strain_rate_indices] = viscosity_derivative;
                  }
              }

            /**
             * Now compute the derivative of the viscosity to the pressure
             */
            const double pressure_difference = in.pressure[i] + (std::fabs(in.pressure[i]) * finite_difference_accuracy);
            in_derivatives.pressure[i] = pressure_difference;
            // Modify the in_derivatives object again to take the original strain rate.
            in_derivatives.strain_rate[i] = in.strain_rate[i];
            const std::vector<double> viscosity_difference =
              calculate_isostrain_viscosities(in_derivatives, i, volume_fractions,
                                              phase_function_values, n_phases_per_composition).composition_viscosities;
            for (unsigned int composition_index = 0; composition_index < viscosity_difference.size(); ++composition_index)
              {
                double viscosity_derivative = viscosity_difference[composition_index] - composition_viscosities[composition_index];
                if (viscosity_difference[composition_index] != 0)
                  {
                    if (in.pressure[i] != 0)
                      {
                        viscosity_derivative /= std::fabs(in.pressure[i]) * finite_difference_accuracy;
                      }
                    else
                      {
                        viscosity_derivative = 0;
                      }
                  }
                composition_dviscosities_dpressure[composition_index] = viscosity_derivative;
              }
            double viscosity_averaging_p = 0; // Geometric
            if (viscosity_averaging == MaterialUtilities::harmonic)
              viscosity_averaging_p = -1;
            if (viscosity_averaging == MaterialUtilities::arithmetic)
              viscosity_averaging_p = 1;
            if (viscosity_averaging == MaterialUtilities::maximum_composition)
              viscosity_averaging_p = 1000;
            derivatives->viscosity_derivative_wrt_strain_rate[i] =
              Utilities::derivative_of_weighted_p_norm_average(out.viscosities[i],
                                                               volume_fractions,
                                                               composition_viscosities,
                                                               composition_viscosities_derivatives,
                                                               viscosity_averaging_p);
            derivatives->viscosity_derivative_wrt_pressure[i] =
              Utilities::derivative_of_weighted_p_norm_average(out.viscosities[i],
                                                               volume_fractions,
                                                               composition_viscosities,
                                                               composition_dviscosities_dpressure,
                                                               viscosity_averaging_p);
          }
      }



      template <int dim>
      ComponentMask
      ViscoPlasticGrainSize<dim>::
      get_volumetric_composition_mask() const
      {
        // Store which components to exclude during the volume fraction computation.
        ComponentMask composition_mask = strain_rheology.get_strain_composition_mask();
        //Providing a mask for a field to hold the sediment age. We don't want this to be included in any averages.
		if (this->introspection().compositional_name_exists("sediment_age"))
            composition_mask.set(this->introspection().compositional_index_for_name("sediment_age"),false);
		//模型新增"rift_grain_size"作为实现存储粒径以及可视化的媒介，不参与相关的模型计算
		if (this->introspection().compositional_name_exists("rift_grain_size"))
            composition_mask.set(this->introspection().compositional_index_for_name("rift_grain_size"), false);
		//新增“viscous_strain_grainsize”观察引入粒度计算时黏性应变的情况。
		//此变量不涉及应变弱化问题，对模型演化过程无影响。
		if (this->introspection().compositional_name_exists("viscous_strain_grainsize"))
			composition_mask.set(this->introspection().compositional_index_for_name("viscous_strain_grainsize"), false);
		//新增"water_content"存储地幔部分的含水状态
	    if (this->introspection().compositional_name_exists("water_content"))
			composition_mask.set(this->introspection().compositional_index_for_name("water_content"), false);
		//新增“melt_fraction”记录熔融状态
		if (this->introspection().compositional_name_exists("melt_fraction"))
			composition_mask.set(this->introspection().compositional_index_for_name("melt_fraction"), false);
		//新增“melt_record”记录每一步的熔融变化
		if (this->introspection().compositional_name_exists("melt_record"))
			composition_mask.set(this->introspection().compositional_index_for_name("melt_record"), false);
		//20240614 新增“max_melt”记录每个区域发生的最大熔融程度，用于计算密度变化
		if (this->introspection().compositional_name_exists("max_melt"))
			composition_mask.set(this->introspection().compositional_index_for_name("max_melt"), false);
        if (use_elasticity)
          {
            for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components ; ++i)
              composition_mask.set(i,false);
          }
		//20251110
		//岩墙相关,strain_dependent内根据弱化模式确定成分场，此处更改为根据成分场存在与否判定
		if(dike_mode)
		{
			if (this->introspection().compositional_name_exists("plastic_strain"))
			    composition_mask.set(this->introspection().compositional_index_for_name("plastic_strain"), false);
		    if (this->introspection().compositional_name_exists("viscous_strain"))
			    composition_mask.set(this->introspection().compositional_index_for_name("viscous_strain"), false);
		}
		//20260121
		//添加蛇纹岩化相关
		//蛇纹岩化程度（最大为1）
		if (this->introspection().compositional_name_exists("serp"))
			composition_mask.set(this->introspection().compositional_index_for_name("serp"), false);
		//蛇纹岩化速率(每年)
		if (this->introspection().compositional_name_exists("d_serp"))
			composition_mask.set(this->introspection().compositional_index_for_name("d_serp"), false);
		
        return composition_mask;
      }


      template <int dim>
      void
      ViscoPlasticGrainSize<dim>::declare_parameters (ParameterHandler &prm)
      {
        Rheology::StrainDependent<dim>::declare_parameters (prm);
        Rheology::Elasticity<dim>::declare_parameters (prm);
        // Reference and minimum/maximum values
        prm.declare_entry ("Minimum strain rate", "1.0e-20", Patterns::Double (0.),
                           "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}.");
        prm.declare_entry ("Reference strain rate","1.0e-15",Patterns::Double (0.),
                           "Reference strain rate for first time step. Units: \\si{\\per\\second}.");
        prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double (0.),
                           "Lower cutoff for effective viscosity. Units: \\si{\\pascal\\second}.");
        prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double (0.),
                           "Upper cutoff for effective viscosity. Units: \\si{\\pascal\\second}.");
        prm.declare_entry ("Reference viscosity", "1e22", Patterns::Double (0.),
                           "Reference viscosity for nondimensionalization. "
                           "To understand how pressure scaling works, take a look at "
                           "\\cite{KHB12}. In particular, the value of this parameter "
                           "would not affect the solution computed by \\aspect{} if "
                           "we could do arithmetic exactly; however, computers do "
                           "arithmetic in finite precision, and consequently we need to "
                           "scale quantities in ways so that their magnitudes are "
                           "roughly the same. As explained in \\cite{KHB12}, we scale "
                           "the pressure during some computations (never visible by "
                           "users) by a factor that involves a reference viscosity. This "
                           "parameter describes this reference viscosity."
                           "\n\n"
                           "For problems with a constant viscosity, you will generally want "
                           "to choose the reference viscosity equal to the actual viscosity. "
                           "For problems with a variable viscosity, the reference viscosity "
                           "should be a value that adequately represents the order of "
                           "magnitude of the viscosities that appear, such as an average "
                           "value or the value one would use to compute a Rayleigh number."
                           "\n\n"
                           "Units: \\si{\\pascal\\second}.");
        // Rheological parameters
        prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                           Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                           "When more than one compositional field is present at a point "
                           "with different viscosities, we need to come up with an average "
                           "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                           "geometric, or maximum composition.");
        prm.declare_entry ("Viscous flow law", "composite",
                           Patterns::Selection("diffusion|dislocation|frank kamenetskii|composite"),
                           "Select what type of viscosity law to use between diffusion, "
                           "dislocation, frank kamenetskii, and composite options. Soon there will be an option "
                           "to select a specific flow law for each assigned composition ");
        prm.declare_entry ("Yield mechanism", "drucker",
                           Patterns::Selection("drucker|limiter"),
                           "Select what type of yield mechanism to use between Drucker Prager "
                           "and stress limiter options.");
        prm.declare_entry ("Allow negative pressures in plasticity", "false",
                           Patterns::Bool (),
                           "Whether to allow negative pressures to be used in the computation "
                           "of plastic yield stresses and viscosities. Setting this parameter "
                           "to true may be advantageous in models without gravity where the "
                           "dynamic stresses are much higher than the lithostatic pressure. "
                           "If false, the minimum pressure in the plasticity formulation will "
                           "be set to zero.");
		//ratio of fractions modification
		prm.declare_entry("Ratio of fractions modification", "0.5",
                          Patterns::Double(0.),
                         "Units: none.");
        // Diffusion creep parameters
        Rheology::DiffusionCreepGrainSize<dim>::declare_parameters(prm);

        // Dislocation creep parameters
        Rheology::DislocationCreepGrainSize<dim>::declare_parameters(prm);

        // Frank-Kamenetskii viscosity parameters
        Rheology::FrankKamenetskii<dim>::declare_parameters(prm);

        // Peierls creep parameters
        Rheology::PeierlsCreep<dim>::declare_parameters(prm);

        prm.declare_entry ("Include Peierls creep", "false",
                           Patterns::Bool (),
                           "Whether to include Peierls creep in the rheological formulation.");

        // Constant viscosity prefactor parameters
        Rheology::ConstantViscosityPrefactors<dim>::declare_parameters(prm);

        // Drucker Prager plasticity parameters
        Rheology::DruckerPrager<dim>::declare_parameters(prm);

        // Stress limiter parameters
        prm.declare_entry ("Stress limiter exponents", "1.0",
                           Patterns::List(Patterns::Double (0.)),
                           "List of stress limiter exponents, $n_{\\text{lim}}$, "
                           "for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "Units: none.");

        // Temperature gradient in viscosity laws to include an adiabat (note units of K/Pa)
        prm.declare_entry ("Adiabat temperature gradient for viscosity", "0.0", Patterns::Double (0.),
                           "Add an adiabatic temperature gradient to the temperature used in the flow law "
                           "so that the activation volume is consistent with what one would use in a "
                           "earth-like (compressible) model. Default is set so this is off. "
                           "Note that this is a linear approximation of the real adiabatic gradient, which "
                           "is okay for the upper mantle, but is not really accurate for the lower mantle. "
                           "Using a pressure gradient of 32436 Pa/m, then a value of "
                           "0.3 K/km = 0.0003 K/m = 9.24e-09 K/Pa gives an earth-like adiabat."
                           "Units: \\si{\\kelvin\\per\\pascal}.");

        prm.declare_entry ("Include viscoelasticity", "false",
                           Patterns::Bool (),
                           "Whether to include elastic effects in the rheological formulation.");
		//20250122
		//新增最小粘度切换，以解决对流和裂谷过程的过渡问题
		prm.declare_entry ("Switch minimum viscosity", "3e19", Patterns::Double (0.),
                           "Lower cutoff for effective viscosity. Units: \\si{\\pascal\\second}.");
		prm.declare_entry ("Switch time", "10000e6", Patterns::Double (0.),
                           "Switch time in year. Units: year.");
	    //20251101
		prm.declare_entry ("Strain weakening mechanism for dike mode", "default",
                           Patterns::Selection("none|plastic weakening with plastic strain only|plastic weakening with plastic strain and viscous weakening with viscous strain|plastic weakening with plastic strain and viscous weakening with grain size|default"),
                           "Strain weakening mechanism for dike mode.");
		
      }



      template <int dim>
      void
      ViscoPlasticGrainSize<dim>::parse_parameters (ParameterHandler &prm,
	                                                const double max_x,
	                                                const double max_y,
													const bool add_cratontic,
													const bool center_cratontic,
													const bool exclude_mantle_lower,
													const bool input_dike_mode,
									                const bool input_constant_M,
													const bool input_constant_dike_length,
									                const double input_M_for_dike,
									                const double input_dike_width,
									                const double input_dike_length,
									                const double input_spreading_rate,
													const std::vector<double> layer_thickness,
													const std::vector<double> cratontic_thicknesses,
													const std::vector<double> cratontic_boundary,
                                                    const std::shared_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {

        //20260320 新增蛇纹岩化相关参数，在visco_plastic_grainsize主函数内声明
		serp_weakening = prm.get_double ("Plastic wekening factor for serpentinization");
		//20241104
		//接收max_x，用于判定是否在俯冲输入区域内
		box_x = max_x;
		//20251104
		//接收岩墙相关参数，用于修正应变速率
		melt_info_path = this->get_output_directory() + "melt_info/";
		dike_mode = input_dike_mode;
		constant_M = input_constant_M;
		//20251201 新增 constant_dike_length
		constant_dike_length = input_constant_dike_length;
		M_for_dike = input_M_for_dike;
		dike_width = input_dike_width;
		dike_length = input_dike_length;
		spreading_rate = input_spreading_rate;
		
		//20251110
		//新增岩墙相关：累积应变类型判定
		if (prm.get ("Strain weakening mechanism for dike mode") == "none")
          dike_weakening_mechanism = dike_none;
        else if (prm.get ("Strain weakening mechanism for dike mode") == "plastic weakening with plastic strain only")
          dike_weakening_mechanism = dike_plastic_weakening_with_plastic_strain_only;
        else if (prm.get ("Strain weakening mechanism for dike mode") == "plastic weakening with plastic strain and viscous weakening with viscous strain")
          dike_weakening_mechanism = dike_plastic_weakening_with_plastic_strain_and_viscous_weakening_with_viscous_strain;
        else if (prm.get ("Strain weakening mechanism for dike mode") == "plastic weakening with plastic strain and viscous weakening with grain size")
          dike_weakening_mechanism = dike_plastic_weakening_with_plastic_strain_and_viscous_weakening_with_grain_size;
        else if (prm.get ("Strain weakening mechanism for dike mode") == "default")
          dike_weakening_mechanism = dike_none;
	    else
          AssertThrow(false, ExcMessage("Not a valid Strain weakening mechanism for dike mode!"));
	    if (dike_weakening_mechanism == dike_plastic_weakening_with_plastic_strain_only
            || dike_weakening_mechanism == dike_plastic_weakening_with_plastic_strain_and_viscous_weakening_with_viscous_strain
	        || dike_weakening_mechanism == dike_plastic_weakening_with_plastic_strain_and_viscous_weakening_with_grain_size)
          {
            AssertThrow(this->introspection().compositional_name_exists("plastic_strain"),
                        ExcMessage("Material model plastic weakening with plastic strain only works if there is a "
                                   "compositional field called plastic_strain."));
          }

        if (dike_weakening_mechanism == dike_plastic_weakening_with_plastic_strain_and_viscous_weakening_with_viscous_strain)
          {
            AssertThrow(this->introspection().compositional_name_exists("viscous_strain"),
                        ExcMessage("Material model visco_plastic with viscous strain weakening only works if there is a "
                                   "compositional field called viscous_strain."));
          }
		if (dike_weakening_mechanism == dike_plastic_weakening_with_plastic_strain_and_viscous_weakening_with_grain_size)
          {
            AssertThrow(this->introspection().compositional_name_exists("viscous_strain_grainsize"),
                        ExcMessage("Material model visco_plastic with grain size only works if there is a "
                                   "compositional field called viscous_strain_grainsize."));
          }
		if (dike_weakening_mechanism != dike_none)
          {
            AssertThrow((this->get_parameters().nonlinear_solver ==
                         Parameters<dim>::NonlinearSolver::single_Advection_single_Stokes
                         ||
                         this->get_parameters().nonlinear_solver ==
                         Parameters<dim>::NonlinearSolver::single_Advection_iterated_Stokes
                         ||
                         this->get_parameters().nonlinear_solver ==
                         Parameters<dim>::NonlinearSolver::single_Advection_iterated_Newton_Stokes
                         ||
                         this->get_parameters().nonlinear_solver ==
                         Parameters<dim>::NonlinearSolver::single_Advection_iterated_defect_correction_Stokes),
                        ExcMessage("The material model will only work with the nonlinear "
                                   "solver schemes 'single Advection, single Stokes', "
                                   "'single Advection, iterated Stokes', "
                                   "'single Advection, iterated Newton Stokes', and "
                                   "'single Advection, iterated defect correction Stokes' "
                                   "when strain weakening is enabled, because more than one nonlinear "
                                   "advection iteration will result in the incorrect value of strain."));
          }
		
        // increment by one for background:
        const unsigned int n_fields = this->n_compositional_fields() + 1;

        strain_rheology.initialize_simulator (this->get_simulator());
        strain_rheology.parse_parameters(prm);

        use_elasticity = prm.get_bool ("Include viscoelasticity");

        if (use_elasticity)
          {
            elastic_rheology.initialize_simulator (this->get_simulator());
            elastic_rheology.parse_parameters(prm);
          }

        // Reference and minimum/maximum values
        min_strain_rate = prm.get_double("Minimum strain rate");
        ref_strain_rate = prm.get_double("Reference strain rate");
        min_visc = prm.get_double ("Minimum viscosity");
        max_visc = prm.get_double ("Maximum viscosity");
        ref_visc = prm.get_double ("Reference viscosity");
		fraction_ratio = prm.get_double("Ratio of fractions modification");
		init_water_content = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Initial water content"))),
                                                                            n_fields,
                                                                            "Initial water content");
		input_from_flank = prm.get_bool ("Put material input channel to flank area");																	

        AssertThrow(max_visc >= min_visc, ExcMessage("Maximum viscosity should be larger or equal to the minimum viscosity. "));

        viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                              prm);
		//20250122
		switch_min_visc = prm.get_double ("Switch minimum viscosity");
		switch_time_in_year = prm.get_double ("Switch time");
		
		//20250916
		input_total_force = prm.get_double ("Total force in one side");

        // Rheological parameters
        if (prm.get ("Viscous flow law") == "composite")
          viscous_flow_law = composite;
        else if (prm.get ("Viscous flow law") == "diffusion")
          viscous_flow_law = diffusion;
        else if (prm.get ("Viscous flow law") == "dislocation")
          viscous_flow_law = dislocation;
        else if (prm.get ("Viscous flow law") == "frank kamenetskii")
          viscous_flow_law = frank_kamenetskii;
        else
          AssertThrow(false, ExcMessage("Not a valid viscous flow law"));

        if (prm.get ("Yield mechanism") == "drucker")
          yield_mechanism = drucker_prager;
        else if (prm.get ("Yield mechanism") == "limiter")
          yield_mechanism = stress_limiter;
        else
          AssertThrow(false, ExcMessage("Not a valid yield mechanism."));

        AssertThrow(use_elasticity == false || yield_mechanism == drucker_prager,
                    ExcMessage("Elastic behavior is only tested with the "
                               "'drucker prager' plasticity option."));

        allow_negative_pressures_in_plasticity = prm.get_bool ("Allow negative pressures in plasticity");

        // Diffusion creep parameters
        diffusion_creep_grainsize.initialize_simulator (this->get_simulator());
        diffusion_creep_grainsize.parse_parameters(prm, 
		                                           max_x,
		                                           max_y,
									               add_cratontic,
									               center_cratontic,
									               exclude_mantle_lower,
									               layer_thickness,
									               cratontic_thicknesses,
									               cratontic_boundary,
		                                           expected_n_phases_per_composition);

        // Dislocation creep parameters
        dislocation_creep_grainsize.initialize_simulator (this->get_simulator());
        dislocation_creep_grainsize.parse_parameters(prm, expected_n_phases_per_composition);

        // Frank Kamenetskii viscosity parameters
        if (viscous_flow_law == frank_kamenetskii)
          {
            frank_kamenetskii_rheology = std_cxx14::make_unique<Rheology::FrankKamenetskii<dim>>();
            frank_kamenetskii_rheology->initialize_simulator (this->get_simulator());
            frank_kamenetskii_rheology->parse_parameters(prm);
          }

        // Peierls creep parameters
        use_peierls_creep = prm.get_bool ("Include Peierls creep");
        if (use_peierls_creep)
          {
            peierls_creep = std_cxx14::make_unique<Rheology::PeierlsCreep<dim>>();
            peierls_creep->initialize_simulator (this->get_simulator());
            peierls_creep->parse_parameters(prm);
          }

        // Constant viscosity prefactor parameters
        constant_viscosity_prefactors.initialize_simulator (this->get_simulator());
        constant_viscosity_prefactors.parse_parameters(prm);

        // Plasticity parameters
        drucker_prager_parameters = drucker_prager_plasticity.parse_parameters(this->n_compositional_fields()+1,
                                                                               prm);

        // Stress limiter parameter
        exponents_stress_limiter  = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Stress limiter exponents"))),
                                                                            n_fields,
                                                                            "Stress limiter exponents");

        // Include an adiabat temperature gradient in flow laws
        adiabatic_temperature_gradient_for_viscosity = prm.get_double("Adiabat temperature gradient for viscosity");
        if (this->get_heating_model_manager().adiabatic_heating_enabled())
          AssertThrow (adiabatic_temperature_gradient_for_viscosity == 0.0,
                       ExcMessage("If adiabatic heating is enabled you should not add another adiabatic gradient"
                                  "to the temperature for computing the viscosity, because the ambient"
                                  "temperature profile already includes the adiabatic gradient."));
      }

      template <int dim>
      void
      ViscoPlasticGrainSize<dim>::
      fill_grainsize_outputs(const MaterialModel::MaterialModelInputs<dim>& in,
                            MaterialModel::MaterialModelOutputs<dim>& out,
                            const unsigned int i,
							const bool plastic_yielding,
                            const std::vector<double>& volume_fractions,
                            const std::vector<double>& phase_function_values,
                            const std::vector<unsigned int>& n_phases_per_composition)const
      {
		  //20260203 可删除rift_grain_size
		  //此函数在上一层的调用部分限制
          const double fill_edot_ii = std::max(std::sqrt(std::fabs(second_invariant(in.strain_rate[i]))),min_strain_rate);
          const std::vector<double>fill_temp_volume = volume_fractions;
		  const double fill_posx = in.position[i][0];
          const double fill_posy = in.position[i][1];
          const double fill_last_grain_size = std::max(std::min(std::abs(in.composition[i][this->introspection().compositional_index_for_name("rift_grain_size")]), diffusion_creep_grainsize.max_grain_size)
		                                               , diffusion_creep_grainsize.minimum_grain_size);
		  const double fill_melt_fraction = this->introspection().compositional_name_exists("melt_fraction") ? std::abs(in.composition[i][this->introspection().compositional_index_for_name("melt_fraction")]) : 0.;
          double fill_output_grainsize = 0;
          double fill_minus_grainsize = 0;
          //输出到conpositional fields中
		  //20260127 蛇纹岩化相关：根据蛇纹岩化程度修改内摩擦角，蛇纹岩化程度大于0.1时，内摩擦角调整为输入值的0.1倍
		  //20260225 蛇纹岩化区域参与塑性应变弱化
		  const double fill_new_crust_phase = this->introspection().compositional_name_exists("new_crust") ? volume_fractions[this->introspection().compositional_index_for_name("new_crust") + 1] : 0.;
		  const double fill_serp_value = this->introspection().compositional_name_exists("serp") ? in.composition[i][this->introspection().compositional_index_for_name("serp")] : 0.;
          for (unsigned int j = 0; j < volume_fractions.size(); ++j)
          {
			  const std::array<double, 3> weakening_factors_pre_grainsize_calculation = dike_mode ? dikemode_compute_strain_weakening_factors(j, in.composition[i]) : strain_rheology.compute_strain_weakening_factors(j, in.composition[i]);
			  const double current_cohesion_pre_grainsize_calculation = drucker_prager_parameters.cohesions[j] * weakening_factors_pre_grainsize_calculation[0];
              const double current_friction_pre_grainsize_calculation = fill_serp_value > 0.1 && fill_new_crust_phase < 0.2 ? drucker_prager_parameters.angles_internal_friction[j] * serp_weakening * weakening_factors_pre_grainsize_calculation[1]: drucker_prager_parameters.angles_internal_friction[j] * weakening_factors_pre_grainsize_calculation[1];
			  const double pressure_for_plasticity_pre_grainsize_calculation = std::max(in.pressure[i],0.0);
			  const double yield_stress_pre_grainsize_calculation = drucker_prager_plasticity.compute_yield_stress(current_cohesion_pre_grainsize_calculation,
                                                                                         current_friction_pre_grainsize_calculation,
                                                                                         pressure_for_plasticity_pre_grainsize_calculation,
                                                                                         drucker_prager_parameters.max_yield_stress);
			  
			  
			  
			  //含水量最多降低到30 ppm
			  const double fill_water_content_for_viscosity = std::max(30. , std::abs(in.composition[i][this->introspection().compositional_index_for_name("water_content")]));
	
              const double fill_temperature_for_viscosity = in.temperature[i] + adiabatic_temperature_gradient_for_viscosity * in.pressure[i];
			  //20240702
              const double fill_current_grain_size = diffusion_creep_grainsize.compute_grainsize(fill_edot_ii, in.pressure[i], fill_temperature_for_viscosity, j,
                                                                                                 fill_last_grain_size, fill_posx, fill_posy, fill_water_content_for_viscosity, init_water_content[j], yield_stress_pre_grainsize_calculation, fill_melt_fraction, phase_function_values,
                                                                                                 n_phases_per_composition);
              fill_output_grainsize += fill_current_grain_size * fill_temp_volume[j];
              fill_minus_grainsize += fill_last_grain_size * fill_temp_volume[j];
          }
          out.reaction_terms[i][this->introspection().compositional_index_for_name("rift_grain_size")] = fill_output_grainsize - in.composition[i][this->introspection().compositional_index_for_name("rift_grain_size")];
		  //输出粘性应变量
		  if (this->introspection().compositional_name_exists("viscous_strain_grainsize") && !plastic_yielding)
		  {
			  const double delta_e_ii = fill_edot_ii * this->get_timestep();
			  //20241104
			  //修正输入区域的应变累积量
			  //20251201 岩墙模式下不在这里更新
			  if(!dike_mode)
				  out.reaction_terms[i][this->introspection().compositional_index_for_name("viscous_strain_grainsize")] = input_from_flank && (fill_posx <= 10e3 || box_x - fill_posx <= 10e3) ?
			                                                                                                          0. -  in.composition[i][this->introspection().compositional_index_for_name("viscous_strain_grainsize")]
			                                                                                                        : std::max(delta_e_ii, -in.composition[i][this->introspection().compositional_index_for_name("viscous_strain_grainsize")]);
		  }
      }

      template <int dim>
      void
      ViscoPlasticGrainSize<dim>::create_plastic_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        if (out.template get_additional_output<PlasticAdditionalOutputsGrainSize<dim> >() == nullptr)
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std_cxx14::make_unique<PlasticAdditionalOutputsGrainSize<dim>> (n_points));
          }
      }

      template <int dim>
      void
      ViscoPlasticGrainSize<dim>::
      fill_plastic_outputs(const unsigned int i,
                           const std::vector<double> &volume_fractions,
                           const bool plastic_yielding,
                           const MaterialModel::MaterialModelInputs<dim> &in,
                                 MaterialModel::MaterialModelOutputs<dim> &out,
						   const std::vector<double>& phase_function_values,
                           const std::vector<unsigned int>& n_phases_per_composition) const
      {
          PlasticAdditionalOutputsGrainSize<dim> *plastic_out = out.template get_additional_output<PlasticAdditionalOutputsGrainSize<dim> >();
        if (plastic_out != nullptr)
          {
			const std::vector<double>out_temp_volume = volume_fractions;
			const double out_edot_ii = std::max(std::sqrt(std::fabs(second_invariant(in.strain_rate[i]))),min_strain_rate);
			const double out_posx = in.position[i][0];
            const double out_posy = in.position[i][1];
			const double out_melt_fraction = this->introspection().compositional_name_exists("melt_fraction") ? std::abs(in.composition[i][this->introspection().compositional_index_for_name("melt_fraction")]) : 0.;
            plastic_out->cohesions[i] = 0;
            plastic_out->friction_angles[i] = 0;
            plastic_out->yielding[i] = plastic_yielding ? 1 : 0;
            plastic_out->creep_visco_ratio[i] = 0;
            plastic_out->volume_0[i] = 0;
			plastic_out->eta_dislocation[i] = 0;
			plastic_out->eta_diffusion[i] = 0;
			plastic_out->eta_for_viscous[i] = 0;
			plastic_out->plastic_stress[i] = 0;
			plastic_out->viscous_stress[i] = 0;
			plastic_out->yield_strength[i] = 0;
			double temp_plastic_stress = 0;
			double temp_viscous_stress = 0;
        // set to weakened values, or unweakened values when strain weakening is not used
		//20260127 蛇纹岩化相关：根据蛇纹岩化程度修改内摩擦角，蛇纹岩化程度大于0.1时，内摩擦角调整为输入值的0.1倍
		//20260225 蛇纹岩化区域参与塑性应变弱化
		  const double out_new_crust_phase = this->introspection().compositional_name_exists("new_crust") ? volume_fractions[this->introspection().compositional_index_for_name("new_crust") + 1] : 0.;
		  const double out_serp_value = this->introspection().compositional_name_exists("serp") ? in.composition[i][this->introspection().compositional_index_for_name("serp")] : 0.;
            for (unsigned int j=0; j < volume_fractions.size(); ++j)
              {
                // Calculate the strain weakening factors and weakened values
                const std::array<double, 3> weakening_factors = dike_mode ? dikemode_compute_strain_weakening_factors(j, in.composition[i]) : strain_rheology.compute_strain_weakening_factors(j, in.composition[i]);
				const double current_cohesions = volume_fractions[j] * (drucker_prager_parameters.cohesions[j] * weakening_factors[0]);
                plastic_out->cohesions[i]   += current_cohesions;
                // Also convert radians to degrees
				const double current_friction_angles = out_serp_value > 0.1 && out_new_crust_phase < 0.2 ? 180.0/numbers::PI * volume_fractions[j] * (drucker_prager_parameters.angles_internal_friction[j] * serp_weakening * weakening_factors[1]) : 180.0/numbers::PI * volume_fractions[j] * (drucker_prager_parameters.angles_internal_friction[j] * weakening_factors[1]);
                plastic_out->friction_angles[i] += current_friction_angles;
				const double current_plastic_stress = volume_fractions[j] * in.pressure[i] * std::cos(current_friction_angles * numbers::PI / 180.0) + current_cohesions * std::sin(current_friction_angles * numbers::PI / 180.0);
				temp_plastic_stress += current_plastic_stress;
				
				//20260203 可删除rift_grain_size
				double yield_stress_pre_grainsize_calculation = 0.;
				if(this->introspection().compositional_name_exists("rift_grain_size"))
				{
					const std::array<double, 3> weakening_factors_pre_grainsize_calculation = dike_mode ? dikemode_compute_strain_weakening_factors(j, in.composition[i]) : strain_rheology.compute_strain_weakening_factors(j, in.composition[i]);
			        const double current_cohesion_pre_grainsize_calculation = drucker_prager_parameters.cohesions[j] * weakening_factors_pre_grainsize_calculation[0];
                    const double current_friction_pre_grainsize_calculation = out_serp_value > 0.1 && out_new_crust_phase < 0.2 ? drucker_prager_parameters.angles_internal_friction[j] * serp_weakening * weakening_factors_pre_grainsize_calculation[1] : drucker_prager_parameters.angles_internal_friction[j] * weakening_factors_pre_grainsize_calculation[1];
			        const double pressure_for_plasticity_pre_grainsize_calculation = std::max(in.pressure[i],0.0);
			        yield_stress_pre_grainsize_calculation = drucker_prager_plasticity.compute_yield_stress(current_cohesion_pre_grainsize_calculation,
                                                                                           current_friction_pre_grainsize_calculation,
                                                                                           pressure_for_plasticity_pre_grainsize_calculation,
                                                                                           drucker_prager_parameters.max_yield_stress);
				}
				
				//20260203 可删除rift_grain_size
                const double out_last_grain_size = this->introspection().compositional_name_exists("rift_grain_size") ? std::max(std::min(std::abs(in.composition[i][this->introspection().compositional_index_for_name("rift_grain_size")]), diffusion_creep_grainsize.max_grain_size)
		                                               , diffusion_creep_grainsize.minimum_grain_size) : diffusion_creep_grainsize.initial_grain_size[j];
				//含水量最多降低到30 ppm
				const double out_water_content_for_viscosity = std::max(30. , std::abs(in.composition[i][this->introspection().compositional_index_for_name("water_content")]));
													   
                const double out_temperature_for_viscosity = in.temperature[i] + adiabatic_temperature_gradient_for_viscosity * in.pressure[i];
				//20240702
                const double out_current_grain_size = diffusion_creep_grainsize.compute_grainsize(out_edot_ii, in.pressure[i], out_temperature_for_viscosity, j,
                                                                                                  out_last_grain_size, out_posx, out_posy, out_water_content_for_viscosity, init_water_content[j], yield_stress_pre_grainsize_calculation, out_melt_fraction, phase_function_values,
                                                                                                  n_phases_per_composition);
			    const double out_viscosity_diffusion = diffusion_creep_grainsize.compute_viscosity(in.pressure[i], out_temperature_for_viscosity, j, out_current_grain_size,
                                                                                                   out_water_content_for_viscosity, init_water_content[j], phase_function_values,
                                                                                                   n_phases_per_composition);
				plastic_out->eta_diffusion[i] += out_viscosity_diffusion * out_temp_volume[j];
				
                const double out_viscosity_dislocation = dislocation_creep_grainsize.compute_viscosity(out_edot_ii, in.pressure[i], out_temperature_for_viscosity, out_water_content_for_viscosity, init_water_content[j], j,
                                                                                              phase_function_values,
                                                                                             n_phases_per_composition);
				plastic_out->eta_dislocation[i] += out_viscosity_dislocation * out_temp_volume[j];
				
				const double current_eta_for_viscous = out_viscosity_diffusion * out_viscosity_dislocation / (out_viscosity_diffusion + out_viscosity_dislocation) * out_temp_volume[j];
				plastic_out->eta_for_viscous[i] += current_eta_for_viscous;
				
				const double current_viscous_stress = 2 * out_edot_ii * current_eta_for_viscous * out_temp_volume[j];
				temp_viscous_stress += current_viscous_stress;
				
				
                plastic_out->creep_visco_ratio[i] += (1 / (1 + out_viscosity_diffusion / out_viscosity_dislocation)) * out_temp_volume[j];
              }
			plastic_out->plastic_stress[i] += temp_plastic_stress;
			plastic_out->viscous_stress[i] += temp_viscous_stress;
			plastic_out->yield_strength[i] = std::min(temp_viscous_stress, temp_plastic_stress);
			plastic_out->volume_0[i] = out_temp_volume[0];
          }
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  namespace Rheology \
  { \
    template class ViscoPlasticGrainSize<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}

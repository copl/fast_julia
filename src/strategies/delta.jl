function ipopt_style_inertia_correction!(newt::abstract_newton_direction, vars::class_variables, settings::class_settings)
	try

		MAX_IT = 50;
		j = 1;
		used_delta = 0.0;

		newt.delta = max(newt.delta, settings.delta_min );

		# try delta = 0.0
		old_delta = newt.delta;
		newt.delta = settings.delta_min;
		inertia = update_newton!(newt, vars, settings);

		if inertia == 1
			used_delta = newt.delta;
		else
			newt.delta = old_delta;

			for j = 2:MAX_IT
				inertia = update_newton!(newt, vars, settings)

				if inertia == 1
					used_delta = newt.delta;

					#toler = 1e-6
					#eig_min, eigen_vector, err = newt.minimum_delta(nlp_eval, newt.delta, toler);

					#println(full(newt.K))
					#suggested_delta = min(1.000,newt.delta * 0.999) #0.5*newt.delta + 0.5*smallest_delta;

					#if err < toler
					#newt.delta = max(-eig_min * 1.001, settings.delta_min)
					#else
					newt.delta = newt.delta * settings.delta_decrease;
					#end

					break
				elseif inertia == 0
					if newt.delta <= settings.delta_min;
						newt.delta = settings.delta_start;
					else
						newt.delta = newt.delta * settings.delta_increase;
					end
				elseif inertia == -1
					error("numerical stability issues when computing KKT system !!!")
				else
					error("inertia_corection")
				end
			end
		end

		if j == MAX_IT
			error("maximum iterations for inertia_corection reached")
		else
			return used_delta, j
		end


	catch e
		println("ERROR in ipopt_style_inertia_correction")
		throw(e)
	end
end

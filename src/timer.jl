type class_algorithm_timer

	start_times::Dict{String,Float64}
	total_times::Dict{String,Float64}

	start::Function
	stop::Function
	reset::Function
	print_stats::Function

	function class_algorithm_timer()
		this = new();
		this.start_times = Dict{String,Float64}()
		this.total_times = Dict{String,Float64}()

		function round_num(num,k=1)
			return round(num*10^k)/10^k
		end

		function show_percent(val,total)
			return round_num(100*val/total)
		end

		this.start = function(str)
			this.start_times[str] = time();
		end

		this.stop = function(str)
			if ~(str in keys(this.total_times))
				this.total_times[str] = 0.0;
			end
			time_passage = (time() - this.start_times[str]);
			this.total_times[str] = this.total_times[str] + time_passage;
		end

		this.reset = function()
			this.total_times = Dict{String,Float64}()
		end

		this.print_stats = function(standardize_by)
			println("======= Time statistics =======")
			if standardize_by != ""
				total = this.total_times[standardize_by]
				println( standardize_by, " = ", round_num(total, 3) )
				for label in keys(this.total_times)
					if label != standardize_by
						println( label, " = ", show_percent(this.total_times[label],total), " %" )
					end
				end
			else
				for label in keys(this.total_times)
					println( label, " = ", round_num(this.total_times[label], 3) )
				end
			end
		end

		return this;
	end
end

global GLOBAL_timer = class_algorithm_timer();

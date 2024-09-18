using Optim
using Distributions
using Plots

steps = 2000
cases = 2000
alphas = collect(0.01:0.01:0.1)

# run the optimization
function run_optimization()

    # initialize vectors
    beta_0 = zeros(length(alphas))
    beta_1 = zeros(length(alphas))

    # new_law = zeros(length(alphas))
    normal_law = zeros(length(alphas))

    for (i, alpha) in enumerate(alphas)
        normseq = rand(Normal(0,1), cases, steps) ./ sqrt(steps)
        for j in 2:steps
            normseq[:, j] =  normseq[:, j] + normseq[:, j-1]
        end
        function squaretoptim(beta)
            bar = max(beta[1],0) .* sqrt.(collect(1:1:steps) ./ steps)  .+ max(beta[2],0)
            success_trial = sum(sum(bar .- normseq .> 0, dims=2) .== steps)
            return 1000*abs(success_trial / cases - (1 - alpha)) + max(beta[2],0) + max(beta[1]*2/3,0)
        end
        res = optimize(squaretoptim, [0.766, 1.271], ParticleSwarm(), Optim.Options(f_calls_limit=1000, iterations=200))
        beta_0[i], beta_1[i] = Optim.minimizer(res)
        # new_law[i] = Optim.minimizer(res)[1]
        # new_law[i] = beta_0
        normal_law[i] = -quantile(Normal(), alpha / 2)
        println("Iteration $i completed")
    end

    # output 
    return beta_0, beta_1, normal_law, alphas
end

# # Run the optimization and time it
# @time beta_0, beta_1, normal_law, alphas = run_optimization()
# # new_law, normal_law, alphas = run_optimization()

# # Plot the results
# p = plot(alphas, [normal_law, beta_0, beta_1, beta_0 ./ normal_law],
# xlabel="α", ylabel="score", 
# label=["Pointwise" "Uniform (beta_0)" "Uniform (beta_1)" "Ratio"], 
# line = [:solid :dot :dash :dashdot], linewidth=2)

# # Display the plot
# display(p)

# # Save the plot
# savefig("optimization_plot_betas.png")
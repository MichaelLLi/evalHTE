using Optim
using Distributions

function run_optimization()
    steps = 2000
    cases = 2000
    alphas = collect(0.01:0.01:0.1)
    new_law = zeros(length(alphas))
    normal_law = zeros(length(alphas))
    
    for (i, alpha) in enumerate(alphas)
        normseq = rand(Normal(0,1), cases, steps) ./ sqrt(steps)
        for j in 2:steps
            normseq[:, j] =  normseq[:, j] + normseq[:, j-1]
        end
        
        function squaretoptim(beta)
            bar = max(beta[1],0) .* sqrt.(collect(1:1:steps) ./ steps)  .+ max(beta[2],0)
            success_trial = 0
            for i in 1:cases
                success_trial += Int(sum(bar .- normseq[i,:] .> 0) == steps)
            end
            return 1000*abs(success_trial / cases - (1 - alpha)) + max(beta[2],0) + max(beta[1]*2/3,0)
        end
        
        res = optimize(squaretoptim, [0.766, 1.271], ParticleSwarm(), Optim.Options(f_calls_limit=1000, iterations=200))
        new_law[i] = Optim.minimizer(res)[1]
        normal_law[i] = -quantile(Normal(), alpha / 2)
        print(i)
    end
    
    return (alphas, new_law, normal_law)
end
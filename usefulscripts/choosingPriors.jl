
using Optim, Distributions

function chooseNormalParams(desired_mean, desired_quantile; makeplot=false)

    function moments(sd)
        current_quantile = quantile(Truncated(Normal(desired_mean, sd), 0, 1), 0.975)
        return(current_quantile)
    end

    function objective(sd)
        current_quantile = moments(sd[1])
        square_error = (current_quantile - desired_quantile)^2
        return(square_error)
    end

    init_guess = [0.01]
    opt = Optim.optimize(objective, [1e-20], [Inf], init_guess, Fminbox(GradientDescent()))

    sd_chosen = opt.minimizer[1]
    moms = moments(sd_chosen)
    
    println("Optimal parameter: sd = ", sd_chosen)
    println("Has upper quantile " * string(moms))

    if makeplot
        x = 0.0001:0.0001:0.02
        y = pdf.(Truncated(Normal(desired_mean, sd_chosen), 0, Inf),x)
        display(plot(x,y))
    end
    
    return(sd_chosen)
    
end


function chooseLogisticParams(desired_mean, desired_quantile; makeplot=false)

    function moments(sd)
        current_quantile = quantile(Truncated(Logistic(desired_mean, sd), 0, 1), 0.975)
        return(current_quantile)
    end

    function objective(sd)
        current_quantile = moments(sd[1])
        square_error = (current_quantile - desired_quantile)^2
        return(square_error)
    end

    init_guess = [0.01]
    opt = Optim.optimize(objective, [1e-20], [Inf], init_guess, Fminbox(GradientDescent()))

    sd_chosen = opt.minimizer[1]
    moms = moments(sd_chosen)
    
    println("Optimal parameter: sd = ", sd_chosen)
    println("Has upper quantile " * string(moms))

    if makeplot
        x = 0.0001:0.0001:0.02
        y = pdf.(Truncated(Logistic(desired_mean, sd_chosen), 0, Inf),x)
        display(plot(x,y))
    end
    
    return(sd_chosen)
    
end



function chooseInvGamParams(desired_mean, desired_quantile; makeplot=false)
    
    # Function to compute current mean and quantile for given parameters
    function moments(params)
        alpha, beta = params
        current_mean = beta / (alpha - 1)
        current_quantile = quantile(InverseGamma(alpha, beta), 0.95)
        return [current_mean, current_quantile]
    end
    
    # Objective function to minimize: squared distance to desired mean and quantile
    function objective(α)


        β = (α[1]-1)*desired_mean

        (current_mean, upper_quant) = moments([α[1], β])
        square_error = (upper_quant - desired_quantile)^2

        return(square_error)

    end
    
    # Initial guess for alpha
    init_guess = [100]
    opt = Optim.optimize(objective, [1.0+1e-20], [Inf], [2.0], Fminbox(GradientDescent()))

    α = opt.minimizer[1]
    β = (α - 1) * desired_mean
    moms = moments([α, β])
    
    println("Optimal parameters: alpha = ", α, ", beta = ", β)
    println("Have mean " * string(moms[1]) * " and upper quant " * string(moms[2]))

    if makeplot
        x = 0.0001:0.0001:0.02
        y = pdf.(InverseGamma(α, β), x)
        display(plot(x,y))
    end
    
    return(α, β)
    
end
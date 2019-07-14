module AutoCorrelation

using Distributions

export pre_kenny, kenny
# Function kenny takes input u: time lag, N: Number of steps(jumps), x_0: Starting Point, alpha: increase rate (constant), 
# beta: slope of linear decrease rate 
# and returns the autocorrelation value for time lag u

"""
`x_0`: Starting Point, `alpha`: increase rate (constant), 
`beta`: slope of linear decrease rate 

External links
* [Gollespie algorithm](https://en.wikipedia.org/wiki/Gillespie_algorithm)
"""
function pre_kenny(N, x_0, alpha, beta)
    x =zeros(Int, N+1)
    t =zeros(N+1)

    # Initialize x and t, first x is x_0 and last t is 0
    x[1] = x_0
    t[N] = 0

    # Gillespie algorithm
    for i in 1:N
        # Specify Rates, rate2 update with x for every step
        rates1 = alpha
        rates2 = beta*x[i]
        kT = rates1+rates2

        # Generate waiting time, uses package distribution, may be improved by using only uniform
        dist = Exponential(1/kT)
        t[i] =rand(dist)

        # Decide which step to take and record
        u_event =rand()
        if (u_event < rates1/kT)
            x[i+1] = x[i]+1
        else
            x[i+1] = x[i]-1
        end
    end
    return x, t
end

"""
    kenny(u, x, t)

`u`: time lag, `x`, `t` are sequences prepared by pre_kenny

Returns the autocorrelation value for time lag u

"""
function kenny(u::Union{Int, Float64}, x, t)
    N = x |> length

    # Calculate the time spend in each step
    min_x = minimum(x)
    max_x = maximum(x)
    t_sum = zeros(max_x+1)

    for num_x in 0:max_x
        for j in 1:N
            if(x[j]==num_x)
                # index starts with one, x is shifted
                t_sum[num_x+1]=t_sum[num_x+1]+t[j]
            end
        end
    end

    #@show t_sum

    # Mean of x given by probability of being in state x times x
    xmean = sum([i*t_sum[i+1] for i in 1:max_x]) / sum(t_sum)

    # Second moment of x
    xsq = sum([(i^2)*t_sum[i+1] for i in 1:max_x]) / sum(t_sum)

    # Variance of x
    xvar=xsq-xmean^2

    #@show xmean
    #@show xsq
    #@show xvar

    # The covariance stable equation for 1 variable
    #@show xvar/xmean/xmean-1/xmean

    # Calculate the Autocorrelation

    # Initialize, expected is expected value of X_t and X_{t+u}
    expected=0
    # Keep track of time
    cum_t = cumsum(t)
    t_sum_2 = zeros(max_x+1,max_x+1)

    # Using Multithreads
    # Probability that x_t=i and x_{t+u}=k
    # time spent that is true for x_t=i and x_{t+u}=k, notice indices are shifted since i,k can be zero
    # Check each step
    @inbounds for i in min_x:max_x
        for k in min_x:max_x
            for j in 1:N
                @simd for p in 1:N
                    # Discuss two end points of each interval
                    time_jp = 0
                    # None edge-case, indicator shows if there time is 0
                    if(j>1 && p>1)
                        indicator=( (cum_t[j-1]+u)<=(cum_t[p-1])<(cum_t[j]+u) || (cum_t[j-1]+u)<(cum_t[p])<=(cum_t[j]+u) || (cum_t[p-1])<=(cum_t[j-1]+u)<(cum_t[p]) || (cum_t[p-1])<(cum_t[j]+u)<=(cum_t[p]) )
                    end

                    if(x[j]==i && x[p]==k && j>1 && p>1 && indicator)
                        # Four cases of overlap of two intervals

                        # 1. t[p] is a subinterval of t[j]
                        if((cum_t[j-1]+u)<=(cum_t[p-1])<(cum_t[j]+u) && (cum_t[j-1]+u)<(cum_t[p])<=(cum_t[j]+u))
                            time_jp=t[p]

                            # 2. t[p] starts in between t[j] but ends after t[j] 
                        elseif((cum_t[j-1]+u)<=(cum_t[p-1])<(cum_t[j]+u))
                            time_jp=(cum_t[j]+u)-(cum_t[p-1])
                            #time_jp=min(t[p],t[j],u)

                            # 3. t[j] is a subinterval of t[p]
                        elseif((cum_t[p-1])<=(cum_t[j-1]+u)<(cum_t[p]) && (cum_t[p-1])<(cum_t[j]+u)<=(cum_t[p]) )
                            time_jp=t[j]
                            # 4. t[p] starts before t[j] but ends in-between t[j]
                        else
                            time_jp=(cum_t[p])-(cum_t[j-1]+u)
                        end
                        #min(t[p],t[j],u)

                        # edge case, first interval equal
                    elseif(j==1 && p==1 && (i==k) && (u<t[1]) && x[j]==i && x[p]==k)
                        time_jp=t[1]-u

                        # edge case,t[j] is first interval
                    elseif(j==1 && p>1 && x[j]==i && x[p]==k)
                        indicator=( (u)<=(cum_t[p-1])<(cum_t[j]+u) || (u)<(cum_t[p])<=(cum_t[j]+u) || (cum_t[p-1])<=(u)<(cum_t[p]) || (cum_t[p-1])<(cum_t[j]+u)<=(cum_t[p]) )
                        if(x[j]==i && x[p]==k && indicator)
                            #
                            if((u)<=(cum_t[p-1])<(cum_t[j]+u) && (u)<(cum_t[p])<=(cum_t[j]+u))
                                time_jp=t[p]
                            elseif((u)<=(cum_t[p-1])<(cum_t[j]+u))
                                time_jp=(cum_t[j]+u)-(cum_t[p-1])
                                #time_jp=min(t[p],t[j],u)
                            elseif((cum_t[p-1])<=(u)<(cum_t[p]) && (cum_t[p-1])<(cum_t[j]+u)<=(cum_t[p]) )
                                time_jp=t[j]
                            else
                                time_jp=(cum_t[p])-(u)
                            end
                        end
                        # If t[j] is not the first interval but t[p] is, time must be zero
                    else
                        #time_jp=0
                    end
                    #t_sum_2[i,k]=t_sum_2[i,k]+time_jp
                    t_sum_2[i+1, k+1] += time_jp
                end
            end
        end
    end


    # expected value, i,k times probability
    for i in 1:max_x+1
        @simd for k in 1:max_x+1
            expected+=(i-1)*(k-1)*t_sum_2[i,k] #sum all elems of matrix
        end
    end
    expected /= sum(t_sum_2)
    # return Autocorrelation (autocovariance if not divided by xvar)
    return (expected-xmean^2)/xvar
end

#end of module
end


using Distributions
#using BenchmarkTools
using PyPlot

# Function kenny takes input u: time lag, N: Number of steps(jumps), x_0: Starting Point, alpha: increase rate (constant), 
# beta: slope of linear decrease rate 
# and returns the autocorrelation value for time lag u

# Function declaration and argument type
function kenny(u::Union{Int, Float64}, N::Int, x_0::Int,
    alpha::Union{Int, Float64}, beta::Union{Int, Float64})

# Initialize x and t, first x is x_0 and last t is 0
    x =zeros(Int, N+1)
    t =zeros(N+1)

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

    # Calculate the time spend in each step
    t_sum = zeros(maximum(x)+1)

    for num_x in 0:maximum(x)
        for j in 1:N
            if(x[j]==num_x)
                # index starts with one, x is shifted
                t_sum[num_x+1]=t_sum[num_x+1]+t[j]
            end
        end
    end

    #@show t_sum

    # Mean of x given by probability of being in state x times x
    xmean = sum([i*t_sum[i+1] for i in 1:maximum(x)]) / sum(t_sum)

    # Second moment of x
    xsq = sum([(i^2)*t_sum[i+1] for i in 1:maximum(x)]) / sum(t_sum)

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
    t_sum_2 = zeros(maximum(x)+1,maximum(x)+1)
#    for i in minimum(x):(maximum(x))
#        for k in minimum(x):(maximum(x))
            #time_jp=np.zeros((N,N))
    
# Using Multithreads
    # Probability that x_t=i and x_{t+u}=k
    Threads.@threads for i in minimum(x):(maximum(x))
        Threads.@threads for k in minimum(x):(maximum(x))
            # time spent that is true for x_t=i and x_{t+u}=k, notice indices are shifted since i,k can be zero
            time_jp_matrix = zeros(maximum(x)+1,maximum(x)+1)
            # Check each step
            for j in 1:N
                for p in 1:N
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
                    time_jp_matrix[i+1, k+1] += time_jp
                end
            end
            t_sum_2+=time_jp_matrix
        end
    end


# expected value, i,k times probability
    for i in 1:maximum(x)+1
        for k in 1:maximum(x)+1
            expected
            expected=expected+(i-1)*(k-1)*t_sum_2[i,k] #sum all elems of matrix
        end
    end
    expected /= sum(t_sum_2)
    # return Autocorrelation (autocovariance if not divided by xvar)
    return (expected-xmean^2)/xvar
end







# Test when u=0, should always be equal to xvar and autocorrelation is 1
@time kenny(0, 200, 10, 1, 0.1)

# Define the scope of u from 0 to 2 with stepsize 0.025
us=[0:0.025:2...];

# Define function that parallel u
function splits(us,N)
    values = zeros(length(us));
    Threads.@threads for i in 1:length(us)
        values[i]=kenny(us[i],N,10,1,0.1)
    end
    values
end

# Run u without parallel
@time result=[kenny(u,2500,10,1,1) for u in us]

# Plot the result using Gadfly
using Gadfly
simu=layer(x=us,y=result,Geom.line,Theme(default_color=colorant"green"))
tru_exp=layer(x->exp(-1x),0,2,Theme(default_color=colorant"orange"))
Gadfly.plot(simu,tru_exp,Guide.manual_color_key("Legend",["Simulation","True Exp"],["green","orange"]),Guide.xlabel("time u"), Guide.ylabel("Autocorrelation"))

# Run u with parallel, need to re-run the previous step
@time result=splits(us,2500)

# Test 10000 steps for 1 u-value
@time kenny(0.5, 10000, 10, 1, 0.1)

# Plot the result using Gadfly
using Gadfly
simu=layer(x=us,y=result,Geom.line,Theme(default_color=colorant"green"))
tru_exp=layer(x->exp(-1x),0,2,Theme(default_color=colorant"orange"))
Gadfly.plot(simu,tru_exp,Guide.manual_color_key("Legend",["Simulation","True Exp"],["green","orange"]),Guide.xlabel("time u"), Guide.ylabel("Autocorrelation"))













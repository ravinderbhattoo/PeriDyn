export minimize!, quasi_static!

function minimize!(env::GeneralEnv, step_size::Float64; max_iter::Int64=50, x_tol::Float64=1.0e-6, f_tol::Float64=1.0e-6)
    for bc in env.boundary_conditions
        apply_bc!(env, bc, :position)
        apply_bc!(env, bc, :velocity)
    end
    update_acc!(env)
    
    dt = env.dt
    ps = env.material_blocks[1].general.particle_size    
    
    k = 1.0e6
    x_tol = ps / k
    f_tol = ps / k / dt^2
    
    println("Bypassing given tolerance values.\nusing x_tol = $(x_tol) and f_tol = $(f_tol).")
    
    opt = Flux.ADAM(step_size)
    mask = false
    for bc in env.boundary_conditions
        if ~(bc.onlyatstart)
            mask = mask .| bc.bool
        end
    end
    mask = .!(mask)
    println("Minimizing...")
    iter = ProgressBar(1:max_iter)
    for i in iter
        f_tol_ = maximum(abs.(env.f[:, mask]))
        
        Δy = Flux.Optimise.apply!(opt, env.y, env.f)
        env.y .+= Δy
        
        x_tol_ = maximum(abs.(Δy[:, mask]))
        
        for bc in env.boundary_conditions
            if ~(bc.onlyatstart)
                apply_bc!(env, bc, :position)
                apply_bc!(env, bc, :velocity)
            end
        end
        update_acc!(env)
        
        if (f_tol_ < f_tol) #|| (x_tol_ < x_tol)
            println("Tolerance reached, iter $i. x_tol $x_tol_ <= $x_tol or f_tol $f_tol_ <= $f_tol")
            break
        end
        if i==max_iter
            println("Maximum iteration, iter $i ($max_iter) reached. x_tol $x_tol_ !<= $x_tol and f_tol $f_tol_ !<= $f_tol")
        end
        set_postfix(iter, F_tol=f_tol_)
    end
    for bc in env.boundary_conditions
        check!(bc, env)
    end
    env.time_step += 1
end


"""
    quasi_static!(envs::Any,N::Int64,step_size::Float64; max_iter::Int64=100, min_step_tol_per::Float64=0.5, filewrite_freq::Int64=10, neigh_update_freq::Int64=50, file_prefix::String="datafile",start_at::Int64=0,write_from::Int=0)

Implements quasi static simulation using minimize for each step.
"""
function quasi_static!(envs::Any,N::Int64,step_size::Float64; max_iter::Int64=100, filewrite_freq::Int64=10, neigh_update_freq::Int64=1, out_dir::String="datafile",start_at::Int64=0,write_from::Int=0, ext::Symbol=:jld)
    mkpath(out_dir)
    _print_data_file!(step) = print_data_file!(envs, out_dir, step; ext=ext)

    _print_data_file!(0)
    update_neighs!(envs)
    
    _print_data_file!(0+write_from)
 
    
    for id in 1:size(envs,1)
        if envs[id].state==2
            apply_bc_at0(envs[id], start_at)
        end
    end
    N = N + start_at
    for i in (1+start_at):N
        for id in 1:size(envs,1)
            if envs[id].state==2
                fill!(envs[id].v,0.0)
                minimize!(envs[id],step_size; max_iter=max_iter)
            end
            if envs[id].Collect! !== nothing
                envs[id].Collect!(envs[id],i)
            end
        end
        if i%filewrite_freq==0.0 || (i==1 && write_from==0)
            _print_data_file!(i+write_from)
        end
        
        if i%neigh_update_freq==0.0
            update_neighs!(envs)
        end
        println(round(i/N*100, digits=3),"%")
    end
end

SOLVERS[:qs] = quasi_static!


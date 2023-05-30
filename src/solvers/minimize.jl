export minimize!, apply_solver!, QSDrag

"""
    QSDrag

Quasi-static solver struct.

# Fields
- `step_size`: Float64, the step size.
- `drag`: Float64, the drag coefficient.
- `max_iter`: Int64, the maximum number of iterations. Default is 1000.
- `x_tol`: Float64, the tolerance for position. Default is 1.0e-3.
- `f_tol`: Float64, the tolerance for force. Default is 1.0e-3.

# Example
```
solver = QSDrag(1.0e-3, 1.0e-3)
```

# See also
- `QSDrag`
- `apply_solver!`
- `minimize!`
"""
struct QSDrag <: QuasiStaticSolver
    step_size
    drag
    @qssolver_gf
end


"""
    QSDrag(step_size, drag; max_iter=1000, x_tol=1.0e-3, f_tol=1.0e-3)

Create a QSDrag solver.

# Arguments
- `step_size`: Float64, the step size.
- `drag`: Float64, the drag coefficient.

# Keyword Arguments
- `max_iter`: Int64, the maximum number of iterations. Default is 1000.
- `x_tol`: Float64, the tolerance for position. Default is 1.0e-3.
- `f_tol`: Float64, the tolerance for force. Default is 1.0e-3.

# Example
```
solver = QSDrag(1.0e-3, 1.0e-3)
```
"""
function QSDrag(step_size, drag; max_iter=1000, x_tol=1.0e-3, f_tol=1.0e-3)
    deviceconvert(QSDrag(step_size, drag, max_iter, x_tol, f_tol))
end


"""
    apply_solver!(env, solver::QSDrag)

Apply the QSDrag solver to the environment.

# Arguments
- `env`: GeneralEnv, the environment.
- `solver`: QSDrag, the QSDrag solver.

# Example
```
solver = QSDrag(1.0e-3, 1.0e-3)
apply_solver!(env, solver)
```
"""
function apply_solver!(env, solver::QSDrag)
    minimize!(env, solver)
end



"""
    minimize!(env, solver::QSDrag)

Minimize the environment using the QSDrag solver.

# Arguments
- `env`: GeneralEnv, the environment.
- `solver`: QSDrag, the QSDrag solver.

# Example
```
solver = QSDrag(1.0e-3, 1.0e-3)
minimize!(env, solver)
```
"""
function minimize!(env::GeneralEnv, solver::QSDrag)
    step_size = solver.step_size
    lambda = solver.drag
    max_iter = solver.max_iter
    x_tol = solver.x_tol
    f_tol = solver.f_tol

    for bc in env.boundary_conditions
        apply_bc!(env, bc, :position)
        apply_bc!(env, bc, :velocity)
    end

    dt = step_size
    ps = env.material_blocks[1].general.particle_size

    x_tol = ps * 1.0e-4
    f_tol = x_tol / dt^2

    log_impinfo("Bypassing given tolerance values.\nusing x_tol = $(x_tol) and f_tol = $(f_tol).")

    function clip(x; a=1.0)
        max.(min.(x, a), -a)
    end

    function clipx(x)
        clip(x; a=ps/10)
    end

    function clipv(v)
        clip(v; a=ps/10/dt)
    end


    mask = false
    for bc in env.boundary_conditions
        if ~(bc.onlyatstart)
            mask = mask .| bc.bool
        end
    end

    mask = .!(mask)
    opt = zero(env.f[:, mask])
    k = lambda / 10

    function apply_!(opt, f, step_size)
        force = f .- lambda * opt .* (1.0 .+ k * abs.(opt))
        Δy = 0.5 * force * step_size^2 .+ opt * step_size
        opt .+= step_size .* force
        Δy, opt = clipx(Δy), clipv(opt)
        return Δy
    end

    log_info("Minimizing...")

    # function obj_fn(y)
    #     env.y[:, mask] .= y
    #     update_acc!(env)
    #     sum(env.f[:, mask].^2)
    # end

    # function grad_fn!(opt, y)
    #     env.y[:, mask] .= y
    #     update_acc!(env)
    #     opt .= env.f[:, mask]
    # end

    # function fg!(F, G, y)
    #     env.y[:, mask] .= y
    #     update_acc!(env)
    #     if G != nothing
    #         G .= env.f[:, mask]
    #     end
    #     if F != nothing
    #         # value = ... code to compute objective function
    #         return sum(env.f[:, mask].^2)
    #     end
    # end

    # result = optimize(Optim.only_fg!(fg!), env.y[:, mask], GradientDescent(),
    #                 Optim.Options(g_tol = f_tol,
    #                 iterations = max_iter,
    #                 store_trace = false,
    #                 show_trace = true,
    #                 show_every = 1)
    #                 )

    # env.y[:, mask] .= Optim.minimizer(result)

    gs = @view env.f[:, mask]
    f_tol_ = maximum(abs.(gs))

    iter = ProgressBar(1:max_iter)
    for i in iter

        update_acc!(env)
        gs = @view env.f[:, mask]
        Δy = apply_!(opt, gs, step_size)
        env.y[:, mask] .+= Δy

        f_tol_, x_tol_ = maximum(abs.(gs)), maximum(abs.(Δy))

        # for bc in env.boundary_conditions
        #     if ~(bc.onlyatstart)
        #         apply_bc!(env, bc, :position)
        #         apply_bc!(env, bc, :velocity)
        #     end
        # end

        if ((f_tol_ < f_tol) && (x_tol_ < x_tol))  #&& i > 100 #|| (x_tol_ < x_tol) ) && i > 100
            log_info("Tolerance reached, iter $i. x_tol $x_tol_ <= $x_tol or f_tol $f_tol_ <= $f_tol")
            break
        end
        if i==max_iter
            log_impinfo("Maximum iteration, iter $i ($max_iter) reached. x_tol $x_tol_ !<= $x_tol and f_tol $f_tol_ !<= $f_tol")
        end
        set_postfix(iter, Iter=i, X_tol=round(x_tol_, digits=6), F_tol=round(f_tol_, digits=6), MaxF=round(maximum(env.f), digits=6))
    end

    for bc in env.boundary_conditions
        check!(bc, env)
    end

    env.time_step += 1
end

SOLVERS[:qs] = QSDrag
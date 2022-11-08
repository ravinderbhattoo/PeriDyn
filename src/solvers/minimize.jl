export minimize!, apply_solver!, QSDrag

struct QSDrag <: QuasiStaticSolver
    step_size
    drag
    @qssolver_gf
end


function QSDrag(step_size, drag; max_iter=1000, x_tol=1.0e-3, f_tol=1.0e-3)
    QSDrag(step_size, drag, max_iter, x_tol, f_tol)
end


function apply_solver!(env, solver::QSDrag)
    minimize!(env, solver)
end

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

    x_tol_ = max(f_tol * dt^2, ps * 10e-3)
    x_tol = min(x_tol, x_tol_)

    function clip(x; a=1.0)
        max.(min.(x, a), -a)
    end

    function clipx(x)
        clip(x; a=ps/10)
    end

    function clipv(v)
        clip(v; a=ps/10/dt)
    end

    println("Bypassing given tolerance values.\nusing x_tol = $(x_tol) and f_tol = $(f_tol).")

    mask = false
    for bc in env.boundary_conditions
        if ~(bc.onlyatstart)
            mask = mask .| bc.bool
        end
    end

    mask = .!(mask)

    opt = zeros(eltype(env.f), size(env.f[:, mask]))

    k = lambda/10

    function apply_!(opt, f)
        force = f .- lambda * opt .* (1.0 .+ k*abs.(opt))
        Δy = opt * step_size .+ 0.5*force*step_size^2
        opt .+= step_size .* force
        Δy = clipx(Δy)
        opt = clipv(opt)
        return Δy
    end



    # function obj_fn(y)
    #     env.y[:, mask] .= y
    #     update_acc!(env)
    #     sum(env.f[:, mask].^2)
    # end

    println("Minimizing...")

    # optimize(obj_fn, env.y[:, mask], LBFGS())


    iter = ProgressBar(1:max_iter)

    for i in iter

        update_acc!(env)

        gs = @view env.f[:, mask]
        Δy = apply_!(opt, gs)
        env.y[:, mask] .+= Δy

        f_tol_ = maximum(abs.(gs))
        x_tol_ = maximum(abs.(Δy))

        # for bc in env.boundary_conditions
        #     if ~(bc.onlyatstart)
        #         apply_bc!(env, bc, :position)
        #         apply_bc!(env, bc, :velocity)
        #     end
        # end

        if ((f_tol_ < f_tol))  && i > 100 #|| (x_tol_ < x_tol) ) && i > 100
            println("Tolerance reached, iter $i. x_tol $x_tol_ <= $x_tol or f_tol $f_tol_ <= $f_tol")
            break
        end
        if i==max_iter
            println("Maximum iteration, iter $i ($max_iter) reached. x_tol $x_tol_ !<= $x_tol and f_tol $f_tol_ !<= $f_tol")
        end
        set_postfix(iter, X_tol=x_tol_, F_tol=f_tol_, MaxF=maximum(env.f))
    end

    for bc in env.boundary_conditions
        check!(bc, env)
    end

    env.time_step += 1
end

SOLVERS[:qs] = QSDrag


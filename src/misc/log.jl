export log_data, log_warning, log_info, log_impinfo

const PRINT_DEFAULT_DATA = Ref{Bool}(true)
const LOGLEVEL = Ref(2)
const GOBBLE_DATA_KEYS = Ref(true)
const WRITE_ONLY_HEADERS = Ref(false)
LOGFILE_HANDLE = ""

"""
    write_headers(envs)

Write the column headers of the log file.

# Arguments
- `envs`: Vector{Environment}, the simulation environments.

"""
function write_headers(envs)
    WRITE_ONLY_HEADERS[] = true
    log_data(Step="")
    for env in envs
        env.time_step += 1
        printdata(env)
        env.time_step -= 1
    end
    println()
    write(LOGFILE_HANDLE, "\n")
    WRITE_ONLY_HEADERS[] = false
end

"""
    log_data(; kwargs...)

Log the data to the log file (and print).

# Arguments
- `kwargs`: Keyword arguments, the data to be logged.

"""
function log_data(; kwargs...)
    if LOGLEVEL[] >= -1
        if ~WRITE_ONLY_HEADERS[]
            if GOBBLE_DATA_KEYS[]
                for (key, val) in kwargs
                    str_ = "$val\t"
                    print(str_)
                    write(LOGFILE_HANDLE, str_)
                end
            else
                for (key, val) in kwargs
                    print( @magenta(@bold "$key:"), " $val,")
                    write(LOGFILE_HANDLE, "$val\t")
                end
            end
        else
            for (key, val) in kwargs
                str_ = "$key\t"
                write(LOGFILE_HANDLE, str_)
                if GOBBLE_DATA_KEYS[]
                    print(@magenta(@bold(str_)))
                end
            end
        end
    end
end

macro swarning(x)
    if LOGLEVEL[] >= 1
        quote
            println(@red @bold @underline "WARNING: " * string($(esc(x))) * "\n")
        end
    else
    end
end

"""
    log_warning(x)

Log a warning message. Prints only if the log level is 1 or higher.

# Arguments
- `x`: String, the warning message.

"""
log_warning(x) = @swarning(x)

macro simpinfo(x)
    if LOGLEVEL[] >= 2
        quote
            print(@yellow @bold "INFO: ")
            println($(esc(x)))
        end
    else
    end
end

"""
    log_impinfo(x)

Log an important information message. Prints only if the log level is 2 or higher.

# Arguments
- `x`: String, the information message.

"""
log_impinfo(x) = @simpinfo(x)

macro sinfo(x)
    if LOGLEVEL[] >= 3
        quote
            println($(esc(x)))
        end
    else
    end
end


"""
    log_info(x)

Log an information message.

# Arguments
- `x`: String, the information message. Prints only if the log level is 3 or higher.

"""
log_info(x) = @sinfo(x)

macro sdetail(x)
    if LOGLEVEL[] >= 4
        quote
            println(@green $(esc(x)))
        end
    else
    end
end

"""
    log_detail(x)

Log a detailed information message. Prints only if the log level is 4 or higher.

# Arguments
- `x`: String, the information message.

"""
log_detail(x) = @sdetail(x)

"""
    set_loglevel(x)

Set the log level.

# Arguments
- `x`: Int64, the log level.

"""
function set_loglevel(x::Int64)
    LOGLEVEL[] = x
    @simpinfo "Log level: $(LOGLEVEL[])"
    refresh()
end



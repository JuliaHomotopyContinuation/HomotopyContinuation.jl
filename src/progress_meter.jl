# This is a slightly modified version of ProgressMeter.jl
#
# The MIT License (MIT)
# Copyright (c) 2013 Timothy E. Holy
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

module ProgressMeter

using Printf: @sprintf

export Progress, ProgressUnknown,next!, update!, cancel, finish!

"""
`ProgressMeter` contains a suite of utilities for displaying progress
in long-running computations. The major functions/types in this module
are:

- `@showprogress`: an easy interface for straightforward situations
- `Progress`: an object for managing progress updates with a predictable number of iterations
- `ProgressThresh`: an object for managing progress updates where termination is governed by a threshold
- `next!` and `update!`: report that progress has been made
- `cancel` and `finish!`: early termination
"""
ProgressMeter

abstract type AbstractProgress end

"""
Holds the five characters that will be used to generate the progress bar.
"""
mutable struct BarGlyphs
    leftend::Char
    fill::Char
    front::Union{Vector{Char}, Char}
    empty::Char
    rightend::Char
end
"""
String constructor for BarGlyphs - will split the string into 5 chars
"""
function BarGlyphs(s::AbstractString)
    glyphs = (s...,)
    if !isa(glyphs, NTuple{5,Char})
        error("""
            Invalid string in BarGlyphs constructor.
            You supplied "$s".
            Note: string argument must be exactly 5 characters long, e.g. "[=> ]".
        """)
    end
    return BarGlyphs(glyphs...)
end

"""
`prog = Progress(n; dt=0.1, desc="Progress: ", color=:green,
output=stdout, barlen=tty_width(desc))` creates a progress meter for a
task with `n` iterations or stages. Output will be generated at
intervals at least `dt` seconds apart, and perhaps longer if each
iteration takes longer than `dt`. `desc` is a description of
the current task.
"""
mutable struct Progress <: AbstractProgress
    n::Int
    dt::Float64
    counter::Int
    tfirst::Float64
    tlast::Float64
    printed::Bool           # true if we have issued at least one status update
    desc::AbstractString    # prefix to the percentage, e.g.  "Computing..."
    barlen::Int             # progress bar size (default is available terminal width)
    barglyphs::BarGlyphs    # the characters to be used in the bar
    color::Symbol           # default to green
    output::IO              # output stream into which the progress is written
    offset::Int             # position offset of progress bar (default is 0)
    numprintedvalues::Int   # num values printed below progress in last iteration
    clear_output_ijulia::Bool  # flush current display in IJulia possibly overwriting other printed things

    function Progress(n::Integer;
                      dt::Real=0.1,
                      desc::AbstractString="Progress: ",
                      color::Symbol=:green,
                      output::IO=stdout,
                      delay::Real=0.0,
                      barlen::Integer=tty_width(desc),
                      barglyphs::BarGlyphs=BarGlyphs('|','█', Sys.iswindows() ? '█' : ['▏','▎','▍','▌','▋','▊','▉'],' ','|',),
                      offset::Int=0,
                      clear_output_ijulia::Bool=false)
        counter = 0
        tfirst = time()
        tlast = tfirst + delay
        printed = false
        new(n, dt, counter, tfirst, tlast, printed, desc, barlen, barglyphs, color, output, offset, 0, clear_output_ijulia)
    end
end

Progress(n::Integer, dt::Real, desc::AbstractString="Progress: ",
         barlen::Integer=tty_width(desc), color::Symbol=:green, output::IO=stdout; kwargs...) =
    Progress(n; dt=dt, desc=desc, barlen=barlen, color=color, output=output, kwargs...)

Progress(n::Integer, desc::AbstractString, offset::Integer=0) = Progress(n, desc=desc, offset=offset)


"""
`prog = ProgressUnknown(; dt=0.1, desc="Progress: ",
color=:green, output=stdout)` creates a progress meter for a task
which has a non-deterministic termination criterion.
Output will be generated at intervals at least `dt` seconds
apart, and perhaps longer if each iteration takes longer than
`dt`. `desc` is a description of the current task.
"""
mutable struct ProgressUnknown <: AbstractProgress
    done::Bool
    dt::Float64
    counter::Int
    triggered::Bool
    tfirst::Float64
    tlast::Float64
    printed::Bool        # true if we have issued at least one status update
    desc::AbstractString # prefix to the percentage, e.g.  "Computing..."
    color::Symbol        # default to green
    output::IO           # output stream into which the progress is written
    numprintedvalues::Int   # num values printed below progress in last iteration
    clear_output_ijulia::Bool      # flush current display in IJulia possibly overwriting other printed things
end

function ProgressUnknown(;dt::Real=0.1, desc::AbstractString="Progress: ",
                          color::Symbol=:green,
                          output::IO=stdout,
                          delay::Real=0.0,
                          clear_output_ijulia::Bool=false)
    tfirst = time()
    tlast = tfirst + delay
    printed = false
    ProgressUnknown(false, dt, 0, false, tfirst, tlast, printed, desc, color, output, 0, clear_output_ijulia)
end

ProgressUnknown(dt::Real, desc::AbstractString="Progress: ",
         color::Symbol=:green, output::IO=stdout; kwargs...) =
    ProgressUnknown(dt=dt, desc=desc, color=color, output=output; kwargs...)

ProgressUnknown(desc::AbstractString; kwargs...) = ProgressUnknown(;desc=desc, kwargs...)

#...length of percentage and ETA string with days is 29 characters
tty_width(desc) = max(0, displaysize(stdout)[2] - (length(desc) + 29))

flush_display(p) = p.clear_output_ijulia && isdefined(Main, :IJulia)

# update progress display
function updateProgress!(p::Progress; showvalues = Any[], valuecolor = :blue, offset::Integer = p.offset, keep = (offset == 0))
    p.offset = offset
    t = time()
    if p.counter >= p.n
        if p.counter == p.n && p.printed
            percentage_complete = 100.0 * p.counter / p.n
            bar = barstring(p.barlen, percentage_complete, barglyphs=p.barglyphs)
            dur = durationstring(t-p.tfirst)
            msg = @sprintf "%s%3u%%%s Time: %s" p.desc round(Int, percentage_complete) bar dur
            !flush_display(p) && print(p.output, "\n" ^ (p.offset + p.numprintedvalues))
            move_cursor_up_while_clearing_lines(p.output, p.numprintedvalues, flush_display(p))
            printover(p.output, msg, p.color, flush_display(p))
            printvalues!(p, showvalues; color = valuecolor)
            if keep
                println(p.output)
            else
                !flush_display(p) && print(p.output, "\r\u1b[A" ^ (p.offset + p.numprintedvalues))
            end
            flush(p.output)
        end
        return nothing
    end

    if t > p.tlast+p.dt
        percentage_complete = 100.0 * p.counter / p.n
        bar = barstring(p.barlen, percentage_complete, barglyphs=p.barglyphs)
        elapsed_time = t - p.tfirst
        est_total_time = 100 * elapsed_time / percentage_complete
        if 0 <= est_total_time <= typemax(Int)
            eta_sec = round(Int, est_total_time - elapsed_time )
            eta = durationstring(eta_sec)
        else
            eta = "N/A"
        end
        msg = @sprintf "%s%3u%%%s  ETA: %s" p.desc round(Int, percentage_complete) bar eta
        !flush_display(p) && print(p.output, "\n" ^ (p.offset + p.numprintedvalues))
        move_cursor_up_while_clearing_lines(p.output, p.numprintedvalues, flush_display(p))
        printover(p.output, msg, p.color, flush_display(p))
        printvalues!(p, showvalues; color = valuecolor)
        !flush_display(p) && print(p.output, "\r\u1b[A" ^ (p.offset + p.numprintedvalues))
        flush(p.output)
        # Compensate for any overhead of printing. This can be
        # especially important if you're running over a slow network
        # connection.
        p.tlast = t + 2*(time()-t)
        p.printed = true
    end
    return nothing
end

function updateProgress!(p::ProgressUnknown; showvalues = Any[], valuecolor = :blue)
    t = time()
    if p.done
        if p.printed
            dur = durationstring(t-p.tfirst)
            msg = @sprintf "%s %d \t Time: %s" p.desc p.counter dur
            move_cursor_up_while_clearing_lines(p.output, p.numprintedvalues, flush_display(p))
            printover(p.output, msg, p.color, flush_display(p))
            printvalues!(p, showvalues; color = valuecolor)
            println(p.output)
            flush(p.output)
        end
        return
    end

    if t > p.tlast+p.dt
        dur = durationstring(t-p.tfirst)
        msg = @sprintf "%s %d \t Time: %s" p.desc p.counter dur
        move_cursor_up_while_clearing_lines(p.output, p.numprintedvalues, flush_display(p))
        printover(p.output, msg, p.color, flush_display(p))
        printvalues!(p, showvalues; color = valuecolor)
        flush(p.output)
        # Compensate for any overhead of printing. This can be
        # especially important if you're running over a slow network
        # connection.
        p.tlast = t + 2*(time()-t)
        p.printed = true
        return
    end
end

# update progress display
"""
`next!(prog, [color])` reports that one unit of progress has been
made. Depending on the time interval since the last update, this may
or may not result in a change to the display.

You may optionally change the color of the display. See also `update!`.
"""
function next!(p::Union{Progress, ProgressUnknown}; options...)
    p.counter += 1
    updateProgress!(p; options...)
end

function next!(p::Union{Progress, ProgressUnknown}, color::Symbol; options...)
    p.color = color
    next!(p; options...)
end

"""
`update!(prog, counter, [color])` sets the progress counter to
`counter`, relative to the `n` units of progress specified when `prog`
was initialized.  Depending on the time interval since the last
update, this may or may not result in a change to the display.

If `prog` is a `ProgressThresh`, `update!(prog, val, [color])` specifies
the current value.

You may optionally change the color of the display. See also `next!`.
"""
function update!(p::Union{Progress, ProgressUnknown}, counter::Int; options...)
    p.counter = counter
    updateProgress!(p; options...)
end

function update!(p::Union{Progress, ProgressUnknown}, counter::Int, color::Symbol; options...)
    p.color = color
    update!(p, counter; options...)
end


"""
`cancel(prog, [msg], [color=:red])` cancels the progress display
before all tasks were completed. Optionally you can specify the
message printed and its color.

See also `finish!`.
"""
function cancel(p::AbstractProgress, msg::AbstractString = "Aborted before all tasks were completed", color = :red; showvalues = Any[], valuecolor = :blue, offset = p.offset, keep = (offset == 0))
    p.offset = offset
    if p.printed
        print(p.output, "\n" ^ (p.offset + p.numprintedvalues))
        move_cursor_up_while_clearing_lines(p.output, p.numprintedvalues, flush_display(p))
        printover(p.output, msg, color, flush_display(p))
        printvalues!(p, showvalues; color = valuecolor)
        if keep
            println(p.output)
        else
            print(p.output, "\r\u1b[A" ^ (p.offset + p.numprintedvalues))
        end
    end
    return
end

"""
`finish!(prog)` indicates that all tasks have been completed.

See also `cancel`.
"""
function finish!(p::Progress; options...)
    while p.counter < p.n
        next!(p; options...)
    end
end

function finish!(p::ProgressUnknown; options...)
    p.done = true
    updateProgress!(p; options...)
end

# Internal method to print additional values below progress bar
function printvalues!(p::AbstractProgress, showvalues; color = false)
    length(showvalues) == 0 && return
    maxwidth = maximum(Int[length(string(name)) for (name, _) in showvalues])
    for (name, value) in showvalues
        msg = "\n  " * rpad(string(name) * ": ", maxwidth+2+1) * string(value)
        (color == false) ? print(p.output, msg) : printstyled(p.output, msg; color=color)
    end
    p.numprintedvalues = length(showvalues)
end

function move_cursor_up_while_clearing_lines(io, numlinesup, clear_output_ijulia)
    if numlinesup > 0 && clear_output_ijulia && isdefined(Main, :IJulia) && Main.IJulia.inited
        Main.IJulia.clear_output(true)
    else
        for _ in 1:numlinesup
            print(io, "\r\u1b[K\u1b[A")
        end
    end
end

function printover(io::IO, s::AbstractString, color::Symbol = :color_normal, clear_output_ijulia=false)
    print(io, "\r")
    printstyled(io, s; color=color)
    if isdefined(Main, :IJulia)
        Main.IJulia.stdio_bytes[] = 0 # issue #76: circumvent IJulia I/O throttling
    elseif isdefined(Main, :ESS) || isdefined(Main, :Atom)
    else
        print(io, "\u1b[K")     # clear the rest of the line
    end
end

function compute_front(barglyphs::BarGlyphs, frac_solid::AbstractFloat)
    barglyphs.front isa Char && return barglyphs.front
    idx = round(Int, frac_solid * (length(barglyphs.front) + 1))
    return idx > length(barglyphs.front) ? barglyphs.fill :
           idx == 0 ? barglyphs.empty :
           barglyphs.front[idx]
end

function barstring(barlen, percentage_complete; barglyphs)
    bar = ""
    if barlen>0
        if percentage_complete == 100 # if we're done, don't use the "front" character
            bar = string(barglyphs.leftend, repeat(string(barglyphs.fill), barlen), barglyphs.rightend)
        else
            n_bars = barlen * percentage_complete / 100
            nsolid = trunc(Int, n_bars)
            frac_solid = n_bars - nsolid
            nempty = barlen - nsolid - 1
            bar = string(barglyphs.leftend,
                         repeat(string(barglyphs.fill), max(0,nsolid)),
                         compute_front(barglyphs, frac_solid),
                         repeat(string(barglyphs.empty), max(0, nempty)),
                         barglyphs.rightend)
        end
    end
    bar
end

function durationstring(nsec)
    days = div(nsec, 60*60*24)
    r = nsec - 60*60*24*days
    hours = div(r,60*60)
    r = r - 60*60*hours
    minutes = div(r, 60)
    seconds = floor(r - 60*minutes)

    hhmmss = @sprintf "%u:%02u:%02u" hours minutes seconds
    if days>9
        return @sprintf "%.2f days" nsec/(60*60*24)
    elseif days>0
        return @sprintf "%u days, %s" days hhmmss
    end
    hhmmss
end

end

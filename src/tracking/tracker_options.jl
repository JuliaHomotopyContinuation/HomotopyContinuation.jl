
###
### Options and Parameters
###

"""
    TrackerParameters

Parameters that control the performance and robustness characteristics of the path tracking
algorithm.
"""
Base.@kwdef mutable struct TrackerParameters
    N::Int = 3
    target_accuracy::Float64 = 1e-14
    predictor_local_order::Int = 5
end
Base.show(io::IO, TP::TrackerParameters) = print_fieldnames(io, TP)


"""
    TrackerOptions(; options...)

The set of options for a [`Tracker`](@ref).

## Options

* `max_steps = 10_000`: The maximal number of steps a tracker attempts
* `max_step_size = Inf`: The maximal size of a step
* `max_initial_step_size = Inf`: The maximal size of the first step
* `min_step_size = 1e-48`: The minimal step size. If a smaller step size would
  be necessary, then the tracking gets terminated.
* `extended_precision = true`: Whether to allow for the use of extended precision,
  if necessary, in some computations. This can greatly improve the ability to track
  numerically difficult paths.
* `terminate_cond = 1e13`: If the relative component-wise condition number
  `cond(H_x, ẋ)` is larger than `terminate_cond` then the path is terminated as too
  ill-conditioned.
* `parameters::TrackerParameters = TrackerParameters()`
"""
Base.@kwdef mutable struct TrackerOptions
    max_steps::Int = 10_000
    max_step_size::Float64 = Inf
    max_initial_step_size::Float64 = Inf
    extended_precision::Bool = true
    min_step_size::Float64 = 1e-48
    min_rel_step_size::Float64 = 0.0
    terminate_cond::Float64 = 1e13
    parameters::TrackerParameters = TrackerParameters()
end

Base.show(io::IO, opts::TrackerOptions) = print_fieldnames(io, opts)

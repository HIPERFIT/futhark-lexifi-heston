import "futlib/math"

import "rand"

module type distance = {
  val distance: []f64 -> []f64 -> f64
}

module absolute_distance: distance = {
  fun distance (quotes: [num_quotes]f64) (prices: [num_quotes]f64): f64 =
    let norm (price: f64) (quote: f64) =
      (let rel = (price - quote) / quote
       in rel * rel)
    in reduce (+) 0.0 (map norm quotes prices)
}

module relative_distance: distance = {
  fun distance (quotes: [num_quotes]f64) (prices: [num_quotes]f64): f64 =
    let norm (price: f64) (quote: f64) =
      (let dif = price - quote
       in dif * dif)
    let a = map norm quotes prices
    let a[0] = a[0] + 1.0
    let a[0] = a[0] - 1.0
    in reduce (+) 0.0 a
}

module type pricer = {
  type pricer_ctx
  val pricer: pricer_ctx -> []f64 -> []f64

  include distance
}

type range = {lower_bound: f64,
              upper_bound: f64,
              initial_value: f64}

-- Pretend this is a sum type with two constructors.
type optimization_variable = ( bool -- fixed?
                             , f64 -- value if fixed
                             , range -- range if not fixed
                             )

fun fixed_value (v: f64): optimization_variable =
  (true, v, {lower_bound=0.0, upper_bound=0.0, initial_value=0.0})

fun optimize_value (r: range): optimization_variable =
  (false, 0.0, r)

module least_squares(P: pricer) = {

  type calibration_result = { parameters: []f64,
                              root_mean_squared_error: f64,
                              quoted_prices: []f64,
                              calibrated_prices: []f64,
                              nb_feval: i32 }

  -- Parameterisation of how the randomised search takes place.
  type mutation = {np: i32, -- Population size
                   cr: f64  -- Crossover probability [0,1]
                   }

  type termination = {max_iterations: i32, max_global: i32, target: f64}

  type status = i32 -- Pretend it's opaque!
  val max_iterations_reached: status = 0
  val max_global_reached: status = 1
  val target_reached: status = 2

  type result = {x0: []f64, f: f64, nb_feval: i32, status: status}

  fun active_vars (vars_to_free_vars: [num_vars]i32)
                  (variables: [num_vars]optimization_variable)
                  (xs: [num_active]f64) =
    map (\fv (fixed,x,_) -> if fixed then x else unsafe xs[fv])
        vars_to_free_vars variables

  fun min_and_idx (a:f64,a_i:i32) (b:f64,b_i:i32) =
    if      a < b     then (a,a_i)
    else if b < a     then (b,b_i)
    else if a_i < b_i then (a, a_i)
    else                   (b, b_i)

  fun optimize (pricer_ctx: P.pricer_ctx)
               (quotes: [num_quotes]f64)
               (vars_to_free_vars: [num_vars]i32)
               (variables: [num_vars]optimization_variable)
               ({np, cr}: mutation)
               (lower_bounds: [num_free_vars]f64)
               (upper_bounds: [num_free_vars]f64)
               ({max_iterations,max_global,target}: termination): result =
    -- The objective function.  This could be factored out into a
    -- function argument (as a parametric module).
    let objective (x: [num_free_vars]f64): f64 =
      P.distance quotes (P.pricer pricer_ctx (active_vars vars_to_free_vars variables x))

    let rng = random_f64.rng_from_seed 0x123
    let rngs = random_f64.split_rng np rng
    let (rngs, rss) = unzip (map (\rng -> random_f64.nrand rng (0.0, 1.0) num_free_vars) rngs)
    let rng = random_f64.join_rng rngs
    let x = (let init_j (lower_bound: f64) (upper_bound: f64) (r: f64) =
               lower_bound + (upper_bound-lower_bound) * r
             let init_i (rs: [num_free_vars]f64) = map init_j lower_bounds upper_bounds rs
             in map init_i rss)
    let fx = map objective x
    let fx[0] = fx[0]+1.0
    let fx[0] = fx[0]-1.0
    let (fx0, best_idx) =
      reduceComm min_and_idx (f64.inf, 0) (zip fx (iota np))

    let mutation (difw: f64) (best_idx: i32) (x: [np][num_free_vars]f64)
                 (rng: random_f64.rng) (i :i32) (x_i: [num_free_vars]f64) =
      (-- We have to draw 'to_draw' distinct elements from 'x', and it
       -- can't be 'i'.  We do this with brute-force looping.
       let (rng,a) = random_i32.rand rng (0,np)
       let (rng,b) = random_i32.rand rng (0,np)
       let (rng,c) = random_i32.rand rng (0,np)
       loop ((rng,a)) = while a == i do random_i32.rand rng (0,np)
       loop ((rng,b)) = while b == i || b == a do random_i32.rand rng (0,np)
       loop ((rng,c)) = while c == i || c == a || c == b do random_i32.rand rng (0,np)
       let (rng,r) = random_f64.rand rng (0.0, 1.0)
       let x_r1 = unsafe if r <= 0.5 then x[best_idx] else x[a]
       let x_r2 = unsafe x[b]
       let x_r3 = unsafe x[c]
       let (rng,j0) = random_i32.rand rng (0,num_free_vars)
       let (rng,rs) = random_f64.nrand rng (0.0, 1.0) num_free_vars
       let auxs = map (+) x_r1 (map (difw*) (map (-) x_r2 x_r3))
       let v_i = map (\j r lower_bound upper_bound aux x_i_j ->
                      if (j == j0 || r <= cr) && lower_bound <= aux && aux <= upper_bound
                      then aux
                      else x_i_j)
                     (iota num_free_vars) rs lower_bounds upper_bounds auxs x_i

       in (rng, v_i))

    let recombination (fx0: f64) (best_idx: i32) (fx: [np]f64) (x: [np][num_free_vars]f64) (v: [np][num_free_vars]f64) =
      (let f_v = map objective v
       let fx' = map f64.min f_v fx
       let x' = map (\f fx_i x_i v_i -> if f < fx_i then v_i else x_i)
                    f_v fx x v
       let f_v[0] = f_v[0] + 1.0
       let f_v[0] = f_v[0] - 1.0
       let (fx0', best_idx') =
         reduceComm min_and_idx (fx0, best_idx) (zip f_v (iota np))
       in (fx0', best_idx', fx', x'))

    -- We are not counting the numer of invocations of the objective
    -- function quite as in LexiFi's code (they use a closure that
    -- increments a counter), but we should be close.
    loop ((rng, ncalls, nb_it,
           (fx0, best_idx, fx, x)) =
          (rng, np, max_iterations,
           (fx0, best_idx, fx, x))) = while nb_it > 0 && max_global > ncalls && fx0 > target do
      (let (rng,differential_weight) = random_f64.rand rng (0.5, 1.0)
       let rngs = random_f64.split_rng np rng
       let (rngs, v) = unzip (map (mutation differential_weight best_idx x) rngs (iota np) x)
       let rng = random_f64.join_rng rngs
       let (fx0, best_idx, fx, x) = recombination fx0 best_idx fx x v
       in (rng, ncalls + np, nb_it - 1,
           (fx0, best_idx, fx, x)))
    let x0 = x[best_idx]
    let status = if      fx0 <= target       then target_reached
                 else if max_global < ncalls then max_global_reached
                 else if nb_it == 0          then max_iterations_reached
                 else 1337 -- never reached
    in {x0=x0, f=fx0, nb_feval=ncalls, status=status}

  fun least_squares
      (pricer_ctx: P.pricer_ctx)
      (max_global: i32)
      (np: i32)
      (variables: [num_vars]optimization_variable)
      (quotes: [num_quotes]f64)
      : calibration_result =
    let (free_vars_to_vars, free_vars) =
      unzip (filter (\(_, (fixed, _, _)) -> !fixed) (zip (iota num_vars) variables))
    let num_free_vars = (shape free_vars)[0]
    let vars_to_free_vars = write free_vars_to_vars (iota num_free_vars)
                                  (replicate num_vars (-1))
    let (x, lower_bounds, upper_bounds) =
      unzip (map (\(_, _, {initial_value, lower_bound, upper_bound}) ->
                  (initial_value, lower_bound, upper_bound)) free_vars)

    let rms_of_error (err: f64) = f64.sqrt(err * (10000.0 / f64 num_quotes))

    let (x, nb_feval) =
      if max_global > 0
      then let res = (optimize pricer_ctx quotes vars_to_free_vars variables
                      {np = np, cr = 0.9} lower_bounds upper_bounds
                      {max_iterations = 0x7FFFFFFF, max_global = max_global, target = 0.0})
           in (#x0 res, #nb_feval res)
      else (x, 0)

    let prices = P.pricer pricer_ctx (active_vars vars_to_free_vars variables x)

    in {parameters = active_vars vars_to_free_vars variables x,
        root_mean_squared_error = rms_of_error (P.distance quotes prices),
        quoted_prices = quotes,
        calibrated_prices = prices,
        nb_feval = nb_feval}
}

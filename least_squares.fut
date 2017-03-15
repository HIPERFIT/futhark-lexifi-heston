import "futlib/math"

import "rand"

module type distance = {
  val distance: []f64 -> []f64 -> f64
}

module absolute_distance: distance = {
  fun distance (quotes: [m]f64) (prices: [m]f64): f64 =
    let norm (price: f64) (quote: f64) =
      (let rel = (price - quote) / quote
       in rel * rel)
    in reduce (+) 0.0 (map norm quotes prices)
}

module relative_distance: distance = {
  fun distance (quotes: [m]f64) (prices: [m]f64): f64 =
    let norm (price: f64) (quote: f64) =
      (let dif = price - quote
       in dif * dif)
    in reduce (+) 0.0 (map norm quotes prices)
}

module type pricer = {
  type parameters
  type pricer_ctx
  val parameters_of_vector: []f64 -> parameters
  val pricer: pricer_ctx -> parameters -> []f64

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

  type calibration_result = { parameters: P.parameters,
                              root_mean_squared_error: f64,
                              quoted_prices: []f64,
                              calibrated_prices: []f64 }

  type t = {np: i32, -- Population size
            cr: f64  -- Crossover probability [0,1]
           }

  type termination = {maxit: i32, maxf: i32, target: f64}

  type status = i32 -- Pretend it's opaque!
  val maxit_reached: status = 0
  val maxf_reached: status = 1
  val target_reached: status = 2

  type result = {x0: []f64, f: f64, nb_feval: i32, status: status}

  fun parameters_of_active_vars (vars_to_free_vars: [num_variables]i32)
                                (variables: [num_variables]optimization_variable)
                                (xs: [num_active]f64) =
    P.parameters_of_vector (map (\fv (fixed,x,_) ->
                                 if fixed then x else unsafe xs[fv])
                            vars_to_free_vars variables)

  fun optimize (pricer_ctx: P.pricer_ctx)
               (quotes: [m]f64)
               (vars_to_free_vars: [num_variables]i32)
               (variables: [num_variables]optimization_variable)
               ({np, cr}: t)
               (lower_bounds: [n]f64)
               (upper_bounds: [n]f64)
               ({maxit,maxf,target}: termination): result =
  -- The optimisation function.  This could be factored out into a
  -- function argument (as a parametric module).
    let f (x: []f64): f64 =
      P.distance quotes (P.pricer pricer_ctx (parameters_of_active_vars vars_to_free_vars variables x))

    let rng = random_f64.rng_from_seed 0x123
    -- number instead.
    let rngs = random_f64.split_rng np rng
    let (rngs, rss) = unzip (map (\rng -> random_f64.nrand rng (0.0, 1.0) n) rngs)
    let rng = random_f64.join_rng rngs
    let x = (let init_j (lower: f64) (upper: f64) (r: f64) =
               (let blo = lower
                let db = upper - blo
                in blo + db * r)
             let init_i (rs: [n]f64) = map init_j lower_bounds upper_bounds rs
             in map init_i rss)
    let fx = map f x
    let (fx0, best_idx) =
      reduce (\(a,a_i) (b,b_i) -> if a < b then (a,a_i) else (b,b_i))
             (f64.inf, 0) (zip fx (iota np))

    let mutation (difw: f64) (best_idx: i32) (x: [np][n]f64)
                 (rng: random_f64.rng) (i :i32) (x_i: [n]f64) =
      (let to_draw = 3
       -- We have to draw 'to_draw' distinct elements from 'x', and
       -- it can't be 'i'.  We do this with a brute-force loop.
       let indices = replicate to_draw i
       let drawn = 0
       loop ((rng,indices,drawn)) = while drawn < to_draw do
         (let (rng, j) = random_i32.rand rng (0,np)
          loop (ok = true) = for l < drawn do ok && unsafe indices[l] != j
          in if ok
             then unsafe let indices[drawn] = j in (rng, indices, drawn+1)
             else (rng, indices, drawn))
       let (rng,r) = random_f64.rand rng (0.0, 1.0)
       let x_r1 = unsafe if r <= 0.5 then x[best_idx] else unsafe x[indices[0]]
       let x_r2 = unsafe x[indices[1]]
       let x_r3 = unsafe x[indices[2]]
       let (rng,j0) = random_i32.rand rng (0,n)
       let (rng,rs) = random_f64.nrand rng (0.0, 1.0) n
       let auxs = map (+) x_r1 (map (difw*) (map (-) x_r2 x_r3))
       let v_i = map (\j r lower_bound upper_bound aux x_i_j ->
                      if (j == j0 || r <= cr) && lower_bound <= aux && aux <= upper_bound
                      then aux
                      else x_i_j)
                     (iota n) rs lower_bounds upper_bounds auxs x_i

       in (rng, v_i))

    let recombination (fx0: f64) (best_idx: i32) (fx: [np]f64) (x: [np][n]f64) (v: [np][n]f64) =
      (let f_v = map f v
       let fx' = map f64.min f_v fx
       let x' = map (\f fx_i x_i v_i -> if f < fx_i then v_i else x_i)
                    f_v fx x v
       let (fx0', best_idx') =
         reduce (\(a,a_i) (b,b_i) -> if a < b then (a,a_i) else (b,b_i))
                (fx0, best_idx) (zip f_v (iota np))
       in (fx0', best_idx', fx', x'))

    -- We are not counting the numer of invocations of the objective
    -- function quite as in LexiFi's code (they use a closure that
    -- increments a counter), but we should be close.
    loop ((rng, ncalls, nb_it,
           (fx0, best_idx, fx, x)) =
          (rng, np*2, maxit,
           (fx0, best_idx, fx, x))) = while nb_it > 0 && maxf > ncalls && fx0 > target do
      (let (rng,differential_weight) = random_f64.rand rng (0.5, 1.0)
       let rngs = random_f64.split_rng np rng
       let (rngs, v) = unzip (map (mutation differential_weight best_idx x) rngs (iota np) x)
       let rng = random_f64.join_rng rngs
       let (fx0, best_idx, fx, x) = recombination fx0 best_idx fx x v
       in (rng, ncalls + np, nb_it - 1,
           (fx0, best_idx, fx, x)))
    let x0 = x[best_idx]
    let status = if      fx0 <= target then target_reached
                 else if maxf < ncalls then maxf_reached
                 else if nb_it == 0    then maxit_reached
                 else 1337 -- never reached
    in {x0=x0, f=fx0, nb_feval=ncalls, status=status}

  fun default_parameters (n: i32) = {np = i32.min (10*n) 40,
                                     cr = 0.9}

  fun least_squares
      (pricer_ctx: P.pricer_ctx)
      (max_global: i32)
      (variables: [num_variables]optimization_variable)
      (quotes: [m]f64)
      : calibration_result =
    let (free_vars_to_vars, free_vars) =
      unzip (filter (\(_, (fixed, _, _)) -> !fixed) (zip (iota num_variables) variables))
    let num_free_vars = (shape free_vars)[0]
    let vars_to_free_vars = write free_vars_to_vars (iota num_free_vars)
                                  (replicate num_variables (-1))
    let (x, lower_bounds, upper_bounds) =
      unzip (map (\(_, _, {initial_value, lower_bound, upper_bound}) ->
                  (initial_value, lower_bound, upper_bound)) free_vars)

    let rms_of_error (err: f64) = f64.sqrt(err * (10000.0 / f64 m))

    let x =
      if max_global > 0
      then #x0 (optimize pricer_ctx quotes vars_to_free_vars variables
                (default_parameters num_free_vars) lower_bounds upper_bounds
                {maxit = 0x7FFFFFFF, maxf = max_global, target = 0.0})
      else x

    let prices = P.pricer pricer_ctx (parameters_of_active_vars vars_to_free_vars variables x)

    in {parameters = parameters_of_active_vars vars_to_free_vars variables x,
        root_mean_squared_error = rms_of_error (P.distance quotes prices),
        quoted_prices = quotes,
        calibrated_prices = prices}
}

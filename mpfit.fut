module mpfit(P: val f: i32 -> i32 -> []f64 -> []f64) = {
  type mp_par = ( bool -- fixed?
                , (bool, f64) -- lower limited? and limit
                , (bool, f64) -- upper limited? and limit
                )

  type config = ( f64 -- ftol
                , f64 -- xtol
                , f64 -- gtol
                , i32 -- maxfev
                )

  -- I think 'm' is the siez of the array produced by mpfit.
  fun mpfit(m: i32) (xall: [n]f64) (pars: [n]mp_par) (config: config): [n]64
}

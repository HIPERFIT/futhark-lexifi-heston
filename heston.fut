-- ==
-- compiled input @ data_1062_quotes
-- compiled input @ data_10000_quotes
-- compiled input @ data_100000_quotes

import "futlib/date"
import "least_squares"
import "price_european_calls"

type calibration_input = { today: date
                         , quotes: []{maturity: date, strike: f64, quote: f64}
                         , max_global: i32
                         , np: i32
                         , strike_weight_bandwidth: f64
                         , maturity_weight_x0: f64
                         , maturity_weight_gamma: f64
                         , integral_iterations: nb_points
                         , variables: []optimization_variable }

fun heston_parameters_from_vector (x: [5]f64) =
  { initial_variance = x[0]
  , long_term_variance = x[1]
  , correlation = x[2]
  , mean_reversion = x[3]
  , variance_volatility = x[4] }

module heston_least_squares = least_squares {
  type pricer_ctx = {day_count_fractions: []f64,
                     quotes: []{maturity: i32, strike: f64, vega: f64, weight: f64},
                     gauss_laguerre_coefficients: ([]f64, []f64)
                     }

  fun pricer ({day_count_fractions, quotes, gauss_laguerre_coefficients}: pricer_ctx) (x: []f64) =
    let heston_parameters = heston_parameters_from_vector x
    let prices = price_european_calls
                 gauss_laguerre_coefficients
                 false 1.0 1.0 1.0
                 heston_parameters
                 day_count_fractions
                 (map (\q -> {maturity=#maturity q, strike=#strike q}) quotes)
    in map (\q p -> #weight q * p / #vega q) quotes prices

  open relative_distance
}

fun distinct_maturities (dates: [n]date): ([]date, [n]i32) =
  let switched (x: date) (i: i32) = i == 0 || unsafe !(same_date x dates[i-1])
  let switches = map switched dates (iota n)
  in (#2 (unzip (filter (\x -> #1 x) (zip switches dates))),
      map (-1) (scan (+) 0 (map i32 switches)))

fun run_calibration({today,
                     quotes,
                     max_global,
                     np,
                     strike_weight_bandwidth,
                     maturity_weight_x0,
                     maturity_weight_gamma,
                     integral_iterations,
                     variables}: calibration_input): heston_least_squares.calibration_result =
  let price_and_vega_of_quote (strike: f64) (maturity: date) (quote: f64) =
    (let (price, vega) = bs_call true today 1. strike maturity quote
     in (price, f64.max 1e-1 vega))
  let strike_weight (p: f64) (x: f64) = f64.exp (p * (f64.log x + 1.0 - x))
  let maturity_weight (x0: f64) (gamma: f64) (x: f64) =
      (let k = 1.0 / (f64.exp(gamma * x0) - 1.0)
       in if x <= x0 then k * (f64.exp(gamma * x) - 1.0) else 1.0)
  let weight (strike: f64) (mat: date) =
    maturity_weight maturity_weight_x0 maturity_weight_gamma (diff_dates today mat) *
    strike_weight strike_weight_bandwidth strike


  let (maturity_dates, quotes_to_maturities) =
    distinct_maturities (map (\q -> #maturity q) quotes)
  let weights = map (\{maturity, strike, quote=_} -> weight strike maturity) quotes
  let prices_and_vegas = map (\{maturity, strike, quote} ->
                              price_and_vega_of_quote strike maturity quote) quotes
  let quotes_for_optimization = map (\(p,v) w -> w * p / v) prices_and_vegas weights
  let quotes_for_ctx =
    map (\{maturity=_, strike, quote=_} w (_,v) i -> { maturity = i
                                                     , strike = strike
                                                     , weight = w
                                                     , vega = v})
        quotes weights prices_and_vegas quotes_to_maturities

  let ctx = { day_count_fractions =
                map (diff_dates today) maturity_dates
            , quotes =
                quotes_for_ctx
            , gauss_laguerre_coefficients =
                gauss_laguerre_coefficients integral_iterations }

  in heston_least_squares.least_squares ctx max_global np variables quotes_for_optimization

fun date_of_int(x: i32) =
  let d = x%100
  let m = (x/100)%100
  let y = x/10000
  in date_of_triple (y, m, d)

let default_variables: []optimization_variable =
    [optimize_value {lower_bound =  1e-6, initial_value = 4e-2, upper_bound = 1.},
     optimize_value {lower_bound =  1e-6, initial_value = 4e-2, upper_bound = 1.},
     optimize_value {lower_bound = -1.  , initial_value = -0.5, upper_bound = 0.},
     optimize_value {lower_bound =  1e-4, initial_value = 1e-2, upper_bound = 4.},
     optimize_value {lower_bound =  1e-4, initial_value = 0.4, upper_bound = 2.}
    ]

fun main (max_global: i32)
         (nb_points: i32)
         (np: i32)
         (today: i32)
         (quotes_maturity: [num_quotes]i32)
         (quotes_strike: [num_quotes]f64)
         (quotes_quote: [num_quotes]f64) =
  let result =
    run_calibration { today = date_of_int today
                    , quotes = map (\m k q -> {maturity = date_of_int m, strike = k, quote = q})
                      quotes_maturity quotes_strike quotes_quote
                    , max_global = max_global
                    , np = np
                    , strike_weight_bandwidth = 0.0
                    , maturity_weight_x0 = 0.0
                    , maturity_weight_gamma = 1.0
                    , integral_iterations = if nb_points == 10 then ten else twenty
                    , variables = default_variables
                      }
  let { initial_variance,
        long_term_variance,
        correlation,
        mean_reversion,
        variance_volatility} = heston_parameters_from_vector (#parameters result)
  in (#root_mean_squared_error result,
      #nb_feval result,
      initial_variance,
      long_term_variance,
      mean_reversion,
      variance_volatility,
      correlation
      )

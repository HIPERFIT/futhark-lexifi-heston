import "least_squares"
import "date"
import "price_european_calls"

type calibration_input = { today: date
                         , quotes: [](date, f64, f64)
                         , max_global: i32
                         , strike_weight_bandwidth: f64
                         , maturity_weight_x0: f64
                         , maturity_weight_gamma: f64
                         , integral_iterations: nb_points
                         , variables: []optimization_variable }

module heston_least_squares = least_squares {
  type parameters = heston_parameters
  type pricer_ctx = {maturities: []{day_count_fraction: f64, vega: f64, weight: f64},
                     quotes: []{maturity: i32, strike: f64},
                     integral_iterations: nb_points
                     }

  fun parameters_of_vector (x: [5]f64) =
    { initial_variance = x[0]
    , long_term_variance = x[1]
    , correlation = x[2]
    , mean_reversion = x[3]
    , variance_volatility = x[4] }

  fun pricer ({maturities, quotes, integral_iterations}: pricer_ctx) (heston_parameters: parameters) =
    let day_count_fractions = map (\q -> #day_count_fraction q) maturities
    let prices = price_european_calls
                 (gauss_laguerre_coefficients integral_iterations)
                 false 1.0 1.0 1.0
                 heston_parameters
                 day_count_fractions quotes
    in prices

  open relative_norm
}

fun distinct_maturities (dates: [n]date): ([]date, [n]i32) =
  let switched (x: date) (i: i32) = i == 0 || !(same_date x dates[i-1])
  let switches = map switched dates (iota n)
  in (#2 (unzip (filter (\x -> #1 x) (zip switches dates))),
      scan (+) 0 (map i32 switches))

entry run_calibration({today,
                       quotes,
                       max_global,
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
    maturity_weight maturity_weight_x0 maturity_weight_gamma (sub_act_365 mat today) *
    strike_weight strike_weight_bandwidth strike

  let quotes_for_optimization =
    map (\(maturity, strike, quote) ->
         let (p, v) = price_and_vega_of_quote strike maturity quote
         in weight strike maturity * p / v)
        quotes

  let (maturities, quotes_to_maturities) =
    distinct_maturities (map (\(m,_,_) -> m) quotes)

  let ctx = { maturities = maturities
            , quotes = quotes
            , integral_iterations = integral_iterations
              }

  in heston_least_squares.least_squares ctx max_global variables quotes_for_optimization

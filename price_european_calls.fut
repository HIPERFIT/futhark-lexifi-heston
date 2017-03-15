import "futlib/math"
import "futlib/complex"
import "date"

module c64 = complex(f64)
type c64 = c64.complex

fun (x: c64) +! (y: c64) = x c64.+ y
fun (x: c64) -! (y: c64) = x c64.- y
fun (x: c64) *! (y: c64) = x c64.* y
fun (x: c64) /! (y: c64) = x c64./ y

val zero: c64 = c64.mk_re 0.0
val one: c64 = c64.mk_re 1.0
val two: c64 = c64.mk_re 2.0

val isqrt2pi = 2.0 * f64.pi ** -0.5

val inv_sqrt2 = 1.0 / f64.sqrt 2.0

fun erfc(x: f64): f64 =
  let a1 =  0.254829592
  let a2 = -0.284496736
  let a3 =  1.421413741
  let a4 = -1.453152027
  let a5 =  1.061405429
  let p  =  0.3275911

  let sign = if x < 0.0 then -1.0 else 1.0
  let x = f64.abs x

  let t = 1.0/(1.0 + p*x)
  let y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*f64.exp(-x*x)

  in sign*y

fun erf (x: f64): f64 =
  let a1 =  0.254829592
  let a2 = -0.284496736
  let a3 =  1.421413741
  let a4 = -1.453152027
  let a5 =  1.061405429
  let p  =  0.3275911

  let t = 1. / (1. + p * x)
  let t2 = t * t
  let t3 = t * t2
  let t4 = t *t3
  let t5 = t * t4
  in 1. - (a1 * t + a2 * t2 + a3 * t3 + a4 * t4 + a5 * t5) * f64.exp ((-x) * x)

fun pnorm (x: f64): f64 =
  let u = x / f64.sqrt 2.
  let erf = if u < 0. then -erf (-u) else erf u
  in 0.5 * (1. + erf)

fun ugaussian_pdf (x: f64) =
  if f64.isinf x then 0.0
  else isqrt2pi * f64.exp (-0.5 * x * x)

fun ugaussian_P (x: f64) =
  if f64.isinf x then (if x > 0.0 then 1.0 else 0.0)
  else 1.0 - 0.5 * erfc (x * inv_sqrt2)

type nb_points = bool -- Pretend it's opaque.
val ten: nb_points = true
val twenty: nb_points = false

fun gauss_laguerre_coefficients (nb: nb_points) =
  if nb == ten then
  ([

   0.1377934705404924298211, 0.7294545495031707904587,
   1.8083429017403143124199, 3.4014336978548351808627,
   5.5524961400642398601235,

   8.3301527467632645596041, 11.8437858379019154142497,
   16.2792578313766149733510, 21.9965858119813617577165,
   29.9206970122737736517138

   ],
   [

   0.3540097386069980256451, 0.8319023010435621090508,
   1.3302885617494675241090, 1.8630639031111377867944,
   2.4502555580731559814467,

   3.1227641551735070279960, 3.9341526954948387029276,
   4.9924148722151180379569, 6.5722024851118048260901,
   9.7846958403678012672344

  ])
  else -- nb == twenty
  ([

    0.070539889691988738596, 0.372126818001613290932,
    0.916582102483245675373, 1.707306531028168317121,
    2.749199255315394108123, 4.048925313808060089116,
    5.615174970938836551682, 7.459017454225468135576,
    9.594392865483230892210, 12.038802560608859337776,
    14.814293416061961039532, 17.948895543122549867121,
    21.478788273773705697067, 25.451702644185488111361,
    29.932554890200286479285, 35.013433968980073984767,
    40.833057239401291838021, 47.619993970772064528774,
    55.810795768051811194255, 66.524416523865880890298

   ],
   [

    0.18108006241898921829, 0.42255676787860191324,
    0.66690954670145563554, 0.91535237278121706073,
    1.16953970728044676086, 1.43135498604326039107,
    1.70298113670862205637, 1.98701590220746626692,
    2.28663572001997383865, 2.60583491316971160856,
    2.94978326190603512558, 3.32539692360021144069,
    3.74225473838238897883, 4.21423782504235067137,
    4.76251619189997654757, 5.42172741864077067930,
    6.25401126576962873571, 7.38731454069754001068,
    9.15132857271978572555, 12.89338863845354232751

   ])

type heston_parameters = { initial_variance: f64
                         , long_term_variance: f64
                         , mean_reversion: f64
                         , variance_volatility: f64
                         , correlation: f64 }

module type pricer_parameter = {
  val normal: bool
  val psi_h: f64 -> heston_parameters -> c64 -> c64
  val psi_bs: c64 -> c64 -> c64 -> c64
  val moneyness_f: f64 -> f64
}

module normal_true: pricer_parameter = {
  val normal = true

  fun psi_h (day_count_fraction: f64) (heston_parameters: heston_parameters) (xi: c64) =
    let {initial_variance = v0,
         long_term_variance = theta,
         mean_reversion = kappa,
         variance_volatility = eta,
         correlation = rho} = heston_parameters
    let kappai = c64.mk_re kappa
    let etai = c64.mk_re eta
    let etai2 = etai *! etai
    let coeff1 = kappai *! c64.mk_re theta /! etai2
    let coeff2 = c64.mk_re v0 /! etai2
    let ti = c64.mk_re day_count_fraction
    let i = c64.mk_im 1.0
    let d0 = kappai -! (c64.mk_im rho) *! etai *! xi
    let d = c64.sqrt (d0 *! d0 +! etai2 *! xi *! xi)
    let a_minus = d0 -! d
    let g = a_minus /! (d0 +! d)
    let e = c64.exp (zero -! d *! ti)
    in c64.exp (xi *! i +!
                coeff1 *! (a_minus *! ti -! two *! c64.log ((one -! g *! e) /! (one -! g))) +!
                coeff2 *! a_minus *! (one -! e) /! (one -! g *! e))

  fun psi_bs (minus_half_sigma2_t: c64) (i: c64) (xi: c64) =
    c64.exp (xi *! i +! minus_half_sigma2_t *! xi *! xi)

  fun moneyness_f (k: f64) = -k
}

module normal_false: pricer_parameter = {
  val normal = false

  fun psi_h (day_count_fraction: f64) (heston_parameters: heston_parameters) (xi: c64) =
    let {initial_variance = v0,
         long_term_variance = theta,
         mean_reversion = kappa,
         variance_volatility = eta,
         correlation = rho} = heston_parameters
    let kappai = c64.mk_re kappa
    let etai = c64.mk_re eta
    let etai2 = etai *! etai
    let coeff1 = kappai *! c64.mk_re theta /! etai2
    let coeff2 = c64.mk_re v0 /! etai2
    let ti = c64.mk_re day_count_fraction
    let i = c64.mk_im 1.0
    let d0 = kappai -! (c64.mk_im rho) *! etai *! xi
    let d = c64.sqrt (d0 *! d0 +! etai2 *! xi *! (i +! xi))
    let a_minus = d0 -! d
    let g = a_minus /! (d0 +! d)
    let e = c64.exp (zero -! d *! ti)
    in c64.exp (coeff1 *! (a_minus *! ti -! two *! c64.log ((one -! g *! e) /! (one -! g))) +!
                coeff2 *! a_minus *! (one -! e) /! (one -! g *! e))

  fun psi_bs (minus_half_sigma2_t: c64) (i: c64) (xi: c64) =
    c64.exp (minus_half_sigma2_t *! xi *! (i *! xi))

  fun moneyness_f (k: f64) = -(f64.log k)
}

fun bs_control (moneyness: f64) (sigma_sqrtt: f64) =
  let d1 = (-f64.log moneyness) / sigma_sqrtt + 0.5 * sigma_sqrtt
  in pnorm d1 - moneyness * pnorm (d1 - sigma_sqrtt)

fun price_european_calls
    (x: [n]f64, w: [n]f64)
    (ap1: bool)
    (spot: f64)
    (df_div: f64)
    (df: f64)
    (heston_parameters: heston_parameters)
    (day_count_fractions: [nmaturities]f64)
    (quotes: [nstrikes]{strike: f64, maturity: i32})
  : [nstrikes]f64 =
       let {initial_variance = v0, long_term_variance = theta, mean_reversion = kappa, correlation = rho, variance_volatility = eta} = heston_parameters
       let maturity_for_quote = map (\q -> #maturity q) quotes
       let strikes = map (\q -> #strike q) quotes
       let f0 = spot * df_div / df
       let kappai = c64.mk_re kappa
       let etai = c64.mk_re eta
       let etai2 = etai *! etai
       let coeff1 = kappai *! c64.mk_re theta /! etai2
       let coeff2 = c64.mk_re v0 /! etai2
       let i = c64.mk_im 1.0
       let psi_h (day_count_fraction: f64) (xi: c64) =
         (let ti = c64.mk_re day_count_fraction
          let d0 = kappai -! (c64.mk_im rho) *! etai *! xi
          let d = c64.sqrt (d0 *! d0 +! etai2 *! (xi *! (i +! xi)))
          let a_minus = d0 -! d
          let g = a_minus /! (d0 +! d)
          let e = c64.exp (zero -! d *! ti)
          in c64.exp (coeff1 *! (a_minus *! ti -! two *! c64.log ((one -! g *! e) /! (one -! g))) +!
                      coeff2 *! a_minus *! (one -! e) /! (one -! g *! e)))
       let sigma2 (day_count_fraction: f64) =
         (if ap1 then v0
          else let eta = -1.
               let eps = 1e-2
               let two_da_time_eps = psi_h day_count_fraction (c64.mk eps eta) -!
                                     psi_h day_count_fraction (c64.mk (-eps) eta)
               let two_db_time_eps = psi_h day_count_fraction (c64.mk_im (eta + eps)) -!
                                     psi_h day_count_fraction (c64.mk_im (eta - eps))
               in 0.5 * (c64.im (two_da_time_eps -! i *! two_db_time_eps)) /
            (day_count_fraction * eps))
       let psi_bs (day_count_fraction: f64) (xi: c64) =
         (let minus_half_sigma2_t = c64.mk_re (-0.5 * day_count_fraction * sigma2 day_count_fraction)
          in c64.exp (minus_half_sigma2_t *! (xi *! (i +! xi))))
       let moneyness = map (/f0) strikes
       let minus_ik = map (\k -> c64.mk_im (- f64.log k)) moneyness

       let iter (j: i32): [nstrikes]f64 =
         (let xj = x[j]
          let wj = w[j]
          let x = c64.mk_re xj
          let mk_w_and_coeff_k (day_count_fraction: f64) =
            (if ap1
             then (let x_minus_half_i = x -! c64.mk_im 0.5
                   in (wj / (0.25 + xj * xj),
                       (psi_bs day_count_fraction x_minus_half_i -! psi_h day_count_fraction x_minus_half_i)))
             else (let x_minus_i = x -! i
                   in (wj,
                       (psi_bs day_count_fraction x_minus_i -! psi_h day_count_fraction x_minus_i) /!
                       (x *! x_minus_i))))
          let (ws, coeff_ks) = unzip (map mk_w_and_coeff_k day_count_fractions)
          in map (\minus_ikk m ->
                  let w = unsafe ws[m]
                  let coeff_k = unsafe coeff_ks[m]
                  in w * c64.re (coeff_k *! c64.exp (x *! minus_ikk)))
                 minus_ik maturity_for_quote)
       -- Writing this as a map-reduce requires way too much memory
       -- for compiler limitation reasons.
       loop (res = replicate nstrikes 0.0) = for j < n do map (+) res (iter j)
       in map (\moneyness resk m ->
               let day_count_fraction = unsafe day_count_fractions[m]
               let sigma_sqrtt = f64.sqrt (sigma2 day_count_fraction * day_count_fraction)
               let bs = bs_control moneyness sigma_sqrtt
               in if moneyness * f0 <= 0.0
                  then df * f64.max 0.0 (f0 *  (1.0 - moneyness))
                  else if moneyness < 0.0
                  then (let scale = if ap1 then f64.sqrt (-moneyness) else 1.0
                        in (-f0) * df * scale * resk / f64.pi +
                         bs + moneyness - 1.0)
                  else (let lb = f64.max 0.0 (1.0 - moneyness)
                        let scale = if ap1 then f64.sqrt moneyness else 1.0
                        in f0 * df * f64.max lb (f64.min 1.0 (scale * resk / f64.pi + bs))))
              moneyness res maturity_for_quote

fun gauss (x: f64) = f64.exp(-0.5 * x * x) / f64.sqrt(2. * f64.pi)

fun bs_call (call: bool) (today: date) (spot: f64) (strike: f64) (maturity: date) (vol: f64) =
  if same_date today maturity || vol <= 1e-15 then
    let forward = spot in
    let p = f64.max 0. (forward - strike) in
    (p, 0.)
  else
    let normal_dist (x: f64) = pnorm x in
    let eps = if call then 1. else -1. in
    let t = diff_dates today maturity in
    let sqrt_t = f64.sqrt t in
    let df_r = 1. in
    let df_d = 1. in
    let fwd = spot * df_d / df_r in
    let (d1, d2) =
      let d (add: f64) = (f64.log(fwd / strike) + add * 0.5 * vol * vol * t) / (vol * sqrt_t)
      in (d 1., d (-1.))
    in
    let (n1, n2) = (normal_dist (eps * d1), normal_dist (eps * d2)) in
    (eps * df_r * (fwd * n1 - strike * n2),
     spot * df_d * sqrt_t * gauss d1)

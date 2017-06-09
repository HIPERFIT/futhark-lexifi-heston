-- ==
-- compiled input @ data_1062_quotes
-- compiled input @ data_10000_quotes
-- compiled input @ data_100000_quotes

import "/futlib/random"
import "/futlib/math"
import "heston"

module heston64 = heston f64 minstd_rand u32_to_f64

let main (max_global: i32)
         (nb_points: i32)
         (np: i32)
         (today: i32)
         (quotes_maturity: [#num_quotes]i32)
         (quotes_strike: [#num_quotes]f64)
         (quotes_quote: [#num_quotes]f64) =
  heston64.heston max_global nb_points np today quotes_maturity quotes_strike quotes_quote

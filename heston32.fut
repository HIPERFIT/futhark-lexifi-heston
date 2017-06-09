-- ==
-- compiled input @ data_1062_quotes
-- compiled input @ data_10000_quotes
-- compiled input @ data_100000_quotes

import "/futlib/random"
import "/futlib/math"
import "heston"

module heston32 = heston f32 minstd_rand

-- We still read the data sets as double precision, and initially
-- convert them to single.  This is included in measurements, but
-- takes a negligible amount of time.
let main (max_global: i32)
         (nb_points: i32)
         (np: i32)
         (today: i32)
         (quotes_maturity: [#num_quotes]i32)
         (quotes_strike: [#num_quotes]f64)
         (quotes_quote: [#num_quotes]f64) =
  heston32.heston max_global nb_points np today quotes_maturity (map f32 quotes_strike) (map f32 quotes_quote)

-- ==
-- compiled input @ data_1062_quotes
-- compiled input @ data_10000_quotes
-- compiled input @ data_100000_quotes

import "futlib/date"
import "rand"
import "heston"

module heston32 = heston f32 random_f32

let main (max_global: i32)
         (nb_points: i32)
         (np: i32)
         (today: i32)
         (quotes_maturity: [#num_quotes]i32)
         (quotes_strike: [#num_quotes]f64)
         (quotes_quote: [#num_quotes]f64) =
  heston32.heston max_global nb_points np today quotes_maturity (map f32 quotes_strike) (map f32 quotes_quote)

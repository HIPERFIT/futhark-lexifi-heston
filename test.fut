fun main(dates: [n]i32) =
  let switched (x: i32) (i: i32) = i == 0 || x != dates[i-1]
  let switches = map switched dates (iota n)
  in (scan (+) 0 (map i32 switches),
      #2 (unzip (filter (\x -> #1 x) (zip switches dates))))

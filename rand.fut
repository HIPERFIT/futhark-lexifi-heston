-- Generate random floats.  Lazy and hackish.  Maybe don't use this for anything important.

module type random = {
  type t
  type rng

  val rng_from_seed: i32 -> rng
  val split_rng: i32 -> rng -> []rng
  val join_rng: []rng -> rng

  val rand: rng -> (t,t) -> (rng,t)
  val nrand: rng -> (t,t) -> i32 -> (rng, []t)
}

module random_i32: (random with t = i32) = {
  type t = i32
  type rng = i32

  -- From http://stackoverflow.com/a/12996028
  let hash(x: i32): i32 =
    let x = ((x >>> 16) ^ x) * 0x45d9f3b
    let x = ((x >>> 16) ^ x) * 0x45d9f3b
    let x = ((x >>> 16) ^ x) in
    x

  let rng_from_seed (x: i32) = x ^ 0b1010101010101

  let split_rng (n: i32) (x: rng): [n]rng =
    map (\i -> x ^ hash i) (iota n)

  let join_rng (xs: [#n]rng): rng =
    xs[0]

  let rand (x: i32) (min: i32, max: i32) =
    let x' = hash x
    in (x', min + x % max)

  let nrand (x: i32) (min: i32, max: i32) (n: i32): (i32, [n]i32) =
    let xs = split_rng n x
    let (xs', ys) = unzip (map (\x -> rand x (min,max)) xs)
    in (join_rng xs', ys)

}

module random_f64: (random with t = f64 with rng = random_i32.rng) = {
  type t = f64
  type rng = random_i32.rng

  let rng_from_seed (x: i32) = random_i32.rng_from_seed x
  let split_rng (n: i32) (x: rng) = random_i32.split_rng n x
  let join_rng (xs: []rng) = random_i32.join_rng xs

  let rand (rng: rng) (min: f64, max: f64) =
    let (rng', x) = random_i32.rand rng (0, 0x7FFFFFFF)
    let y = f64 x / f64 0x7FFFFFFF -- range [0,1)
    in (rng', min + y * (max-min))

  let nrand (rng: rng) (min: f64, max: f64) (n: i32): (rng, [n]f64) =
    let (rng', xs) = random_i32.nrand rng (0, 0x7FFFFFFF) n
    let xs' = map (\x -> let y = f64 x / f64 0x7FFFFFFF
                         in min + y * (max-min))
                  xs
    in (rng', xs')
}

module random_f32: (random with t = f32 with rng = random_i32.rng) = {
  type t = f32
  type rng = random_i32.rng

  let rng_from_seed (x: i32) = random_i32.rng_from_seed x
  let split_rng (n: i32) (x: rng) = random_i32.split_rng n x
  let join_rng (xs: []rng) = random_i32.join_rng xs

  let rand (rng: rng) (min: f32, max: f32) =
    let (rng', x) = random_i32.rand rng (0, 0x7FFFFFFF)
    let y = f32 x / f32 0x7FFFFFFF -- range [0,1)
    in (rng', min + y * (max-min))

  let nrand (rng: rng) (min: f32, max: f32) (n: i32): (rng, [n]f32) =
    let (rng', xs) = random_i32.nrand rng (0, 0x7FFFFFFF) n
    let xs' = map (\x -> let y = f32 x / f32 0x7FFFFFFF
                         in min + y * (max-min))
                  xs
    in (rng', xs')
}

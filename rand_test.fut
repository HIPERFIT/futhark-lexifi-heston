import "rand"

fun main(x: i32, n: i32, min: f64, max: f64) =
  let seqout = replicate n 0.0
  let rng = random_f64.rng_from_seed x
  let (_, parout) = random_f64.nrand rng (min,max) n
  loop ((seqout,rng)) = for i < n do (let (rng', x) = random_f64.rand rng (min,max)
                                   let seqout[i] = x
                                   in (seqout,rng'))
  in (reduce (+) 0.0 seqout / f64 n,
      reduce (+) 0.0 parout / f64 n)

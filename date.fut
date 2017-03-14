-- A simple date library

module type date = {
  type date
  val sub_act_365: date -> date -> f64
  val add_act_365: date -> f64 -> date
  val add_days: date -> i32 -> date
  val sub_days: date -> i32 -> date
  val days_between: date -> date -> f64
  val diff_dates: date -> date -> f64

  
  val triple_of_date: date -> (i32,i32,i32)
  val date_of_triple: (i32,i32,i32) -> date

  val same_date: date -> date -> bool
}

open ({
  type date = i32

  type gregorian = { year: i32,
                     month: i32,
                     day: i32,
                     hour: i32,
                     minute: i32 }

  val hours_in_day = 24
  val minutes_in_day = hours_in_day * 60
  val fminutes_in_day = f64 minutes_in_day
  val minutes_to_noon = (hours_in_day / 2) * 60

  fun date_of_gregorian ({year = y, month = m, day = d, hour = hr, minute = mn}: gregorian) =
    ((if m == 1 || m == 2 then
       ( 1461 * ( y + 4800 - 1 ) ) / 4 +
         ( 367 * ( m + 10 ) ) / 12 -
         ( 3 * ( ( y + 4900 - 1 ) / 100 ) ) / 4
      else
       ( 1461 * ( y + 4800 ) ) / 4 +
         ( 367 * ( m - 2 ) ) / 12 -
         ( 3 * ( ( y + 4900 ) / 100 ) ) / 4) + d - 32075 - 2444238) * minutes_in_day
  + hr * 60 + mn

  fun gregorian_of_date (minutes_since_epoch: i32) =
    let jul = minutes_since_epoch / minutes_in_day in
    let l = jul + 68569 + 2444238 in
    let n = ( 4 * l ) / 146097 in
    let l = l - ( 146097 * n + 3 ) / 4 in
    let i = ( 4000 * ( l + 1 ) ) / 1461001 in
    let l = l - ( 1461 * i ) / 4 + 31 in
    let j = ( 80 * l ) / 2447 in
    let d = l - ( 2447 * j ) / 80 in
    let l = j / 11 in
    let m = j + 2 - ( 12 * l ) in
    let y = 100 * ( n - 49 ) + i + l in
    let daytime = minutes_since_epoch % minutes_in_day in
    if daytime == minutes_to_noon
    then {year = y, month = m, day = d, hour = 12, minute = 0}
    else {year = y, month = m, day = d, hour = daytime / 60, minute = daytime % 60}

  fun check_date (year: i32) (month: i32) (day: i32) =
    1 <= day &&
    1 <= month && month <= 12 &&
    1980 <= year && year <= 2299 &&
    (day <= 28 ||
     if month == 2 then
     day == 29 && year % 4 == 0 && (year == 2000 || (year % 100 != 0))
     else if month == 4 || month == 6 || month == 9 || month == 11
     then day <= 30
     else day <= 1)

  fun date_of_triple (year: i32, month: i32, day: i32) =
    date_of_gregorian {year=year, month=month, day=day, hour=12, minute=0}

  fun triple_of_date (x: date) =
    let {year, month, day, hour = _, minute = _} = gregorian_of_date x
    in (year, month, day)

  fun int_of_date (x: date) = x
  fun date_of_int (x: i32) = x

  val fminutes_in_365 = f64 (minutes_in_day * 365)
  val inv_fminutes_in_365 = 1.0 / fminutes_in_365
  val inv_fminutes_in_day = 1.0 / fminutes_in_day

  fun sub_act_365 (t1: date) (t2: date) =
    f64 (int_of_date t1 - int_of_date t2) * inv_fminutes_in_365
  fun add_act_365 (t: date) (dt: f64) =
    date_of_int (i32 (f64 (int_of_date t) + fminutes_in_365 * dt))
  fun add_days (t1: date) (displ: i32) =
    date_of_int(int_of_date t1 + displ * minutes_in_day)
  fun sub_days (t1: date) (displ: i32) =
    date_of_int(int_of_date t1 - displ * minutes_in_day)
  fun days_between (t1: date) (t2: date) =
    (f64 (int_of_date t2 - int_of_date t1)) * inv_fminutes_in_day
  fun diff_dates (t1: date) (t2: date) =
    (days_between t1 t2) / 365.0

  fun same_date (x: date) (y: date) = x == y

} : date)

test_that("change_476_type_ids_to_open_intervals replaces :R(5,9) suffix", {
  input <- c("Del(C):Ins(C):R(5,9)", "Ins(C):R(5,9)")
  out <- change_476_type_ids_to_open_intervals(input)
  expect_equal(out, c("Del(C):Ins(C):R(5,)", "Ins(C):R(5,)"))
})

test_that("change_476_type_ids_to_open_intervals leaves non-matching strings unchanged", {
  input <- c("Del(T):R(3,5)", "Ins(C):R(5,9):extra", "Del(C):R(5,9)x")
  out <- change_476_type_ids_to_open_intervals(input)
  expect_equal(out, input)
})

test_that("change_476_type_ids_to_open_intervals handles empty input", {
  expect_equal(change_476_type_ids_to_open_intervals(character(0)), character(0))
})

test_that("change_89_type_ids_to_open_intervals replaces all bounded repeat suffixes", {
  input <- c(
    "Ins(2,):R(5,9)",
    "Ins(C):R(7,9)",
    "[Del(T):R(8,9)]",
    "[Ins(T):R(8,9)]",
    "Del(2,):U(1,2):R(5,9)",
    "Del(3,):U(3,):R(3,9)"
  )
  expected <- c(
    "Ins(2,):R(5,)",
    "Ins(C):R(7,)",
    "[Del(T):R(8,)]",
    "[Ins(T):R(8,)]",
    "Del(2,):U(1,2):R(5,)",
    "Del(3,):U(3,):R(3,)"
  )
  expect_equal(change_89_type_ids_to_open_intervals(input), expected)
})

test_that("change_89_type_ids_to_open_intervals leaves non-matching strings unchanged", {
  input <- c("Del(T):R(3,5)", "Ins(C):R(5,8)", "Del(3,):U(3,):R(3,)")
  expect_equal(change_89_type_ids_to_open_intervals(input), input)
})

test_that("change_89_type_ids_to_open_intervals handles empty input", {
  expect_equal(change_89_type_ids_to_open_intervals(character(0)), character(0))
})

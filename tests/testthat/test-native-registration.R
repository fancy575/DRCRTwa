test_that("the compiled fitting routine is loaded and registered", {
  info <- getNativeSymbolInfo(
    "_DRCRTwa_drcrtwa_fit_cpp",
    PACKAGE = "DRCRTwa",
    withRegistrationInfo = TRUE
  )
  expect_s3_class(info, "NativeSymbolInfo")
  expect_identical(info$numParameters, 20L)
})

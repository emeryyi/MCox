context("Projections.f90")

test_that("ProjectL1 gives correct results",{
    expect_equal(ProjectL1(10, 1), 1)
    expect_equal(ProjectL1(-10, 1), -1)
    expect_equal(ProjectL1(0, 1), 0)
    expect_equal(ProjectL1(0:10, 1), c(rep(0,10),1))
    expect_equal(ProjectL1(0:4-2, 1), c(-0.5,0,0,0,0.5))
})
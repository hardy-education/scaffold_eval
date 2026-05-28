#' Tidy variance components from an lme4 fit.
#'
#' @param fit lmerMod or glmerMod object.
#' @param add_glmer_residual If TRUE, append binomial residual variance pi^2/3.
varcorr_summary <- function(fit, add_glmer_residual = FALSE) {
  vc <- fit |>
    lme4::VarCorr() |>
    janitor::clean_names() |>
    tibble::as_tibble() |>
    dplyr::filter(is.na(var2))

  if (add_glmer_residual && !any(vc$grp == "Residual", na.rm = TRUE)) {
    vc <- dplyr::bind_rows(
      vc,
      tibble::tibble(
        grp = "Residual",
        var1 = NA_character_,
        var2 = NA_character_,
        vcov = pi^2 / 3,
        sdcor = NA_real_
      )
    )
  }

  vc |>
    dplyr::mutate(
      pct = round(vcov / sum(vcov) * 100, 2),
      agent = stringr::str_detect(grp, "agent"),
      bench = stringr::str_detect(grp, "benchmark"),
      model = stringr::str_detect(grp, "model"),
      task = stringr::str_detect(grp, "task"),
      run = stringr::str_detect(grp, "run"),
      resid = grp == "Residual",
      model_agent = model & agent,
      dplyr::across(agent:model_agent, \(x) vcov * as.numeric(x))
    )
}

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(splines)
})

# ============================================================
# Glioblastoma, IDH-wildtype Fixed-Horizon Mortality Estimator
# Shiny app for the corrected 2024 NCDB Brain PUF model
#
# Expected deployment structure:
#   shiny/
#   ├── app.R
#   ├── manifest.json
#   ├── data/
#   │   └── processed/
#   │       └── gbm_model_objects.rds
#   └── www/
#       └── ohsu_logo.png
# ============================================================

# ----------------------------
# 1. Load model object
# ----------------------------
model_object_candidates <- c(
  file.path("data", "processed", "gbm_model_objects.rds"),
  "gbm_model_objects.rds",
  file.path("..", "data", "processed", "gbm_model_objects.rds")
)

model_object_path <- model_object_candidates[file.exists(model_object_candidates)][1]

if (is.na(model_object_path)) {
  stop(
    paste0(
      "Could not find gbm_model_objects.rds.\n\n",
      "Searched these paths:\n",
      paste(
        normalizePath(model_object_candidates, winslash = "/", mustWork = FALSE),
        collapse = "\n"
      ),
      "\n\nPlace the current model object at shiny/data/processed/gbm_model_objects.rds."
    ),
    call. = FALSE
  )
}

obj <- readRDS(model_object_path)
models <- obj$models

if (is.null(models) || length(models) == 0) {
  stop("gbm_model_objects.rds does not contain obj$models.", call. = FALSE)
}

if (is.null(obj$df_ref) || is.null(obj$modal_profile)) {
  stop(
    "gbm_model_objects.rds must contain obj$df_ref and obj$modal_profile.",
    call. = FALSE
  )
}

df_ref <- as.data.frame(obj$df_ref)
modal_profile <- as.data.frame(obj$modal_profile)

horizons <- suppressWarnings(as.integer(names(models)))
horizons <- horizons[is.finite(horizons)]
horizons <- sort(horizons)
models <- models[as.character(horizons)]

expected_horizons <- c(6L, 12L, 18L, 24L, 36L)
if (!identical(horizons, expected_horizons)) {
  stop(
    paste0(
      "The model object must contain 6-, 12-, 18-, 24-, and 36-month models. ",
      "Found: ", paste(horizons, collapse = ", "), "."
    ),
    call. = FALSE
  )
}

model_n <- if (!is.null(obj$model_n)) as.integer(obj$model_n) else nrow(df_ref)

# Current manuscript/analysis constants from the corrected 2024 PUF pipeline.
analytic_n <- 24976L
analytic_deaths <- 20336L
training_n <- 17483L
validation_n <- 7493L

validation_auc <- c(
  `6` = 0.849,
  `12` = 0.796,
  `18` = 0.782,
  `24` = 0.781,
  `36` = 0.809
)
validation_brier <- c(
  `6` = 0.128,
  `12` = 0.183,
  `18` = 0.176,
  `24` = 0.149,
  `36` = 0.104
)
calibration_slope <- c(
  `6` = 0.971,
  `12` = 0.983,
  `18` = 1.031,
  `24` = 1.035,
  `36` = 0.949
)

# ----------------------------
# 2. Utility functions
# ----------------------------
clamp_num <- function(x, lower, upper, default) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) x <- default
  min(max(x, lower), upper)
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  tab <- sort(table(as.character(x)), decreasing = TRUE)
  names(tab)[1]
}

safe_levels <- function(var) {
  if (!var %in% names(df_ref)) return(character(0))
  x <- df_ref[[var]]
  if (is.factor(x)) return(levels(x))
  sort(unique(as.character(x[!is.na(x)])))
}

select_default <- function(var, fallback = NULL) {
  if (!var %in% names(df_ref)) return(fallback)
  mv <- mode_value(df_ref[[var]])
  lv <- safe_levels(var)
  if (!is.na(mv) && mv %in% lv) return(mv)
  if (length(lv) > 0) return(lv[1])
  fallback
}

make_choices <- function(var, labels = NULL) {
  lv <- safe_levels(var)
  if (length(lv) == 0) return(character(0))
  if (is.null(labels)) return(stats::setNames(lv, lv))

  labels <- labels[lv]
  missing_labels <- is.na(labels)
  labels[missing_labels] <- lv[missing_labels]
  stats::setNames(lv, labels)
}

apply_ref_classes <- function(nd) {
  for (nm in intersect(names(nd), names(df_ref))) {
    if (is.factor(df_ref[[nm]])) {
      nd[[nm]] <- factor(as.character(nd[[nm]]), levels = levels(df_ref[[nm]]))
    } else if (is.numeric(df_ref[[nm]]) || is.integer(df_ref[[nm]])) {
      nd[[nm]] <- suppressWarnings(as.numeric(nd[[nm]]))
    } else {
      nd[[nm]] <- as.character(nd[[nm]])
    }
  }
  nd
}

fmt_pct <- function(x) {
  ifelse(is.na(x), "—", sprintf("%.1f%%", 100 * x))
}

model_predictors <- unique(unlist(lapply(models, function(fit) {
  tryCatch(
    all.vars(stats::delete.response(stats::terms(fit))),
    error = function(e) character(0)
  )
})))

has_removed_facility_predictors <- any(
  c("facility_type_cd", "facility_location_cd") %in% model_predictors
)

model_object_is_current <- isTRUE(model_n == training_n) &&
  "2023" %in% safe_levels("dx_year") &&
  !has_removed_facility_predictors

# ----------------------------
# 3. Defaults and choices
# ----------------------------
age_min <- 18
age_max <- 90
age_default <- clamp_num(
  round(median(df_ref$age_years, na.rm = TRUE)),
  age_min,
  age_max,
  65
)

tumor_min <- 1
tumor_max <- 200
tumor_default <- clamp_num(
  round(median(df_ref$tumor_size_harmonized_mm, na.rm = TRUE)),
  tumor_min,
  tumor_max,
  44
)

sex_choices <- make_choices("sex_cat")
race_choices <- make_choices(
  "race_cat",
  c("White" = "White", "Black" = "Black", "Other" = "Other race")
)
ethnicity_choices <- make_choices("ethnicity_cat")
insurance_choices <- make_choices(
  "insurance_cat",
  c(
    "Private" = "Private insurance",
    "Medicare" = "Medicare",
    "Medicaid" = "Medicaid",
    "Not insured" = "Not insured"
  )
)
income_choices <- make_choices(
  "income_quartile",
  c(
    "Q4 highest income" = "Q4, highest area-level income",
    "Q3" = "Q3",
    "Q2" = "Q2",
    "Q1 lowest income" = "Q1, lowest area-level income"
  )
)
education_choices <- make_choices(
  "education_quartile",
  c(
    "Q4 highest education" = "Q4, highest area-level education",
    "Q3" = "Q3",
    "Q2" = "Q2",
    "Q1 lowest education" = "Q1, lowest area-level education"
  )
)
charlson_choices <- make_choices(
  "charlson_deyo_cat",
  c("0" = "0", "1" = "1", "2+" = "2 or more")
)
mgmt_choices <- make_choices(
  "mgmt_status",
  c("Unmethylated" = "Unmethylated", "Methylated" = "Methylated")
)
site_choices <- make_choices("primary_site_group")
surgery_choices <- make_choices(
  "surgery_extent",
  c(
    "No surgery" = "No surgery",
    "STR" = "Subtotal resection (STR)",
    "GTR" = "Gross total resection (GTR)"
  )
)
radiation_choices <- make_choices("radiation_status")
chemo_choices <- make_choices("chemo_status")
year_choices <- make_choices("dx_year")
year_default <- if ("2023" %in% unname(year_choices)) {
  "2023"
} else {
  select_default("dx_year")
}

logo_ui <- if (file.exists(file.path("www", "ohsu_logo.png"))) {
  img(src = "ohsu_logo.png", class = "ohsu-logo")
} else {
  div("OHSU", class = "ohsu-logo-fallback")
}

model_status_ui <- if (!model_object_is_current) {
  div(
    class = "model-warning",
    tags$strong("Model-object update required. "),
    paste0(
      "The interface reflects the corrected 2024-PUF analysis, but the loaded ",
      "gbm_model_objects.rds does not match the final derivation cohort. Replace ",
      "shiny/data/processed/gbm_model_objects.rds before clinical demonstration."
    )
  )
} else {
  NULL
}

# ----------------------------
# 4. UI
# ----------------------------
ui <- page_fluid(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter"),
    primary = "#1f4e79",
    bg = "#f4f7fb",
    fg = "#243447"
  ),

  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML("
      :root {
        --page-max: 1320px;
        --card-radius: 24px;
        --shadow-soft: 0 8px 28px rgba(31, 52, 73, 0.07);
        --border-soft: #e7edf5;
        --text-main: #243447;
        --text-muted: #5b6b7f;
        --bg-soft: #f4f7fb;
        --accent: #1f4e79;
      }
      body { background: var(--bg-soft); }
      .app-container { max-width: var(--page-max); margin: 0 auto; padding: 24px 22px 36px 22px; }
      .app-header { background: #ffffff; border-radius: 28px; padding: clamp(18px, 2.2vw, 30px); margin-bottom: 18px; box-shadow: var(--shadow-soft); border: 1px solid var(--border-soft); }
      .header-grid { display: grid; grid-template-columns: minmax(70px, 96px) 1fr; gap: 20px; align-items: center; }
      .logo-wrap { display: flex; align-items: center; justify-content: center; }
      .ohsu-logo { width: clamp(58px, 6vw, 92px); height: auto; display: block; }
      .ohsu-logo-fallback { width: 88px; height: 88px; border-radius: 22px; display: flex; align-items: center; justify-content: center; background: #1f4e79; color: white; font-weight: 900; letter-spacing: 0.06em; }
      .header-title { margin: 0 0 8px 0; font-weight: 800; line-height: 1.04; font-size: clamp(1.9rem, 3.4vw, 3.2rem); color: var(--text-main); max-width: 1000px; }
      .ohsu-subtitle { color: var(--text-muted); margin: 0 0 3px 0; font-size: 1.05rem; }
      .ohsu-dept { color: #738396; margin: 0; font-size: 0.98rem; }
      .model-warning { max-width: var(--page-max); margin: 0 auto 18px auto; padding: 14px 18px; border-radius: 14px; border: 1px solid #e6b800; background: #fff8d8; color: #594600; }
      .input-card, .metric-card, .plot-card, .detail-card { background: #ffffff; border: 1px solid var(--border-soft) !important; border-radius: var(--card-radius) !important; box-shadow: var(--shadow-soft); }
      .metric-card .card-body, .plot-card .card-body, .detail-card .card-body { padding: 22px; }
      .input-card .card-body { padding: 18px 18px 16px 18px; }
      .sticky-panel { position: sticky; top: 24px; max-height: calc(100vh - 48px); overflow-y: auto; padding-right: 4px; scrollbar-width: thin; }
      .section-title { font-weight: 800; color: var(--text-main); margin-bottom: 14px; line-height: 1.06; font-size: clamp(1.55rem, 2vw, 2rem); }
      .plot-title { font-weight: 800; color: var(--text-main); margin-bottom: 10px; font-size: 1.15rem; }
      .form-label { font-weight: 650; color: #2f4257; margin-bottom: 5px; font-size: 0.97rem; }
      .shiny-input-container { margin-bottom: 10px; }
      .form-control, .form-select { border-radius: 14px !important; border: 1px solid #d4dde8 !important; min-height: 44px; box-shadow: none !important; }
      .btn-primary { background-color: #245789 !important; border-color: #245789 !important; border-radius: 14px !important; font-weight: 750; min-height: 46px; margin-top: 6px; }
      .input-note { margin: 12px 2px 2px 2px; color: #65758a; font-size: 0.86rem; line-height: 1.4; }
      .metric-grid { display: grid; grid-template-columns: repeat(5, minmax(0, 1fr)); gap: 14px; margin-bottom: 18px; }
      .metric-card { min-height: 118px; }
      .metric-value { font-size: clamp(1.45rem, 2vw, 2.05rem); line-height: 1; font-weight: 800; color: var(--accent); margin-bottom: 10px; }
      .metric-label { font-size: 0.92rem; color: var(--text-muted); line-height: 1.35; }
      .detail-card h3 { font-size: 1.08rem; font-weight: 750; color: var(--text-main); margin-top: 0; margin-bottom: 0.8rem; }
      .detail-card ul { margin-bottom: 0; padding-left: 1.15rem; }
      .detail-card li { color: #425466; margin-bottom: 0.48rem; line-height: 1.5; }
      .block-gap { height: 18px; }
      .detail-grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 26px 38px; }
      @media (max-width: 1199px) {
        .metric-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); }
        .sticky-panel { position: static; max-height: none; overflow-y: visible; padding-right: 0; }
      }
      @media (max-width: 767px) {
        .app-container { padding: 18px 14px 28px 14px; }
        .header-grid, .metric-grid, .detail-grid { grid-template-columns: 1fr; }
        .header-grid { text-align: center; }
        .plot-card .shiny-plot-output { height: 420px !important; }
      }
    "))
  ),

  div(
    class = "app-container",

    div(
      class = "app-header",
      div(
        class = "header-grid",
        div(class = "logo-wrap", logo_ui),
        div(
          h1("Glioblastoma, IDH-wildtype Mortality Risk Estimator", class = "header-title"),
          p("Oregon Health & Science University", class = "ohsu-subtitle"),
          p("Department of Neurological Surgery", class = "ohsu-dept")
        )
      )
    ),

    model_status_ui,

    layout_columns(
      col_widths = c(4, 8),

      div(
        class = "sticky-panel",
        card(
          class = "input-card",
          card_body(
            h2("Patient and treatment characteristics", class = "section-title"),
            numericInput(
              "age",
              "Age at diagnosis (years)",
              value = age_default,
              min = age_min,
              max = age_max,
              step = 1
            ),
            selectInput("sex", "Sex", choices = sex_choices, selected = select_default("sex_cat"), selectize = FALSE),
            selectInput("race", "Race", choices = race_choices, selected = select_default("race_cat"), selectize = FALSE),
            selectInput("ethnicity", "Ethnicity", choices = ethnicity_choices, selected = select_default("ethnicity_cat"), selectize = FALSE),
            selectInput("insurance", "Insurance", choices = insurance_choices, selected = select_default("insurance_cat"), selectize = FALSE),
            selectInput("income", "Area-level income quartile", choices = income_choices, selected = select_default("income_quartile"), selectize = FALSE),
            selectInput("education", "Area-level education quartile", choices = education_choices, selected = select_default("education_quartile"), selectize = FALSE),
            selectInput("cdcc", "Charlson-Deyo comorbidity score", choices = charlson_choices, selected = select_default("charlson_deyo_cat"), selectize = FALSE),
            numericInput(
              "tsize_mm",
              "Tumor size (mm)",
              value = tumor_default,
              min = tumor_min,
              max = tumor_max,
              step = 1
            ),
            selectInput("mgmt", "MGMT promoter methylation", choices = mgmt_choices, selected = select_default("mgmt_status"), selectize = FALSE),
            selectInput("site", "Primary site group", choices = site_choices, selected = select_default("primary_site_group"), selectize = FALSE),
            selectInput("surgery", "Surgery extent", choices = surgery_choices, selected = select_default("surgery_extent"), selectize = FALSE),
            selectInput("radiation", "Radiation therapy", choices = radiation_choices, selected = select_default("radiation_status"), selectize = FALSE),
            selectInput("chemo", "Chemotherapy", choices = chemo_choices, selected = select_default("chemo_status"), selectize = FALSE),
            if (length(year_choices) > 0) {
              selectInput(
                "dx_year",
                "Diagnosis year represented in model",
                choices = year_choices,
                selected = year_default,
                selectize = FALSE
              )
            },
            actionButton("calc", "Estimate mortality risk", class = "btn-primary w-100"),
            div(
              class = "input-note",
              "Treatment entries describe recorded or planned first-course care. Changing them does not estimate the causal benefit of a treatment."
            )
          )
        )
      ),

      div(
        div(
          class = "metric-grid",
          lapply(horizons, function(h) {
            card(
              class = "metric-card",
              card_body(
                div(textOutput(paste0("risk_", h)), class = "metric-value"),
                div(paste0(h, "-month mortality risk"), class = "metric-label")
              )
            )
          })
        ),
        div(class = "block-gap"),
        card(
          class = "plot-card",
          card_body(
            h2("Estimated fixed-horizon mortality risk", class = "plot-title"),
            plotOutput("riskplot", height = "560px")
          )
        )
      )
    ),

    div(class = "block-gap"),

    card(
      class = "detail-card",
      card_body(
        h2("Model details, validation, and intended use", class = "section-title"),
        div(
          class = "detail-grid",
          div(
            h3("Cohort and intended use"),
            tags$ul(
              tags$li(
                paste0(
                  "The complete-case analytic cohort included ",
                  format(analytic_n, big.mark = ","),
                  " adults with molecularly defined glioblastoma, IDH-wildtype; ",
                  format(analytic_deaths, big.mark = ","),
                  " deaths were observed."
                )
              ),
              tags$li(
                paste0(
                  "Models were developed in ", format(training_n, big.mark = ","),
                  " patients and internally evaluated in ",
                  format(validation_n, big.mark = ","), " patients."
                )
              ),
              tags$li("The calculator estimates all-cause mortality risk at 6, 12, 18, 24, and 36 months after diagnosis."),
              tags$li("It is intended to support prognostic counseling and risk communication, not to replace multidisciplinary clinical judgment.")
            )
          ),
          div(
            h3("Data and predictors"),
            tags$ul(
              tags$li("Data source: 2024 National Cancer Database Brain Participant User File."),
              tags$li("Study population: adults diagnosed from 2018 through 2023 with primary brain-site glioblastoma histology, microscopic confirmation, and Brain Molecular Marker 05 consistent with IDH-wildtype glioblastoma."),
              tags$li("Predictors include age, sex, race, ethnicity, insurance, area-level income and education, Charlson-Deyo score, tumor size, MGMT promoter methylation, primary site group, surgery extent, radiation, chemotherapy, and diagnosis year."),
              tags$li("Diagnosis years after 2023 were not represented during model development and require cautious interpretation.")
            )
          ),
          div(
            h3("Internal validation"),
            tags$ul(
              tags$li("Separate nonlinear logistic models were fit for each horizon using inverse probability of censoring weighting."),
              tags$li("Age and tumor size were modeled with natural splines, with clinically selected treatment, MGMT, surgery, and age interactions."),
              tags$li(
                paste0(
                  "Validation AUCs at 6, 12, 18, 24, and 36 months were ",
                  paste(sprintf("%.3f", validation_auc), collapse = ", "), ", respectively."
                )
              ),
              tags$li(
                paste0(
                  "Calibration slopes were ",
                  paste(sprintf("%.3f", calibration_slope), collapse = ", "),
                  "; Brier scores were ",
                  paste(sprintf("%.3f", validation_brier), collapse = ", "), "."
                )
              ),
              tags$li("Confidence intervals were estimated from 1,000 bootstrap resamples of the full training and validation pipeline.")
            )
          ),
          div(
            h3("Interpretation and limitations"),
            tags$ul(
              tags$li("Treatment variables describe first-course treatment and must not be interpreted as causal or counterfactual treatment effects."),
              tags$li("The registry does not include Karnofsky Performance Status, neurologic deficits, postoperative residual tumor volume, recurrence, treatment completion, tumor treating fields, or longitudinal treatment changes."),
              tags$li("The complete-case design may limit transportability, and the NCDB represents Commission on Cancer-accredited facilities rather than the entire United States population."),
              tags$li("Independent external validation and, where necessary, recalibration are required before routine clinical implementation.")
            )
          )
        )
      )
    )
  )
)

# ----------------------------
# 5. Server
# ----------------------------
server <- function(input, output, session) {
  observe({
    current_age <- suppressWarnings(as.numeric(input$age))
    if (!is.na(current_age) && current_age > age_max) {
      updateNumericInput(session, "age", value = age_max)
    }
    if (!is.na(current_age) && current_age < age_min) {
      updateNumericInput(session, "age", value = age_min)
    }
  })

  observe({
    current_size <- suppressWarnings(as.numeric(input$tsize_mm))
    if (!is.na(current_size) && current_size > tumor_max) {
      updateNumericInput(session, "tsize_mm", value = tumor_max)
    }
    if (!is.na(current_size) && current_size < tumor_min) {
      updateNumericInput(session, "tsize_mm", value = tumor_min)
    }
  })

  newdata <- eventReactive(input$calc, {
    nd <- modal_profile[1, , drop = FALSE]

    if ("age_years" %in% names(nd)) {
      nd$age_years <- clamp_num(input$age, age_min, age_max, age_default)
    }
    if ("tumor_size_harmonized_mm" %in% names(nd)) {
      nd$tumor_size_harmonized_mm <- clamp_num(
        input$tsize_mm,
        tumor_min,
        tumor_max,
        tumor_default
      )
    }
    if ("sex_cat" %in% names(nd)) nd$sex_cat <- input$sex
    if ("race_cat" %in% names(nd)) nd$race_cat <- input$race
    if ("ethnicity_cat" %in% names(nd)) nd$ethnicity_cat <- input$ethnicity
    if ("insurance_cat" %in% names(nd)) nd$insurance_cat <- input$insurance
    if ("income_quartile" %in% names(nd)) nd$income_quartile <- input$income
    if ("education_quartile" %in% names(nd)) nd$education_quartile <- input$education
    if ("charlson_deyo_cat" %in% names(nd)) nd$charlson_deyo_cat <- input$cdcc
    if ("mgmt_status" %in% names(nd)) nd$mgmt_status <- input$mgmt
    if ("primary_site_group" %in% names(nd)) nd$primary_site_group <- input$site
    if ("surgery_extent" %in% names(nd)) nd$surgery_extent <- input$surgery
    if ("radiation_status" %in% names(nd)) nd$radiation_status <- input$radiation
    if ("chemo_status" %in% names(nd)) nd$chemo_status <- input$chemo
    if ("dx_year" %in% names(nd) && !is.null(input$dx_year)) {
      nd$dx_year <- input$dx_year
    }

    apply_ref_classes(nd)
  }, ignoreNULL = FALSE)

  risk_data <- eventReactive(input$calc, {
    req(newdata())

    risks <- vapply(names(models), function(h) {
      p <- tryCatch(
        as.numeric(predict(models[[h]], newdata = newdata(), type = "response")),
        error = function(e) NA_real_
      )
      pmin(pmax(p, 0), 1)
    }, numeric(1))

    data.frame(
      horizon_months = horizons,
      risk = as.numeric(risks)
    )
  }, ignoreNULL = FALSE)

  lapply(horizons, function(h) {
    output[[paste0("risk_", h)]] <- renderText({
      d <- risk_data()
      fmt_pct(d$risk[d$horizon_months == h][1])
    })
  })

  output$riskplot <- renderPlot({
    d <- risk_data()

    ggplot(d, aes(x = horizon_months, y = risk)) +
      geom_line(linewidth = 1.4, color = "#1f6feb") +
      geom_point(size = 3.2, color = "#1f6feb") +
      scale_x_continuous(
        breaks = horizons,
        limits = c(min(horizons), max(horizons))
      ) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        labels = function(x) paste0(round(100 * x), "%")
      ) +
      labs(
        x = "Months from diagnosis",
        y = "Predicted mortality risk"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#e7edf5", linewidth = 0.7),
        axis.title = element_text(color = "#2f4257", face = "bold"),
        axis.text = element_text(color = "#425466"),
        plot.background = element_rect(fill = "#ffffff", color = NA),
        panel.background = element_rect(fill = "#ffffff", color = NA),
        plot.margin = margin(10, 10, 8, 8)
      )
  }, res = 120)
}

shinyApp(ui, server)

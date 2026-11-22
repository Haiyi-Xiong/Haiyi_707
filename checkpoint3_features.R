suppressPackageStartupMessages({library(data.table);library(dplyr);library(tidyr)})
cp_build_table <- function(set_name) {
# 2.1 Build the modeling table
# Uses set-a only
set.seed(707)

# 1. Define relative paths
data_dir <- "data"
out_dir <- "output"

dir.create(out_dir, showWarnings = FALSE)

# 2. Read the ICU records
files <- list.files(
  file.path(data_dir, paste0("set-", set_name)),
  pattern = "\\.txt$",
  full.names = TRUE
)

stopifnot(length(files) == 4000)

raw <- rbindlist(lapply(files, function(f) {
  x <- fread(f)
  
  stopifnot(
    all(c("Time", "Parameter", "Value") %in% names(x))
  )
  
  x[, RecordID := as.integer(
    tools::file_path_sans_ext(basename(f))
  )]
  
  x
}), fill = TRUE)

# 3. Convert measurement values to numeric
raw[, Value := as.numeric(Value)]

# Treat -1 as missing
raw[Value == -1, Value := NA_real_]

# 4. Convert time into ICU hours
time_parts <- tstrsplit(raw$Time, ":", fixed = TRUE)

raw[, ICU_hour :=
      as.numeric(time_parts[[1]]) +
      as.numeric(time_parts[[2]]) / 60]

# Restrict measurements to hours 0–48
raw <- raw[
  !is.na(ICU_hour) &
    ICU_hour >= 0 &
    ICU_hour <= 48
]

# 5. Extract static patient characteristics
static_vars <- c(
  "Age", "Gender", "Height", "ICUType", "Weight"
)

static <- raw[
  Time == "00:00" & Parameter %chin% static_vars,
  .(Value = Value[1]),
  by = .(RecordID, Parameter)
]

static <- dcast(
  static,
  RecordID ~ Parameter,
  value.var = "Value"
)

# Height and weight cannot be zero
static[!is.na(Height) & Height <= 0,
       Height := NA_real_]

static[!is.na(Weight) & Weight <= 0,
       Weight := NA_real_]

# 6. Extract time-varying measurements
measurements <- raw[
  !Parameter %chin% c(static_vars, "RecordID") &
    !is.na(Value)
]

setorder(
  measurements,
  RecordID,
  Parameter,
  ICU_hour
)

# 7. Summarize each measurement
summaries <- measurements[
  ,
  .(
    count = .N,
    first = Value[1],
    last = Value[.N],
    min = min(Value),
    max = max(Value),
    mean = mean(Value)
  ),
  by = .(RecordID, Parameter)
]

# 8. Convert long data to wide data
wide <- as_tibble(summaries) |>
  pivot_wider(
    id_cols = RecordID,
    names_from = Parameter,
    values_from = c(
      count, first, last, min, max, mean
    ),
    names_glue = "{Parameter}_{.value}"
  )

# 9. Replace missing measurement counts with zero
count_cols <- grep(
  "_count$",
  names(wide),
  value = TRUE
)

wide <- wide |>
  mutate(
    across(
      all_of(count_cols),
      ~ replace_na(.x, 0L)
    )
  )

# 10. Read outcome data
outcomes <- fread(
  file.path(data_dir, paste0("Outcomes-", set_name, ".txt"))
)

stopifnot(
  nrow(outcomes) == 4000,
  anyDuplicated(outcomes$RecordID) == 0
)

# 11. Join predictors and outcomes
model_table <- as_tibble(static) |>
  left_join(wide, by = "RecordID") |>
  left_join(
    as_tibble(outcomes) |>
      select(
        RecordID,
        `In-hospital_death`,
        `SAPS-I`,
        SOFA
      ),
    by = "RecordID"
  ) |>
  rename(
    death = `In-hospital_death`
  ) |>
  mutate(
    death = as.integer(death),
    Gender = factor(Gender),
    ICUType = factor(ICUType)
  )

# 11-2

# Fix missing measurement counts after joining

model_table <- model_table |>
  mutate(
    across(
      ends_with("_count"),
      ~ tidyr::replace_na(.x, 0L)
    )
  )

# Verify that measurement counts contain no NA
stopifnot(
  sum(is.na(
    model_table |> select(ends_with("_count"))
  )) == 0
)

# 12. Keep benchmark scores separately
benchmarks <- model_table |>
  select(
    RecordID,
    death,
    `SAPS-I`,
    SOFA
  )

# SAPS-I and SOFA are not predictors
model_table <- model_table |>
  select(
    -`SAPS-I`,
    -SOFA
  )

# 13. Validate the dataset
stopifnot(
  nrow(model_table) == 4000,
  anyDuplicated(model_table$RecordID) == 0,
  !anyNA(model_table$death),
  all(model_table$death %in% c(0, 1))
)


list(predictors=model_table, benchmarks=benchmarks)
}


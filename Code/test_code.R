# Ensure Date is Date
group_data <- group_data %>%
  mutate(Date = as.Date(Date))

# First acquisition
first_acq <- group_data %>%
  filter(DiffHI == HC) %>%
  mutate(Date = as.Date(Date)) %>%
  arrange(Date) %>%
  group_by(Code) %>%
  slice(1) %>%             
  ungroup()

# Collapse by month
first_acq <- first_acq %>%
  mutate(month_date = as.Date(format(Date, "%Y-%m-01")))

# Create one event per month with ≥1 new learner
event_map <- first_acq %>%
  group_by(month_date) %>%
  summarise(n_new = n(), .groups = "drop") %>%
  arrange(month_date) %>%
  mutate(event_period = row_number())

# Attach back to individuals
first_acq <- first_acq %>%
  left_join(event_map[, c("month_date", "event_period")],
            by = "month_date")

# Attach to individuals in ILV_all
ILV_all <- ILV_all %>%
  left_join(
    first_acq[, c("Code", "event_period")],
    by = c("Alias" = "Code")
  )

ILV_all$event_period <- ifelse(is.na(ILV_all$event_period), 
                               max(na.omit(ILV_all$event_period)) + 1,
                               ILV_all$event_period) 

# Attach to observations
group_data <- group_data %>%
  mutate(month_date = as.Date(format(Date, "%Y-%m-01"))) %>%
  left_join(event_map[, c("month_date", "event_period")],
            by = "month_date")
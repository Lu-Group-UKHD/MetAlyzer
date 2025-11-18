# Helper function to calculate conditional mean based on significance

Calculates mean of a metric vector based on q-values. Prioritizes values
where qval \<= q_thresh. If none exist, uses all non-NA values.

## Usage

``` r
calculate_conditional_mean(metric_vec, qval_vec, q_thresh)
```

## Arguments

- metric_vec:

  Numeric vector of the metric to average (e.g., log2FC).

- qval_vec:

  Numeric vector of corresponding q-values.

- q_thresh:

  Significance threshold for q-values.

## Value

The calculated conditional mean, or NA.

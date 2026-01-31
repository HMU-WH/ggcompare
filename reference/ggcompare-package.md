# ggcompare: Mean Comparison in 'ggplot2'

Add mean comparison annotations to a 'ggplot'. This package provides an
easy way to indicate if two or more groups are significantly different
in a 'ggplot'. Usually you do not need to specify the test method, you
only need to tell stat_compare() whether you want to perform a
parametric test or a nonparametric test, and stat_compare() will
automatically choose the appropriate test method based on your data. For
comparisons between two groups, the p-value is calculated by t-test
(parametric) or Wilcoxon rank sum test (nonparametric). For comparisons
among more than two groups, the p-value is calculated by One-way ANOVA
(parametric) or Kruskal-Wallis test (nonparametric).

## See also

Useful links:

- <https://hmu-wh.github.io/ggcompare/>

- <https://github.com/HMU-WH/ggcompare/>

- Report bugs at <https://github.com/HMU-WH/ggcompare/issues/>

## Author

**Maintainer**: Hao Wang <wanghao8772@gmail.com>
([ORCID](https://orcid.org/0009-0000-9477-9585))

matstat
-------

Erlang statistics module. The current stdlib in OTP is lacking a statistics module.
I hope to add one sooner or later, in the meanwhile I will add stuff here.

In other words: the API will be subject to change until I feel it could be fixated.

I have two aspects of this:

 1. A somewhat complete set of functions that can be used as building blocks for functions to query measurements
 2. Sane return values that can be used as arguments directly with other functions

### Feature Areas ###

 * Statistics:
   * frequencies - `itemfreq/1`
   * histogram - `histogram_new/3`, `histogram_add/2`, `histogram_counts/2`, `histogram_property/1`
   * median - `cmedian/1`
   * mean - `tmean/1,2`, `gmean/1`, `hmean/1`
   * standard error - `tsem/1,2`
   * attributes - `tvar/1,2`, `tstd/1,2`, `tmin/1,2`, `tmax/1,2`
   * skew - TODO
   * kurtosis - TODO
 * Correlations:
   * Pearson’s product-moment correlation coefficient (Pearson's r) - `pearsonsr/1`
   * Point biserial correlation coefficient - TODO
   * Spearman’s rank correlation - TODO
   * Kendall rau rank correlation coefficient: tau a, tau b  - TODO
 * Regression:
   * Simple - `linregress/1`
   * Multiple - TODO
   * Polynomial regression - TODO
 * Factorial Analysis:
   * Extraction (PCA and Principal Axis) - TODO
   * Rotation (Varimax, Equimax, Quartimax) - TODO
   * Velicer’s MAP test - TODO
 * Variance Analysis: 
   * One-way ANOVA - TODO
 * Tests: F, T, Levene, U-Mannwhitney. - TODO


### Naming ###

 * Drop the leading 't' for all trimmed functions? Makes sense if they are the default ones.
 * Rename limit option 'inf' to 'none'? 'inf' works great for positive limit but we don't want two atoms, i.e 'inf', '-inf'.
 * `linregress/1` could get a prettier name.
 * `itemfreq/1` could get a prettier name.

### Done ###

#### Part of continuous stats function set ####

The continuous stats function set keeps a state for values added which may be queried.

 * `matstat:new() -> stats()` - Create a new stats state.
 * `matstat:new([{min | max, number()}]) -> stats()` - Create a new stats state.
 * `matstat:add([number()] | number(), stats()) -> stats()`.

The following functions may also query the stats state.

 *  Compute the trimmed mean
   * `tmean(stats() | [number()]) -> Mean :: float()`
   * `tmean([number()], {'inf' | L :: number(), 'inf' | U :: number()}) -> Mean :: float()`
 *  Compute the trimmed minimum
   * `tmin(stats() | [number()]) -> Minimum :: number()`
   * `tmin([number()], 'inf' | number()) -> Minimum :: number()`
 * Compute the trimmed maximum
   * `tmax(stats() | [number()]) -> Maximum :: number()` 
   * `tmax([number()], 'inf' | number()) -> Maximum :: number()` 
 * Compute the trimmed variance
   * `tvar(stats() | [number()]) -> Variance :: float()`
   * `tvar([number()], {'inf' | number(), 'inf' | number()}) -> Variance :: float()`
 * Compute the trimmed sample standard deviation 
   * `tstd(stats() | [number()]) -> StdDev :: float()`
   * `tstd([number()], {'inf' | number(), 'inf' | number()}) -> StdDev :: float()`
 * Compute the trimmed standard error of the mean
   * `tsem(stats() | [number()]) -> StdErr :: float()`
   * `tsem([number()], {'inf' | number(), 'inf' | number()}) -> StdErr :: float()`

#### Standalone for now ####

 * `cmedian([number()]) -> number()` - Returns the computed median value from a list of numbers
 * `gmean([number()]) -> Mean :: float()` - Compute the geometric mean along the specified axis.
 * `hmean([number()]) -> Mean :: float()` - Calculates the harmonic mean along the specified axis.
 * `linregress([{ X :: number(), Y :: number()}) -> {{Slope :: number(), Intercept :: number()}, RSq :: float()}` - Calculate a regression line
 * `itemfreq([term()]) -> [{term(), integer()}]` - Returns a 2D list of item frequencies. Highest frequency first.
 * `pearsonr([{number(),number()}) -> float()` - Calculates a Pearson correlation coefficient (and the p-value for testing *not yet impl.*)
 * `histogram([number()], Nbins :: integer()) -> [{number(), integer()}]` - Separates the range into several bins and returns the number of instances of a in each bin.
   * `histogram_new/3`, `histogram_add/2`, `histogram_counts/2`, `histogram_property/1`
 * `msn/1 -> {Mean :: float(), StdDev :: float()}` - calculate mean and *sampled* standard deviation,

Function list from SciPy.stats Statistical Functions.

### ToDo ###

Will Implement in Prio order:

 * `chisquare(f_obs[, f_exp, ddof])` - Calculates a one-way chi square test.
 * skew(a[, axis, bias]) - Computes the skewness of a data set.
 * kurtosis(a[, axis, fisher, bias]) - Computes the kurtosis (Fisher or Pearson) of a dataset.
 * zmap(scores, compare[, axis, ddof]) - Calculates the relative z-scores.
 * zscore(a[, axis, ddof]) - Calculates the z score of each value in the sample, relative to the sample mean and standard deviation.
 * moment(a[, moment, axis]) - Calculates the nth moment about the mean for a sample.
 * spearmanr(a[, b, axis]) - Calculates a Spearman rank-order correlation coefficient and the p-value

Still flaky about:

 * mode(a[, axis]) - Returns an array of the modal (most common) value in the passed array.
 * variation(a[, axis]) - Computes the coefficient of variation, the ratio of the biased standard deviation to the mean.
 * describe(a[, axis]) - Computes several descriptive statistics of the passed array.
 * skewtest(a[, axis]) - Tests whether the skew is different from the normal distribution.
 * kurtosistest(a[, axis]) - Tests whether a dataset has normal kurtosis
 * normaltest(a[, axis]) - Tests whether a sample differs from a normal distribution.
 * scoreatpercentile(a, per[, limit, ...]) - Calculate the score at the given per percentile of the sequence a.
 * percentileofscore(a, score[, kind]) - The percentile rank of a score relative to a list of scores.
 * cumfreq(a[, numbins, defaultreallimits, weights]) - Returns a cumulative frequency histogram, using the histogram function.
 * relfreq(a[, numbins, defaultreallimits, weights]) - Returns a relative frequency histogram, using the histogram function.
 * `obrientransform(*args)` - Computes a transform on input data (any number of columns).
 * signaltonoise(a[, axis, ddof]) - The signal-to-noise ratio of the input data.
 * `bayes_mvs(data[, alpha])` - Bayesian confidence intervals for the mean, var, and std.
 * threshold(a[, threshmin, threshmax, newval]) - Clip array to a given value.
 * trimboth(a, proportiontocut) - Slices off a proportion of items from both ends of an array.
 * trim1(a, proportiontocut[, tail]) - Slices off a proportion of items from ONE end of the passed array
 * `f_oneway(*args)` - Performs a 1-way ANOVA.
 * pointbiserialr(x, y) - Calculates a point biserial correlation coefficient and the associated p-value.
 * `kendalltau(x, y[, initial_lexsort])` - Calculates Kendall’s tau, a correlation measure for ordinal data.
 * `ttest_1samp(a, popmean[, axis])` - Calculates the T-test for the mean of ONE group of scores a.
 * `ttest_ind(a, b[, axis, equal_var])` - Calculates the T-test for the means of TWO INDEPENDENT samples of scores.
 * `ttest_rel(a, b[, axis])` - Calculates the T-test on TWO RELATED samples of scores, a and b.
 * kstest(rvs, cdf[, args, N, alternative, mode]) - Perform the Kolmogorov-Smirnov test for goodness of fit
 * `ks_2samp(data1, data2)` - Computes the Kolmogorov-Smirnof statistic on 2 samples.
 * `mannwhitneyu(x, y[, use_continuity])` - Computes the Mann-Whitney rank test on samples x and y.
 * tiecorrect(rankvals) - Tie correction factor for ties in the Mann-Whitney U and
 * ranksums(x, y) - Compute the Wilcoxon rank-sum statistic for two samples.
 * wilcoxon(x[, y]) - Calculate the Wilcoxon signed-rank test.
 * `kruskal(*args)` - Compute the Kruskal-Wallis H-test for independent samples
 * `friedmanchisquare(*args)` - Computes the Friedman test for repeated measurements
 * ansari(x, y) - Perform the Ansari-Bradley test for equal scale parameters
 * `bartlett(*args)` - Perform Bartlett’s test for equal variances
 * `levene(*args, **kwds)` - Perform Levene test for equal variances.
 * shapiro(x[, a, reta]) - Perform the Shapiro-Wilk test for normality.
 * anderson(x[, dist]) - Anderson-Darling test for data coming from a particular distribution
 * `binom_test(x[, n, p])` - Perform a test that the probability of success is p.
 * `fligner(*args, **kwds)` - Perform Fligner’s test for equal variances.
 * mood(x, y) - Perform Mood’s test for equal scale parameters.
 * `oneway(*args, **kwds)` - Test for equal means in two or more samples from the normal distribution.

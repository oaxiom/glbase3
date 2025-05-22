"""

stats sub_class for expression objects

"""

import statistics
import numpy
import scipy.stats

class stats:
    def __init__(self, parent):
        self.parent = parent

    def print_stats(self, log=None):
        """
        **Purpose**
            Print out the min, max, mean, median and stddev for each condition in the microarray data

        **Arguments**
            log (Optional)
                if true log2 transform the data.

                (Data with an expression of zero will be excluded from the log calculation and
                mean scores)

        **Returns**
            a dictionary containing min, max, mean, median, stdev for each condition. In the order
            of "conditions".
        """
        means = []
        medians = []
        stds = []
        mins = []
        maxs = []
        sums = []

        print("All Array:")
        print("\tmin", self.parent.stats.min())
        print("\tmax", self.parent.stats.max())
        print("\tmean", self.parent.stats.mean())
        print("\tstddev", self.parent.stats.stdev())

        for i, k in enumerate(self.parent.getConditionNames()):
            expn_values = self.parent.numpy_array_all_data[:,i]
            if log:
                expn_values = [v for v in expn_values if int(v*10000) > 0]
                expn_values = numpy.log2(expn_values)

            aaverage = numpy.average(expn_values)
            amedian = numpy.median(expn_values)
            astd = numpy.std(expn_values)
            amin = numpy.min(expn_values)
            amax = numpy.max(expn_values)
            asum = numpy.sum(expn_values)

            print("Condition: %s" % k)
            print("\tmin", amin)
            print("\tmax", amax)
            print("\tsum", asum)
            print("\tmean", aaverage)
            print("\tmedian", amedian)
            print("\tstddev", astd)

            means.append(aaverage)
            medians.append(amedian)
            stds.append(astd)
            mins.append(amin)
            maxs.append(amax)
            sums.append(asum)

        return {"mean": means, "median": medians, "stdev": stds, "min": mins, "max": maxs, "sums": sums,
            "conditions": self.parent.getConditionNames()}

    def all(self, log=None):
        """
        **Purpose**
            Get the min, max, mean, median and stddev for each condition in the expression data

            TODO: Add test for normality

        **Arguments**
            log (Optional)
                if true log2 transform the data.

                (Data with an expression of zero will be excluded from the log calculation and
                mean scores)

        **Returns**
            a dictionary containing min, max, mean, median, stdev for each condition. In the order
            of "conditions".
        """
        means = []
        medians = []
        stds = []
        mins = []
        maxs = []

        for k in self._conditions:
            expn_values = self.parent.numpy_array_all_data
            if log:
                expn_values = [v for v in expn_values if int(v*10000) > 0]
                expn_values = numpy.log2(expn_values)

            aaverage = numpy.mean(expn_values)
            amedian = numpy.median(expn_values)
            astd = numpy.std(expn_values)
            amin = self.min()
            amax = self.max()

            means.append(aaverage)
            medians.append(amedian)
            stds.append(astd)
            mins.append(amin)
            maxs.append(amax)

        return {"mean": means, "median": medians, "stdev": stds, "min": mins, "max": maxs,
            "conditions": self.parent._conditions}

    def min(self):
        """
        **Purpose**
            get the minimum value in the expression data

        """
        return self.parent.numpy_array_all_data.min()

    def max(self):
        """
        **Purpose**
            get the maximum value in the expression data

        """
        return self.parent.numpy_array_all_data.max()

    def mean(self):
        """
        **Purpose**
            get the mean value in the expression data

        """
        return self.parent.numpy_array_all_data.mean()

    def stdev(self):
        """
        **Purpose**
            get the standard deviation of the expression data
            (Does not calculate if your data is actually normal)

        """
        return self.parent.numpy_array_all_data.std()

    def ttest(self, condition1, condition2, equal_var=False):
        """
        **Purpose**
            perform a t-test between condition1 and condition2.

            Assumes the data is normally distributed. The t-test is only valid if the
            data is normally distributed.

            This performs and independent t-test of the samples. And assumes the
            two samples are independent. An appropriate use of this test would be
            to compare between a set of genes across two different biological samples.

            See Scipy.stats.ttest_ind() for more details

        **Arguments**
            condition1
                the name of condition 1 to use

            condition2
                the name of condition 2 to use

            equal_var (Optional, default=False)
                Generally the variation in samples is not expected to be the same for a typical
                set of expression data.

                So, set this to False and perform a 'Welch's t-test instead'. You can set it to True
                and perform a Student's independent t-test which assumes sample variances are the same.

                See:
                http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html

                For more details.

        **Returns**
            The t-statistic and the p-value
        """
        assert condition1 in self.parent._conditions, "condition1 '%s' not on this array" % condition1
        assert condition2 in self.parent._conditions, "condition2 '%s' not on this array" % condition2

        data1 = self.parent.getDataForCondition(condition1)
        data2 = self.parent.getDataForCondition(condition2)

        return scipy.stats.ttest_ind(data1, data2, equal_var=equal_var)

    def pearsonr(self, condition1, condition2):
        """
        **Purpose**
            perform a Pearson correlation test between condition1 and condition2.

            See Scipy.stats.pearsonr() for more details

        **Arguments**
            condition1
                the name of condition 1 to use

            condition2
                the name of condition 2 to use

        **Returns**
            The pearson correlation co-efficient and the p-value
        """
        return self.__unified_stats("pearsonr", condition1, condition2)

    def spearmanr(self, condition1, condition2):
        """
        **Purpose**
            perform a Spearman correlation test between condition1 and condition2.

            See Scipy.stats.pearsonr() for more details

        **Arguments**
            condition1
                the name of condition 1 to use

            condition2
                the name of condition 2 to use

        **Returns**
            The pearson correlation co-efficient and the p-value
        """
        return self.__unified_stats("spearmanr", condition1, condition2)

    def kruskal(self, list_of_conditions):
        """
        **Purpose**
            perform a Kruskal-Wallis H-test for independent samples for several conditions.

            See Scipy.stats.kruskal() for more details

        **Arguments**
            list_of_conditions (Required)
                a list of condition names in this expression object to test.

        **Returns**
            The kruskal correlation co-efficient and the p-value
        """
        datas = [self.parent.getDataForCondition(c) for c in list_of_conditions]

        t, p = scipy.stats.kruskal(*datas)
        return (t, p)

    def mannwhitneyu(self, condition1, condition2):
        """
        **Purpose**
            perform a Mann Whitney U-test for independent samples for two conditions.

            See Scipy.stats.mannwhitneyu() for more details

        **Arguments**
            condition1
                the name of condition 1 to use

            condition2
                the name of condition 2 to use

        **Returns**
            The u-statistic and two-sided p-value
        """
        data1 = self.parent.getDataForCondition(condition1)
        data2 = self.parent.getDataForCondition(condition2)

        return scipy.stats.mannwhitneyu(data1, data2, alternative='two-sided')

    def wilcoxon(self, condition1, condition2):
        """
        **Purpose**
            perform a 'Wilcoxon signed-rank' test between condition1 and condition2.

            See Scipy.stats.wilcoxon() for more details

        **Arguments**
            condition1
                the name of condition 1 to use

            condition2
                the name of condition 2 to use

        **Returns**
            The wilcoxon z-statistic and p-value
        """
        return self.__unified_stats("wilcoxon", condition1, condition2)

    def __unified_stats(self, test, condition1, condition2):
        assert condition1 in self.parent._conditions, "condition1 '%s' not on this array" % condition1
        assert condition2 in self.parent._conditions, "condition2 '%s' not on this array" % condition2

        data1 = self.parent.getDataForCondition(condition1)
        data2 = self.parent.getDataForCondition(condition2)

        test_funcs = {"wilcoxon": scipy.stats.wilcoxon,
            "spearmanr": scipy.stats.spearmanr,
            "pearsonr": scipy.stats.pearsonr
            }

        return test_funcs[test](data1, data2)

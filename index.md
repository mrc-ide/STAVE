
# STAVE

STAVE is an R package designed for storing and managing aggregate
spatial-temporal genetic data. By “aggregate” we mean that data consist solely
of counts (numerator and denominator), rather than individual-level details. By
“spatial-temporal” we mean that each observation is tied to a specific
geographic location (latitude and longitude), and a precise collection time (a
single day). Once data are loaded into STAVE it becomes straightforward to
calculate the prevalence of different genetic variants at each sampling
location, making this a useful format for data harmonization prior to spatial
analysis.

At the core of STAVE is a single class (an R6 object) that serves as the main data
container. This class enables users to efficiently import, store, and manipulate
genetic data. All operations, such as calculating prevalence, are performed
through dedicated member functions.

This site describes the [STAVE data format](TODO). It also contains
[installation instructions](TODO) and simple [tutorials](TODO) to help you get
started.


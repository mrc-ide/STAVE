
# STAVE

STAVE stands for **S**patial **T**emporal **A**ggregate **V**ariant
**E**ncoding. It is an R package designed for storing and managing aggregate
spatial-temporal genetic data. By “aggregate” we mean that data consist solely
of counts, i.e. numerators and denominators, rather than individual-level
details. By “spatial-temporal” we mean that each observation is tied to a
specific geographic location (latitude and longitude), and a precise collection
time (a single day). This format is ideally suited for downstream analysis
methods that assume a point process, but could also be used for areal methods.
Although STAVE has been developed with *Plasmodium* genetic data in mind, it is
general enough to accommodate data from other organisms.

At the core of STAVE is a single class (an R6 object) that serves as the main data
container. This class enables users to efficiently import, store, and manipulate
genetic data. All operations, such as calculating prevalence, are performed
through dedicated member functions.

This site describes the [STAVE data format](TODO). It also contains
[installation instructions](TODO) and simple [tutorials](TODO) to help you get
going.

## Motivations for developing STAVE

### 1. Avoiding duplication of information

STAVE achieves a relatively efficient encoding through a relational database
structure. There are three tables that make up this database:

1. **Studies**: This is a record of where the data came from. This could be
academic publications, reports, or other sources
of information.
2. **Surveys**: A survey is defined here as a discrete instance of data
collection. Survey-level data include information on location (latitude and
longitude) and time (day) of collection. There may be multiple surveys within a
given study, for example different sampling sites or different collection
periods may come from the same academic publication.
3. **Counts**: This is where the aggregate genetic data are stored. Data include
the name of the variant that was observed (i.e. the gene, locus, and mutation),
how many times it was observed, and how many times it was successfully tested
for.

Separating the data into tables helps eliminate duplication and improve
efficiency. For instance, if the data originate from an academic paper, the same
list of authors will typically apply to all surveys reported within that paper.
Repeating this information for each survey is not only inefficient but also
increases the risk of data entry errors. To address this, we store shared
information, such as authorship, at the study level and link it to survey-level
data using a relational key. This approach ensures that the information is
stored only once while maintaining clear connections between studies and their
associated surveys.

### 2. Nailing down the input format

STAVE performs close to 100 different checks on the data when it is read in.
This includes everything from character formatting to checking that prevalence
cannot exceed 100%. Having this extremely pedantic series of checks ensures that
we know exactly what we have once loaded, which makes life easier when it comes
to downstream analysis.

### 3. Flexible encoding genetic variants

Encoding *Plasmodium* genetic data can be challenging due to the complexity and
variability of the information involved. Sequence data may be recorded at the
single-codon level or across multiple codons, and heterozygous calls — where
more than one allele is observed — can occur at some codons but not others.
These heterozygous calls may also be phased or unphased, adding another layer of
detail. Furthermore, patterns of missingness can vary between individuals,
causing the denominator to vary as we look along the genome. In many cases,
particularly when working with aggregate counts extracted from published
studies, some of this detailed information may be incomplete or entirely
unavailable. Therefore, we require an encoding system that is both flexible and
expressive: capable of capturing this complexity when the data are available,
but not overly restrictive in situations where certain aspects are missing.

A drawback of this flexible encoding is that calculating prevalence becomes less
straightforward. This trade-off is somewhat unavoidable, as flexibility in data
encoding and simplicity in prevalence calculation are somewhat opposing goals.
In STAVE, we prioritize achieving a clean and flexible data structure, ensuring
that all relevant information is captured. Prevalence calculations are then
performed using dedicated member functions that operate on the encoded data.



## Example workflow

Input data take the form of three tables. We will use the pre-loaded example
data here, if you are using your own data then you need to match this input format:
```{r}
library(STAVE)

data("example_input")
View(example_input$studies)
View(example_input$surveys)
View(example_input$counts)
```

Create a STAVE object and append the data:
```{r}
# create new object
s <- STAVE_object$new()

# append data using a member function
s$append_data(studies_dataframe = example_input$studies,
              surveys_dataframe = example_input$surveys,
              counts_dataframe = example_input$counts)

# check how many studies are now loaded
s
```

Once data are loaded, we can always view the different tables using *get* functions. However, we cannot alter the values directly.
```{r}
s$get_studies() |> head()
s$get_surveys() |> head()
s$get_counts() |> head()
```

We can calculate the prevalence of any variant using `get_prevalence()`:
```{r}
p <- s$get_prevalence("mdr1:184:F")
View(p)
```

We can also return a list of all variants in the data:
```{r}
s$get_variants()
```

Finally, we can selectively drop studies from the data using their study ID:
```{r}
s$drop_study(drop_study_ID = "wwarn_10297_Nelson")
s
```


<img src="reference/figures/hex.png" alt="STAVE Logo" style="max-height: 250px; margin-top:50px; margin-bottom:20px;">

STAVE is a lightweight R package focused on a specific use case: working with aggregated genetic data for molecular surveillance. It provides three core functionalities:

1. **Storage of aggregate genetic data using [variant
string format](https://github.com/mrc-ide/variantstring)**. STAVE encodes numerator/denominator-style data in a compact string representation tailored for amino acid-level variation. This format supports both single codons and multi-locus haplotypes, making it especially suitable for drug resistance markers.
2. **Linkage of genetic data to precise space-time coordinates**. Rather than relying on administrative names (which may be ambiguous or change over time), STAVE links each data point to specific geographic coordinates and timestamps, ensuring consistency across datasets.
3. **Computation of prevalence from encoded data**. Calculating prevalence isn’t always straightforward. For instance, the variant of interest may represent only a subset of a broader haplotype, or the data may include mixed (heterozygous) calls. STAVE handles these complexities with dedicated logic.

STAVE is deliberately kept simple and focused. It’s designed to do just these three things—and to do them reliably. While it doesn’t cover everything (for example, it doesn’t deal with individual-level data or generate maps), it helps streamline one tricky part of working with genetic surveillance data. It’s intended to fit neatly into a broader workflow, making life a little easier for researchers working in this space.

The [How it works](articles/relational_structure.html) section describes the
STAVE data format, which you will need to understand before importing data. Once
you are ready, checkout the [Installation](articles/installation.html) and
[Tutorials](articles/reading_in_data.html) sections for practical examples.


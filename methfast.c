#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "cgranges.h"

#define MAX_LINE_LENGTH 1024

typedef struct {
    float fraction;
    int coverage;
} MethInterval;

typedef struct {
    cgranges_t *cr;
    MethInterval *meth_intervals;
    int num_intervals;
    int capacity;
} MethRanges;

// Define a struct to hold the result arrays and sums
typedef struct {
    float *meth_coverage;         // Array of methylated coverages
    float *unmeth_coverage;       // Array of unmethylated coverages
    int *total_coverage;          // Array of total coverages
    float *fraction_methylation;  // Array of methylation fractions
    int num_positions;            // Total number of positions (size of arrays)

    // Sums of each coverage type
    float sum_meth_coverage;
    float sum_unmeth_coverage;
    int sum_total_coverage;
} MethStats;

// Function to free MethStats struct
void free_meth_stats(MethStats *stats) {
    free(stats->meth_coverage);
    free(stats->unmeth_coverage);
    free(stats->total_coverage);
    free(stats->fraction_methylation);
    free(stats);
}

// Initialize MethRanges structure
MethRanges *init_meth_ranges() {
    MethRanges *ranges = (MethRanges *)malloc(sizeof(MethRanges));
    ranges->cr = cr_init();
    ranges->num_intervals = 0;
    ranges->capacity = 1024;
    ranges->meth_intervals = (MethInterval *)malloc(ranges->capacity * sizeof(MethInterval));
    return ranges;
}

// Free MethRanges structure
void free_meth_ranges(MethRanges *ranges) {
    cr_destroy(ranges->cr);
    free(ranges->meth_intervals);
    free(ranges);
}

// Add an interval to MethRanges
void add_interval(MethRanges *ranges, const char *chrom, int start, int end, float fraction, int coverage) {
    if (ranges->num_intervals == ranges->capacity) {
        ranges->capacity *= 2;
        ranges->meth_intervals = (MethInterval *)realloc(ranges->meth_intervals, ranges->capacity * sizeof(MethInterval));
    }

    ranges->meth_intervals[ranges->num_intervals].fraction = fraction;
    ranges->meth_intervals[ranges->num_intervals].coverage = coverage;

    // This previously added 1 to the start to make it 1-based, but it seems cgranges expected 0-based coordinates.
    cr_add(ranges->cr, chrom, start, end, ranges->num_intervals);
    ranges->num_intervals++;
}

// Index intervals in MethRanges
void index_meth_ranges(MethRanges *ranges) {
    cr_index(ranges->cr);
}


// collect_meth_stats function
MethStats *collect_meth_stats(MethRanges *ranges, const char *chrom, int start, int end) {
    int64_t *overlap_indices = NULL;
    int64_t m_overlap = 0;
    int64_t num_overlaps = cr_overlap(ranges->cr, chrom, start, end, &overlap_indices, &m_overlap);

    // Allocate MethStats struct and arrays for each metric
    MethStats *stats = (MethStats *)malloc(sizeof(MethStats));
    stats->meth_coverage = (float *)malloc(num_overlaps * sizeof(float));
    stats->unmeth_coverage = (float *)malloc(num_overlaps * sizeof(float));
    stats->total_coverage = (int *)malloc(num_overlaps * sizeof(int));
    stats->fraction_methylation = (float *)malloc(num_overlaps * sizeof(float));
    stats->num_positions = num_overlaps;

    // Initialize sums to zero
    stats->sum_meth_coverage = 0.0;
    stats->sum_unmeth_coverage = 0.0;
    stats->sum_total_coverage = 0;

    // Populate the arrays in the MethStats struct and calculate sums
    for (int64_t i = 0; i < num_overlaps; i++) {
        int idx = overlap_indices[i];
        int meth_start = cr_start(ranges->cr, idx);
        int meth_end = cr_end(ranges->cr, idx);

        if (meth_start < start) meth_start = start;
        if (meth_end > end) meth_end = end;

        float fraction = ranges->meth_intervals[idx].fraction;
        int coverage = ranges->meth_intervals[idx].coverage;
        
        // Compute meth and unmeth coverages
        float meth_coverage = fraction * coverage;
        float unmeth_coverage = (1.0f - fraction) * coverage;

        // Store values in the arrays
        stats->meth_coverage[i] = meth_coverage;
        stats->unmeth_coverage[i] = unmeth_coverage;
        stats->total_coverage[i] = coverage;
        stats->fraction_methylation[i] = fraction;

        // Update sums
        stats->sum_meth_coverage += meth_coverage;
        stats->sum_unmeth_coverage += unmeth_coverage;
        stats->sum_total_coverage += coverage;
    }

    free(overlap_indices);

    return stats;
}

// Check if a file is gzipped or bzipped
int is_gzipped(const char *filepath) {
    FILE *fp = fopen(filepath, "rb");
    if (!fp) {
        perror("fopen");
        return 0;
    }

    // Read the first three bytes from the header
    unsigned char header[3];
    size_t n = fread(header, 1, sizeof(header), fp);
    fclose(fp);

    // If we can't read enough bytes, it's not a valid gzipped file
    if (n < sizeof(header))
        return 0;

    // Check for gzip magic numbers: 0x1F, 0x8B and the deflate method (0x08)
    if (header[0] != 0x1F || header[1] != 0x8B || header[2] != 0x08)
        return 0;

    // File is gzipped (or BGZipped, as BGZF is a subtype of gzip)
    return 1;
}

MethRanges *parse_meth_bed(const char *filepath, int frac_col, int cov_col, int meth_col, int unmeth_col) {
    MethRanges *ranges = init_meth_ranges();

    char chrom[100], line[1024];
    int start, end;
    float fraction;
    int coverage, methylated, unmethylated;

    FILE *plain_file = NULL;
    gzFile gz_file = NULL;

    int is_gz = is_gzipped(filepath);
    if (is_gz) {
        gz_file = gzopen(filepath, "r");
        if (!gz_file) {
            fprintf(stderr, "Error opening gzipped file: %s\n", filepath);
            exit(EXIT_FAILURE);
        }
    } else {
        plain_file = fopen(filepath, "r");
        if (!plain_file) {
            fprintf(stderr, "Error opening file: %s\n", filepath);
            exit(EXIT_FAILURE);
        }
    }

    while ((is_gz ? gzgets(gz_file, line, sizeof(line)) : fgets(line, sizeof(line), plain_file)) != NULL) {
        char *fields[20];
        int field_count = 0;
        char *token = strtok(line, "\t\n");
        while (token && field_count < 20) {
            fields[field_count++] = token;
            token = strtok(NULL, "\t\n");
        }

        if (field_count < 4) continue;  // Skip if there aren't enough columns

        strncpy(chrom, fields[0], sizeof(chrom) - 1);
        chrom[sizeof(chrom) - 1] = '\0'; // Ensure null termination
        start = atoi(fields[1]);
        end = atoi(fields[2]);

        if (meth_col > 0 && meth_col <= field_count && unmeth_col > 0 && unmeth_col <= field_count) {
            methylated = atoi(fields[meth_col - 1]);
            unmethylated = atoi(fields[unmeth_col - 1]);
            coverage = methylated + unmethylated;
            fraction = coverage > 0 ? (float)methylated / coverage : 0.0;
        } else if (meth_col > 0 && meth_col <= field_count && cov_col > 0 && cov_col <= field_count) {
            methylated = atoi(fields[meth_col - 1]);
            coverage = atoi(fields[cov_col - 1]);
            fraction = coverage > 0 ? (float)methylated / coverage : 0.0;
        } else if (cov_col > 0 && cov_col <= field_count && frac_col > 0 && frac_col <= field_count) {
            fraction = atof(fields[frac_col - 1]);
            coverage = atoi(fields[cov_col - 1]);
        } else {
            fprintf(stderr, "Error: invalid column indices\n");
            exit(EXIT_FAILURE);
        }

        add_interval(ranges, chrom, start, end, fraction, coverage);
    }

    if (is_gz) gzclose(gz_file);
    else fclose(plain_file);

    index_meth_ranges(ranges);
    return ranges;
}

// Process target BED entries and compute weighted methylation fraction for each
void process_targets(MethRanges *ranges, const char *target_filepath, FILE *outfp) {
    FILE *file = fopen(target_filepath, "r");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", target_filepath);
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE_LENGTH];
    char chrom[100];
    int start, end;

    while (fgets(line, sizeof(line), file)) {
        // Split the line into tokens using tab as the delimiter
        char *token = strtok(line, "\t");
        if (!token) continue;

        // Get chrom, start, and end from the tokens
        strncpy(chrom, token, sizeof(chrom) - 1);
        chrom[sizeof(chrom) - 1] = '\0'; // Ensure null termination

        token = strtok(NULL, "\t");
        if (!token) continue;
        start = atoi(token);

        token = strtok(NULL, "\t");
        if (!token) continue;
        end = atoi(token);

        MethStats *stats = collect_meth_stats(ranges, chrom, start, end);

        // calculate weighted fraction
        float weighted_fraction = stats->sum_total_coverage > 0
            ? stats->sum_meth_coverage / stats->sum_total_coverage
            : 0.0;

        // Print interval, weighted fraction, and number of positions to outfp
        fprintf(outfp, "%s\t%d\t%d\t%d\t%d\t%.4f\n",
                chrom, start, end, stats->num_positions, stats->sum_total_coverage, weighted_fraction);

        free_meth_stats(stats);
    }

    fclose(file);
}


// Main function
int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr,
            "Usage: %s <methylation_bed(.gz)> <target_bed> "
            "[-f <frac_col>] [-c <cov_col>] [-m <meth_col>] [-u <unmeth_col>] "
            "[-o <output_file>]\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    int frac_col = 4, cov_col = 5, meth_col = 0, unmeth_col = 0;
    FILE *outfp = stdout;  // Default to stdout
    int i;

    for (i = 3; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
            frac_col = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc) {
            cov_col = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
            meth_col = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-u") == 0 && i + 1 < argc) {
            unmeth_col = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            outfp = fopen(argv[++i], "w");
            if (!outfp) {
                fprintf(stderr, "Error: Unable to open output file.\n");
                return EXIT_FAILURE;
            }
        }
    }

    MethRanges *ranges = parse_meth_bed(argv[1], frac_col, cov_col, meth_col, unmeth_col);
    process_targets(ranges, argv[2], outfp);

    free_meth_ranges(ranges);

    // Close outfp if it's not stdout
    if (outfp && outfp != stdout) {
        fclose(outfp);
    }

    return EXIT_SUCCESS;
}
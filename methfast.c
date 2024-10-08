#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "cgranges.h"

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

typedef struct {
    float mean_fraction;
    int total_coverage;
    int total_meth_counts;
    int total_unmeth_counts;
    int num_cpgs;

} MethInfo;

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

    cr_add(ranges->cr, chrom, start+1, end, ranges->num_intervals);
    ranges->num_intervals++;
}

// Index intervals in MethRanges
void index_meth_ranges(MethRanges *ranges) {
    cr_index(ranges->cr);
}

// Calculate weighted methylation fraction for a given target interval
float compute_weighted_fraction(MethRanges *ranges, const char *chrom, int start, int end) {
    int64_t *overlap_indices = NULL;
    int64_t m_overlap = 0;
    int64_t num_overlaps = cr_overlap(ranges->cr, chrom, start, end, &overlap_indices, &m_overlap);

    float meth_coverage = 0.0;
    int total_coverage = 0;

    for (int64_t i = 0; i < num_overlaps; i++) {
        int idx = overlap_indices[i];
        int meth_start = cr_start(ranges->cr, idx);
        int meth_end = cr_end(ranges->cr, idx);

        if (meth_start < start) meth_start = start;
        if (meth_end > end) meth_end = end;

        float fraction = ranges->meth_intervals[idx].fraction;
        int coverage = ranges->meth_intervals[idx].coverage;

        meth_coverage += fraction * coverage;
        total_coverage += coverage;
    }

    free(overlap_indices);

    return total_coverage > 0 ? meth_coverage / total_coverage : 0;
}

// Check if file is gzipped by extension
int is_gzipped(const char *filepath) {
    size_t len = strlen(filepath);
    return (len > 3 && strcmp(filepath + len - 3, ".gz") == 0);
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

        if (meth_col > 0 && unmeth_col > 0) {
            methylated = atoi(fields[meth_col - 1]);
            unmethylated = atoi(fields[unmeth_col - 1]);
            coverage = methylated + unmethylated;
            fraction = coverage > 0 ? (float)methylated / coverage : 0.0;
        } else {
            fraction = atof(fields[frac_col - 1]);
            coverage = atoi(fields[cov_col - 1]);
        }

        add_interval(ranges, chrom, start, end, fraction, coverage);
    }

    if (is_gz) gzclose(gz_file);
    else fclose(plain_file);

    index_meth_ranges(ranges);
    return ranges;
}

// Process target BED entries and compute weighted methylation fraction for each
void process_targets(MethRanges *ranges, const char *target_filepath) {
    FILE *file = fopen(target_filepath, "r");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", target_filepath);
        exit(EXIT_FAILURE);
    }

    char chrom[100];
    int start, end;

    while (fscanf(file, "%s\t%d\t%d\n", chrom, &start, &end) == 3) {
        float weighted_fraction = compute_weighted_fraction(ranges, chrom, start, end);
        printf("%s\t%d\t%d\t%.4f\n", chrom, start, end, weighted_fraction);
    }

    fclose(file);
}

// Main function
int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <methylation_bed(.gz)> <target_bed> [-f <frac_col>] [-c <cov_col>] [-m <meth_col>] [-u <unmeth_col>]\n", argv[0]);
        return EXIT_FAILURE;
    }

    int frac_col = 4, cov_col = 5, meth_col = 0, unmeth_col = 0;
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) frac_col = atoi(argv[++i]);
        else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc) cov_col = atoi(argv[++i]);
        else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) meth_col = atoi(argv[++i]);
        else if (strcmp(argv[i], "-u") == 0 && i + 1 < argc) unmeth_col = atoi(argv[++i]);
    }

    MethRanges *ranges = parse_meth_bed(argv[1], frac_col, cov_col, meth_col, unmeth_col);
    process_targets(ranges, argv[2]);

    free_meth_ranges(ranges);
    return EXIT_SUCCESS;
}
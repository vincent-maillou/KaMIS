/**
 * kamis.h — C API for the KaMIS Maximum Independent Set library.
 *
 * This header exposes the KaMIS evolutionary MIS solver through a
 * stable, C-compatible interface suitable for ctypes bindings.
 *
 * The graph is represented using CSR (Compressed Sparse Row) format,
 * identical to the METIS convention:
 *
 *   - xadj[n+1] : index pointer array
 *   - adjncy[m] : adjacency array
 *
 * Both arrays are 0‑based (C‑indexing).
 *
 * === Typical usage ===
 *
 *   // Configure the solver
 *   KamisConfig config;
 *   kamis_default_config(&config);
 *   config.time_limit = 60.0;
 *
 *   // Solve
 *   int solution_size = 0;
 *   int* solution = NULL;
 *   int ret = kamis_solve_mis(n, xadj, adjncy, &solution, &solution_size, &config);
 *   if (ret == 0) {
 *       // solution[0..solution_size-1] holds the MIS vertex indices
 *   }
 *   kamis_free_solution(solution);
 */

#ifndef KAMIS_H
#define KAMIS_H

#ifdef __cplusplus
extern "C" {
#endif

/** Return codes */
#define KAMIS_OK             0
#define KAMIS_ERROR_INPUT   -1
#define KAMIS_ERROR_MEMORY  -2

/** Kernelization modes */
#define KAMIS_KERNEL_FASTKER 0   /**< FastKer (default, faster) */
#define KAMIS_KERNEL_FULL    1   /**< Full branch-and-reduce kernelization */

/** Preset configuration identifiers */
#define KAMIS_CONFIG_STANDARD      0
#define KAMIS_CONFIG_SOCIAL        1
#define KAMIS_CONFIG_FULL_STANDARD 2
#define KAMIS_CONFIG_FULL_SOCIAL   3

/* ---------------------------------------------------------------------------
 * Configuration structure
 * --------------------------------------------------------------------------- */
typedef struct {
    /** Seed for the random number generator (0 = use time). */
    int seed;

    /** Time limit in seconds for the evolutionary algorithm. */
    double time_limit;

    /** Size of the population used by the evolutionary algorithm. */
    int population_size;

    /** Imbalance parameter for KaHIP partitioning calls (e.g., 3.0). */
    double imbalance;

    /** KaHIP partitioning mode (0 = FAST, 1 = ECO, 2 = STRONG,
        3 = FASTSOCIAL, 4 = ECOSOCIAL, 5 = STRONGSOCIAL). */
    int kahip_mode;

    /** Kernelization style: KAMIS_KERNEL_FASTKER or KAMIS_KERNEL_FULL. */
    int kernel_mode;

    /** Number of repetitions per round. */
    int repetitions;

    /** Insert a solution if no new solution has been inserted
        for this number of operations. */
    int insert_threshold;

    /** Update the pool of node separators after this many rounds. */
    int pool_threshold;

    /** Factor for timing-based renewal of the separator pool. */
    double pool_renewal_factor;

    /** Whether to randomize the imbalance on each KaHIP call. */
    int randomize_imbalance;

    /** Whether to diversify initial solutions using different RNG seeds. */
    int diversify;

    /** Whether to use tournament selection instead of random. */
    int enable_tournament_selection;

    /** Use the vertex-cover approach for the multiway combine operator. */
    int use_multiway_vc;

    /** Number of blocks for multiway combine operators. */
    int multiway_blocks;

    /** Number of individuals in a tournament. */
    int tournament_size;

    /** Percentage for the mutation of a solution (0–100). */
    int flip_coin;

    /** Force parameter for the mutation. */
    int force_k;

    /** Number of candidates for forced insertion. */
    int force_cand;

    /** Number of initial node separators to construct. */
    int number_of_separators;

    /** Number of initial partitions to construct. */
    int number_of_partitions;

    /** Number of initial k-separators to construct. */
    int number_of_k_separators;

    /** Number of initial k-partitions to construct. */
    int number_of_k_partitions;

    /** Number of iterations for the inner ILS solver. */
    int ils_iterations;

    /** Use the Hopcroft-Karp bipartite vertex cover for combine. */
    int use_hopcroft;

    /** Fraction of IS nodes to remove before reduction (0.0–1.0). */
    double remove_fraction;

    /** Whether to optimize candidates for the ILS. */
    int optimize_candidates;

    /** Remove IS nodes from best individual before recursive reduction. */
    int extract_best_nodes;

    /** Use the adaptive greedy starting solution. */
    int start_greedy_adaptive;

    /** Whether to write progress to stdout. */
    int console_log;
} KamisConfig;

/* ---------------------------------------------------------------------------
 * API functions
 * --------------------------------------------------------------------------- */

/**
 * Fill `config` with sensible default values.
 */
void kamis_default_config(KamisConfig* config);

/**
 * Apply a named preset to `config`.
 *
 * @param preset  One of KAMIS_CONFIG_* constants.
 * @param config  Config to overwrite.
 */
void kamis_apply_preset(int preset, KamisConfig* config);

/**
 * Solve the Maximum Independent Set problem.
 *
 * @param n              Number of vertices.
 * @param xadj           Index pointer array of length n+1 (CSR, 0‑based).
 * @param adjncy         Adjacency array of length xadj[n] (CSR, 0‑based).
 * @param solution       Output – pointer to the array of MIS vertex indices.
 *                       Caller must free with kamis_free_solution().
 * @param solution_size  Output – number of vertices in the solution.
 * @param config         Solver configuration.
 *
 * @return KAMIS_OK on success, or a negative error code.
 */
int kamis_solve_mis(
    int          n,
    const int*   xadj,
    const int*   adjncy,
    int**        solution,
    int*         solution_size,
    const KamisConfig* config
);

/**
 * Free memory allocated by kamis_solve_mis().
 */
void kamis_free_solution(int* solution);

#ifdef __cplusplus
}
#endif

#endif /* KAMIS_H */

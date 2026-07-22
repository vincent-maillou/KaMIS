/**
 * kamis_api.cpp
 * C-compatible wrapper around the KaMIS evolutionary MIS solver.
 *
 * See include/kamis.h for the public API documentation.
 */

#include "kamis.h"

#include <cstring>
#include <memory>
#include <vector>

#include "definitions.h"
#include "data_structure/graph_access.h"
#include "graph_io.h"
#include "mis_config.h"
#include "mis_log.h"
#include "reduction_evolution.h"
#include "population_mis.h"
#include "timer.h"
#include "mis/kernel/branch_and_reduce_algorithm.h"
#include "mis/kernel/ParFastKer/fast_reductions/src/full_reductions.h"

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

namespace {

void copy_config_to_misconfig(const KamisConfig* src, MISConfig& dst) {
    dst.seed                        = src->seed;
    dst.time_limit                  = src->time_limit;
    dst.population_size             = static_cast<unsigned int>(src->population_size);
    dst.imbalance                   = src->imbalance;
    dst.kahip_mode                  = static_cast<unsigned int>(src->kahip_mode);
    dst.fullKernelization           = (src->kernel_mode == KAMIS_KERNEL_FULL);
    dst.repetitions                 = static_cast<unsigned int>(src->repetitions);
    dst.insert_threshold            = static_cast<unsigned int>(src->insert_threshold);
    dst.pool_threshold              = static_cast<unsigned int>(src->pool_threshold);
    dst.pool_renewal_factor         = src->pool_renewal_factor;
    dst.randomize_imbalance         = (src->randomize_imbalance != 0);
    dst.diversify                   = (src->diversify != 0);
    dst.enable_tournament_selection = (src->enable_tournament_selection != 0);
    dst.use_multiway_vc             = (src->use_multiway_vc != 0);
    dst.multiway_blocks             = static_cast<unsigned int>(src->multiway_blocks);
    dst.tournament_size             = static_cast<unsigned int>(src->tournament_size);
    dst.flip_coin                   = src->flip_coin;
    dst.force_k                     = static_cast<unsigned int>(src->force_k);
    dst.force_cand                  = static_cast<unsigned int>(src->force_cand);
    dst.number_of_separators        = static_cast<unsigned int>(src->number_of_separators);
    dst.number_of_partitions        = static_cast<unsigned int>(src->number_of_partitions);
    dst.number_of_k_separators      = static_cast<unsigned int>(src->number_of_k_separators);
    dst.number_of_k_partitions      = static_cast<unsigned int>(src->number_of_k_partitions);
    dst.ils_iterations              = static_cast<unsigned int>(src->ils_iterations);
    dst.use_hopcroft                = (src->use_hopcroft != 0);
    dst.remove_fraction             = src->remove_fraction;
    dst.optimize_candidates         = (src->optimize_candidates != 0);
    dst.extract_best_nodes          = (src->extract_best_nodes != 0);
    dst.start_greedy_adaptive       = (src->start_greedy_adaptive != 0);
    dst.console_log                 = (src->console_log != 0);

    // Output/IO defaults (not exposed in KamisConfig, kept sensible)
    dst.print_repetition            = false;
    dst.print_population            = false;
    dst.print_log                   = false;
    dst.write_graph                 = false;
    dst.check_sorted                = false;
    dst.all_reductions              = true;

    // Convergence defaults
    dst.reduction_threshold         = 350;
    dst.best_limit                  = 0.0;
}

template <typename Reducer>
int run_mis(MISConfig& mis_config, graph_access& G,
            std::vector<int>& solution_out) {
    reduction_evolution<Reducer> evo;
    std::vector<bool> independent_set(G.number_of_nodes(), false);
    std::vector<NodeID> best_nodes;

    evo.perform_mis_search(mis_config, G, independent_set, best_nodes);

    // Collect solution
    solution_out.clear();
    solution_out.reserve(G.number_of_nodes());
    for (NodeID node = 0; node < G.number_of_nodes(); node++) {
        if (independent_set[node]) {
            solution_out.push_back(static_cast<int>(node));
        }
    }

    return KAMIS_OK;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Public C API
// ---------------------------------------------------------------------------

void kamis_default_config(KamisConfig* config) {
    std::memset(config, 0, sizeof(*config));

    config->seed                    = 0;
    config->time_limit              = 1000.0;
    config->population_size         = 50;
    config->imbalance               = 0.03;
    config->kahip_mode              = 0;  // FAST
    config->kernel_mode             = KAMIS_KERNEL_FASTKER;
    config->repetitions             = 50;
    config->insert_threshold        = 150;
    config->pool_threshold          = 250;
    config->pool_renewal_factor     = 10.0;
    config->randomize_imbalance     = 1;
    config->diversify               = 1;
    config->enable_tournament_selection = 1;
    config->use_multiway_vc         = 0;
    config->multiway_blocks         = 64;
    config->tournament_size         = 2;
    config->flip_coin               = 1;
    config->force_k                 = 1;
    config->force_cand              = 4;
    config->number_of_separators    = 10;
    config->number_of_partitions    = 10;
    config->number_of_k_separators  = 10;
    config->number_of_k_partitions  = 10;
    config->ils_iterations          = 15000;
    config->use_hopcroft            = 0;
    config->remove_fraction         = 0.10;
    config->optimize_candidates     = 1;
    config->extract_best_nodes      = 1;
    config->start_greedy_adaptive   = 0;
    config->console_log             = 0;

    // Default kernel style for social graphs
}

void kamis_apply_preset(int preset, KamisConfig* config) {
    kamis_default_config(config);

    switch (preset) {
    case KAMIS_CONFIG_STANDARD:
        // Already set by defaults
        break;
    case KAMIS_CONFIG_SOCIAL:
        config->kahip_mode = 3;  // FASTSOCIAL
        break;
    case KAMIS_CONFIG_FULL_STANDARD:
        config->population_size          = 250;
        config->time_limit               = 36000.0;
        config->number_of_separators     = 30;
        config->number_of_partitions     = 30;
        config->number_of_k_separators   = 30;
        config->number_of_k_partitions   = 30;
        config->flip_coin                = 10;
        config->pool_threshold           = 200;
        break;
    case KAMIS_CONFIG_FULL_SOCIAL:
        config->population_size          = 250;
        config->time_limit               = 36000.0;
        config->number_of_separators     = 30;
        config->number_of_partitions     = 30;
        config->number_of_k_separators   = 30;
        config->number_of_k_partitions   = 30;
        config->flip_coin                = 10;
        config->pool_threshold           = 200;
        config->kahip_mode               = 3;  // FASTSOCIAL
        break;
    default:
        break;
    }
}

int kamis_solve_mis(
    int n,
    const int* xadj,
    const int* adjncy,
    int** solution,
    int* solution_size,
    const KamisConfig* config)
{
    if (n <= 0 || xadj == NULL || adjncy == NULL ||
        solution == NULL || solution_size == NULL || config == NULL) {
        return KAMIS_ERROR_INPUT;
    }

    // Validate CSR structure
    if (xadj[0] != 0) {
        return KAMIS_ERROR_INPUT;
    }

    try {
        // Build the graph
        graph_access G;
        // build_from_metis expects non-const int* but doesn't modify them
        G.build_from_metis(n, const_cast<int*>(xadj), const_cast<int*>(adjncy));

        // Translate configuration
        MISConfig mis_config;
        copy_config_to_misconfig(config, mis_config);

        // Set up logging
        mis_log::instance()->restart_total_timer();
        mis_log::instance()->set_config(mis_config);

        // Run the solver
        std::vector<int> sol;
        int ret = KAMIS_OK;

        if (mis_config.fullKernelization) {
            ret = run_mis<branch_and_reduce_algorithm>(mis_config, G, sol);
        } else {
            ret = run_mis<full_reductions>(mis_config, G, sol);
        }

        if (ret != KAMIS_OK) {
            return ret;
        }

        // Copy solution to output array
        int sz = static_cast<int>(sol.size());
        int* out = static_cast<int*>(std::malloc(sz * sizeof(int)));
        if (out == NULL) {
            return KAMIS_ERROR_MEMORY;
        }

        std::memcpy(out, sol.data(), sz * sizeof(int));
        *solution      = out;
        *solution_size = sz;

        return KAMIS_OK;
    } catch (const std::bad_alloc&) {
        return KAMIS_ERROR_MEMORY;
    } catch (...) {
        return KAMIS_ERROR_INPUT;
    }
}

void kamis_free_solution(int* solution) {
    std::free(solution);
}

/**
 * This is a single-file C++ 11 implementation of the MAX-MIN Ant System
 * algorithm for solving the TSP and ATSP, as described in:
 *
 * Stützle, Thomas, and Holger H. Hoos. "MAX–MIN Ant System." Future generation
 * computer systems 16.8 (2000): 889-914.
 *
 * It is intended mainly for educational puroposes and may not offer the best
 * possible performance.
 *
 * Licensed under terms of MIT license (see LICENSE)
 *
 * Copyright (c) 2018 Rafał Skinderowicz, rafal.skinderowicz@us.edu.pl
 */
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

using namespace std;

std::default_random_engine &get_rng() {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static default_random_engine instance(seed);

    return instance;
}

uint32_t get_random_uint32(uint32_t min, uint32_t max_inclusive) {
    uniform_int_distribution<uint32_t> distribution(min, max_inclusive);
    return distribution(get_rng());
}

double get_random_double(double from = 0.0, uint32_t to_exclusive = 1.0) {
    uniform_real_distribution<double> distribution(from, to_exclusive);
    return distribution(get_rng());
}

struct ProblemInstance {
    uint32_t dimension_;
    bool is_symmetric_ = true;
    vector<double> distance_matrix_;
    vector<vector<uint32_t>> nearest_neighbor_lists_;

    ProblemInstance(uint32_t dimension, const vector<double> &distance_matrix,
                    bool is_symmetric)
        : dimension_(dimension), is_symmetric_(is_symmetric),
          distance_matrix_(distance_matrix) {

        assert(dimension >= 2);
    }

    void initialize_nn_lists(uint32_t nn_list_size) {
        assert(dimension_ > 1);

        nn_list_size = min(dimension_ - 1, nn_list_size);

        nearest_neighbor_lists_.resize(dimension_);
        vector<uint32_t> neighbors(dimension_);
        for (uint32_t i = 0; i < dimension_; ++i) {
            neighbors[i] = i;
        }

        for (uint32_t node = 0; node < dimension_; ++node) {
            sort(neighbors.begin(), neighbors.end(),
                 [this, node](uint32_t a, uint32_t b) {
                     return get_distance(node, a) < get_distance(node, b);
                 });
            assert(get_distance(node, neighbors.at(0)) <=
                   get_distance(node, neighbors.at(1)));

            auto &nn_list = nearest_neighbor_lists_.at(node);
            nn_list.clear();
            nn_list.reserve(nn_list_size);
            uint32_t count = 0;
            for (uint32_t i = 0; count < nn_list_size; ++i) {
                if (neighbors[i] != node) { // node is not its own neighbor
                    nn_list.push_back(neighbors[i]);
                    ++count;
                }
            }
        }
    }

    double get_distance(uint32_t from, uint32_t to) const {
        assert((from < dimension_) && (to < dimension_));

        return distance_matrix_[from * dimension_ + to];
    }

    const vector<uint32_t> &get_nearest_neighbors(uint32_t node) const {
        assert(node < nearest_neighbor_lists_.size());

        return nearest_neighbor_lists_[node];
    }

    double calculate_route_length(const vector<uint32_t> &route) const {
        double distance = 0;
        if (!route.empty()) {
            auto prev_node = route.back();
            for (auto node : route) {
                distance += get_distance(prev_node, node);
                prev_node = node;
            }
        }
        return distance;
    }
};

/**
 * Tries to load a Traveling Salesman Problem (or ATSP) instance in TSPLIB
 * format from file at 'path'. Only the instances with 'EDGE_WEIGHT_TYPE:
 * EUC_2D' or 'EXPLICIT' are supported.
 *
 * Throws runtime_error if the file is in unsupported format or if an error was
 * encountered.
 *
 * Returns the loaded problem instance.
 */
ProblemInstance load_tsplib_instance(const char *path) {
    enum EdgeWeightType { EUC_2D, EXPLICIT };

    ifstream in(path);

    if (!in.is_open()) {
        throw runtime_error(string("Cannot open TSP instance file: ") + path);
    }

    string line;

    uint32_t dimension = 0;
    vector<double> distances;
    EdgeWeightType edge_weight_type{EUC_2D};
    bool is_symmetric = true;

    while (getline(in, line)) {
        cout << "Read line: " << line << endl;
        if (line.find("TYPE") == 0) {
            if (line.find(" TSP") != string::npos) {
                is_symmetric = true;
            } else if (line.find(" ATSP") != string::npos) {
                is_symmetric = false;
            } else {
                throw runtime_error("Unknown problem type");
            }
        } else if (line.find("DIMENSION") != string::npos) {
            istringstream line_in(line.substr(line.find(':') + 1));
            if (!(line_in >> dimension)) {
                throw runtime_error(string("Cannot read instance dimension"));
            }
        } else if (line.find("EDGE_WEIGHT_TYPE") != string::npos) {
            if (line.find(" EUC_2D") != string::npos) {
                edge_weight_type = EUC_2D;
            } else if (line.find(" EXPLICIT") != string::npos) {
                edge_weight_type = EXPLICIT;
            } else {
                throw runtime_error(string("Unsupported edge weight type"));
            }
        } else if (line.find("NODE_COORD_SECTION") != string::npos) {
            vector<pair<double, double>> coords;

            while (getline(in, line)) {
                if (line.find("EOF") == string::npos) {
                    istringstream line_in(line);
                    uint32_t id;
                    pair<double, double> point;
                    line_in >> id >> point.first >> point.second;
                    if (line_in.bad()) {
                        cerr << "Error while reading coordinates";
                    }
                    coords.push_back(point);
                } else {
                    break;
                }
            }

            distances.resize(dimension * dimension, 0);

            for (uint32_t i = 0; i < dimension; ++i) {
                auto from = coords.at(i);

                for (uint32_t j = 0; j < dimension; ++j) {
                    if (i != j) {
                        auto to = coords.at(j);
                        auto dx = to.first - from.first;
                        auto dy = to.second - from.second;
                        double distance = int(sqrt(dx * dx + dy * dy) + 0.5);
                        distances.at(i * dimension + j) = distance;
                    }
                }
            }
        } else if (line.find("EDGE_WEIGHT_SECTION") != string::npos) {
            assert(dimension > 0);
            if (edge_weight_type != EXPLICIT) {
                throw runtime_error("Expected EXPLICIT edge weight type");
            }

            distances.reserve(dimension * dimension);
            while (getline(in, line)) {
                if (line.find("EOF") != string::npos) {
                    break;
                }
                istringstream line_in(line);
                double distance;
                while (line_in >> distance) {
                    distances.push_back(distance);
                }
            }
            assert(distances.size() == dimension * dimension);
        }
    }
    in.close();

    assert(dimension > 2);

    return ProblemInstance(dimension, distances, is_symmetric);
}

struct Ant {
    vector<uint32_t> visited_; // A list of visited nodes, i.e. a route
    vector<uint8_t> is_visited_;
    double cost_ = std::numeric_limits<double>::max();

    void initialize(uint32_t dimension) {
        visited_.clear();
        visited_.reserve(dimension);
        is_visited_.clear();
        is_visited_.resize(dimension, false);
    }

    void visit(uint32_t node) {
        assert(!is_visited_.at(node));

        visited_.push_back(node);
        is_visited_.at(node) = true;
    }

    bool is_visited(uint32_t node) const {
        assert(node < is_visited_.size());

        return is_visited_[node];
    }

    bool all_visited() const {
        return find(is_visited_.begin(), is_visited_.end(), false) ==
               is_visited_.end();
    }
};

struct PheromoneMemory {
    uint32_t dimension_;
    vector<double> pheromone_values_; // For every edge (a,b),
                                      // where 0 <= a, b < dimension_
    double min_pheromone_value_;

    PheromoneMemory(uint32_t dimension, double min_pheromone_value = 0)
        : dimension_(dimension), min_pheromone_value_(min_pheromone_value) {
        pheromone_values_.resize(dimension * dimension, min_pheromone_value);
    }

    double get(uint32_t from, uint32_t to) const {
        assert((from < dimension_) && (to < dimension_));

        return pheromone_values_[from * dimension_ + to];
    }

    void evaporate_from_all(double evaporation_rate,
                            double min_pheromone_value) {

        for (auto &value : pheromone_values_) {
            value = max(min_pheromone_value, value * (1 - evaporation_rate));
        }
    }

    void increase(uint32_t from, uint32_t to, double deposit,
                  double max_pheromone_value, bool is_symmetric) {

        assert((from < dimension_) && (to < dimension_));

        auto &value = pheromone_values_[from * dimension_ + to];
        value = min(max_pheromone_value, value + deposit);

        if (is_symmetric) {
            pheromone_values_[to * dimension_ + from] = value;
        }
    }
};

/**
 * This creates a solution using nearest neighbor heuristic that always selects
 * a clostest of the (yet) unvisited nodes (cities).
 */
Ant create_solution_nn(const ProblemInstance &instance,
                       uint32_t start_node = 0) {
    Ant ant;
    ant.initialize(instance.dimension_);

    uint32_t current_node = start_node;
    ant.visit(current_node);

    for (uint32_t i = 1; i < instance.dimension_; ++i) {
        uint32_t next_node = current_node;
        const auto &candidates = instance.get_nearest_neighbors(current_node);

        for (auto node : candidates) {
            if (!ant.is_visited(node)) {
                next_node = node;
                break;
            }
        }

        if (next_node == current_node) { // All closest nodes were visited,
                                         // we have to check the rest
            double min_distance = numeric_limits<double>::max();

            for (uint32_t node = 0; node < instance.dimension_; ++node) {
                if (!ant.is_visited(node)) {
                    auto distance = instance.get_distance(current_node, node);
                    if (distance < min_distance) {
                        min_distance = distance;
                        next_node = node;
                    }
                }
            }
        }

        assert(next_node != current_node);

        ant.visit(next_node);
        current_node = next_node;
    }
    return ant;
}

/**
 * This are based on the article mentioned.
 */
struct MMASParameters {
    double rho_ = 0.98;
    uint32_t ants_count_ = 10;
    double beta_ = 2;
    uint32_t cand_list_size_ = 15;
    double p_best_ = 0.05; // Prob. that the constructed sol. will contain only
                           // the edges with the highest pheromone values

    double get_evaporation_rate() const { return 1 - rho_; }
};

const uint32_t MaxCandListSize = 64;

/*
 * Moves 'ant' from its current node to a next one chosen according to the MMAS
 * rules. Returns the selected node.
 */
uint32_t move_ant_mmas(const ProblemInstance &instance,
                       const PheromoneMemory &pheromone,
                       const vector<double> &heuristic, Ant &ant) {
    assert(!ant.visited_.empty());

    const auto dimension = instance.dimension_;
    const auto current_node = ant.visited_.back();
    const uint32_t offset = current_node * dimension;

    // A list of the nearest unvisited neighbors of 'current_node':
    uint32_t cand_list[MaxCandListSize];
    uint32_t cand_list_size = 0;
    for (auto node : instance.get_nearest_neighbors(current_node)) {
        if (!ant.is_visited(node)) {
            cand_list[cand_list_size] = node;
            ++cand_list_size;
        }
    }

    uint32_t chosen_node = current_node;

    if (cand_list_size > 0) { // Select from the closest nodes
        double products_prefix_sum[MaxCandListSize] = {0};
        double total = 0;
        for (uint32_t i = 0; i < cand_list_size; ++i) {
            const auto node = cand_list[i];
            const auto product =
                pheromone.get(current_node, node) * heuristic[offset + node];
            total += product;
            products_prefix_sum[i] = total;
        }

        chosen_node = cand_list[cand_list_size - 1];
        const auto r = get_random_double() * total;
        for (uint32_t i = 0; i < cand_list_size; ++i) {
            if (r < products_prefix_sum[i]) {
                chosen_node = cand_list[i];
                break;
            }
        }
    } else { // Select from the rest of the unvisited nodes the one with the
             // maximum product of pheromone and heuristic
        double max_product = 0;

        for (uint32_t node = 0u; node < dimension; ++node) {
            if (!ant.is_visited(node)) {
                const auto product = pheromone.get(current_node, node) *
                                     heuristic[offset + node];
                if (product > max_product) {
                    max_product = product;
                    chosen_node = node;
                }
            }
        }
    }
    assert(chosen_node != current_node);

    ant.visit(chosen_node);
    return chosen_node;
}

/**
 * This is based on Eq. 11 from the original MAX-MIN paper:
 *
 * Stützle, Thomas, and Holger H. Hoos. "MAX–MIN ant system." Future generation
 * computer systems 16.8 (2000): 889-914.
 */
pair<double, double> calc_trail_limits_mmas(const MMASParameters &params,
                                            uint32_t instance_dimension,
                                            double solution_cost) {
    const auto tau_max = 1 / (solution_cost * (1. - params.rho_));
    const auto avg = instance_dimension / 2.;
    const auto p = pow(params.p_best_, 1. / instance_dimension);
    const auto tau_min = min(tau_max, tau_max * (1 - p) / ((avg - 1) * p));
    return make_pair(tau_min, tau_max);
}

/**
 * Runs the MMAS for the given number of iterations.
 * Returns the best solution (ant) found.
 */
Ant run_mmas(const ProblemInstance &instance, const MMASParameters &params,
             uint32_t iterations) {

    const auto greedy_sol = create_solution_nn(instance);
    const auto greedy_cost =
        instance.calculate_route_length(greedy_sol.visited_);
    const auto initial_limits =
        calc_trail_limits_mmas(params, instance.dimension_, greedy_cost);
    auto min_pheromone = initial_limits.first;
    auto max_pheromone = initial_limits.second;

    PheromoneMemory pheromone(instance.dimension_, max_pheromone);

    vector<double> heuristic;
    heuristic.reserve(instance.dimension_ * instance.dimension_);

    for (auto distance : instance.distance_matrix_) {
        heuristic.push_back(1 / pow(distance, params.beta_));
    }

    vector<Ant> ants(params.ants_count_);
    Ant best_ant;

    for (uint32_t iteration = 0; iteration < iterations; ++iteration) {
        for (auto &ant : ants) {
            ant.initialize(instance.dimension_);
            auto start_node = get_random_uint32(0, instance.dimension_ - 1);
            ant.visit(start_node);
            for (uint32_t j = 1; j < instance.dimension_; ++j) {
                move_ant_mmas(instance, pheromone, heuristic, ant);
            }
            ant.cost_ = instance.calculate_route_length(ant.visited_);
        }

        auto &iteration_best = ants.front();
        bool new_best_found = false;

        // Have we found an improved solution?
        for (auto &ant : ants) {
            if (ant.cost_ < best_ant.cost_) {
                best_ant = ant;
                new_best_found = true;

                cout << "New best solution found with the cost: "
                     << best_ant.cost_ << " at iteration " << iteration << endl;
            }
            if (ant.cost_ < iteration_best.cost_) {
                iteration_best = ant;
            }
        }

        if (new_best_found) {
            auto limits = calc_trail_limits_mmas(params, instance.dimension_,
                                                 best_ant.cost_);
            min_pheromone = limits.first;
            max_pheromone = limits.second;
        }

        pheromone.evaporate_from_all(params.get_evaporation_rate(),
                                     min_pheromone);

        // Deposit pheromone on the edges belonging to the iteration best ant
        const auto &update_ant = iteration_best;
        const double deposit = 1.0 / update_ant.cost_;
        auto prev_node = update_ant.visited_.back();
        for (auto node : update_ant.visited_) {
            // The global update of the pheromone trails
            pheromone.increase(prev_node, node, deposit, max_pheromone,
                               instance.is_symmetric_);
            prev_node = node;
        }
    }
    return best_ant;
}

int main(int argc, char *argv[]) {
    string path = "kroA100.tsp";
    if (argc >= 2) {
        path = argv[1];
    }
    try {
        MMASParameters params;
        auto instance = load_tsplib_instance(path.c_str());
        instance.initialize_nn_lists(params.cand_list_size_);
        params.ants_count_ = instance.dimension_;

        // Create 1 million solutions
        run_mmas(instance, params, 1000000 / params.ants_count_);
    } catch (runtime_error e) {
        cout << "An error has occurred: " << e.what() << endl;
    }
    return 0;
}
